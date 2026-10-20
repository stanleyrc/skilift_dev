"""Unit tests for gos_sc_signature_versions.py (synthetic fixtures only; no patient data)."""
import gzip
import json
import os
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gos_sc_signature_versions as V  # noqa: E402


def synth_sigs(path):
    """Three synthetic signatures: A = C>T at CpG, B = flat, C = T>C."""
    A = np.array([5.0 if (c[2:5] == "C>T" and c[6] == "G") else 0.01 for c in V.SBS96])
    B = np.ones(96)
    C = np.array([3.0 if c[2:5] == "T>C" else 0.01 for c in V.SBS96])
    with open(path, "w") as fh:
        fh.write("MutationType\tSIGA\tSBS5\tSIGC\n")
        for i, c in enumerate(V.SBS96):
            fh.write(f"{c}\t{A[i] / A.sum()}\t{B[i] / 96}\t{C[i] / C.sum()}\n")
    return np.column_stack([A / A.sum(), B / 96, C / C.sum()])


class Core(unittest.TestCase):
    def test_channel_order_and_context(self):
        self.assertEqual(len(V.SBS96), 96)
        self.assertEqual(V.SBS96[0], "A[C>A]A")
        self.assertEqual(V.sbs96_context("ACG", "C", "T"), "A[C>T]G")
        # purine reference is reverse-complemented
        self.assertEqual(V.sbs96_context("CGT", "G", "A"), "A[C>T]G")
        self.assertIsNone(V.sbs96_context("ANG", "N", "T"))
        self.assertIsNone(V.sbs96_context("BCG", "C", "T"))  # IUPAC code in the reference

    def test_selection_recovers_two_signatures(self):
        with tempfile.TemporaryDirectory() as d:
            S = synth_sigs(os.path.join(d, "s.tsv"))
            names, S2 = V.read_signature_table(os.path.join(d, "s.tsv"))
            rng = np.random.default_rng(0)
            ct = rng.multinomial(3000, 0.6 * S[:, 0] + 0.4 * S[:, 2]).astype(float)
            cat = V.choose_catalogue(ct, S2, names, n_boot_select=30)
            self.assertEqual(sorted(cat["signatures"]), ["SIGA", "SIGC"])
            fit = V.refit_bootstrap(ct, S2, names, cat["signatures"], n_boot=100)
            fr = {a["signature"]: a["fraction"] for a in fit["activities"]}
            self.assertAlmostEqual(fr["SIGA"], 0.6, delta=0.05)
            for a in fit["activities"]:
                self.assertLessEqual(a["lo"], a["fraction"]); self.assertGreaterEqual(a["hi"], a["fraction"])

    def test_too_few(self):
        with tempfile.TemporaryDirectory() as d:
            synth_sigs(os.path.join(d, "s.tsv"))
            names, S = V.read_signature_table(os.path.join(d, "s.tsv"))
            cat = V.choose_catalogue(np.ones(96) * 0 + np.eye(96)[0] * 10, S, names)
            self.assertEqual(cat["status"], "not_available")

    def test_posterior_rows_and_argmax(self):
        with tempfile.TemporaryDirectory() as d:
            synth_sigs(os.path.join(d, "s.tsv"))
            names, S = V.read_signature_table(os.path.join(d, "s.tsv"))
            post = V.channel_posterior(S, [0, 2], np.array([100.0, 100.0]))
            self.assertTrue(np.allclose(post.sum(axis=1), 1))
            self.assertGreater(post[V.CH_INDEX["A[C>T]G"], 0], 0.9)
            self.assertGreater(post[V.CH_INDEX["A[T>C]A"], 1], 0.9)

    def test_callable_singleton(self):
        r = {"alt": 4, "depth": 10, "adf_alt": 2, "adr_alt": 2, "n_cov_tumor": 30}
        self.assertTrue(V.is_callable_singleton(r))
        self.assertFalse(V.is_callable_singleton(dict(r, adr_alt=0)))
        self.assertFalse(V.is_callable_singleton(dict(r, alt=2)))
        self.assertFalse(V.is_callable_singleton(dict(r, depth=40)))


def write_fixture(d):
    """Synthetic patient SYN: tree sites in snv_matrix.sites.json, VCF evidence with singletons, a tiny FASTA."""
    rng = np.random.default_rng(1)
    os.makedirs(os.path.join(d, "data", "SYN"))
    S = synth_sigs(os.path.join(d, "s.tsv"))
    # FASTA: one contig of random sequence; singletons placed at C or T positions
    seq = "".join(rng.choice(list("ACGT"), 20000))
    with open(os.path.join(d, "ref.fa"), "w") as fh:
        fh.write(">chrS\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")
    with open(os.path.join(d, "ref.fa.fai"), "w") as fh:
        fh.write(f"chrS\t{len(seq)}\t6\t60\t61\n")
    ids, cat, ctx, alt_cells, cs, mc = [], [], [], [], [], []
    p_tr = 0.7 * S[:, 0] + 0.3 * S[:, 1]
    p_sub = 0.3 * S[:, 0] + 0.7 * S[:, 2]
    for k in range(1200):
        truncal = k < 600
        ch = V.SBS96[rng.choice(96, p=p_tr if truncal else p_sub)]
        ids.append(f"chrT_{k + 1}_C_T"); ctx.append(ch)
        cat.append("truncal" if truncal else ("subclonal" if k < 1000 else "private"))
        alt_cells.append(10 if k < 1000 else 2); cs.append(0.9); mc.append(2.0 if k % 2 else 0)
    cols = {"id": {"type": "str", "values": ids},
            "category": {"type": "dict", "dict": ["truncal", "subclonal", "private"], "codes": [["truncal", "subclonal", "private"].index(c) for c in cat]},
            "clade_score": {"type": "num", "values": cs}, "alt_cells": {"type": "num", "values": alt_cells},
            "map_confidence": {"type": "num", "values": mc}, "context": {"type": "dict", "dict": V.SBS96, "codes": [V.SBS96.index(c) for c in ctx]},
            "cellphy_input": {"type": "bool", "values": [1] * len(ids)}}
    json.dump({"n_sites": len(ids), "cells": ["c1", "c2", "n1"], "columns": cols}, open(os.path.join(d, "data", "SYN", "snv_matrix.sites.json"), "w"))
    rows = []
    for i, m in enumerate(ids):
        rows.append([m, "tree", ".", "30", "0.1", "somatic", 10, 0, 0, 40, 2, "", "", 20, 60, 10, 10, 20, 20, 0.33])
    n_s = 0
    for pos in range(2, len(seq) - 2):
        if n_s >= 400:
            break
        b = seq[pos - 1]
        if b in "CT" and rng.random() < 0.2:
            alt = "T" if b == "C" else "C"
            role = "normal" if n_s % 10 == 0 else "tree_tumor"
            rows.append([f"chrS_{pos}_{b}_{alt}", "singleton", ".", "5", "1", "somatic", 0 if role == "normal" else 1, 1 if role == "normal" else 0, 0,
                         30, 2, "n1" if role == "normal" else "c1", role, 4, 10, 2, 2, 3, 3, 0.4])
            n_s += 1
    with gzip.open(os.path.join(d, "vcf_sites.tsv.gz"), "wt") as fh:
        fh.write("\t".join(["mutation", "kind", "filter", "qual", "germline_q", "branch_filter", "n_alt_tumor", "n_alt_normal", "n_alt_other",
                            "n_cov_tumor", "n_cov_normal", "carrier", "carrier_role", "alt", "depth", "adf_alt", "adr_alt", "adf_ref", "adr_ref", "med_vaf"]) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    with open(os.path.join(d, "art.tsv"), "w") as fh:
        fh.write("MutationType\tART\n")
        for c in V.SBS96:
            fh.write(f"{c}\t{1.0 if c[2:5] == 'C>A' else 0.01}\n")


class Patient(unittest.TestCase):
    def test_end_to_end(self):
        with tempfile.TemporaryDirectory() as d:
            write_fixture(d)
            out = V.run_patient("SYN", os.path.join(d, "data"), os.path.join(d, "vcf_sites.tsv.gz"), os.path.join(d, "s.tsv"),
                                [f"ART={os.path.join(d, 'art.tsv')}"], os.path.join(d, "ref.fa"), n_boot=50, n_boot_select=20)
            vs = {v["id"]: v for v in out["versions"]}
            self.assertEqual(len(vs), len(V.VERSIONS))
            self.assertEqual(vs["v2_current"]["n"], 1200)
            self.assertEqual(vs["v6_subclonal_no_singletons"]["n"], 600)
            self.assertEqual(vs["v7_ge3cells"]["n"], 1000)          # private sites have 2 alt cells
            self.assertEqual(vs["v8_qc_pass"]["n"], 600)            # map_confidence >= 1 on odd sites
            self.assertEqual(out["group_counts"]["singleton"] + out["group_counts"]["normal_singleton"], 400)
            self.assertGreater(vs["v4_tree_cells"]["n_singleton"], 300)
            self.assertIn("ART", vs["v9_artifact_cols"]["sets"][0]["signatures"])
            self.assertEqual(len(vs["v1_legacy"]["sets"]), 3)
            # flows conserve SNV counts: v2 -> v6 table counts every v2 SNV once
            f = [x for x in out["flows"]["all"] if x["from"] == "v2_current" and x["to"] == "v6_subclonal_no_singletons"][0]
            self.assertEqual(sum(r[2] for r in f["rows"]), 1200)
            self.assertEqual(len(out["tree_sites"]["id"]), 1200)
            self.assertEqual(len(out["tree_sites"]["assign"]["v2_current"]), 1200)
            g = out["diagnostics"]["groups"]["singleton"]
            self.assertIn("SBS5", g["cosine"])
            self.assertIn("read_end_distance", g["reads"])
            json.dumps(out)  # serialisable
            st = os.path.join(d, "sites.tsv.gz")
            V.write_site_table(st, *V.build_snvs("SYN", os.path.join(d, "data"), os.path.join(d, "vcf_sites.tsv.gz"), os.path.join(d, "ref.fa"))[:1],
                               {v["id"]: {} for v in V.VERSIONS}, [v["id"] for v in V.VERSIONS])
            self.assertEqual(sum(1 for _ in gzip.open(st, "rt")), 1)  # header only when nothing is assigned
            coh = V.cohort({"SYN": out, "SYN2": out})
            self.assertEqual(coh["patients"], ["SYN", "SYN2"])

    def test_no_vcf_evidence(self):
        with tempfile.TemporaryDirectory() as d:
            write_fixture(d)
            out = V.run_patient("SYN", os.path.join(d, "data"), "", os.path.join(d, "s.tsv"), [], "", n_boot=20, n_boot_select=10)
            vs = {v["id"]: v for v in out["versions"]}
            self.assertIn("not available", out["sources"]["singletons"])
            self.assertEqual(vs["v9_artifact_cols"]["status"], "not_available")
            self.assertEqual(vs["v3_ge2cells"]["n"], 0)


if __name__ == "__main__":
    unittest.main()
