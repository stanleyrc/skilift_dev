"""Unit tests for gos_sc_splice_confidence.py (synthetic fixtures only)."""
import csv
import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gos_sc_splice_confidence as C  # noqa: E402

COLS = ["gene", "cluster_id", "verdict", "cis_dna", "protein_consequence", "annotation_class", "effect", "frame", "nmd_50nt",
        "protein_len", "ref_protein_len", "pfam_affected", "terminal_intron", "patient_qb_q", "patient_given_batch_q",
        "bwh75_batch_dpsi", "tn_k", "tn_dpsi", "tn_lo", "tn_hi", "tn_q", "tn_perm_p", "tn_I2", "tn_same_dir", "glmm_p_tn",
        "st_best", "st_k", "st_dpsi", "st_lo", "st_hi", "st_q", "st_perm_p", "st_I2", "st_same_dir", "glmm_p_st"]


def row(**kw):
    r = {c: "" for c in COLS}
    r.update({k: str(v) for k, v in kw.items()})
    return r


ROWS = [
    row(gene="GA", cluster_id="c_tn", verdict="keep: within-patient differential, malignant vs non-malignant cells, consistent across patients",
        tn_k=5, tn_dpsi=-0.6, tn_lo=-0.7, tn_hi=-0.5, tn_q=1e-20, tn_same_dir=5, terminal_intron="first",
        protein_consequence="known alternative protein-coding isoform"),
    row(gene="GB", cluster_id="c_st", verdict="keep: within-patient differential between malignant cell states, consistent across patients",
        st_best="NPC", st_k=7, st_dpsi=0.3, st_lo=0.15, st_hi=0.45, st_q=0.008, st_same_dir=7, terminal_intron="internal"),
    row(gene="GC", cluster_id="c_cis", verdict="keep: patient-specific junction with a somatic splice-site variant in cis (P1 GC c.1-1G>C)",
        cis_dna="P1 GC c.1-1G>C", protein_consequence="novel junction: frameshift / PTC, predicted NMD (loss of productive mRNA)", nmd_50nt="TRUE"),
    row(gene="GD", cluster_id="c_desc", verdict="descriptive only: tumour-to-tumour difference (survives cell-level and batch tests; no within-patient contrast)"),
    row(gene="GE", cluster_id="c_tent", verdict="tentative: within-patient contrast in one patient only", tn_k=1, tn_dpsi=-0.8, tn_same_dir=1),
    row(gene="GF", cluster_id="c_rescue", verdict="drop: differs between library batches within BWH75", tn_k=4, tn_same_dir=4, tn_q=0.07,
        tn_dpsi=-0.35, bwh75_batch_dpsi=-0.33),
    row(gene="GG", cluster_id="c_b75", verdict="drop: differs between library batches within BWH75", tn_k=4, tn_same_dir=2, tn_q=1, bwh75_batch_dpsi=0.7),
    row(gene="GH", cluster_id="c_ns", verdict="drop: not significant once cells (not reads) are the unit"),
    row(gene="GI", cluster_id="c_batch", verdict="drop: between-patient difference explained by library batch"),
    row(gene="GJ", cluster_id="c_par", verdict="drop: paralog / repeat gene family (mapping artefact)"),
]


def write_table(path, rows=ROWS):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


class Classify(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.tsv = os.path.join(self.d, "events_annotated.tsv")
        write_table(self.tsv)
        self.t = C.read_table(self.tsv)

    def test_tiers(self):
        tiers = {k: v["tier"] for k, v in self.t.items()}
        self.assertEqual(tiers, {"c_tn": "consistent", "c_st": "consistent", "c_cis": "consistent", "c_desc": "suggestive",
                                 "c_tent": "suggestive", "c_rescue": "suggestive", "c_b75": "low", "c_ns": "low",
                                 "c_batch": "low", "c_par": "low"})

    def test_contrast_labels_and_detail(self):
        self.assertIn("non-malignant", self.t["c_tn"]["detail"]["contrast"])
        self.assertIn("oligodendrocytes", self.t["c_tn"]["reason"])
        self.assertIn("5/5 patients", self.t["c_tn"]["reason"])
        self.assertTrue(self.t["c_tn"]["detail"]["first_intron"])
        self.assertEqual(self.t["c_tn"]["detail"]["consequence"], "known isoform switch")
        self.assertEqual(self.t["c_st"]["detail"]["contrast"], "NPC-like vs other malignant cells")
        self.assertIn("DNA", self.t["c_cis"]["reason"])
        self.assertEqual(self.t["c_cis"]["detail"]["consequence"], "novel, NMD-type (frameshift / PTC)")
        self.assertIs(self.t["c_cis"]["detail"]["nmd"], True)

    def test_drop_reasons(self):
        self.assertIn("cells (not reads)", self.t["c_ns"]["reason"])
        self.assertIn("library batch", self.t["c_batch"]["reason"])
        self.assertIn("BWH75", self.t["c_b75"]["reason"])
        self.assertIn("paralog", self.t["c_par"]["reason"])
        self.assertIn("batch-sensitive", self.t["c_rescue"]["reason"])

    def test_rescue_threshold(self):
        t = C.read_table(self.tsv, rescue_min_patients=5)
        self.assertEqual(t["c_rescue"]["tier"], "low")


class Files(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.tsv = os.path.join(self.d, "t.tsv")
        write_table(self.tsv)
        self.data = os.path.join(self.d, "data")
        self.out = os.path.join(self.d, "out")
        p = {"format": "gos-sc-splicing/1", "patient": "P1",
             "findings": [{"id": "patient_specific:c_cis:P1", "kind": "patient_specific", "gene": "GC", "cluster_id": "c_cis", "text": "GC cryptic"},
                          {"id": "known_variant:V3:P1", "kind": "known_variant", "gene": "GV", "cluster_id": None, "variant_id": "V3",
                           "text": "V3 carriers", "numbers": {"n_cells_alt": 7, "n_cells": 9}},
                          {"id": "clone_differential:c_ns:P1", "kind": "clone_differential", "gene": "GH", "cluster_id": "c_ns", "text": "GH clone"}],
             "clusters": [{"id": "c_cis", "gene": "GC"}, {"id": "c_new", "gene": "GZ"}],
             "events": [{"id": "e1", "cluster_id": "c_cis", "gene": "GC", "summary": "GC event", "finding_id": "patient_specific:c_cis:P1"},
                        {"id": "e2", "cluster_id": "c_new", "gene": "GZ", "summary": "GZ event"},
                        {"id": "e3", "cluster_id": None, "gene": "GV", "summary": "V3 event", "finding_id": "known_variant:V3:P1"}]}
        coh = {"format": "gos-sc-splicing-cohort/1", "findings": p["findings"], "clusters": [{"id": c} for c in ("c_tn", "c_desc", "c_ns")]}
        sup = {"patients": {"P1": {"GC": {"dna": [{}], "rna": [{"gene": "GC", "text": "GC cryptic"}], "support": "DNA+RNA"},
                                   "GZ": {"dna": [], "rna": [{"gene": "GZ", "text": "GZ event"}], "support": "RNA"},
                                   "GQ": {"dna": [{}], "rna": [], "support": "DNA"}}}}
        vig = {"patients": {"P1": {"items": [{"kind": "splice", "evidence": [{"kind": "splice", "finding": "clone_differential:c_ns:P1"}]},
                                             {"kind": "snv", "evidence": []}]}}}
        os.makedirs(os.path.join(self.data, "P1", "rna"))
        os.makedirs(os.path.join(self.data, "_cohort", "rna"))
        for path, doc in ((("P1", "rna", "splicing.json"), p), (("_cohort", "rna", "splicing.json"), coh),
                          (("_cohort", "splice_support.json"), sup), (("_cohort", "vignettes.json"), vig)):
            with open(os.path.join(self.data, *path), "w") as f:
                json.dump(doc, f)
        self.before = {q: open(os.path.join(self.data, q)).read() for q in ("P1/rna/splicing.json", "_cohort/vignettes.json")}

    def test_staged_run(self):
        C.main(["--data", self.data, "--deepdive", self.tsv, "--out", self.out])
        # the data dir is only read
        for q, s in self.before.items():
            self.assertEqual(open(os.path.join(self.data, q)).read(), s)
        p = json.load(open(os.path.join(self.out, "P1", "rna", "splicing.json")))
        f = {x["id"]: x for x in p["findings"]}
        self.assertEqual(f["patient_specific:c_cis:P1"]["confidence"], "consistent")
        self.assertEqual(f["known_variant:V3:P1"]["confidence"], "consistent")
        self.assertIn("7/9", f["known_variant:V3:P1"]["confidence_reason"])
        self.assertEqual(f["clone_differential:c_ns:P1"]["confidence"], "low")
        c = {x["id"]: x["confidence"] for x in p["clusters"]}
        self.assertEqual(c, {"c_cis": "consistent", "c_new": "low"})
        e = {x["id"]: x for x in p["events"]}
        self.assertEqual(e["e1"]["confidence"], "consistent")
        self.assertNotIn("confidence_detail", e["e1"])
        self.assertEqual(e["e2"]["confidence_reason"], C.NOT_RETESTED)
        self.assertEqual(e["e3"]["confidence"], "consistent")
        self.assertEqual(p["confidence_meta"]["counts"]["findings"], {"consistent": 2, "suggestive": 0, "low": 1})
        coh = json.load(open(os.path.join(self.out, "_cohort", "rna", "splicing.json")))
        self.assertEqual([x["confidence"] for x in coh["clusters"]], ["consistent", "suggestive", "low"])
        sup = json.load(open(os.path.join(self.out, "_cohort", "splice_support.json")))["patients"]["P1"]
        self.assertEqual(sup["GC"]["rna_confidence"], "consistent")
        self.assertEqual(sup["GZ"]["rna"][0]["confidence"], "low")
        self.assertIsNone(sup["GQ"]["rna_confidence"])
        vig = json.load(open(os.path.join(self.out, "_cohort", "vignettes.json")))["patients"]["P1"]["items"]
        self.assertEqual(vig[0]["confidence"], "low")
        self.assertNotIn("confidence", vig[1])

    def test_idempotent(self):
        C.main(["--data", self.data, "--deepdive", self.tsv, "--out", self.out])
        first = open(os.path.join(self.out, "P1", "rna", "splicing.json")).read()
        C.main(["--data", self.out, "--deepdive", self.tsv, "--out", self.out, "--patients", "P1"])
        a, b = json.loads(first), json.load(open(os.path.join(self.out, "P1", "rna", "splicing.json")))
        for d in (a, b):
            d["confidence_meta"].pop("created")
        self.assertEqual(a, b)

    def test_requires_target(self):
        with self.assertRaises(SystemExit):
            C.main(["--data", self.data, "--deepdive", self.tsv])


if __name__ == "__main__":
    unittest.main()
