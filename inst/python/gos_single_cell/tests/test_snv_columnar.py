"""Synthetic tests of the columnar SNV site sets and the CN-mapped driver rows (no patient data).
python3 -m unittest test_snv_columnar"""
import gzip
import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gos_sc_snv_columnar as sc  # noqa: E402


def snv_fixture():
    return {"schemaVersion": 1, "format": "compact",
            "variants": [{"id": "chr1_100_A_G", "category": "truncal", "anchor": "c1|c3", "node": 5, "driver": False, "gene": "G1", "clade_score": 0.9},
                         {"id": "chr2_200_C_T", "category": "subclonal", "anchor": "c2|c3", "node": 7, "driver": True, "gene": "G2", "clade_score": None},
                         {"id": "chr3_300_G_A", "category": "private", "anchor": "c1", "node": 1, "driver": False, "gene": None, "clade_score": 0.2}],
            "cells": {"c1": [[0, 3, 2, 1], [2, 0, 4, None]], "c2": [[1, 5, 0, 0]], "c3": [[0, 1, 1, None], [1, 2, 70000, 1]]}}


class Columns(unittest.TestCase):
    def test_roundtrip_types(self):
        for vals in (["a", "b", "a", None, "a"], [1, 2.5, None], [True, False, None], [{"x": 1}, None]):
            col = sc.encode_column(vals)
            self.assertEqual(sc.decode_column(col), vals)
        self.assertEqual(sc.encode_column(["a", "a", "b", "a"])["type"], "dict")
        self.assertEqual(sc.encode_column(["a", "b", "c"])["type"], "str")

    def test_matrix_roundtrip_and_alignment(self):
        entries = [[(0, 3, 2, 1), (2, 1, 1, -1)], [], [(1, 5, 0, 0)]]
        blob, offsets, head = sc.encode_matrix(entries)
        self.assertEqual(head % 8, 0)
        self.assertEqual(len(blob), head + 8 * 3)
        self.assertEqual(sc.decode_matrix(blob, 3), [[(0, 3, 2, 1), (2, 1, 1, -1)], [], [(1, 5, 0, 0)]])

    def test_pl_byte(self):
        self.assertEqual(sc.pl_byte(0, 60, 60), 1)            # confident REF -> P(ALT) ~ 0
        self.assertEqual(sc.pl_byte(60, 30, 0), 250)          # confident ALT (hom) -> ~1
        self.assertEqual(sc.pl_byte(60, 0, 60), 250)          # confident het -> ~1
        self.assertEqual(sc.pl_byte(0, 0, 0), sc.PL_FLAT)     # flat PL = no reads: no evidence, not 2/3
        self.assertEqual(sc.pl_byte(None, 3, 3), sc.PL_FLAT)
        self.assertEqual(sc.pl_byte(float("nan"), 3, 3), sc.PL_FLAT)
        mid = sc.pl_byte(3, 0, 20)                              # uncertain: strictly between
        self.assertTrue(100 < mid < 250)
        # per-site prior: one REF read (PL 0,3,30) is not grey at a site almost nobody carries, but is at a truncal site
        self.assertLess(sc.pl_byte(0, 3, 30, q=0.02), 10)
        self.assertGreater(sc.pl_byte(0, 3, 30), 80)            # flat prior: visible grey haze
        self.assertGreater(sc.pl_byte(0, 3, 30, q=0.95), sc.pl_byte(0, 3, 30, q=0.02))
        self.assertEqual(sc.pl_byte(0, 0, 0, q=0.5), sc.PL_FLAT)
        self.assertEqual(sc.pl_byte(10, 10, 0) , 1 + round(249 * (0.1 + 1) / (0.1 + 0.1 + 1)))

    def test_pl_roundtrip_and_tree_set(self):
        entries = [[(0, 3, 2, 1, 250), (2, 1, 1, -1)], [(1, 5, 0, 0, 1)]]
        blob, _, _ = sc.encode_matrix(entries)
        self.assertEqual(sc.decode_matrix(blob, 2, with_pl=True), [[(0, 3, 2, 1, 250), (2, 1, 1, -1, 0)], [(1, 5, 0, 0, 1)]])
        self.assertEqual(sc.decode_matrix(blob, 2)[0][0], (0, 3, 2, 1))
        snv = {"variants": [{"id": "chr1_1_A_C"}], "cells": {"c1": [[0, 2, 0, 0, 1]], "c2": [[0, None, None, None, 255]], "c3": [[0, 1, 1, 1]]}}
        v, e, cells = sc.tree_set(snv)
        self.assertEqual(e[0], [(0, 2, 0, 0, 1), (1, 65535, 65535, -1, 255), (2, 1, 1, 1)])

    def test_clipping(self):
        blob, _, _ = sc.encode_matrix([[(0, 1, 70000, 1)]])
        self.assertEqual(sc.decode_matrix(blob, 1)[0][0][2], 65534)

    def test_call_without_counts(self):
        snv = {"variants": [{"id": "chr1_1_A_C"}], "cells": {"c1": [[0, None, None, 1]], "c2": [[0, None, None, None]]}}
        v, e, cells = sc.tree_set(snv)
        self.assertEqual(e[0], [(0, 65535, 65535, 1)])
        blob, _, _ = sc.encode_matrix(e)
        self.assertEqual(sc.decode_matrix(blob, 1)[0], [(0, 65535, 65535, 1)])


class TreeSet(unittest.TestCase):
    def test_tree_set_groups_and_ranges(self):
        variants, entries, cells = sc.tree_set(snv_fixture())
        self.assertEqual(cells, ["c1", "c2", "c3"])
        self.assertEqual(entries[0], [(0, 3, 2, 1), (2, 1, 1, -1)])
        with tempfile.TemporaryDirectory() as d:
            t = sc.write_set(d, "snv_matrix", {"patient": "SYN"}, variants, entries, cells)
            blob = open(os.path.join(d, "snv_matrix.cells.bin"), "rb").read()
            tab = json.load(open(os.path.join(d, "snv_matrix.sites.json")))
            self.assertEqual(tab["format"], sc.FORMAT)
            ids = sc.decode_column(tab["columns"]["id"])
            # sorted by anchor: c1 < c1|c3 < c2|c3
            self.assertEqual(ids, ["chr3_300_G_A", "chr1_100_A_G", "chr2_200_C_T"])
            dec = sc.decode_matrix(blob, 3)
            self.assertEqual(dec[1], [(0, 3, 2, 1), (2, 1, 1, -1)])
            # each group's byte range holds exactly its sites' entries
            for g in t["matrix"]["groups"]:
                n = sum(len(dec[i]) for i in range(g["first_site"], g["first_site"] + g["n_sites"]))
                self.assertEqual(g["byte_end"] - g["byte_start"], 8 * n)
            self.assertEqual(sc.decode_column(tab["columns"]["clade_score"]), [0.2, 0.9, None])


def cn_row(mid, gene, ef, status="subclonal", **kw):
    r = {"mutation": mid, "gene": gene, "func": "exonic", "exonic_func": ef, "protein": "p.X1Y", "status": status, "node": "n9", "anchor": "c1|c2",
         "clade_cells": "2", "posterior": "0.95", "confidence": "high", "multiplicity": "one copy", "lost_in": "", "lost_cells": "0", "modal_cn": "3",
         "t_alt_cells": "2", "alt_in_clade": "2", "alt_out_clade": "0", "d_iid": "12.5"}
    r.update(kw)
    return r


class CnMapped(unittest.TestCase):
    def test_driver_flag_needs_cancer_gene_and_protein_change(self):
        genes = {"TP53"}
        self.assertTrue(sc.cn_variant(cn_row("chr17_1_C_T", "TP53", "stopgain"), genes)["driver"])
        self.assertFalse(sc.cn_variant(cn_row("chr17_2_C_T", "TP53", "synonymous_SNV"), genes)["driver"])
        self.assertFalse(sc.cn_variant(cn_row("chr15_3_C_A", "HERC2", "stopgain"), genes)["driver"])     # not a cancer gene
        v = sc.cn_variant(cn_row("chr17_1_C_T", "TP53", "stopgain"), genes)
        self.assertTrue(v["cn_mapped"])
        self.assertEqual(v["category"], "cn_mapped")

    def test_cn_set_and_driver_rows(self):
        with tempfile.TemporaryDirectory() as w:
            os.makedirs(os.path.join(w, "SYN", "final"))
            os.makedirs(os.path.join(w, "SYN", "map"))
            rows = [cn_row("chr17_100_C_T", "TP53", "stopgain", status="trunk", lost_in="n12", lost_cells="3"), cn_row("chr15_5_C_A", "HERC2", "stopgain")]
            with gzip.open(os.path.join(w, "SYN", "final", "cn_mapped.tsv.gz"), "wt") as fh:
                fh.write("\t".join(rows[0].keys()) + "\n")
                for r in rows:
                    fh.write("\t".join(r[k] for k in rows[0].keys()) + "\n")
            with gzip.open(os.path.join(w, "SYN", "final", "cn_mapped_cells.tsv.gz"), "wt") as fh:
                fh.write("mutation\tcells\n")
                fh.write("chr17_100_C_T\t0:2:3,1:4:0,2:0:0\n")
                fh.write("chr15_5_C_A\t1:1:2\n")
            json.dump({"cells": ["c1", "c2", "n1"]}, open(os.path.join(w, "SYN", "map", "cell_index.json"), "w"))
            variants, entries, cells, _ = sc.cn_set(w, "SYN", {"TP53"})
            self.assertEqual(len(variants), 2)
            d = sc.driver_rows("SYN", variants, entries, cells, {"c1", "c2"})
            self.assertEqual(len(d), 1)
            self.assertEqual(d[0]["evidence"], sc.EVIDENCE)
            self.assertEqual(d[0]["clonality"], "clonal")
            self.assertEqual(d[0]["cells"], "1/2")
            self.assertEqual(d[0]["Variant_g"], "17:100-100 C>T")
            self.assertIn("lost in n12", d[0]["variant_summary"])
            # an existing event at the same site is annotated, not duplicated; reruns replace earlier CN-mapped rows
            existing = [{"gene": "TP53", "seqnames": "17", "start": "100", "Variant_g": "17:100-100 C>T", "evidence": "known driver (rescued)"},
                        {"gene": "X", "seqnames": "1", "start": "5", "source": sc.SOURCE}]
            ev, n_add, n_ann = sc.merge_driver_rows(existing, d)
            self.assertEqual((n_add, n_ann), (0, 1))
            self.assertEqual(len(ev), 1)
            self.assertEqual(ev[0]["cn_mapped_branch"]["node"], "n9")
            ev2, n_add2, _ = sc.merge_driver_rows([], d)
            self.assertEqual(n_add2, 1)


if __name__ == "__main__":
    unittest.main()
