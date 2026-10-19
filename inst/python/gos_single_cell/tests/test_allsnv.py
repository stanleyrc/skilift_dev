"""Synthetic fixtures only (no patient data): python3 -m pytest test_allsnv.py  (or python3 test_allsnv.py)"""
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gos_sc_allsnv as a  # noqa: E402


def preset(af, ar, rf, rr, ca, cc):
    return dict(alt_f=af, alt_r=ar, ref_f=rf, ref_r=rr, other=0, **{"del": 0}, cells_alt=ca, cells_cov=cc)


class AllSnv(unittest.TestCase):
    def test_rna_like_needs_spliced_alt_excess(self):
        self.assertTrue(a.rna_like(dict(alt=20, alt_spl=8, ref=100, ref_spl=3))[0])      # 40% of alt vs 3% of ref
        self.assertFalse(a.rna_like(dict(alt=20, alt_spl=1, ref=100, ref_spl=0))[0])     # a single spliced read
        self.assertFalse(a.rna_like(dict(alt=20, alt_spl=8, ref=20, ref_spl=8))[0])      # spliced alignments everywhere (a spliced gene region)
        self.assertFalse(a.rna_like(dict(alt=0, alt_spl=0, ref=5, ref_spl=0))[0])

    def test_row_clonality_and_cells(self):
        cells = [f"c{i}" for i in range(10)]
        clone = ["Normal", "Normal"] + ["A"] * 8
        calls = dict(cells=cells, clone=clone)
        o = dict(cells=dict(alt=[0, 0, 2, 1, 1, 3, 1, 0, 0, 0], dp=[4, 4, 5, 5, 6, 9, 3, 6, 0, 0]),
                 presets=dict(caller=dict(tumor=preset(5, 3, 14, 13, 5, 6), normal=preset(0, 0, 4, 4, 0, 2))), splice=dict(alt=8, alt_spl=0, ref=27, ref_spl=0))
        c = dict(gene="NF1", protein="p.R1*", type="Trunc", tier=2, ref="C", alt="T", chrom="chr17", pos=31200000, stats={}, role="Tumor Suppressor")
        v = dict(level="supported", label="Supported", summary="s", checks=[])
        row = a.make_row("P", c, calls, o, v, (False, None))
        self.assertEqual(row["cell_ids"], "c2,c3,c4,c5,c6")           # normals excluded, >= 1 alt read
        self.assertEqual(row["covered_ref_cell_ids"], "c7")
        self.assertEqual(row["cells"], "5/6")
        self.assertEqual(row["clonality"], "subclonal")               # 5/6 = 0.83 < 0.85
        self.assertEqual(row["source"], a.SOURCE)
        self.assertTrue(row["rescued"])
        self.assertEqual(row["coverage_fraction"], 0.75)


if __name__ == "__main__":
    unittest.main()
