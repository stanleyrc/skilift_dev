"""Unit tests for gos_sc_homdel_qc.py (synthetic fixtures only)."""
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gos_sc_homdel_qc as H  # noqa: E402

FOCAL = [{"chromosome": "9", "startPoint": 21_900_000, "endPoint": 22_100_000, "y": 0}, {"chromosome": "9", "startPoint": 1, "endPoint": 21_899_999, "y": 2}]
WHOLE = [{"chromosome": "22", "startPoint": 1, "endPoint": 50_800_000, "y": 0}]


class T(unittest.TestCase):
    def test_focal_not_broad(self):
        e = H.annotate({"Genome_Location": "9:21967753-21995301", "cell_ids": "a,b"}, lambda c: FOCAL)
        self.assertEqual((e["cn0_segment_mb_median"], e["broad_cn0"], e["cn0_cells_measured"]), (0.2, False, 2))

    def test_whole_chromosome_is_broad(self):
        bc = []
        e = H.annotate({"Genome_Location": "22:29603556-29698598", "cell_ids": "a,b,c"}, lambda c: WHOLE, bc)
        self.assertTrue(e["broad_cn0"])
        self.assertEqual(bc, ["a", "b", "c"])
        self.assertAlmostEqual(e["cn0_segment_mb_median"], 50.8, places=1)

    def test_missing_cells(self):
        e = H.annotate({"Genome_Location": "22:1-2", "cell_ids": "a"}, lambda c: None)
        self.assertEqual((e["cn0_segment_mb_median"], e["broad_cn0"], e["cn0_cells_measured"]), (None, False, 0))


if __name__ == "__main__":
    unittest.main()
