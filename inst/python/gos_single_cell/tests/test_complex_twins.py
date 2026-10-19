"""Unit tests for gos_sc_complex_twins.py (synthetic fixtures only)."""
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gos_sc_complex_twins as T  # noqa: E402


def cpx(variant, loc, cells):
    return {"type": "Complex SV", "vartype": variant.split(" ")[0], "Variant": variant, "Genome_Location": loc, "cell_ids": ",".join(cells), "n_cells": len(cells)}


class Twins(unittest.TestCase):
    def test_same_cells_shifted_footprint_are_twins(self):
        ev = [cpx("TIC (6,8)", "6:100-2000,8:5000-9000", ["a", "b", "c"]), cpx("rDelDup (6,8)", "6:100-2000,8:5001-8999", ["a", "b", "c"]),
              cpx("DM (7)", "7:1-100", ["a", "b", "c"]), cpx("TIC (1,2)", "1:1-10,2:1-10", ["a", "b"])]
        groups = T.annotate(ev, "ZZ1")
        self.assertEqual(groups, [[0, 1]])
        self.assertEqual(ev[0]["twin_group"], "ZZ1:1")
        self.assertEqual(ev[1]["twin_names"], ["TIC (6,8)", "rDelDup (6,8)"])
        self.assertNotIn("twin_group", ev[2])  # same cells, different footprint
        self.assertNotIn("twin_group", ev[3])

    def test_far_footprint_not_twin(self):
        ev = [cpx("TIC (6,8)", "6:100-2000,8:5000-9000", ["a", "b"]), cpx("rDelDup (6,8)", "6:100-2000,8:50000-90000", ["a", "b"])]
        self.assertEqual(T.annotate(ev, "ZZ1"), [])

    def test_rerun_clears_stale_marks(self):
        ev = [cpx("TIC (6,8)", "6:1-2", ["a"])]
        ev[0]["twin_group"] = "old"
        T.annotate(ev, "ZZ1")
        self.assertNotIn("twin_group", ev[0])


if __name__ == "__main__":
    unittest.main()
