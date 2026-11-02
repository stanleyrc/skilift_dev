"""Synthetic test of gos_sc_snv_pl_patch (no patient data). python3 -m unittest test_snv_pl_patch"""
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gos_sc_snv_columnar as sc  # noqa: E402
import gos_sc_snv_pl_patch as pp  # noqa: E402


class Patch(unittest.TestCase):
    def test_patch_sets_pl_byte_per_site_and_cell(self):
        entries = [[(0, 3, 2, 1), (2, 1, 1, -1)], [], [(1, 5, 0, 0)]]
        blob, _, head = sc.encode_matrix(entries)
        pl = np.zeros(3 * 3, dtype=np.uint8)
        pl[0 * 3 + 0] = 250
        pl[0 * 3 + 2] = 255
        pl[2 * 3 + 1] = 1
        new, n_pl, n = pp.patch(blob, 3, 3, pl, head)
        self.assertEqual((n_pl, n), (3, 3))
        self.assertEqual(sc.decode_matrix(new, 3, with_pl=True), [[(0, 3, 2, 1, 250), (2, 1, 1, -1, 255)], [], [(1, 5, 0, 0, 1)]])
        self.assertEqual(len(new), len(blob))


if __name__ == "__main__":
    unittest.main()
