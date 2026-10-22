"""Synthetic tests: chrX CN of dryclean-artefact cells is 'not assessed' (no patient data).
python3 -m unittest test_chrx_not_assessed"""
import json
import os
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import chrx_flags as cf  # noqa: E402
import clonality as cl  # noqa: E402
import gos_sc_chrx_not_assessed as na  # noqa: E402


def make_cn(d, cells, x_cn, a_cn=2):
    """cn_binned.json + u16 matrix with one chr1 bin and one chrX bin (100 kb)."""
    os.makedirs(d, exist_ok=True)
    meta = {"cells": cells, "levels": {"100kb": {"bin_size": 100000, "n_bins": 2, "missing_value": 65535,
                                                  "chromosomes": [{"chr": "1", "first_bin": 0, "n_bins": 1}, {"chr": "X", "first_bin": 1, "n_bins": 1}],
                                                  "files": {"cn": "cn.u16"}}}}
    json.dump(meta, open(os.path.join(d, "cn_binned.json"), "w"))
    m = np.array([[a_cn, x_cn[c]] for c in cells], dtype="<u2")
    m.tofile(os.path.join(d, "cn.u16"))


class T(unittest.TestCase):
    def test_not_assess_events(self):
        ev = [{"gene": "ATRX", "type": "SCNA", "vartype": "HOMDEL", "Variant": "HOMDEL", "Genome_Location": "X:1-2", "segment_cn": 0},
              {"gene": "TP53", "type": "SCNA", "vartype": "HOMDEL", "Variant": "HOMDEL", "Genome_Location": "17:1-2", "segment_cn": 0},
              {"gene": "STAG2", "type": "Missense", "Variant": "p.R259*", "Genome_Location": "X:5"}]
        new, n = cf.not_assess_events(ev)
        self.assertEqual(n, 1)
        self.assertEqual(new[0]["vartype"], "NOT_ASSESSED")
        self.assertEqual(new[0]["Variant"], cf.NOT_ASSESSED)
        self.assertIsNone(new[0]["segment_cn"])
        self.assertEqual(new[1]["vartype"], "HOMDEL")   # autosome homdel untouched
        self.assertIs(new[2], ev[2])                    # mutation untouched
        self.assertEqual(ev[0]["vartype"], "HOMDEL")    # input not mutated

    def test_cellcn_masks_flagged_chrx_only(self):
        with tempfile.TemporaryDirectory() as t:
            make_cn(os.path.join(t, "P"), ["a", "b", "c"], {"a": 0, "b": 0, "c": 1})
            cn = cl.CellCN(t, "P", chrx_flagged={"a", "b"})
            self.assertEqual(cn.at("chrX", 5), {"a": None, "b": None, "c": 1})
            self.assertEqual(cn.at("1", 5), {"a": 2, "b": 2, "c": 2})
            st = cn.state("X", 5)
            self.assertEqual(st["modal_cn"], 1)
            self.assertEqual(st["chrx_not_assessed"]["cells"], 2)
            cn_all = cl.CellCN(t, "P", chrx_flagged={"a", "b", "c"})
            st = cn_all.state("X", 5)
            self.assertIsNone(st["modal_cn"])
            self.assertIn("not assessed", st["chrx_not_assessed"]["reason"])
            self.assertEqual(cl.CellCN(t, "P").at("X", 5)["a"], 0)   # no flags: unchanged

    def test_cli_stages_only_flagged_homdel_rows(self):
        with tempfile.TemporaryDirectory() as t:
            data, out = os.path.join(t, "d"), os.path.join(t, "o")
            ev = [{"gene": "ATRX", "type": "SCNA", "vartype": "HOMDEL", "Variant": "HOMDEL", "Genome_Location": "X:1-2", "segment_cn": 0}]
            for c in ("flag1", "ok1"):
                os.makedirs(os.path.join(data, c))
                json.dump(ev, open(os.path.join(data, c, "filtered.events.json"), "w"))
            summ = na.run(data, out, {"P": {"flag1"}})
            self.assertEqual(summ["P"]["rows_changed"], 1)
            self.assertTrue(os.path.exists(os.path.join(out, "flag1", "filtered.events.json")))
            self.assertFalse(os.path.exists(os.path.join(out, "ok1")))
            self.assertEqual(json.load(open(os.path.join(data, "flag1", "filtered.events.json")))[0]["vartype"], "HOMDEL")


if __name__ == "__main__":
    unittest.main()
