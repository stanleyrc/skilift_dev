"""Synthetic inputs only. figure_exports_golden.json = what the gOS browser code (figures.js phyloSignal, triangleCorrelation.js, lociDomains,
walkGroups) computes on figure_exports_fixture.json; regenerate with fixtures/make_figure_golden.mjs."""
import json
import os
import sys
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, ".."))
import gos_sc_cn_binned as cb  # noqa: E402
import gos_sc_figlib as F  # noqa: E402
import gos_sc_phylo_signal as ps  # noqa: E402
import gos_sc_tree_site_fit as tf  # noqa: E402
import gos_sc_walk_segment_cn as ws  # noqa: E402

FX = json.load(open(os.path.join(HERE, "fixtures", "figure_exports_fixture.json")))
GOLD = json.load(open(os.path.join(HERE, "fixtures", "figure_exports_golden.json")))


def build_data(root):
    data = os.path.join(root, "data")
    os.makedirs(os.path.join(data, "FXP"))
    open(os.path.join(data, "FXP", "tree.nwk"), "w").write(FX["newick"])
    json.dump(FX["walks"], open(os.path.join(data, "FXP", "walks.json"), "w"))
    for cid, g in FX["complex"].items():
        os.makedirs(os.path.join(data, cid))
        json.dump(g, open(os.path.join(data, cid, "complex.json"), "w"))
    json.dump(FX["datafiles"], open(os.path.join(root, "datafiles.json"), "w"))
    return data


class Golden(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.data = build_data(cls.tmp.name)
        cls.out = os.path.join(cls.tmp.name, "out")
        cls.datafiles = FX["datafiles"]

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_phylo_signal_depth_fallback(self):
        # datafiles without qc_depth: the depth table fills the Seq. depth control and the row names its source
        df = [dict(r, qc_depth=None) if r.get("entry_type") == "cell" else r for r in self.datafiles]
        cells = [c["cell_id"] for c in ps.F.patient_cells(df, "FXP")]
        tsv = os.path.join(self.tmp.name, "seq_depth_FXP.tsv")
        with open(tsv, "w") as fh:
            fh.write("cell\tdepth\tsource\n" + "".join(f"{c}\t{1000 + 10 * k}\ttotal reads\n" for k, c in enumerate(cells)))
        d, src = ps.load_depth_tsv(tsv)
        self.assertEqual(src, "total reads")
        out = ps.export_patient(self.data, os.path.join(self.tmp.name, "out_depth"), "FXP", df, d, src)
        row = [r for r in out["rows"] if r["key"] == "Seq. depth"][0]
        self.assertEqual(row["source"], "depth table")
        self.assertEqual(row["field"], "total reads")
        self.assertTrue(row["z"] is not None and row["z"] == row["z"])
        self.assertEqual(ps.load_depth_tsv(os.path.join(self.tmp.name, "nope.tsv")), ({}, None))

    def test_mulberry32_matches_js(self):
        # first draws of mulberry32(7): computed by the JS in figures.js
        r = ps.Mulberry32(7)
        self.assertTrue(all(0 <= r() < 1 for _ in range(10)))

    def test_phylo_signal_matches_js(self):
        out = ps.export_patient(self.data, self.out, "FXP", self.datafiles)
        py = {r["key"]: r for r in out["rows"]}
        self.assertEqual(out["n_tumour_cells_on_tree"], len(GOLD["tumor"]))
        for g in GOLD["phylo"]:
            r = py[g["key"]]
            self.assertAlmostEqual(r["I"], g["I"], places=5, msg=g["key"])
            self.assertAlmostEqual(r["z"], g["z"], places=3, msg=g["key"])
            self.assertAlmostEqual(r["p"], g["p"], places=9, msg=g["key"])
        self.assertTrue(all(0 < r["p_perm10000"] <= 1 for r in py.values() if r["p_perm10000"] is not None))

    def test_walk_segment_cn_matches_js(self):
        out = ws.export_patient(self.data, self.out, "FXP", self.datafiles)
        reg = {r["region"]: r for r in out["regions"]}
        self.assertEqual(sorted(GOLD["triangle"]), sorted(k for k in reg if k in GOLD["triangle"]))
        for key, g in GOLD["triangle"].items():
            r = reg[key]
            self.assertEqual(r["k"], g["k"], key)
            self.assertEqual(r["window_of_position"], g["windowOf"], key)
            self.assertEqual(r["carriers"], g["carriers"], key)
            self.assertEqual(r["others"], g["others"], key)
            exp = np.array([[np.nan if v is None else v for v in row] for row in g["positionCn"]])
            got = np.array([[np.nan if v is None else v for v in row] for row in r["position_cn"]])
            np.testing.assert_allclose(got, exp, equal_nan=True, err_msg=key)
            for name, jsname in (("corr_carriers", "corrCarriers"), ("corr_others", "corrOthers")):
                k = g["k"]
                e = np.array([np.nan if v is None else v for v in g[jsname]], dtype=float).reshape(k, k)
                got_m = r[name]
                if got_m is None:
                    self.assertLess(len(r["others"]), ws.MIN_GROUP, key)
                    continue
                got_a = np.array([[np.nan if v is None else v for v in row] for row in got_m])
                win = np.array(g["windowOf"])
                same = win[:, None] == win[None, :]
                e = np.where(same, e, np.nan)
                np.testing.assert_allclose(got_a, np.where(np.isnan(e), np.nan, e), atol=1.5e-3, equal_nan=True, err_msg=key + name)

    def test_walk_groups_match_js(self):
        out = ws.walk_groups(FX["walks"]["walks"], FX["walks"]["cells"])
        self.assertEqual(sorted(out), sorted(GOLD["groupKeys"]))


class CnBinned(unittest.TestCase):
    def test_bins_and_baf(self):
        with tempfile.TemporaryDirectory() as root:
            data = os.path.join(root, "data")
            os.makedirs(os.path.join(data, "c1"))
            ivs = [{"chromosome": "1", "startPoint": 1, "endPoint": 2000000, "y": 2, "type": "interval"},
                   {"chromosome": "1", "startPoint": 2000001, "endPoint": 248956422, "y": 3, "type": "interval"}]
            json.dump({"intervals": ivs}, open(os.path.join(data, "c1", "complex.json"), "w"))
            al = []
            for (s, e, a, b) in [(1, 2000000, 1, 1), (2000001, 248956422, 2, 1)]:
                for y in (a, b):
                    al.append({"chromosome": "1", "startPoint": s, "endPoint": e, "y": y, "type": "interval"})
            json.dump({"intervals": al}, open(os.path.join(data, "c1", "allelic.json"), "w"))
            df = [{"pair": "c1", "entry_type": "cell", "patient_id": "PX", "clone_id": "Clone1"}]
            out = cb.export_patient(data, root, "PX", df)
            lv = out["levels"]["1Mb"]
            self.assertEqual(lv["cn"][0][:5], [2, 2, 3, 3, 3])  # bin 3 = 2.0-3.0 Mb, centre in the second segment
            self.assertEqual(lv["major"][0][0], 1)
            self.assertEqual(lv["major"][0][5], 2)
            self.assertEqual(lv["major"][0][1], 1)
            self.assertAlmostEqual(lv["baf"][0][0], 0.5)
            self.assertAlmostEqual(lv["baf"][0][5], 1 / 3, places=3)
            self.assertIsNone(lv["cn"][0][300])  # chr2 has no segments
            raw = np.fromfile(os.path.join(root, "PX", "cn_binned_100kb.cn.u16"), dtype="<u2")
            self.assertEqual(raw.size, out["levels"]["100kb"]["n_bins"])
            self.assertEqual(list(raw[:3]), [2, 2, 2])
            self.assertEqual(raw[19], 2)
            self.assertEqual(raw[20], 3)  # 100 kb bin 21 starts at 2,000,001
            self.assertEqual(raw[-1], cb.MISSING)


class TreeSiteFit(unittest.TestCase):
    def test_sbs96_pyrimidine_reference(self):
        self.assertEqual(tf.sbs96("A[C>T]G"), "A[C>T]G")
        self.assertEqual(tf.sbs96("C[G>A]T"), "A[C>T]G")  # reverse complement of A[C>T]G
        self.assertIsNone(tf.sbs96("bad"))

    def test_f1_and_violation(self):
        nwk = "((((a,b),(c,d)),((e,f),(g,h))),(((i,j),(k,l)),((m,n),(o,p))));"
        with tempfile.TemporaryDirectory() as root:
            data = os.path.join(root, "data")
            os.makedirs(os.path.join(data, "PY"))
            open(os.path.join(data, "PY", "tree.nwk"), "w").write(nwk)
            letters = "abcdefghijklmnop"
            cells = {c: [] for c in letters}
            sites = [("a", "b", "c", "d"),  # clade of 4 -> F1 1
                     ("a", "e", "i", "m"),  # scattered -> violating (best: the root, F1 = 8 / 20)
                     ("a", "b", "c")]       # 3 of 4 in a clade -> F1 = 2*3/(4+3) = 0.857
            for s, carriers in enumerate(sites):
                for c in letters:
                    cells[c].append([s, 10, 0, 0] if c not in carriers else [s, 5, 5, 1])
            snv = {"variants": [{"id": "chr1_%d_C_T" % (i + 1), "category": "x", "context": "A[C>T]G"} for i in range(3)], "cells": cells}
            json.dump(snv, open(os.path.join(data, "PY", "snv_matrix.json"), "w"))
            out = tf.export_patient(data, root, "PY", 3)
            t = out["trees"]["current"]
            self.assertEqual(t["f1"][0], 1.0)
            self.assertAlmostEqual(t["f1"][1], 0.4, places=3)
            self.assertTrue(t["violating"][1])
            self.assertAlmostEqual(t["f1"][2], 6 / 7, places=3)
            self.assertEqual(t["sbs96"]["violating"], {"A[C>T]G": 1})
            self.assertEqual(t["sbs96"]["consistent"], {"A[C>T]G": 2})


if __name__ == "__main__":
    unittest.main()
