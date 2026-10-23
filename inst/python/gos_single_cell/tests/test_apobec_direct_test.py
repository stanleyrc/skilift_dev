"""Synthetic fixtures only: python3 test_apobec_direct_test.py"""
import os, sys, unittest
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gos_sc_apobec_direct_test as t
import gos_sc_portrait as pr


def clade_fixture(n_cg_clade, n_site=100, n_cell=12, seed=0):
    rng = np.random.default_rng(seed)
    ctxs = ["T[C>G]A"] * n_cg_clade + ["A[T>C]G"] * (n_site - n_cg_clade)
    ctxs += ["A[T>C]G"] * 300 + ["T[C>G]T"] * 6          # truncal sites
    ctxs += ["A[T>C]G"] * 200 + ["T[C>G]A"] * 8          # other clade private
    variants = [{"id": f"v{i}", "context": c, "category": "private" if i < n_site else "truncal"} for i, c in enumerate(ctxs)]
    cells = {}
    cell_names = [f"P_MR1_pl{1 + k % 2}_{k}a" for k in range(n_cell)]
    for c in cell_names:
        cells[c] = [[i, 10, int(rng.integers(4, 8)), 1] for i in range(n_site)]
    cells["P_MR1_pl1_nrm"] = [[i, 10, 0, 0] for i in range(n_site)]
    clade = {"id": "P_n1", "n_cells": n_cell, "cells": cell_names, "private_snvs": [f"v{i}" for i in range(n_site)]}
    vpos = {v["id"]: n for n, v in enumerate(variants)}
    truncal = [f"v{i}" for i in range(n_site, n_site + 306)]
    other = [f"v{i}" for i in range(n_site + 306, len(variants))]
    return variants, vpos, cells, clade, truncal, other


class Direct(unittest.TestCase):
    def test_classify(self):
        self.assertEqual(t.classify("T[C>G]A"), "cg_tcw")
        self.assertEqual(t.classify("T[C>T]T"), "ct_tcw")
        self.assertIsNone(t.classify("T[C>G]C"))
        self.assertIsNone(t.classify("A[C>G]A"))
        self.assertIsNone(t.classify(None))

    def test_plate_of(self):
        self.assertEqual(t.plate_of("BWH70_MR_1_pl1_10a", "BWH70"), "MR_1|pl1")
        self.assertEqual(t.plate_of("pool-Q-BWH70-MR-B2-CD45pos-pl2-4c", "BWH70"), "MR-B2-CD45pos|pl2")
        self.assertEqual(t.plate_of("weird", "P"), "unparsed")

    def test_enriched_clean_clade_supported(self):
        v, vpos, cells, clade, tr, ot = clade_fixture(12)
        r = t.test_clade("P", clade, v, vpos, cells, ["P_MR1_pl1_nrm"], list(cells), tr, ot, 0.02)
        r["available"] = True
        t.finalize_verdicts([r])
        self.assertEqual(r["cg_tcw"]["k"], 12)
        self.assertEqual(r["verdict"]["status"], "supported", r["verdict"])

    def test_unenriched_clade_not_supported(self):
        v, vpos, cells, clade, tr, ot = clade_fixture(2)
        r = t.test_clade("P", clade, v, vpos, cells, ["P_MR1_pl1_nrm"], list(cells), tr, ot, 0.02)
        r["available"] = True
        t.finalize_verdicts([r])
        self.assertEqual(r["verdict"]["status"], "not supported")

    def test_normal_cell_alt_reads_fail_the_clade(self):
        v, vpos, cells, clade, tr, ot = clade_fixture(12)
        cells["P_MR1_pl1_nrm"] = [[i, 6, 4, 0] if v[i]["context"] == "T[C>G]A" else [i, 10, 0, 0] for i in range(100)]
        r = t.test_clade("P", clade, v, vpos, cells, ["P_MR1_pl1_nrm"], list(cells), tr, ot, 0.02)
        r["available"] = True
        t.finalize_verdicts([r])
        self.assertEqual(r["verdict"]["status"], "not supported")
        self.assertTrue(any("normal" in f for f in r["verdict"]["failed"]))

    def test_bh(self):
        q = t.bh([0.01, 0.04, 0.03, 0.5])
        self.assertAlmostEqual(q[0], 0.04)
        self.assertTrue(all(0 <= x <= 1 for x in q))


class Cards(unittest.TestCase):
    def test_apobec_cards_downgrade(self):
        def res(clade, k, status, q):
            return {"available": True, "patient": "P", "clade": clade, "n_private_snvs": 100, "verdict": {"status": status, "failed": [], "passed": []},
                    "cg_tcw": {"k": k, "frac": k / 100, "truncal": {"frac": 0.01, "p_greater": 0.001}, "other_clades": {"frac": 0.02}, "q_enrichment": q}}
        class Inp: pass
        inp = Inp()
        inp.sub = {"clade_signatures": {"P": {"clades": [{"clade": "P_n2", "fitted": True, "fit": {"flag_poor_fit": True}}]}}}
        inp.apobec = {"results": [res("P_n1", 12, "not supported", 0.4), res("P_n2", 15, "supported", 0.01)]}
        cards = pr.apobec_cards(inp)
        self.assertEqual(len(cards), 1)
        self.assertEqual(cards[0]["status"], "candidate")
        self.assertIn("Not supported after correction: clade n1", cards[0]["text"])
        inp.apobec = {"results": [res("P_n1", 12, "not supported", 0.4)]}
        self.assertEqual(pr.apobec_cards(inp)[0]["status"], "not supported")
        inp.apobec = None
        self.assertEqual(pr.apobec_cards(inp), [])


if __name__ == "__main__":
    unittest.main()
