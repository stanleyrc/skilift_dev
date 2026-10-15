"""Unit tests on synthetic data for gos_sc_rna_covariate_de.py (run: python -m unittest scripts/test_sc_rna_covariate_de.py)."""
import os
import sys
import unittest

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gos_sc_rna_covariate_de as m  # noqa: E402


class FakePatient:
    def __init__(self, n=40, seed=0):
        rng = np.random.default_rng(seed)
        self.n = n
        self.genes = np.array(["G0", "G1", "G2"])
        self.X = np.vstack([rng.normal(0, 1, n), rng.normal(0, 1, n), rng.normal(0, 1, n)]).astype(np.float32)
        self.X[0, : n // 2] += 3  # G0 higher in the first half
        self.X = np.clip(self.X, 0, None)
        self.E = np.expm1(self.X)
        self.nz = self.X > 0
        self.plate = np.array(["p1", "p2"] * (n // 2))
        self.cycling = np.array([False] * n)
        self.qc = np.array(["lib0mt0"] * n)
        self.mtv = np.ones(n)
        self.lib = np.ones(n)
        self.keep = np.arange(n)
        self.rank_cache = {}


class T(unittest.TestCase):
    def test_bh_monotone_and_bounded(self):
        q = m.bh([0.001, 0.04, 0.5, 0.2])
        self.assertTrue(np.all(q <= 1) and np.all(q >= 0))
        self.assertAlmostEqual(q[0], 0.004)

    def test_event_gene_class(self):
        self.assertEqual(m.event_gene_class("EGFR amp"), ("EGFR", "amp"))
        self.assertEqual(m.event_gene_class("DM (any)"), (None, "complex"))

    def test_stratum_z_matches_mannwhitney_sign_and_size(self):
        P = FakePatient()
        idx = np.arange(P.n)
        a = np.arange(P.n) < P.n // 2
        z = m.stratum_z(P, idx, a)
        ref = stats.mannwhitneyu(P.X[0][a], P.X[0][~a], alternative="two-sided", method="asymptotic").pvalue
        self.assertGreater(z[0], 3)
        self.assertAlmostEqual(2 * stats.norm.sf(abs(z[0])), ref, delta=ref * 0.5 + 1e-6)

    def test_de_patient_and_combine(self):
        P = FakePatient()
        a = np.arange(P.n) < P.n // 2
        r = m.de_patient(P, a, ~a)
        self.assertIsNotNone(r)
        self.assertEqual(r["nA"] + r["nB"], P.n)
        self.assertGreater(r["lfc"][0], 1)
        rows = m.combine([("P1", r, P.genes), ("P2", r, P.genes)])
        g0 = next(x for x in rows if x["gene"] == "G0")
        self.assertEqual((g0["n_patients"], g0["n_concordant"]), (2, 2))
        self.assertTrue(m.is_robust({**g0, "q_val": 0.01, "avg_log2FC": 1.0}))
        self.assertFalse(m.is_robust({**g0, "n_patients": 1, "n_concordant": 1}))

    def test_strata_fall_back_when_groups_confounded_with_plate(self):
        P = FakePatient()
        a = P.plate == "p1"
        plan, adjust = m.strata_plan(P, a, ~a)
        self.assertNotIn("plate", adjust)

    def test_too_small_groups_return_none(self):
        P = FakePatient()
        a = np.zeros(P.n, dtype=bool)
        a[:2] = True
        self.assertIsNone(m.de_patient(P, a, ~a))


if __name__ == "__main__":
    unittest.main()
