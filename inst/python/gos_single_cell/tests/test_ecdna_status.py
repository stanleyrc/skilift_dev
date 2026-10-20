"""Synthetic fixtures only (no patient data): python3 -m unittest test_ecdna_status"""
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gos_sc_ecdna_status as e  # noqa: E402


def simulate(n_cells, mean, ecdna, rng, cv_noise=0.05, gens=12):
    """Per-cell CN: chromosomal = constant; ecDNA = random binomial segregation for `gens` divisions."""
    if not ecdna:
        true = np.full(n_cells, float(mean))
    else:
        cn = np.array([mean])
        for _ in range(gens):
            dau = []
            for c in cn:
                k = rng.binomial(2 * c, 0.5)
                dau += [k, 2 * c - k]
            cn = np.array(dau)
            cn = cn[cn > 0]
            cn = np.round(cn * mean / cn.mean()).astype(int)  # selection keeps the mean
        true = rng.choice(cn, n_cells).astype(float)
    scale = rng.normal(1, cv_noise, n_cells)
    vt = np.full(n_cells, 0.02 * mean)
    x = true * scale + rng.normal(0, np.sqrt(vt))
    return x, vt


class Dispersion(unittest.TestCase):
    def test_stable_amplicon_near_zero(self):
        rng = np.random.default_rng(1)
        x, vt = simulate(60, 40, False, rng)
        D, m = e.dispersion(x, vt, 0.05 ** 2)
        self.assertLess(D, 0.5)
        self.assertAlmostEqual(m, 40, delta=2)

    def test_random_segregation_overdispersed(self):
        rng = np.random.default_rng(2)
        x, vt = simulate(60, 40, True, rng)
        D, _ = e.dispersion(x, vt, 0.05 ** 2)
        self.assertGreater(D, 1)

    def test_bootstrap_axis(self):
        rng = np.random.default_rng(3)
        x = rng.normal(10, 1, (5, 30))
        D, m = e.dispersion(x, np.zeros_like(x), 0.0)
        self.assertEqual(D.shape, (5,))

    def test_cv2_ref_multiplicative(self):
        rng = np.random.default_rng(4)
        scale = rng.normal(1, 0.1, 200)
        W = 2 * scale[:, None] * np.ones((200, 50))
        c0, c90 = e.cv2_ref(W, np.zeros_like(W), np.arange(200))
        self.assertAlmostEqual(float(c0), 0.01, delta=0.003)

    def test_group_stats_calls(self):
        rng = np.random.default_rng(5)
        n = 40
        W = 2 * rng.normal(1, 0.05, n)[:, None] * np.ones((n, 30))
        for ecdna, want in ((False, "stable"), (True, "unstable")):
            x, vt = simulate(n, 40, ecdna, rng)
            ch = np.array([(i, i + 1) for i in range(0, n, 2)])
            r = e.group_stats(x, vt, W, np.zeros_like(W), np.arange(n), np.arange(n), ch, x, vt, rng)
            self.assertEqual(r["call"], want, r)


class Rules(unittest.TestCase):
    def test_call_group(self):
        self.assertEqual(e.call_group(9, 50, 3, 9, 6, 5), "not enough")
        self.assertEqual(e.call_group(20, 4, 3, 9, 6, 5), "not enough")
        self.assertEqual(e.call_group(20, 50, 3, 9, 6, 5), "unstable")
        self.assertEqual(e.call_group(20, 50, 3, 9, 6, 0.2), "indeterminate")
        self.assertEqual(e.call_group(20, 50, 3, 9, 2, np.nan), "unstable (sister untested)")
        self.assertEqual(e.call_group(20, 50, -0.2, 0.3, 6, 0), "stable")
        self.assertEqual(e.call_group(20, 50, -0.2, 0.8, 6, 0), "indeterminate")

    def test_verdict(self):
        cells = {"A": {"a", "b"}, "B": {"c", "d"}, "AB": {"a", "b", "c", "d"}}
        self.assertEqual(e.verdict({"ALL": "unstable", "A": "unstable", "B": "stable"}, cells)[0], "mixed")
        self.assertEqual(e.verdict({"AB": "unstable", "A": "stable"}, cells)[0], "nested")
        self.assertEqual(e.verdict({"ALL": "unstable", "A": "unstable"}, cells)[0], "ecDNA-like")
        self.assertEqual(e.verdict({"A": "stable"}, cells)[0], "chromosomal / integrated")
        self.assertEqual(e.verdict({"A": "indeterminate", "B": "not enough"}, cells)[0], "indeterminate")
        self.assertEqual(e.verdict({"A": "not enough"}, cells)[0], "not enough")

    def test_cherries(self):
        nwk = "(((a:1,b:1):1,(c:1,(d:1,e:1):1):1):1,f:2);"
        self.assertEqual(sorted(e.cherries_from_newick(nwk)), [("a", "b"), ("d", "e")])


class Amplicons(unittest.TestCase):
    def test_detect_focal_not_arm_gain(self):
        chroms = [{"chr": "7", "first_bin": 0, "n_bins": 1594}]
        n = 20
        cn = np.full((n, 1594), 2.0)
        cn[:, 700:1594] = 5.0              # whole-arm gain to 5 (q arm): not focal
        cn[:10, 550:556] = 40.0            # focal amplicon on 7p in 10 cells
        amps, _ = e.detect_amplicons(cn, np.full(n, 2.0), np.ones(n, bool), chroms)
        self.assertEqual(len(amps), 1)
        self.assertEqual(amps[0]["start"], 550 * 100000 + 1)
        self.assertEqual(amps[0]["n_amp_cells"], 10)


if __name__ == "__main__":
    unittest.main()
