"""Synthetic fixtures only: python3 test_signature_assignment.py"""
import os, sys, unittest
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gos_sc_signature_assignment as sa


def toy_profiles(rng, K=96, m=3):
    S = rng.random((K, m)) ** 4
    return S / S.sum(axis=0, keepdims=True)


class Assignment(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(3)
        self.S = toy_profiles(self.rng)
        self.act = np.array([600.0, 300.0, 100.0])

    def test_posterior_rows_sum_to_one_and_follow_bayes(self):
        P = sa.channel_posterior(self.S, self.act)
        self.assertTrue(np.allclose(P.sum(axis=1), 1))
        k = 5
        w = self.S[k] * self.act
        self.assertTrue(np.allclose(P[k], w / w.sum()))

    def test_zero_channel_falls_back_without_nan(self):
        S = self.S.copy()
        S[0, :] = 0
        P = sa.channel_posterior(S, self.act)
        self.assertFalse(np.isnan(P).any())
        self.assertTrue(np.allclose(P[0], self.act / self.act.sum()))

    def test_patient_level_fractions_recover_activities_in_expectation(self):
        mix = self.S @ self.act
        mix = mix / mix.sum()
        cnt = self.rng.multinomial(200000, mix).astype(float)
        P = sa.channel_posterior(self.S, self.act)
        f = sa.group_fractions(cnt, P)
        self.assertTrue(np.allclose(f, self.act / self.act.sum(), atol=0.01))

    def test_absent_signature_cannot_appear(self):
        # a group made purely of a channel profile of signature 3 (not in the set) gets only set-signature mass and a poor fit
        S_set = self.S[:, :2]
        act = np.array([1.0, 1.0])
        P = sa.channel_posterior(S_set, act)
        cnt = self.rng.multinomial(400, self.S[:, 2]).astype(float)
        f = sa.group_fractions(cnt, P)
        self.assertEqual(len(f), 2)
        self.assertAlmostEqual(f.sum(), 1.0)
        q = sa.fit_quality(cnt, S_set, act, 300, np.random.default_rng(1))
        self.assertLess(q["p_fit_refit"], 0.01)
        pooled = self.rng.multinomial(4000, (self.S[:, :2] @ act) / (self.S[:, :2] @ act).sum()).astype(float)
        q2 = sa.fit_quality(cnt, S_set, act, 300, np.random.default_rng(1), pooled=pooled)
        self.assertLess(q2["p_vs_pooled"], 0.01)
        same = self.rng.multinomial(150, pooled / pooled.sum()).astype(float)
        self.assertGreater(sa.fit_quality(same, S_set, act, 300, np.random.default_rng(1), pooled=pooled)["p_vs_pooled"], 0.01)

    def test_group_from_the_patient_mixture_is_not_flagged(self):
        mix = self.S @ self.act
        mix = mix / mix.sum()
        flagged = 0
        for s in range(10):
            cnt = np.random.default_rng(100 + s).multinomial(150, mix).astype(float)
            q = sa.fit_quality(cnt, self.S, self.act, 300, np.random.default_rng(s))
            flagged += q["p_fit_refit"] < 0.01
        self.assertLessEqual(flagged, 1)

    def test_bootstrap_ci_covers_point(self):
        P = sa.channel_posterior(self.S, self.act)
        idx = self.rng.integers(0, 96, 120)
        pt, B = sa.bootstrap_fractions(idx, P, 200, np.random.default_rng(0))
        self.assertTrue(np.all(np.quantile(B, 0.025, axis=0) <= pt + 1e-9))
        self.assertTrue(np.all(np.quantile(B, 0.975, axis=0) >= pt - 1e-9))

    def test_patient_assignment_end_to_end(self):
        channels = [f"ch{k}" for k in range(96)]
        cos = {"A": self.S[:, 0], "B": self.S[:, 1], "C": self.S[:, 2]}
        ci = {c: k for k, c in enumerate(channels)}
        mix = self.S @ self.act
        mix = mix / mix.sum()
        ch = self.rng.choice(96, size=400, p=mix)
        variants = [{"id": f"v{i}", "context": channels[k], "category": "truncal" if i < 150 else ("private" if i < 330 else "subclonal"), "cellphy_input": True} for i, k in enumerate(ch)]
        clades = {"clades": [{"id": "root", "parent": None, "is_trunk": True, "private_snvs": []},
                             {"id": "c1", "parent": "root", "private_snvs": [f"v{i}" for i in range(150, 330)]},
                             {"id": "tiny", "parent": "root", "private_snvs": ["v1", "v2"]}]}
        sj = {"sets": [{"name": "all", "activities": [{"signature": "A", "activity": 600}, {"signature": "B", "activity": 300}, {"signature": "C", "activity": 100}, {"signature": "X", "activity": 5}]}]}
        out, tests = sa.patient_assignment("P", {"variants": variants}, clades, sj, cos, ci, 100, 100, np.random.default_rng(0))
        self.assertTrue(out["available"])
        self.assertEqual(out["dropped_not_in_cosmic"], ["X"])
        self.assertEqual(out["signatures"], ["A", "B", "C"])
        g = {x["group"]: x for x in out["groups"]}
        self.assertFalse(g["tiny"]["fitted"])
        self.assertTrue(g["c1"]["fitted"])
        self.assertAlmostEqual(sum(f["frac"] for f in g["c1"]["fractions"]), 1.0, places=2)
        self.assertEqual(len(out["mutations"]["ids"]), 400)
        self.assertTrue(all(abs(sum(r) - 1) < 0.01 for r in out["mutations"]["probs"]))
        blocks, shifts = sa.cohort_blocks({"P": out}, {"P": tests})
        self.assertEqual(blocks["P"]["added"], [])
        self.assertTrue(all("q" in t for t in tests))

    def test_missing_inputs_reported_not_faked(self):
        out, tests = sa.patient_assignment("P", None, None, None, {}, {}, 10, 10, np.random.default_rng(0))
        self.assertFalse(out["available"])
        self.assertEqual(tests, [])


if __name__ == "__main__":
    unittest.main()
