"""Synthetic tests of clonality.py (no patient data).  python3 -m unittest test_clonality"""
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import clonality as cl  # noqa: E402


def reference(rate_by_depth):
    """Synthetic snv_matrix: 50 truncal sites seen in 200 cells at the given depths with the given detection rates."""
    rng = np.random.default_rng(1)
    variants = [{"category": "truncal"} for _ in range(50)] + [{"category": "subclonal"}]
    cells = {}
    for j in range(200):
        entries = []
        for i in range(50):
            d = list(rate_by_depth)[(i + j) % len(rate_by_depth)]
            alt = int(rng.random() < rate_by_depth[d])
            entries.append([i, d - alt, alt, 1 if alt else None])
        entries.append([50, 0, 5, 1])          # subclonal site: ignored
        cells[f"c{j}"] = entries
    return cl.truncal_detection_by_depth({"variants": variants, "cells": cells})


def tree(n):
    """Two clades A (first half) and B (second half), plus a small sub-clade A1 of A."""
    out = {}
    for i in range(n):
        out[f"t{i}"] = ["A"] if i < n // 2 else ["B"]
        if i < n // 4:
            out[f"t{i}"].append("A1")
    return out


class Reference(unittest.TestCase):
    def test_rates_by_depth_bin(self):
        ref = reference({1: 0.4, 2: 0.6})
        self.assertEqual(ref["n_sites"], 50)
        self.assertAlmostEqual(ref["rate"][cl.depth_bin(1)], 0.4, delta=0.05)
        self.assertAlmostEqual(ref["rate"][cl.depth_bin(2)], 0.6, delta=0.05)
        # sparse bin borrows the nearest one
        self.assertAlmostEqual(cl.expected_rate(ref, 3), ref["rate"][cl.depth_bin(2)])

    def test_no_truncal_sites(self):
        self.assertIsNone(cl.truncal_detection_by_depth({"variants": [{"category": "subclonal"}], "cells": {}}))

    def test_poisson_binomial(self):
        self.assertAlmostEqual(cl.poisson_binomial_cdf(0, [0.5, 0.5]), 0.25)
        self.assertAlmostEqual(cl.poisson_binomial_cdf(1, [0.5, 0.5]), 0.75)


class Clonality(unittest.TestCase):
    def setUp(self):
        self.ref = reference({1: 0.4, 2: 0.6})
        self.cells = [f"t{i}" for i in range(80)]
        self.clades = tree(80)

    def test_clonal_with_dropout(self):
        # 80 covered cells at depth 1-2, carriers ~ truncal detection (~0.5) spread over both clades: clonal, although cell_fraction ~ 0.5
        rng = np.random.default_rng(7)
        depth = {c: 1 + (i % 2) for i, c in enumerate(self.cells)}
        car = {c: bool(rng.random() < (0.4 if depth[c] == 1 else 0.6)) for c in self.cells}
        r = cl.dropout_aware_clonality(car, depth, self.clades, self.ref, n_perm=500)
        self.assertLess(r["observed_detection"], 0.85)
        self.assertEqual(r["status"], "clonal")
        self.assertGreaterEqual(r["p_detection"], cl.ALPHA_DETECTION)
        self.assertEqual(r["n_covered"], 80)

    def test_clade_restricted_low_detection_is_subclonal(self):
        # carriers only in sub-clade A1 (20 of 80 cells), detection 0.25 vs expected ~0.5
        depth = {c: 1 + (i % 2) for i, c in enumerate(self.cells)}
        car = {c: "A1" in self.clades[c] and i % 4 != 3 for i, c in enumerate(self.cells)}
        r = cl.dropout_aware_clonality(car, depth, self.clades, self.ref, n_perm=500)
        self.assertEqual(r["status"], "subclonal")
        self.assertLess(r["p_detection"], cl.ALPHA_DETECTION)
        self.assertLess(r["p_clade_restriction"], cl.ALPHA_CLADE)
        self.assertEqual(r["best_clade"], "A1")

    def test_low_detection_not_clade_restricted_is_indeterminate(self):
        # detection 0.1 spread evenly over both clades: below expectation but not a clade -> indeterminate (locus dropout / artefact)
        depth = {c: 1 + (i % 2) for i, c in enumerate(self.cells)}
        car = {c: i % 10 == 0 for i, c in enumerate(self.cells)}
        r = cl.dropout_aware_clonality(car, depth, self.clades, self.ref, n_perm=500)
        self.assertEqual(r["status"], "indeterminate")

    def test_low_coverage_is_indeterminate(self):
        depth = {c: 1 for c in self.cells[:6]}
        car = {c: True for c in self.cells[:3]}
        r = cl.dropout_aware_clonality(car, depth, self.clades, self.ref)
        self.assertEqual(r["status"], "indeterminate")
        self.assertEqual(r["n_covered"], 6)
        self.assertIsNone(r["p_detection"])

    def test_no_reference(self):
        depth = {c: 1 for c in self.cells}
        lo = cl.dropout_aware_clonality({c: i % 2 == 0 for i, c in enumerate(self.cells)}, depth, self.clades, None)
        hi = cl.dropout_aware_clonality({c: i % 20 != 0 for i, c in enumerate(self.cells)}, depth, self.clades, None)
        self.assertEqual(lo["status"], "indeterminate")
        self.assertEqual(hi["status"], "clonal")

    def test_apply_keeps_raw_fraction(self):
        row = {"cell_fraction": 0.5, "clonality": "subclonal", "source": "called SNV (not a tree site)", "Variant_g": "5:1295113-1295113 G>A"}
        cl.apply(row, {"status": "clonal"})
        self.assertEqual(row["cell_fraction"], 0.5)
        self.assertEqual(row["clonality"], "clonal")
        self.assertEqual(row["clonality_fraction_rule"], "subclonal")
        self.assertEqual(row["clonality_method"], "dropout-aware")
        self.assertEqual(cl.event_clonality(row), "clonal")
        self.assertIsNone(cl.event_clonality({"clonality": "clonal"}))
        self.assertEqual(cl.site_id_of(row), "chr5:1295113:G>A")

    def test_cohort_reference_pools_counts(self):
        a, b = reference({1: 0.4, 2: 0.6}), reference({1: 0.5, 2: 0.7})
        c = cl.cohort_reference({"P1": a, "P2": None, "P3": b})
        self.assertEqual(c["n"][0], a["n"][0] + b["n"][0])
        self.assertEqual(sorted(c["patients"]), ["P1", "P3"])
        r = cl.dropout_aware_clonality({c_: i % 2 == 0 for i, c_ in enumerate(self.cells)}, {c_: 1 for c_ in self.cells}, self.clades, dict(c, scope="cohort"), n_perm=200)
        self.assertEqual(r["reference_scope"], "cohort (no truncal SNVs in this patient)")
        self.assertIsNone(cl.cohort_reference({"P": None}))

    def test_clinical_row_inherits(self):
        import gos_sc_clinical_events as ce
        ev = {"gene": "TERT", "vartype": "Hotspot", "Genome_Location": "5:1295113-1295113", "clonality": "clonal", "clonality_method": "dropout-aware", "clonality_stats": {"status": "clonal"}}
        row = ce.inherit_clonality({}, ev)
        self.assertEqual(row["clonality"], "clonal")
        self.assertEqual(ce.inherit_clonality({}, {"clonality": "subclonal"}), {})

    def test_newick(self):
        cl_ = cl.parse_newick_clades("((a:1,b:1)n1:1,(c,d)n2)root;")
        self.assertIn({"a", "b"}, cl_)
        self.assertIn({"a", "b", "c", "d"}, cl_)


if __name__ == "__main__":
    unittest.main()
