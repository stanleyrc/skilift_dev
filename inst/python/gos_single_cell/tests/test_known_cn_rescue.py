"""Synthetic tests of the CN-aware clonality extension and the Tier K known-driver rule (no patient data).
python3 -m unittest test_known_cn_rescue"""
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import clonality as cl  # noqa: E402
import gos_sc_known_cn_rescue as kr  # noqa: E402


def ref(rate=0.6):
    """Reference with the same detection rate at every depth bin (rates are per bin, plenty of observations)."""
    n = [1000] * len(cl.DEPTH_BINS)
    det = [int(rate * 1000)] * len(cl.DEPTH_BINS)
    return {"bins": [list(b) for b in cl.DEPTH_BINS], "labels": [cl.bin_label(i) for i in range(len(cl.DEPTH_BINS))], "n": n, "detected": det,
            "rate": [rate] * len(cl.DEPTH_BINS), "n_sites": 50, "n_cells": 100, "source": "synthetic"}


def simulate(cn_by_cell, m_by_cell, depth=4, seed=0, r=ref()):
    """Carriers / alt reads drawn from the CN-aware model itself (per-cell VAF m / CN)."""
    rng = np.random.default_rng(seed)
    car, alt, dp = {}, {}, {}
    for c, k in cn_by_cell.items():
        p = cl.detection_at_vaf(r, depth, m_by_cell[c] / k)
        hit = rng.random() < p
        car[c] = bool(hit)
        alt[c] = int(hit) * max(1, int(round(depth * m_by_cell[c] / k)))
        dp[c] = depth
    return car, alt, dp


class DetectionAtVaf(unittest.TestCase):
    def test_half_vaf_equals_reference(self):
        r = ref(0.6)
        self.assertAlmostEqual(cl.detection_at_vaf(r, 4, 0.5), 0.6, places=6)

    def test_monotone_in_vaf(self):
        r = ref(0.6)
        v = [cl.detection_at_vaf(r, 4, f) for f in (1 / 3, 0.5, 2 / 3, 1.0)]
        self.assertTrue(all(a < b for a, b in zip(v, v[1:])))
        self.assertEqual(cl.detection_at_vaf(r, 4, 0), 0.0)


class CnAware(unittest.TestCase):
    def test_clonal_on_gain(self):
        cells = {f"t{i}": 3 for i in range(60)}           # clonal gain, mutation on the single (ungained) allele: VAF 1/3
        car, alt, dp = simulate(cells, {c: 1 for c in cells}, seed=1)
        res = cl.cn_aware_clonality(car, dp, cells, {}, ref(), alt)
        self.assertEqual(res["status"], "clonal")
        self.assertTrue(res["cn_aware"])
        self.assertEqual(res["modal_cn"], 3)

    def test_allele_loss_in_cn1_subset_is_clonal(self):
        cells = {f"t{i}": (2 if i < 60 else 1) for i in range(80)}
        m = {c: (1 if k == 2 else 0) for c, k in cells.items()}   # the CN1 cells lost the mutant allele
        car, alt, dp = simulate(cells, m, seed=2)              # m = 0: detection 0 in the CN1 cells
        clades = {c: (["L"] if cells[c] == 1 else ["K"]) for c in cells}
        res = cl.cn_aware_clonality(car, dp, cells, clades, ref(), alt)
        self.assertEqual(res["status"], "clonal")
        self.assertEqual([x["cn"] for x in res["allele_loss"]], [1])

    def test_subclonal_clade(self):
        cells = {f"t{i}": 2 for i in range(80)}
        clades = {c: (["A"] if int(c[1:]) < 20 else ["B"]) for c in cells}
        m = {c: (1 if clades[c] == ["A"] else 0) for c in cells}
        car, alt, dp = simulate(cells, m, depth=8, seed=3)
        res = cl.cn_aware_clonality(car, dp, cells, clades, ref(), alt)
        self.assertEqual(res["status"], "subclonal")
        self.assertEqual(res["best_clade"], "A")

    def test_scattered_low_detection_is_indeterminate(self):
        cells = {f"t{i}": 2 for i in range(80)}
        clades = {c: (["A"] if int(c[1:]) < 40 else ["B"]) for c in cells}
        car = {c: int(c[1:]) % 10 == 0 for c in cells}       # 8 carriers spread over both clades
        dp = {c: 8 for c in cells}
        res = cl.cn_aware_clonality(car, dp, cells, clades, ref(), {c: int(v) for c, v in car.items()})
        self.assertEqual(res["status"], "indeterminate")

    def test_unknown_cn_falls_back(self):
        cells = {f"t{i}": None for i in range(30)}
        dp = {c: 4 for c in cells}
        res = cl.cn_aware_clonality({c: True for c in cells}, dp, cells, {}, ref(), None)
        self.assertIn("not applied", res["cn_reason"])
        self.assertEqual(res["n_cn_unknown"], 30)


    def test_amplicon_falls_back_to_dropout_aware(self):
        cells = {f"t{i}": 30 for i in range(40)}           # amplicon: CN above MAX_CN_MODEL
        res = cl.cn_aware_clonality({c: True for c in cells}, {c: 6 for c in cells}, cells, {}, ref(), None)
        self.assertEqual(res["status"], "clonal")
        self.assertEqual(res["n_cn_excluded"], 40)
        self.assertIn("not applied", res["cn_reason"])


class KnownRule(unittest.TestCase):
    known = kr.norm_known({"genes": ["TP53", "PIK3CA"], "hotspot_residues": {"PIK3CA": ["E545", "H1047"]},
                           "tert_promoter": [{"name": "TERT C228T", "chrom": "chr5", "pos": 1295113, "ref": "G", "alt": "A"}]})

    def test_hotspot_residue(self):
        info = {"Gene.refGene": "PIK3CA", "Func.refGene": "exonic", "ExonicFunc.refGene": "nonsynonymous_SNV",
                "AAChange.refGene": "PIK3CA:NM_X:exon10:c.A1634C:p.E545A"}
        b = kr.known_basis("chr3", 100, "A", "C", info, self.known)
        self.assertEqual(b[1], 1)
        self.assertEqual(b[3], "p.E545A")

    def test_tert(self):
        self.assertEqual(kr.known_basis("chr5", 1295113, "G", "A", {}, self.known)[4], "TERT")

    def test_cosmic_and_novel(self):
        base = {"Gene.refGene": "TP53", "Func.refGene": "exonic", "ExonicFunc.refGene": "nonsynonymous_SNV", "AAChange.refGene": "TP53:NM_X:exon5:c.C1T:p.R175H"}
        self.assertIsNone(kr.known_basis("chr17", 1, "C", "T", dict(base, cosmic70="ID\\x3dCOSM1\\x3bOCCURENCE\\x3d1(lung)"), self.known))
        b = kr.known_basis("chr17", 1, "C", "T", dict(base, cosmic70="ID\\x3dCOSM1\\x3bOCCURENCE\\x3d2(lung),3(skin)"), self.known)
        self.assertEqual(b[1], 2)
        self.assertIsNone(kr.known_basis("chr1", 1, "C", "T", dict(base, **{"Gene.refGene": "OTHER"}), self.known))

    def test_read_qc(self):
        site = {"counts": {"t1": (3, 2, 1, 1), "t2": (2, 2, 1, 1), "t3": (4, 1, 1, 0), "n1": (6, 0, 0, 0)}, "info": {"gnomad40_genome_AF": "."}}
        q = kr.read_qc(site, ["t1", "t2", "t3"], ["n1"])
        self.assertTrue(q["pass"], q["checks"])
        site["counts"]["n1"] = (5, 1, 1, 0)
        self.assertFalse(kr.read_qc(site, ["t1", "t2", "t3"], ["n1"])["pass"])


if __name__ == "__main__":
    unittest.main()
