"""Synthetic fixtures only (no patient data): python3 -m pytest test_tert_sv.py  (or python3 test_tert_sv.py)"""
import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gos_sc_tert_sv as t  # noqa: E402


def fusion(genes, variant, tier="3"):
    return {"gene": genes, "fusion_genes": genes, "vartype": "fusion", "type": "Fusion", "Variant": variant, "Tier": tier, "cells": "10/12"}


class TertSv(unittest.TestCase):
    def test_detects_whole_cds_tert_3prime(self):
        self.assertTrue(t.is_tert_promoter_sv(fusion("GENEA::TERT", "In-Frame Fusion Exon 1 (p.22)::Exon 1 (p.1)")))
        self.assertTrue(t.is_tert_promoter_sv(fusion("GENEB::TERT", "In-Frame Fusion Exon -1 (p.5'UTR)::Exon 1 (p.1)")))
        self.assertTrue(t.is_tert_promoter_sv(fusion("GENEB::TERT", "Fusion Exon 2 (p.40)::Exon -1 (p.5'UTR)")))

    def test_rejects_other_shapes(self):
        self.assertFalse(t.is_tert_promoter_sv(fusion("TERT::GENEA", "In-Frame Fusion Exon 1 (p.1)::Exon 3 (p.80)")))  # TERT 5'
        self.assertFalse(t.is_tert_promoter_sv(fusion("GENEA::TERT", "In-Frame Fusion Exon 1 (p.22)::Exon 5 (p.600)")))  # TERT truncated
        self.assertFalse(t.is_tert_promoter_sv({"gene": "TERT", "vartype": "HOMDEL", "type": "SCNA", "Variant": "HOMDEL"}))
        self.assertFalse(t.is_tert_promoter_sv(fusion("GENEA::TERT", "")))

    def test_process_marks_backs_up_and_is_idempotent(self):
        with tempfile.TemporaryDirectory() as d:
            f = os.path.join(d, "filtered.events.json")
            ev = [fusion("GENEA::TERT", "In-Frame Fusion Exon 1 (p.22)::Exon 1 (p.1)"), fusion("GENEC::GENED", "Exon 1 (p.1)::Exon 2 (p.5)")]
            json.dump(ev, open(f, "w"))
            n, labels = t.process(f)
            self.assertEqual(n, 1)
            self.assertEqual(len(labels), 1)
            out = json.load(open(f))
            self.assertEqual(out[0]["Tier"], "1")
            self.assertEqual(out[0]["Tier_original"], "3")
            self.assertTrue(out[0]["tert_promoter_sv"])
            self.assertNotIn("tert_promoter_sv", out[1])
            self.assertEqual(len([x for x in os.listdir(d) if "_tertsv" in x]), 1)
            n2, _ = t.process(f)
            self.assertEqual(n2, 0)
            self.assertEqual(json.load(open(f))[0]["Tier_original"], "3")


if __name__ == "__main__":
    unittest.main()
