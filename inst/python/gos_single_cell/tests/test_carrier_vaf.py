"""Synthetic fixtures only: python3 test_carrier_vaf.py"""
import json, os, sys, tempfile, unittest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gos_sc_subclonal_drivers as s


class CarrierVaf(unittest.TestCase):
    def test_flags_low_single_read_carriers_not_real_ones(self):
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(os.path.join(d, "P"))
            variants = [{"id": "good", "category": "subclonal", "clade_score": 0.9, "alt_cells": 5},
                        {"id": "err"}, {"id": "real"}]
            cells = {}
            for k in range(8):
                cells[f"c{k}"] = [[0, 10, 7, 1], [1, 20, 1, 1], [2, 10, 6, 1]]   # good 0.7; err 1/20 per cell; real 0.6
            json.dump({"variants": variants, "cells": cells}, open(os.path.join(d, "P", "snv_matrix.json"), "w"))
            hits = [{"patient": "P", "snv_id": "err"}, {"patient": "P", "snv_id": "real"}, {"patient": "P", "snv_id": "absent"}]
            s.carrier_vaf_flags(d, hits)
            self.assertEqual(len(hits[0]["carrier_vaf_check"]["flags"]), 2)
            self.assertEqual(hits[1]["carrier_vaf_check"]["flags"], [])
            self.assertFalse(hits[2]["carrier_vaf_check"]["available"])
            self.assertAlmostEqual(hits[0]["carrier_vaf_check"]["baseline_vaf"], 0.7)


if __name__ == "__main__":
    unittest.main()
