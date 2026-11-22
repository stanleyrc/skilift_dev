"""Unit tests for gos_sc_sex_mask.py (synthetic fixtures only)."""
import array
import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gos_sc_sex_mask as M  # noqa: E402


def graph():
    iv = [{"chromosome": "1", "startPoint": 1, "endPoint": 100, "iid": 1, "y": 2, "type": "interval"},
          {"chromosome": "X", "startPoint": 1, "endPoint": 100, "iid": 2, "y": 2, "type": "interval"},
          {"chromosome": "Y", "startPoint": 1, "endPoint": 100, "iid": 3, "y": 997, "type": "interval"}]
    con = [{"cid": 1, "source": 1, "sink": -2, "type": "REF"}, {"cid": 2, "source": 2, "sink": -3, "type": "ALT"},
           {"cid": 3, "source": None, "sink": 1, "type": "LOOSE"}]
    return {"settings": {}, "intervals": iv, "connections": con}


CLIN = {"clinical_meta": {"patients": {"ZZ1": {"sex": "Female"}, "ZZ2": {"sex": "Male"}}}}


class SexMask(unittest.TestCase):
    def test_mask_graph_drops_y_and_its_connections(self):
        g, a, b = M.mask_graph(graph())
        self.assertEqual([i["chromosome"] for i in g["intervals"]], ["1", "X"])
        self.assertEqual([c["cid"] for c in g["connections"]], [1, 3])
        self.assertEqual((a, b), (1, 1))
        self.assertIn("not available", g["masked_chromosomes"]["Y"])
        self.assertEqual(M.mask_graph(g)[1], None)  # idempotent

    def test_sex_lookup(self):
        self.assertEqual(M.sex_of(CLIN, M.patient_of("ZZ1_PRI_REC_clean_rr3")), "female")
        self.assertEqual(M.sex_of(CLIN, "ZZ2"), "male")
        self.assertEqual(M.sex_of(CLIN, "ZZ9"), "unknown")

    def test_run_female_only_with_backup_and_cn_binned(self):
        with tempfile.TemporaryDirectory() as d:
            for rep, cell in (("ZZ1_rr", "ZZ1_rr__c1"), ("ZZ2", "ZZ2_c1")):
                os.makedirs(os.path.join(d, cell)); os.makedirs(os.path.join(d, rep), exist_ok=True)
                json.dump(graph(), open(os.path.join(d, cell, "complex.json"), "w"))
            lv = {"n_bins": 4, "chromosomes": [{"chr": "1", "first_bin": 0, "n_bins": 2}, {"chr": "Y", "first_bin": 2, "n_bins": 2}],
                  "files": {"cn": "cn.u16"}, "missing_value": 65535}
            json.dump({"levels": {"1Mb": lv}}, open(os.path.join(d, "ZZ1_rr", "cn_binned.json"), "w"))
            open(os.path.join(d, "ZZ1_rr", "cn.u16"), "wb").write(array.array("H", [2, 2, 5, 9, 2, 3, 7, 7]).tobytes())
            df = [{"pair": "ZZ1_rr", "entry_type": "patient"}, {"pair": "ZZ1_rr__c1", "entry_type": "cell", "patient_id": "ZZ1_rr"},
                  {"pair": "ZZ2", "entry_type": "patient"}, {"pair": "ZZ2_c1", "entry_type": "cell", "patient_id": "ZZ2"}]
            out = M.run(d, df, ["ZZ1_rr", "ZZ2"], CLIN, ".bak_t", log=lambda *a: None)
            self.assertEqual(out["ZZ1_rr"]["cells_masked"], 1)
            self.assertEqual(out["ZZ2"]["cells_masked"] if "cells_masked" in out["ZZ2"] else out["ZZ2"]["cells"], 0)
            g = json.load(open(os.path.join(d, "ZZ1_rr__c1", "complex.json")))
            self.assertNotIn("Y", [i["chromosome"] for i in g["intervals"]])
            self.assertTrue(os.path.exists(os.path.join(d, "ZZ1_rr__c1", "complex.json.bak_t")))
            self.assertIn("Y", [i["chromosome"] for i in json.load(open(os.path.join(d, "ZZ2_c1", "complex.json")))["intervals"]])
            a = array.array("H"); a.frombytes(open(os.path.join(d, "ZZ1_rr", "cn.u16"), "rb").read())
            self.assertEqual(list(a), [2, 2, 65535, 65535, 2, 3, 65535, 65535])
            self.assertTrue(os.path.exists(os.path.join(d, "ZZ1_rr", "sex_mask.json")))
            out2 = M.run(d, df, ["ZZ1_rr"], CLIN, ".bak_t", log=lambda *a: None)   # second run changes nothing
            self.assertEqual(out2["ZZ1_rr"]["cells_masked"], 0)


if __name__ == "__main__":
    unittest.main()
