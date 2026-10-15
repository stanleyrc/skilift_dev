"""Synthetic test of gos_sc_hotspot_rescue.build_rows (no patient data).  python3 -m unittest test_hotspot_rescue"""
import unittest

import gos_sc_hotspot_rescue as hr

CELLS = [f"P_c{i}" for i in range(10)]
CLONE = ["Clone 1"] * 8 + ["Normal"] * 2
SID = "chr5:1295113:G>A"


def make(rows, tum, verdict_level="supported", capped=False):
    calls = {"cells": CELLS, "clone": CLONE, "dna": {SID: {"alt": [0] * 10, "dp": [0] * 10}}}
    ev = {"sites": {SID: {"sample_every": 1, "presets": {"all": {"tumor": tum, "normal": {"alt_f": 0, "alt_r": 0, "ref_f": 1, "ref_r": 0, "other": 0, "del": 0, "cells_alt": 0, "cells_cov": 1}}},
                          "reads": {"capped": capped, "cols": [], "rows": rows}}}}
    regions = {"hotspots": [{"name": "TERT C228T", "gene": "TERT", "chrom": "chr5", "pos": 1295113, "ref": "G", "alt": "A", "kind": "hotspot"},
                            {"name": "IDH1 R132H", "gene": "IDH1", "chrom": "chr2", "pos": 208248388, "ref": "C", "alt": "T", "kind": "hotspot"}]}
    verdicts = {SID: {"level": verdict_level, "label": verdict_level, "summary": "s", "checks": [{"key": "signal", "status": "pass", "text": "t"}, {"key": "bq", "status": "warn", "text": "low bq"}]}}
    return calls, ev, regions, verdicts


class Rescue(unittest.TestCase):
    def test_called_with_three_states(self):
        # cells 0-3 carry (alt), 4-5 covered without alt, 6-7 no reads, 8 = normal cell with a ref read
        rows = [[0, "A"], [0, "A"], [1, "A"], [2, "A"], [3, "A"], [4, "R"], [5, "R"], [8, "R"]]
        tum = {"alt_f": 3, "alt_r": 2, "ref_f": 1, "ref_r": 1, "other": 0, "del": 0, "cells_alt": 4, "cells_cov": 6}
        calls, ev, regions, verdicts = make(rows, tum)
        out, report = hr.build_rows("P", calls, ev, regions, verdicts)
        self.assertEqual(len(out), 1)
        r = out[0]
        self.assertEqual(r["cell_ids"], "P_c0,P_c1,P_c2,P_c3")
        self.assertEqual(r["covered_ref_cell_ids"], "P_c4,P_c5")
        self.assertEqual((r["cells"], r["cell_fraction"], r["n_cells"]), ("4/6", 0.667, 4))
        self.assertEqual((r["tumour_cells"], r["covered_tumour_cells"], r["low_coverage"]), (8, 6, False))
        self.assertEqual((r["type"], r["Tier"], r["source"], r["rescued"]), ("Promoter", 1, hr.SOURCE, True))
        self.assertEqual(r["hotspot_evidence"]["verdict"]["reasons"], ["low bq"])
        self.assertEqual(r["hotspot_evidence"]["tumor"]["alt_f"], 3)
        self.assertEqual([x["name"] for x in report if x["called"]], ["TERT C228T"])
        self.assertIn("not genotyped", [x for x in report if x["name"] == "IDH1 R132H"][0]["reason"])

    def test_low_coverage_flag(self):
        rows = [[0, "A"], [1, "A"], [2, "A"]]
        tum = {"alt_f": 2, "alt_r": 1, "ref_f": 0, "ref_r": 0, "other": 0, "del": 0, "cells_alt": 3, "cells_cov": 3}
        out, _ = hr.build_rows("P", *make(rows, tum)[:3], make(rows, tum)[3])
        self.assertTrue(out[0]["low_coverage"])  # 3 of 8 tumour cells

    def test_not_called(self):
        rows = [[0, "A"], [1, "A"], [2, "A"]]
        tum = {"alt_f": 2, "alt_r": 1, "ref_f": 0, "ref_r": 0, "other": 0, "del": 0, "cells_alt": 3, "cells_cov": 3}
        for level in ("artefact", "absent", "germline"):
            out, rep = hr.build_rows("P", *make(rows, tum, level)[:3], make(rows, tum, level)[3])
            self.assertEqual(out, [])
            self.assertEqual(rep[0]["reason"], f"verdict {level}")
        out, rep = hr.build_rows("P", *make(rows[:2], {**tum, "cells_alt": 2, "cells_cov": 2})[:3], make(rows[:2], tum)[3])
        self.assertEqual(rep[0]["reason"], "only 2 carrier cells")

    def test_falls_back_to_region_calls_when_rows_do_not_match(self):
        rows = [[0, "A"]]
        tum = {"alt_f": 2, "alt_r": 1, "ref_f": 0, "ref_r": 0, "other": 0, "del": 0, "cells_alt": 3, "cells_cov": 3}
        calls, ev, regions, verdicts = make(rows, tum)
        calls["dna"][SID] = {"alt": [1, 1, 2, 0, 0, 0, 0, 0, 1, 0], "dp": [1, 2, 2, 1, 0, 0, 0, 0, 1, 0]}
        out, _ = hr.build_rows("P", calls, ev, regions, verdicts)
        self.assertEqual(out[0]["cell_ids"], "P_c0,P_c1,P_c2")          # the normal cell 8 is not a carrier
        self.assertEqual(out[0]["covered_ref_cell_ids"], "P_c3")
        self.assertIn("region_calls", out[0]["hotspot_evidence"]["cell_source"])

    def test_not_duplicated_when_already_called(self):
        rows = [[0, "A"], [1, "A"], [2, "A"]]
        tum = {"alt_f": 2, "alt_r": 1, "ref_f": 0, "ref_r": 0, "other": 0, "del": 0, "cells_alt": 3, "cells_cov": 3}
        out, rep = hr.build_rows("P", *make(rows, tum)[:3], make(rows, tum)[3], existing_events=[{"vartype": "SNV", "seqnames": "5", "start": "1295113"}])
        self.assertEqual(out, [])
        self.assertEqual(rep[0]["reason"], "already in the SNV calls")


if __name__ == "__main__":
    unittest.main()
