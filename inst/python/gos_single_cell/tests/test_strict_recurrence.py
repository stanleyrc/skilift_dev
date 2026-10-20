"""Synthetic fixtures only: strict vs all-hits filter assembly in the exporter's helper path."""
import os, sys, unittest
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import driver_filters as F


class StrictAssembly(unittest.TestCase):
    def test_excluded_reason_keeps_row_and_informational_does_not_exclude(self):
        codes = ["read_end", "same_allele_other", "segdup"]
        ex, info = F.split_codes(codes)
        self.assertEqual(ex, ["read_end"])
        self.assertEqual(sorted(info), ["same_allele_other", "segdup"])
        self.assertEqual(F.split_codes(["segdup", "alt_mapq_lower"])[0], [])

    def test_all_excluding_codes_listed(self):
        for c in ("exon_split", "split_only", "read_end", "low_carrier_vaf", "single_read", "paralog", "recurrent_allele"):
            self.assertIn(c, F.EXCLUDING)
        for c in ("segdup", "low_mappability", "alt_mapq_lower"):
            self.assertIn(c, F.INFORMATIONAL)


if __name__ == "__main__":
    unittest.main()
