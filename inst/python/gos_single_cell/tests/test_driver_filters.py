"""Synthetic fixtures only: python3 test_driver_filters.py"""
import os, sys, unittest
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import driver_filters as F


class Codes(unittest.TestCase):
    def test_exon_split_and_read_end(self):
        c = {"alt_reads": 10, "ref_reads": 100, "alt_exon_split": 8, "alt_sa": 8, "ref_sa": 5, "alt_end10": 10, "ref_end10": 30}
        self.assertEqual(F.read_qc_codes(c), ["exon_split", "read_end"])

    def test_read_end_ref_threshold_is_40_percent(self):
        c = {"alt_reads": 10, "ref_reads": 100, "alt_end10": 9, "ref_end10": 39}
        self.assertEqual(F.read_qc_codes(c), ["read_end"])
        c["ref_end10"] = 41
        self.assertEqual(F.read_qc_codes(c), [])

    def test_split_only_when_ref_not_split(self):
        c = {"alt_reads": 6, "ref_reads": 60, "alt_sa": 5, "ref_sa": 3}
        self.assertEqual(F.read_qc_codes(c), ["split_only"])

    def test_too_few_alt_reads_no_flag(self):
        self.assertEqual(F.read_qc_codes({"alt_reads": 2, "ref_reads": 50, "alt_end10": 2, "alt_exon_split": 2}), [])

    def test_recurrent_allele_needs_three_patients_with_two_cells(self):
        op = {"A": {"cells_with_alt": 3}, "B": {"cells_with_alt": 2}, "C": {"cells_with_alt": 1}}
        self.assertNotIn("recurrent_allele", F.read_qc_codes({}, op))
        op["C"]["cells_with_alt"] = 2
        self.assertIn("recurrent_allele", F.read_qc_codes({}, op))


class Paralog(unittest.TestCase):
    def test_cigar_walk(self):
        cig = F.parse_cigar("10M2I20M3D30M")
        self.assertEqual(F.query_to_ref(cig, 100, 5), 105)
        self.assertIsNone(F.query_to_ref(cig, 100, 10))          # insertion
        self.assertEqual(F.query_to_ref(cig, 100, 12), 110)      # first base after the insertion
        self.assertEqual(F.query_to_ref(cig, 100, 32), 110 + 20 + 3)  # after deletion
        self.assertIsNone(F.query_to_ref(F.parse_cigar("5S20M"), 0, 2))

    def test_identity(self):
        ident, cov = F.aln_identity(F.parse_cigar("90M"), 9, 100)
        self.assertAlmostEqual(ident, 0.9)
        self.assertAlmostEqual(cov, 0.9)

    def test_paralog_call_alt_matches_homolog(self):
        homs = [{"chrom": "chrX", "pos": 5, "base": "G", "identity": 0.95, "coverage": 1.0}]
        self.assertTrue(F.paralog_call("G", "A", homs)["paralog"])
        self.assertFalse(F.paralog_call("T", "A", homs)["paralog"])        # alt differs from the paralog base
        homs[0]["identity"] = 0.6
        self.assertFalse(F.paralog_call("G", "A", homs)["paralog"])        # distant homology does not count
        homs[0].update(identity=0.95, base="A")
        self.assertFalse(F.paralog_call("G", "A", homs)["paralog"])        # paralog carries ref, not alt


class Intervals(unittest.TestCase):
    def test_point_query_1based(self):
        ix = F.IntervalIndex([("chr1", 100, 200), ("chr1", 150, 300), ("chr2", 0, 10)])
        self.assertTrue(ix.contains("chr1", 101))
        self.assertTrue(ix.contains("chr1", 300))
        self.assertFalse(ix.contains("chr1", 301))
        self.assertFalse(ix.contains("chr1", 100))   # 0-based start 100 = 1-based 101
        self.assertFalse(ix.contains("chr3", 5))


class SiteCodes(unittest.TestCase):
    def test_site_codes_and_split(self):
        q = {"paralog": {"paralog": True}, "segdup": True, "umap_k100": False,
             "mapq": {"alt_reads": 10, "ref_reads": 50, "mean_alt": 50.0, "mean_ref": 59.0}}
        codes = F.site_codes(q)
        self.assertEqual(codes, ["paralog", "segdup", "low_mappability", "alt_mapq_lower"])
        ex, info = F.split_codes(codes)
        self.assertEqual(ex, ["paralog"])
        self.assertEqual(info, ["segdup", "low_mappability", "alt_mapq_lower"])

    def test_mapq_gap_small_is_not_flagged(self):
        q = {"mapq": {"alt_reads": 10, "ref_reads": 50, "mean_alt": 58.0, "mean_ref": 59.0}, "umap_k100": True}
        self.assertEqual(F.site_codes(q), [])
        self.assertEqual(F.site_codes(None), [])


if __name__ == "__main__":
    unittest.main()
