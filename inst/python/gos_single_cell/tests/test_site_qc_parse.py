"""Synthetic fixtures only: python3 test_site_qc_parse.py"""
import os, sys, tempfile, unittest
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gos_sc_driver_site_qc as Q
import driver_filters as F


class Parse(unittest.TestCase):
    def test_load_transcripts_only_requested(self):
        with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as f:
            f.write(">0 ENST1 chrX+ 100-109,200-209\nACGTACGTAC\nGGGGGGGGGG\n>1 ENST2 chr2- 5-14\nTTTTTTTTTT\n")
        seq, meta = Q.load_transcripts(f.name, {"0"})
        os.unlink(f.name)
        self.assertEqual(seq, {"0": "ACGTACGTACGGGGGGGGGG"})
        self.assertEqual(meta["0"], ("ENST1", "chrX", 100, 209))
        self.assertNotIn("1", meta)

    def test_parse_sam_reverse_strand_query_index(self):
        # 201-base query, site at index 100; reverse alignment: index in SAM orientation is qlen-1-100
        line = "0\t16\tchrX\t1000\t60\t201M\t*\t0\t0\tA\t*\tNM:i:8"
        r = Q.parse_sam(line + "\n", 100)["0"][0]
        self.assertEqual((r["qidx"], r["pos0"], r["rev"], r["nm"]), (100, 999, True, 8))
        line2 = "0\t16\tchrX\t1000\t60\t5S196M\t*\t0\t0\tA\t*\tNM:i:2"
        r2 = Q.parse_sam(line2 + "\n", 100)["0"][0]
        self.assertEqual(r2["qidx"], 100)
        self.assertEqual(F.query_to_ref(r2["cigar"], r2["pos0"], r2["qidx"]), 999 + 95)

    def test_unmapped_skipped_and_comp(self):
        self.assertEqual(Q.parse_sam("0\t4\t*\t0\t0\t*\t*\t0\t0\tA\t*\n", 50), {})
        self.assertEqual(Q.comp("A") + Q.comp("c") + Q.comp("N"), "TGN")


if __name__ == "__main__":
    unittest.main()
