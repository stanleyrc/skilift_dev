"""Synthetic inputs only."""
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import gos_sc_clinical_recovery as r  # noqa: E402
import gos_sc_clinical_sites as s  # noqa: E402


class Recovery(unittest.TestCase):
    def test_snv_status(self):
        self.assertEqual(r.snv_status(50, 60, 500, vaf=0.45)[0], "recovered")
        self.assertEqual(r.snv_status(1, 60, 500, vaf=0.001)[0], "partial")
        self.assertEqual(r.snv_status(0, 60, 500, vaf=0)[0], "not_recovered")
        self.assertEqual(r.snv_status(0, 2, 500, vaf=0)[0], "not_covered")
        self.assertEqual(r.snv_status(0, 60, 5, vaf=0)[0], "not_covered")

    def test_germline_like(self):
        self.assertTrue(r.germline_like(30, 60))
        self.assertFalse(r.germline_like(1, 60))
        self.assertFalse(r.germline_like(5, 10))

    def test_cell_state_and_classes(self):
        self.assertEqual(r.cell_state(6, 2), {"gain", "amp"})
        self.assertEqual(r.cell_state(3, 2), {"gain"})
        self.assertEqual(r.cell_state(0, 2), {"loss", "deep"})
        self.assertEqual(r.cell_state(4, 4), set())
        self.assertEqual([r.clinical_class(t) for t in ("Amplification", "Gain", "Loss", "Two copy deletion", "Homozygous loss")],
                         ["amp", "gain", "loss", "deep", "deep"])
        self.assertEqual([r.cnv_status(f) for f in (0.9, 0.3, 0.02, None)], ["recovered", "partial", "not_recovered", "not_covered"])

    def test_weighted_median_and_junction(self):
        self.assertEqual(r.length_weighted_median([(1, 10), (2, 100), (3, 10)]), 2)
        j = [("2", 1000, "2", 90000), ("3", 5, "4", 5)]
        self.assertEqual(r.junction_match(j, ("2", 1200), ("2", 90500)), 500)
        self.assertEqual(r.junction_match(j, ("2", 90500), ("2", 1200)), 500)
        self.assertIsNone(r.junction_match(j, ("2", 1200), ("2", 900000)))

    def test_gene_tokens_exact(self):
        self.assertEqual(r.gene_tokens("upstream of METTL9::SLX1B"), {"upstream", "of", "METTL9", "SLX1B"})
        self.assertNotIn("MET", r.gene_tokens("METTL9::CDR2"))


class Sites(unittest.TestCase):
    def tx(self, strand):
        # two exons, CDS starts at the first base of exon 1 (plus) / last base of the last-listed exon in transcript order
        if strand == "+":
            return {"id": "T", "chrom": "chr1", "strand": "+", "exons": [(100, 109), (200, 209)], "start_g": 103, "canonical": True}
        return {"id": "T", "chrom": "chr1", "strand": "-", "exons": [(200, 209), (100, 109)], "start_g": 207, "canonical": True}

    def test_plus_strand_positions(self):
        t = self.tx("+")
        self.assertEqual(s.pos_g(t, "1", None), 103)
        self.assertEqual(s.pos_g(t, "7", None), 109)
        self.assertEqual(s.pos_g(t, "8", None), 200)  # next exon
        self.assertEqual(s.pos_g(t, "-2", None), 101)
        self.assertEqual(s.pos_g(t, "7", "+2"), 111)  # intronic offset

    def test_minus_strand_positions(self):
        t = self.tx("-")
        self.assertEqual(s.pos_g(t, "1", None), 207)
        self.assertEqual(s.pos_g(t, "9", None), 109)
        self.assertEqual(s.pos_g(t, "9", "+1"), 108)

    def test_codon_table(self):
        self.assertEqual(s.CODON["CGA"], "R")
        self.assertEqual(s.CODON["TGA"], "*")
        self.assertEqual(s.CODON["GAG"], "E")


if __name__ == "__main__":
    unittest.main()
