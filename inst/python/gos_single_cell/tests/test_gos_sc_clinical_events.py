"""Synthetic inputs only."""
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import gos_sc_clinical_events as ce  # noqa: E402

EVENTS = [
    {"gene": "TP53", "type": "Missense", "vartype": "SNV", "Variant": "p.R175H", "Genome_Location": "17:100-100", "seqnames": "17", "start": "100", "end": "100", "cell_ids": "c1,c2"},
    {"gene": "EGFR", "type": "SCNA", "vartype": "AMP", "Genome_Location": "7:5-50", "seqnames": "7", "cell_ids": "c1"},
    {"gene": "A::B", "fusion_genes": "A::B", "type": "Fusion", "vartype": "fusion", "Genome_Location": "1:1-2", "cell_ids": "c3"},
]


class ClinicalEvents(unittest.TestCase):
    def test_event_key_matches_gos_driver_key(self):
        self.assertEqual(ce.event_key(EVENTS[0]), "TP53 SNV|17:100-100")
        self.assertEqual(ce.event_key(EVENTS[1]), "EGFR AMP|7:5-50")
        self.assertEqual(ce.event_key(EVENTS[2]), "A::B|1:1-2")
        self.assertEqual(ce.event_key({"gene": "X", "seqnames": "2", "start": "3", "end": "4"}), "X|2:3-4")

    def test_snv_type(self):
        self.assertEqual(ce.snv_type("p.C379*", "c.1135T>A"), ("Trunc", "SNV"))
        self.assertEqual(ce.snv_type("p.R2P", "c.1G>C"), ("Missense", "SNV"))
        self.assertEqual(ce.snv_type("p.V301del", "c.901_903del"), ("Missense", "INDEL"))
        self.assertEqual(ce.snv_type(None, "c.894+1G>A"), ("Splice", "SNV"))
        self.assertEqual(ce.snv_type(None, "c.379+40G>T"), ("Noncoding", "SNV"))
        self.assertEqual(ce.snv_type(None, "c.-65G>C"), ("Noncoding", "SNV"))
        self.assertEqual(ce.snv_type(None, "c.1799_1814+4del"), ("Splice", "INDEL"))

    def test_find_snv_event(self):
        self.assertIs(ce.find_snv_event(EVENTS, "TP53", "p.R175H", "chr17", 100, 100), EVENTS[0])
        self.assertIs(ce.find_snv_event(EVENTS, "TP53", None, "chr17", 90, 110), EVENTS[0])
        self.assertIsNone(ce.find_snv_event(EVENTS, "TP53", "p.R213Q", "chr17", 200, 200))
        self.assertIsNone(ce.find_snv_event(EVENTS, "EGFR", None, "chr7", 5, 50))  # SCNA never matches an SNV

    def test_find_cn_event(self):
        self.assertIs(ce.find_cn_event(EVENTS, "EGFR", "amp"), EVENTS[1])
        self.assertIs(ce.find_cn_event(EVENTS, "EGFR", "gain"), EVENTS[1])
        self.assertIsNone(ce.find_cn_event(EVENTS, "EGFR", "loss"))
        self.assertIsNone(ce.find_cn_event(EVENTS, "EGFR", "deep"))

    def test_base_row_and_summary(self):
        clin = {"status": "recovered", "germline_like": True}
        r = ce.base_row("P", "snv", "G", "SNV", "Missense", "p.A1T", "1:5-5", "1", 5, 5, ["b", "a"], 10, 12, clin, "k", vaf=0.5)
        self.assertEqual((r["n_cells"], r["cell_ids"], r["cells"], r["cell_fraction"], r["source"], r["germline_like"]), (2, "a,b", "2/10", 0.2, "clinical assay", True))
        self.assertNotIn("Tier", r)
        r2 = ce.base_row("P", "snv", "H", "SNV", "Missense", "p.A1T", None, None, None, None, [], None, 12, {"status": "not_assessable"}, "k2")
        r2["matches_event"] = "x"
        s = ce.summarize([r, r2])
        self.assertEqual((s["findings"], s["germline_like"], s["recovered"], s["not_assessable"], s["linked_to_filtered_event"]), (2, 1, 1, 1, 1))


if __name__ == "__main__":
    unittest.main()


class Annot(unittest.TestCase):
    def test_clinvar_and_gnomad_parsing(self):
        import gos_sc_clinical_annot as a
        line = "17\t43094567\t12345\tC\tT\t.\t.\tCLNSIG=Pathogenic/Likely_pathogenic;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts;CLNDN=Fanconi_anemia|not_provided"
        c = a.clinvar_from_line(line, "2026-10-11")
        self.assertEqual((c["variation_id"], c["stars"], c["plp"], c["condition"]), ("12345", 2, True, "Fanconi_anemia; not_provided"))
        self.assertTrue(c["link"].endswith("/12345/"))
        self.assertFalse(a.is_plp("Conflicting_classifications_of_pathogenicity"))
        self.assertFalse(a.is_plp("Benign"))
        self.assertEqual(a.clinvar_none("X")["text"], "not in ClinVar (release X)")
        g = a.gnomad_from_line("1\t100\trs1\tA\tC,G\t.\tPASS\tAC=3,5;AN=100,100;AF=0.03,0.05;nhomalt=0,1;popmax=nfe,afr;AF_popmax=0.04,0.09", "G")
        self.assertEqual((g["ac"], g["af"], g["nhomalt"], g["popmax"], g["af_popmax"]), (5, 0.05, 1, "afr", 0.09))
        self.assertIsNone(a.gnomad_from_line("1\t100\t.\tA\tC\t.\tPASS\tAC=3", "T"))
        self.assertIn("variant/17-100-A-G?dataset=gnomad_r2_1", a.gnomad_link("chr17", 100, "A", "G"))

    def test_lof_tags_and_note(self):
        import gos_sc_clinical_annot as a
        self.assertTrue(a.is_lof("Trunc", None))
        self.assertTrue(a.is_lof("Splice", "c.894+1G>A"))
        self.assertFalse(a.is_lof("Splice", "c.379+40G>T"))
        self.assertFalse(a.is_lof("Missense", "c.1G>C"))
        cv = {"found": True, "plp": True, "condition": "Fanconi_anemia"}
        self.assertEqual(a.highlight_tags(cv, True, True), ["ClinVar P/LP", "LoF, germline"])
        self.assertEqual(a.highlight_tags(None, True, False), ["LoF"])
        self.assertIn("recessive", a.recessive_note(cv, True))
        self.assertIsNone(a.recessive_note({"found": True, "condition": "not_provided"}, True))
        self.assertIsNone(a.recessive_note(cv, False))
