"""Synthetic report text only (no patient data)."""
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import gos_sc_clinical as g  # noqa: E402

TEXT = (
    "ZZZ1 Oncopanel\r\r ==== ONCOPANEL ==== \r\rTumor Mutational Burden/Megabase: 5.5 \r\rMutations: \r"
    "Tier 1 variants: None identified. \r\rTier 2 variants: \rGENEA c.100G>A (p.V34M), exon 2 - in 40% of 300 reads# \r\r"
    "Structural Variants: \rTier 1 variants: None identified. \r\rCopy Number Variants: \r"
    "Cytoband/Size            Type of Alteration  Genes \r2p24.3                    Amplification       GENEB (estimated 50 copies) \r"
    "7 Chromosomal level      Gain                G1, G2, \rG3 \r9p Arm level              Loss                CDKN2A \r\r"
    "Investigational Variants \rMutations: \rTier 4 variants: \rGENEC c.12+1G>A () - in 54% of 336 reads### \r\r"
    "Structural Variants: \rTier 3 variants: \rRearrangement - GENED exon 3 (chr2:1000) :: X1 intron 2 \r(chr2:2000) \r"
    "Copy Number Variants: \rCytoband/Size            Type of Alteration  Genes \r4q12                      Gain                KDR \r"
    "ZZZ2 – failed oncopanel\r"
    "ZZZ3\r==GENEXUS==\rNothing was detected from Genexus\r"
    "ZZZ4\r==GENEXUS==\rSUMMARY:\r Findings with Potential Clinical Significance: TP53 \rRESULTS: \rSequence Variants: \r"
    "Findings with Potential Clinical Significance: \rGene Variant VAF(%) (Transcript) \rTP53 p.P190L 50.3%  (NM_000546.5) \r"
    "Fusions and Translocations: \r A rearrangement was identified involving AAA Exon 4 and BBB Exon 2 \r"
    "Copy Number Variants: \r EGFR gain \r No PDGFRA amp\rTERT C228T (c.-124C>T): present \rTERT C250T (c.-146C>T): absent \r"
    "Homozygous loss of CDKN2A\r"
)


class ClinicalParse(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p = g.parse_text(g.decode(TEXT.encode("utf-8")))

    def test_patients_and_status(self):
        self.assertEqual(list(self.p), ["ZZZ1", "ZZZ2", "ZZZ3", "ZZZ4"])
        self.assertEqual(self.p["ZZZ2"]["assays"][0]["status"], "failed")
        self.assertEqual(self.p["ZZZ3"]["assays"][0]["status"], "nothing detected")
        self.assertEqual(self.p["ZZZ1"]["assays"][0]["status"], "ok")

    def test_oncopanel_variants(self):
        a = self.p["ZZZ1"]["assays"][0]
        self.assertEqual(a["tmb_per_mb"], 5.5)
        v = a["variants"]
        self.assertEqual([x["gene"] for x in v], ["GENEA", "GENEC"])
        self.assertEqual((v[0]["protein"], v[0]["tier"], v[0]["vaf"], v[0]["reads"], v[0]["exon"]), ("p.V34M", 2, 0.4, 300, 2))
        self.assertIsNone(v[1]["protein"])
        self.assertEqual((v[1]["cdna"], v[1]["tier"], v[1]["investigational"]), ("c.12+1G>A", 4, True))
        self.assertIn("GENEA c.100G>A", v[0]["raw"])

    def test_cnv_rows(self):
        c = self.p["ZZZ1"]["assays"][0]["cnvs"]
        self.assertEqual((c[0]["genes"], c[0]["copies"], c[0]["scope"]), (["GENEB"], 50, "band"))
        self.assertEqual((c[1]["genes"], c[1]["scope"], c[1]["type"]), (["G1", "G2", "G3"], "chromosome", "Gain"))
        self.assertEqual(c[2]["scope"], "arm")
        self.assertTrue(c[3]["investigational"])

    def test_sv_joined_over_lines(self):
        s = self.p["ZZZ1"]["assays"][0]["svs"][0]
        self.assertEqual((s["kind"], s["left"]["pos"], s["right"]["pos"], s["tier"]), ("Rearrangement", 1000, 2000, 3))

    def test_genexus(self):
        a = self.p["ZZZ4"]["assays"][0]
        self.assertEqual((a["variants"][0]["gene"], a["variants"][0]["vaf"]), ("TP53", 0.503))
        self.assertEqual((a["svs"][0]["left"]["feature"], a["svs"][0]["right"]["exon"]), ("AAA", 2))
        self.assertEqual([(c["genes"][0], c["type"]) for c in a["cnvs"]], [("EGFR", "Gain"), ("CDKN2A", "Homozygous loss")])
        self.assertEqual(a["negatives"][0]["gene"], "PDGFRA")
        self.assertEqual(a["tert_promoter"]["C228T"]["status"], "present")
        self.assertEqual(a["tert_promoter"]["C250T"]["status"], "absent")

    def test_meta(self):
        m = g.parse_meta("Location\tMGMT_MET\tAge\tSex\tMRCA_mut\nZZZ1\tLeft Frontal\tN/A\t70\tFemale\t12\n")
        self.assertEqual(m["ZZZ1"], {"location": "Left Frontal", "mgmt": None, "age": 70, "sex": "Female", "old_truncal_mutations": 12})

    def test_status(self):
        st = g.parse_status("# comment\npatient\tdisease_status\nZZZ1\tPrimary\nZZZ2\tRecurrent\nZZZ3\t\n")
        self.assertEqual(st, {"ZZZ1": "Primary", "ZZZ2": "Recurrent"})

    def test_missing_vaf_is_none(self):
        r = g.parse_oncopanel("Tier 1 variants: \nTERT  () - in % of  reads \n")
        self.assertIsNone(r["variants"][0]["vaf"])
        self.assertIsNone(r["variants"][0]["reads"])


if __name__ == "__main__":
    unittest.main()
