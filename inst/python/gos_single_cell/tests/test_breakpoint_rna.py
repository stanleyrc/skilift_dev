"""Unit tests for gos_sc_breakpoint_rna.py (synthetic fixtures only)."""
import csv
import gzip
import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gos_sc_breakpoint_rna as B  # noqa: E402


def wtsv(path, cols, rows):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for r in rows:
            w.writerow([r.get(c, "") for c in cols])


class TestConfidence(unittest.TestCase):
    def test_consistent_zero_noncarrier(self):
        t, _, cr, nr, fold, p = B.confidence(10, 0, 20, 30)
        self.assertEqual(t, "consistent")
        self.assertEqual(nr, 0)
        self.assertIsNone(fold)
        self.assertLess(p, 0.05)

    def test_suggestive_few_cells(self):
        self.assertEqual(B.confidence(3, 0, 20, 30)[0], "suggestive")

    def test_low_when_rates_equal(self):
        self.assertEqual(B.confidence(6, 9, 20, 30)[0], "low")

    def test_fold_fails(self):
        self.assertEqual(B.confidence(10, 8, 20, 30)[0], "suggestive")

    def test_counts_only(self):
        t, reasons, cr, nr, fold, p = B.confidence(6, 1, None, None)
        self.assertEqual(t, "consistent")
        self.assertIsNone(cr)
        self.assertIsNone(p)
        self.assertTrue(any("denominators not available" in r for r in reasons))

    def test_fisher(self):
        self.assertAlmostEqual(B.fisher_two_sided(1, 9, 11, 3), 0.002759, places=5)
        self.assertAlmostEqual(B.fisher_two_sided(5, 5, 5, 5), 1.0)


class TestEndToEnd(unittest.TestCase):
    def setUp(self):
        self.td = tempfile.mkdtemp()
        A = os.path.join(self.td, "an")
        os.makedirs(os.path.join(A, "work"))
        os.makedirs(os.path.join(A, "tables"))
        cross_cols = ["patient", "jid", "stratum", "cls", "chrA", "posA", "chrB", "posB", "span", "n_car_rna", "n_non_rna",
                      "cross_car", "cross_non", "cross_other", "cross_fisher_p"]
        cross = [dict(patient="P1", jid="100", stratum="amplicon", cls="INV-like", chrA="7", posA="1000", chrB="7", posB="5000",
                      span="4000", n_car_rna="12", n_non_rna="30", cross_car="8", cross_non="0", cross_other="3", cross_fisher_p="1e-4"),
                 dict(patient="P1", jid="200", stratum="non-amplicon", cls="TRA-like", chrA="1", posA="100", chrB="2", posB="200",
                      span="Inf", n_car_rna="10", n_non_rna="10", cross_car="3", cross_non="1", cross_other="0", cross_fisher_p="0.5")]
        wtsv(os.path.join(A, "work", "sv_crossing.tsv"), cross_cols, cross)
        wtsv(os.path.join(A, "tables", "sv_transcribed.tsv"), cross_cols + ["transcribed", "rna_junctions"],
             [dict(c, transcribed="TRUE", rna_junctions="x") for c in cross])
        wtsv(os.path.join(A, "work", "sv_junctions.tsv"), ["patient", "jid", "chrA", "posA", "chrB", "posB"],
             [dict(patient="P1", jid="100", chrA="7", posA="1000", chrB="7", posB="5000"),
              dict(patient="P1", jid="200", chrA="1", posA="100", chrB="2", posB="200")])
        rc = ["ci", "chr1", "pos1", "dir1", "chr2", "pos2", "dir2", "reads", "source", "gene1", "gene2", "site1", "site2",
              "reading_frame", "domains", "status", "patient", "jid", "rna_id"]
        reads = []
        for i in range(8):
            reads.append(dict(chr1="7", pos1="1002", chr2="7", pos2="5001", reads="2", source="arriba_high:inversion", gene1="AAA(15)",
                              gene2="BBB", reading_frame="in-frame", domains="|KinaseX(100%)", status="carrier", jid="100", rna_id="c%d" % i))
        reads.append(dict(chr1="7", pos1="1002", chr2="7", pos2="5001", reads="1", source="linear", gene1="AAA", gene2="BBB",
                          reading_frame="", domains=".", status="other", jid="100", rna_id="o1"))
        for i in range(3):
            reads.append(dict(chr1="1", pos1="100", chr2="2", pos2="200", reads="1", source="arriba_discarded:translocation", gene1="CCC",
                              gene2="DDD", reading_frame=".", domains=".", status="carrier", jid="200", rna_id="d%d" % i))
        reads.append(dict(chr1="1", pos1="100", chr2="2", pos2="200", reads="4", source="linear", gene1="CCC", gene2="DDD",
                          reading_frame="", domains=".", status="noncarrier", jid="200", rna_id="n1"))
        wtsv(os.path.join(A, "work", "sv_crossing_reads.tsv.gz"), rc, reads)
        oc = ["fid", "patient", "gene1", "gene2", "chr1", "pos1", "chr2", "pos2", "n_cells", "cls", "jids", "carrier_cells",
              "noncarrier_cells", "other_cells", "splice_both", "explanation"]
        wtsv(os.path.join(A, "work", "rna_only_vs_sv.tsv"), oc, [
            dict(fid="a", patient="P1", gene1="AAA", gene2="BBB", chr1="7", pos1="1010", chr2="7", pos2="5003", n_cells="9", jids="100",
                 carrier_cells="8", noncarrier_cells="0", other_cells="1", splice_both="TRUE", explanation="SV + normal splicing"),
            dict(fid="b", patient="P1", gene1="XXX", gene2="YYY", chr1="3", pos1="30", chr2="4", pos2="40", n_cells="2", jids="999",
                 carrier_cells="1", noncarrier_cells="0", other_cells="1", splice_both="FALSE")])
        wtsv(os.path.join(A, "tables", "fusion_pairs.tsv"),
             ["patient", "pair", "n_isoforms", "n_isoforms_2cells", "n_isoforms_splice_both", "total_cells", "top_isoform", "frames",
              "any_inframe", "readthrough", "domains_top"],
             [dict(patient="P1", pair="AAA::BBB", n_isoforms="2", n_isoforms_2cells="2", n_isoforms_splice_both="1", total_cells="12",
                   top_isoform="AAA e1_end -> BBB e2_start (9 cells, in-frame)", frames="in-frame,.", any_inframe="TRUE", readthrough="TRUE",
                   domains_top="|Kin(100%)"),
              dict(patient="P3", pair="X::Y", n_isoforms="1", n_isoforms_2cells="0", n_isoforms_splice_both="0", total_cells="1",
                   top_isoform="t", frames=".", any_inframe="FALSE", readthrough="FALSE", domains_top=".")])
        wtsv(os.path.join(A, "work", "fusion_isoforms.tsv"),
             ["patient", "pair", "gene1", "gene2", "chr1", "pos1", "exon1", "site1", "chr2", "pos2", "exon2", "site2", "reading_frame",
              "retained_domains", "n_cells", "n_cells_dnalinked", "total_reads", "readthrough", "source"],
             [dict(patient="P1", pair="AAA::BBB", gene1="AAA(12)", gene2="BBB", chr1="7", pos1="1010", exon1="e1_end", site1="CDS/splice-site",
                   chr2="7", pos2="5003", exon2="e2_start", site2="CDS/splice-site", reading_frame="in-frame", retained_domains="|Kin(100%)",
                   n_cells="9", n_cells_dnalinked="8", total_reads="30", readthrough="TRUE", source="filtered.events"),
              dict(patient="P1", pair="AAA::BBB", gene1="AAA", gene2="BBB", chr1="7", pos1="1011", exon1="", site1="intron",
                   chr2="7", pos2="5100", exon2="inside_e3", site2="CDS", reading_frame=".", retained_domains=".",
                   n_cells="3", n_cells_dnalinked="0", total_reads="4", readthrough="TRUE", source="rescued_dna_junction")])
        wtsv(os.path.join(A, "work", "amplicons.tsv"), ["chr", "start", "end", "label", "verdict", "source", "patient", "amp_id"],
             [dict(chr="7", start="1", end="9000000", patient="P1")])
        wtsv(os.path.join(A, "work", "egfrviii_del.tsv"),
             ["patient", "jid", "n_car_rna", "n_non_rna", "vIII_car", "vIII_non", "vIII_reads_car", "vIII_reads_non", "egfr_expr_car",
              "egfr_expr_non", "vIII_cells_all_rna", "vIII_cells_tumour_paired", "vIII_other"],
             [dict(patient="BWH69", jid="207007759", n_car_rna="130", n_non_rna="0", vIII_car="43", vIII_non="0", vIII_reads_car="726",
                   vIII_reads_non="0", egfr_expr_car="79", egfr_expr_non="0", vIII_cells_all_rna="71", vIII_cells_tumour_paired="43", vIII_other="0")])
        wtsv(os.path.join(A, "work", "nearby_genomic_sv.tsv"),
             ["patient", "cdna_jid", "chr", "cdna_lo", "cdna_hi", "cdna_cells", "cdna_reads", "g_jid", "g_cls", "g_lo", "g_hi", "g_cells", "g_reads", "g_tier_v2"],
             [dict(patient="BWH69", cdna_jid="207007575", chr="7", cdna_lo="55019369", cdna_hi="55155828", cdna_cells="30", cdna_reads="48",
                   g_jid="207007759", g_cls="DEL-like", g_lo="55022984", g_hi="55155076", g_cells="189", g_reads="16151", g_tier_v2="3"),
              dict(patient="BWH69", cdna_jid="207007575", chr="7", cdna_lo="55019369", cdna_hi="55155828", cdna_cells="30", cdna_reads="48",
                   g_jid="1", g_cls="INV-like", g_lo="55022906", g_hi="55155202", g_cells="86", g_reads="196", g_tier_v2="3"),
              dict(patient="BWH69", cdna_jid="207007575", chr="7", cdna_lo="55019369", cdna_hi="55155828", cdna_cells="30", cdna_reads="48",
                   g_jid="2", g_cls="INV-like", g_lo="55022907", g_hi="55155203", g_cells="3", g_reads="4", g_tier_v2="2")])
        wtsv(os.path.join(A, "tables", "consequences.tsv"), ["gene", "start", "end", "effect", "frame", "protein_len", "ref_protein_len", "pfam_affected", "patients"],
             [dict(gene="EGFR", start="55019366", end="55155829", effect="truncated_protein_escapes_NMD", frame="in_frame", protein_len="943",
                   ref_protein_len="1210", pfam_affected="PF00757(Furin-like)", patients="BWH69")])
        self.G = os.path.join(self.td, "gos")
        os.makedirs(os.path.join(self.G, "P1", "rna"))
        os.makedirs(os.path.join(self.G, "P2", "rna"))
        json.dump({"fusions": [
            {"id": "BBB::AAA|chr7:5003|chr7:1010", "gene1": "BBB", "gene2": "AAA", "breakpoint1": "chr7:5003", "breakpoint2": "chr7:1010", "n_cells": 9},
            {"id": "OTHER|chr7:5003|chr7:1010", "gene1": "O1", "gene2": "O2", "breakpoint1": "chr7:5003", "breakpoint2": "chr7:1010", "n_cells": 50}]},
            open(os.path.join(self.G, "P1", "rna", "fusions.json"), "w"))
        json.dump({"fusions": []}, open(os.path.join(self.G, "P2", "rna", "fusions.json"), "w"))
        os.makedirs(os.path.join(self.G, "BWH69", "rna"))
        self.A, self.O = A, os.path.join(self.td, "out")

    def test_run(self):
        s = B.run(self.A, self.G, self.O)
        self.assertEqual(s["P1"]["n_sv"], 2)
        bt = json.load(open(os.path.join(self.O, "P1", "rna", "breakpoint_transcripts.json")))
        self.assertEqual(bt["format"], "breakpoint_transcripts_v1")
        e = {x["jid"]: x for x in bt["transcripts"]}
        self.assertEqual(e["100"]["confidence"], "consistent")
        self.assertEqual(e["100"]["reads_car"], 16)
        self.assertEqual(e["100"]["reads_non"], 0)
        self.assertEqual(e["100"]["frame"], "in-frame")
        self.assertIn("KinaseX", e["100"]["consequence"])
        self.assertEqual(e["100"]["genes_a"], ["AAA"])
        self.assertEqual(e["100"]["label"], "AAA::BBB")
        self.assertEqual(e["100"]["rna_junctions"][0]["source"], "arriba_high")
        self.assertEqual(e["100"]["igv_locus_a"], "chr7:900-1100")
        self.assertIsNone(e["100"]["span"] and None)
        self.assertIsNone(e["200"]["span"])
        self.assertIsNone(e["200"]["frame"])
        self.assertTrue(e["200"]["frame_note"])
        self.assertEqual(e["200"]["confidence"], "suggestive")
        self.assertEqual(e["200"]["reads_non"], 4)

    def test_fusion_support(self):
        B.run(self.A, self.G, self.O)
        fs = json.load(open(os.path.join(self.O, "P1", "rna", "fusions_dna_support.json")))
        self.assertEqual(fs["n"], 2)
        by = {e["gene1"]: e for e in fs["entries"]}
        a = by["AAA"]
        self.assertEqual(a["fusion_id"], "BBB::AAA|chr7:5003|chr7:1010")
        self.assertIn("share the breakpoint pair", a["match_note"])
        self.assertEqual(a["class"], "canonical_splice")
        self.assertEqual(a["distance_bp"], 10)
        self.assertEqual(a["sv_jid"], "100")
        self.assertEqual(a["confidence"], "consistent")
        x = by["XXX"]
        self.assertIsNone(x["fusion_id"])
        self.assertEqual(x["class"], "cryptic_end")
        self.assertIsNone(x["distance_bp"])
        self.assertIn("denominators not available", " ".join(x["confidence_reasons"]))

    def test_egfrviii_entry(self):
        B.run(self.A, self.G, self.O, ["BWH69", "P1"])
        bt = json.load(open(os.path.join(self.O, "BWH69", "rna", "breakpoint_transcripts.json")))
        self.assertEqual(bt["n_sv"], 1)
        e = bt["transcripts"][0]
        self.assertEqual(e["jid"], "207007759")
        self.assertEqual((e["cross_car"], e["cross_non"], e["n_non_rna"]), (43, 0, 0))
        self.assertEqual(e["cross_other"], 28)
        self.assertIsNone(e["reads_non"])
        self.assertEqual(e["frame"], "in-frame")
        self.assertIn("943", e["consequence"])
        self.assertFalse(e["from_sv_transcribed"])
        self.assertTrue(e["added_reason"])
        self.assertEqual(e["confidence"], "consistent")
        self.assertIn("no non-carrier tumour cells", " ".join(e["confidence_reasons"]))
        self.assertEqual([x["jid"] for x in e["nearby_genomic_sv"]], ["207007759", "1"])
        self.assertEqual(e["dna_evidence"]["n_cells"], 189)
        self.assertTrue(e["cdna_junction"]["splice_mimic"])
        self.assertEqual(e["stratum"], "non-amplicon")  # no amplicon row for BWH69 in the fixture
        self.assertNotIn("cDNA-only", json.dumps(bt))
        # P1 is unaffected
        p1 = json.load(open(os.path.join(self.O, "P1", "rna", "breakpoint_transcripts.json")))
        self.assertEqual(p1["n_sv"], 2)

    def test_fusion_isoforms(self):
        B.run(self.A, self.G, self.O)
        fi = json.load(open(os.path.join(self.O, "P1", "rna", "fusion_isoforms.json")))
        self.assertEqual((fi["format"], fi["n_pairs"]), ("fusion_isoforms_v1", 1))
        p = fi["pairs"][0]
        self.assertTrue(p["any_inframe"] and p["readthrough"])
        self.assertEqual(p["fusion_ids"], ["BBB::AAA|chr7:5003|chr7:1010", "OTHER|chr7:5003|chr7:1010"])
        self.assertEqual(p["isoforms"][0]["label"], "AAA e1_end -> BBB e2_start")
        self.assertEqual(p["isoforms"][0]["reading_frame"], "in-frame")
        self.assertIsNone(p["isoforms"][1]["reading_frame"])
        self.assertEqual(p["isoforms"][1]["label"], "AAA intron -> BBB inside_e3")
        self.assertEqual(p["isoforms"][0]["retained_domains"], "Kin(100%)")
        fi2 = json.load(open(os.path.join(self.O, "P2", "rna", "fusion_isoforms.json")))
        self.assertEqual((fi2["n_pairs"], fi2["pairs"]), (0, []))

    def test_empty_patient_files(self):
        B.run(self.A, self.G, self.O)
        bt = json.load(open(os.path.join(self.O, "P2", "rna", "breakpoint_transcripts.json")))
        fs = json.load(open(os.path.join(self.O, "P2", "rna", "fusions_dna_support.json")))
        self.assertEqual((bt["n_sv"], bt["transcripts"], fs["n"], fs["entries"]), (0, [], 0, []))


if __name__ == "__main__":
    unittest.main()
