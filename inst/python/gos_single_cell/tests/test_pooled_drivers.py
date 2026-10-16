import csv, gzip, json, os, sys, tempfile, unittest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gos_sc_pooled_drivers as G

HDR = "patient gene chrom pos ref alt type consequence protein oncokb_tier oncokb_effect tumour_alt tumour_dp normal_alt normal_dp carriers covered carrier_vaf_median carriers_2plus_frac clade_best clade_f1 rna_alt rna_dp rna_alt_cells other_patients_alt_cells verdict why".split()


def cand(**k):
    d = dict(patient="P1", gene="NF1", chrom="chr17", pos="100", ref="C", alt="A", type="SNV", consequence="coding", protein="p.C379*", oncokb_tier="2", oncokb_effect="LoF",
             tumour_alt="30", tumour_dp="40", normal_alt="0", normal_dp="10", carriers="2", covered="3", carrier_vaf_median="1.0", carriers_2plus_frac="0.9",
             clade_best="P1_n1", clade_f1="1.0", rna_alt="1", rna_dp="2", rna_alt_cells="1", other_patients_alt_cells="0", verdict="high-confidence", why="w")
    d.update(k)
    return d


def make(tmp, cands, events):
    scan, data = os.path.join(tmp, "scan"), os.path.join(tmp, "data")
    os.makedirs(os.path.join(scan, "P1", "s0")); os.makedirs(os.path.join(scan, "annot")); os.makedirs(os.path.join(data, "P1"))
    with gzip.open(os.path.join(scan, "P1", "s0", "obs.tsv.gz"), "wt") as fh:
        fh.write("cell\tnormal\tchrom\tpos\ttype\tref\talt\n")
        for cell, nrm in [("c1", "False"), ("c2", "False"), ("c1", "False"), ("n1", "True")]:
            fh.write(f"{cell}\t{nrm}\tchr17\t100\tSNV\tC\tA\n")
        fh.write("c3\tFalse\tchr17\t100\tSNV\tC\tG\n")
    with open(os.path.join(scan, "P1", "s0", "cells.tsv"), "w") as fh:
        fh.write("cell\tnormal\nc1\t0\nc2\t0\nc3\t0\nn1\t1\n")
    with open(os.path.join(scan, "annot", "annotation.tsv"), "w") as fh:
        fh.write("patient\tmutation\tgene\tclassification\tprotein\ttype\ttier\trole\teffect\teffect_description\tgene_summary\tvariant_summary\ttherapeutics\tprognoses\tdiagnoses\n")
        fh.write('"P1"\t"chr17_100_C_A"\t"NF1"\t"x"\t"p.C379*"\t"Trunc"\t2\t"Tumor Suppressor"\t"Likely Loss-of-function"\t"d"\t"g"\t"v"\tNA\tNA\tNA\n')
    cp = os.path.join(tmp, "c.tsv")
    with open(cp, "w") as fh:
        w = csv.DictWriter(fh, HDR, delimiter="\t"); w.writeheader(); [w.writerow(c) for c in cands]
    ep = os.path.join(data, "P1", "filtered.events.json")
    json.dump(events, open(ep, "w"))
    return data, scan, cp, ep


class T(unittest.TestCase):
    def test_field_mapping(self):
        with tempfile.TemporaryDirectory() as t:
            data, scan, cp, ep = make(t, [cand()], [])
            res = G.run(data, cp, scan=scan, tag="bak_x", log=lambda *a: None)
            r = json.load(open(ep))[0]
            self.assertEqual((r["gene"], r["vartype"], r["type"], r["Variant"], r["Tier"], r["role"]), ("NF1", "SNV", "Trunc", "p.C379*", 2, "Tumor Suppressor"))
            self.assertEqual((r["Genome_Location"], r["Variant_g"], r["seqnames"], r["start"], r["end"]), ("17:100-100", "17:100-100 C>A", "17", "100", "100"))
            self.assertEqual((r["alt"], r["ref"], r["n_cells"], r["n_normal_cells"], r["cell_ids"], r["source"]), (30, 10, 2, 1, "c1,c2", "pooled driver scan"))
            self.assertIsNone(r["therapeutics"]); self.assertEqual(r["evidence"]["scan_sha"], G.SCAN_SHA); self.assertEqual(r["evidence"]["clade_best"], "P1_n1")
            self.assertTrue(res[0]["row"]["_check"]["cell_ids_in_patient"])

    def test_no_annotation_is_null(self):
        row = G.build_row(cand(gene="X", protein=""), ["c1"], [])
        self.assertIsNone(row["Tier"]); self.assertIsNone(row["role"]); self.assertIsNone(row["effect_description"])

    def test_truncating_tsg_rule_tier(self):
        # no OncoKB entry: a truncating mutation in a listed tumour suppressor gets tier 2 by rule, labelled; other genes stay null
        d = G.build_row(cand(gene="TP53", type="DEL", ref="TTC", alt="T", consequence="frameshift", protein=""), ["c1"], [])
        self.assertEqual((d["Tier"], d["effect"], d["role"]), (2, "Likely Loss-of-function", "Tumor Suppressor"))
        self.assertIn("OncoKB rule", d["annotation_rule"])
        x = G.build_row(cand(gene="X", type="DEL", ref="TTC", alt="T", consequence="frameshift", protein=""), ["c1"], [])
        self.assertIsNone(x["Tier"]); self.assertIsNone(x["annotation_rule"])

    def test_tert_and_indel(self):
        t = G.build_row(cand(gene="TERT", chrom="chr5", pos="1295113", ref="G", alt="A", consequence="flank/unknown", protein=""), ["c1"], [])
        self.assertEqual((t["vartype"], t["type"], t["Variant"], t["Tier"]), ("Hotspot", "Promoter", "C228T (promoter, -124 C>T)", 1))
        self.assertTrue(t["rescued"])
        d = G.build_row(cand(type="DEL", ref="TTC", alt="T", consequence="frameshift", protein=""), ["c1"], [], snp=dict(consequence="frameshift_variant", hgvs_c="c.626_627delGA", hgvs_p="p.Arg209fs"))
        self.assertEqual((d["vartype"], d["type"], d["consequence"], d["Variant"], d["end"]), ("INDEL", "Trunc", "Frame_Shift_Del", "p.R209fs", "102"))
        s = G.build_row(cand(consequence="splice acceptor", protein=""), ["c1"], [], snp=dict(consequence="splice_acceptor_variant", hgvs_c="c.994-1G>C", hgvs_p=None))
        self.assertEqual((s["type"], s["Variant"]), ("Splice", "c.994-1G>C"))

    def test_idempotent_and_skip_existing(self):
        with tempfile.TemporaryDirectory() as t:
            other = {"gene": "NF1", "Genome_Location": "17:100-100", "Variant_g": "17:100-100 C>A", "source": None}
            data, scan, cp, ep = make(t, [cand()], [other, {"gene": "Z", "source": "x"}])
            G.run(data, cp, scan=scan, tag="bak_x", log=lambda *a: None)
            self.assertEqual(len(json.load(open(ep))), 2)       # existing row -> skipped
            data, scan, cp, ep = make(t + "/b", [cand()], [{"gene": "Z", "source": "x"}]) if os.makedirs(t + "/b") is None else None
            for _ in range(3):
                G.run(data, cp, scan=scan, tag="bak_x", log=lambda *a: None)
            ev = json.load(open(ep))
            self.assertEqual([e["source"] for e in ev], ["x", "pooled driver scan"])
            self.assertTrue(os.path.exists(ep + ".bak_x"))
            self.assertEqual(json.load(open(ep + ".bak_x")), [{"gene": "Z", "source": "x"}])

    def test_dry_run_and_filter(self):
        with tempfile.TemporaryDirectory() as t:
            data, scan, cp, ep = make(t, [cand(), cand(gene="Y", pos="5", verdict="review")], [])
            before = open(ep).read()
            res = G.run(data, cp, dry=True, scan=scan, log=lambda *a: None)
            self.assertEqual(open(ep).read(), before); self.assertFalse(os.path.exists(ep + ".bak_x"))
            self.assertEqual(len(res), 1)


if __name__ == "__main__":
    unittest.main()
