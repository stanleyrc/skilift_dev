"""Synthetic fixtures only (no patient data): python3 -m unittest test_amp_events"""
import json
import os
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gos_sc_amp_events as a  # noqa: E402

LV = dict(bin_size=100000, n_bins=40, missing_value=65535, files={"cn": "cn.u16"},
          chromosomes=[dict(chr="1", first_bin=0, n_bins=20), dict(chr="2", first_bin=20, n_bins=20)])


def fixture(n=12, amp_cells=8, cn_amp=9):
    cn = np.full((n, 40), 2.0)
    cn[:amp_cells, 5:8] = cn_amp          # chr1 0.5-0.8 Mb amplicon
    cn[:, 30:] = 4.0                      # chr2 half gain (broad): must not count
    return cn


GENES = [dict(gene="AMPG", chrom="1", start=550001, end=700000), dict(gene="FLAT", chrom="1", start=1500001, end=1600000),
         dict(gene="BROAD", chrom="2", start=3100001, end=3300000)]


class T(unittest.TestCase):
    def test_amplified_genes(self):
        r = a.amplified_genes(fixture(), LV, [False] * 12, GENES)
        self.assertEqual(list(r), ["AMPG"])
        self.assertEqual(len(r["AMPG"]["cells"]), 8)
        self.assertEqual(r["AMPG"]["segment_cn"], 9)

    def test_min_cells_and_normal(self):
        self.assertEqual(a.amplified_genes(fixture(amp_cells=4), LV, [False] * 12, GENES), {})
        normal = [True] * 4 + [False] * 8      # 4 of the 8 carriers are normal cells
        self.assertEqual(a.amplified_genes(fixture(), LV, normal, GENES), {})

    def test_low_cn_not_amp(self):
        self.assertEqual(a.amplified_genes(fixture(cn_amp=4), LV, [False] * 12, GENES), {})

    def _row(self, n_cells=8):
        info = a.amplified_genes(fixture(amp_cells=n_cells), LV, [False] * 12, GENES)["AMPG"]
        return a.make_row("PX", "AMPG", info, ["c%d" % i for i in range(12)], 12)

    def test_row_fields(self):
        r = self._row()
        self.assertEqual((r["vartype"], r["type"], r["n_cells"], r["cells"], r["segment_cn"], r["estimated_altered_copies"]), ("AMP", "SCNA", 8, "8/12", 9, 7))
        self.assertEqual(r["Genome_Location"], "1:550001-700000")

    def test_merge_adds_and_is_idempotent_and_keeps_other_rows(self):
        other = dict(gene="TP53", type="Missense", vartype="Missense", n_cells=3)
        new, s = a.merge_events([other], [self._row()])
        self.assertEqual(new[0], other)
        self.assertEqual(s, [("AMPG", "added", 8)])
        again, s2 = a.merge_events(new, [self._row()])
        self.assertEqual(again, new)
        self.assertEqual(s2, [("AMPG", "extended", 8)] if False else [])  # same source and cells: no-op

    def test_merge_extends_existing_amp_row(self):
        old = dict(gene="AMPG", type="SCNA", vartype="AMP", cell_ids="c0,c1,c11", n_cells=3, cells="3/12", segment_cn=5,
                   estimated_altered_copies=3, role="Oncogene", Tier=2)
        new, s = a.merge_events([old], [self._row()])
        self.assertEqual(len(new), 1)
        self.assertEqual(new[0]["n_cells"], 9)          # c0,c1 shared; c11 kept; +6 new
        self.assertEqual(new[0]["role"], "Oncogene")    # annotation of the live row kept
        self.assertEqual(new[0]["segment_cn"], 9)

    def test_export_roundtrip_and_base_hash(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "data", "PX"); os.makedirs(p); os.makedirs(os.path.join(d, "data", "_cohort"))
            fixture().astype("<u2").tofile(os.path.join(p, "cn.u16"))
            json.dump(dict(cells=["c%d" % i for i in range(12)], normal=[False] * 12, levels={"100kb": LV}), open(os.path.join(p, "cn_binned.json"), "w"))
            json.dump(dict(genes=GENES), open(os.path.join(d, "data", "_cohort", "gene_spans.json"), "w"))
            json.dump(dict(amplicons=[dict(label="AMPG,FLAT")]), open(os.path.join(p, "ecdna_status.json"), "w"))
            json.dump([dict(gene="TP53", type="Missense")], open(os.path.join(p, "filtered.events.json"), "w"))
            r = a.export_patient(os.path.join(d, "data"), os.path.join(d, "out"), "PX")
            self.assertEqual([x["gene"] for x in r["rows"]], ["AMPG"])
            self.assertEqual(r["base_sha256"], a.sha256(os.path.join(p, "filtered.events.json")))
            out = os.path.join(d, "merged.json")
            s = a.merge_cli(os.path.join(p, "filtered.events.json"), os.path.join(d, "out", "PX", "amp_events_added.json"), out)
            self.assertEqual(s, [("AMPG", "added", 8)])
            self.assertEqual(len(json.load(open(out))), 2)

    def test_export_skips_gene_already_called(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "data", "PX"); os.makedirs(p); os.makedirs(os.path.join(d, "data", "_cohort"))
            fixture().astype("<u2").tofile(os.path.join(p, "cn.u16"))
            json.dump(dict(cells=["c%d" % i for i in range(12)], normal=[False] * 12, levels={"100kb": LV}), open(os.path.join(p, "cn_binned.json"), "w"))
            json.dump(dict(genes=GENES), open(os.path.join(d, "data", "_cohort", "gene_spans.json"), "w"))
            json.dump(dict(amplicons=[dict(label="AMPG")]), open(os.path.join(p, "ecdna_status.json"), "w"))
            json.dump([dict(gene="AMPG", type="SCNA", vartype="AMP", n_cells=9, cell_ids="c0")], open(os.path.join(p, "filtered.events.json"), "w"))
            r = a.export_patient(os.path.join(d, "data"), os.path.join(d, "out"), "PX")
            self.assertEqual(r["rows"], [])


if __name__ == "__main__":
    unittest.main()
