"""Unit tests for gos_sc_arm_events.py and gos_sc_junctions_table.py (synthetic fixtures only)."""
import io
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gos_sc_arm_events as A  # noqa: E402
import gos_sc_junctions_table as J  # noqa: E402

CYTO = "chromosome\tstartPoint\tendPoint\ttitle\tstain\n1\t0\t5000000\tp12\tgneg\n1\t5000000\t10000000\tq11\tgneg\n2\t0\t4000000\tp12\tgneg\n2\t4000000\t10000000\tq11\tgneg\nX\t0\t9\tp1\tgneg\n"


def cn_fixture(rows, normal=None, clones=None):
    # chromosome 1 and 2, 10 Mb each at 1 Mb bins
    starts = [i * 1_000_000 + 1 for i in range(10)] * 2
    return {"cells": [f"c{i}" for i in range(len(rows))], "normal": normal or [False] * len(rows), "clone_id": clones or ["Clone 1"] * len(rows),
            "levels": {"1Mb": {"bin_chr": ["1"] * 10 + ["2"] * 10, "bin_start": starts, "bin_end": [s + 999_999 for s in starts], "cn": rows}}}


def arms():
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
        fh.write(CYTO)
    try:
        return A.read_arms(fh.name)
    finally:
        os.unlink(fh.name)


class ArmTests(unittest.TestCase):
    def test_arms_skip_sex_chromosomes(self):
        self.assertEqual(sorted(arms()), ["1", "2"])
        self.assertEqual(arms()["1"]["p"], (0, 5_000_000))

    def test_gain_loss_and_chromosome_level(self):
        base = [2] * 20
        gain1p = [3] * 5 + [2] * 5 + [2] * 10
        loss_chr2 = [2] * 10 + [1] * 10
        ev, n, _ = A.arm_events(cn_fixture([base, gain1p, loss_chr2, loss_chr2]), arms())
        self.assertEqual(n, 4)
        by = {e["id"]: e for e in ev}
        self.assertEqual(by["1p_gain"]["cells"], ["c1"])
        self.assertEqual(by["2p_loss"]["n_cells"], 2)
        self.assertEqual(by["2_loss"]["level"], "chromosome")
        self.assertEqual(by["2_loss"]["cell_fraction"], 0.5)
        self.assertNotIn("1q_gain", by)

    def test_normal_cells_excluded_and_baseline_shifts_with_ploidy(self):
        tri = [3] * 20
        ev, n, _ = A.arm_events(cn_fixture([tri, [2] * 20, tri], normal=[False, True, False]), arms())
        self.assertEqual((n, ev), (2, []))  # a triploid cell is not a gain of every arm
        # arm 2p one copy above a cell that is otherwise CN 3
        ev2, _, _ = A.arm_events(cn_fixture([[3] * 10 + [4] * 4 + [3] * 6]), arms())
        self.assertEqual([e['id'] for e in ev2], ['2p_gain'])

    def test_clone_normal_label_excluded(self):
        ev, n, _ = A.arm_events(cn_fixture([[2] * 10 + [1] * 10, [2] * 20], clones=["Normal", "Clone 1"]), arms())
        self.assertEqual((n, ev), (1, []))


V3 = "jid\tchrA\tposA\tstrA\tchrB\tposB\tstrB\tclass\tspan\tn_reads_normal\tn_cells_normal\tn_cells_tumour\tn_reads_tumour\tsplice_mimic\tpath\ttierlike\tcalled_v3\n" \
     "1\tchr7\t55000100\t+\tchr7\t55300100\t-\tDEL-like\t300000\t0\t0\t12\t40\tFALSE\toriginal\tTRUE\tTRUE\n" \
     "2\tchr1\t10\t+\tchr3\t20\t+\tTRA\tNA\t3\t1\t1\t3\tFALSE\tsingle_cell_k3\tFALSE\tTRUE\n"


class JunctionTests(unittest.TestCase):
    idx = J.gene_index([{"gene": "EGFR", "chrom": "7", "start": 55_000_000, "end": 55_200_000}, {"gene": "LANCL2", "chrom": "7", "start": 55_250_000, "end": 55_350_000}])

    def test_backbone_scope_keeps_backbone_only_and_annotates_genes(self):
        r = J.rows(io.StringIO(V3), self.idx)
        self.assertEqual(len(r), 1)
        self.assertEqual((r[0]["genes"], r[0]["n_cells_tumour"], r[0]["n_reads_normal"], r[0]["tier"]), ("EGFR,LANCL2", 12, 0, "original"))
        self.assertEqual(r[0]["id"], "chr7:55000100+|chr7:55300100-")

    def test_called_v3_scope_and_missing_span(self):
        r = J.rows(io.StringIO(V3), self.idx, "called_v3")
        self.assertEqual(len(r), 2)
        self.assertEqual((r[1]["span"], r[1]["genes"], r[1]["n_reads_normal"]), (None, "", 3))


if __name__ == "__main__":
    unittest.main()
