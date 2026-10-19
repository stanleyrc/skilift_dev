"""Unit tests for gos_sc_event_families.py (synthetic fixtures only)."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gos_sc_event_families as F  # noqa: E402


def graph(intervals, conns):
    return {"intervals": [dict(iid=i, chromosome=c, startPoint=s, endPoint=e) for i, c, s, e in intervals],
            "connections": [dict(cid=k, source=a, sink=b, type="ALT", title="X", weight=1) for k, a, b in conns]}


# master strand letters are the opposite of the graph's: a graph DEL (+-) is -+ in the master
MASTER = {"1:100-|1:5000+": dict(chrom1="chr1", pos1=100, strand1="-", chrom2="chr1", pos2=5000, strand2="+"),
          "1:100+|2:70-": dict(chrom1="chr1", pos1=100, strand1="+", chrom2="chr2", pos2=70, strand2="-")}


class Convention(unittest.TestCase):
    def test_both_ends_use_sign_of_id(self):
        g = graph([(1, "1", 1, 100), (2, "1", 5000, 9000), (3, "1", 200, 4999)], [(1, 1, 2), (2, -2, -3), (3, 1, 3)])
        a = F.cell_alts(g)
        self.assertEqual(a[0][1:3], (("1", 100, "+"), ("1", 9000, "+")))   # source>0 end/+, sink>0 end/+
        self.assertEqual(a[1][1:3], (("1", 5000, "-"), ("1", 200, "-")))   # negative ids: start/-
        self.assertEqual(a[2][1:3], (("1", 100, "+"), ("1", 4999, "+")))


class Matching(unittest.TestCase):
    def setUp(self):
        self.idx = F.build_index(MASTER)

    def test_exact_and_within_tol_and_swapped_order(self):
        self.assertEqual(F.match_junction(("1", 100, "+"), ("1", 5000, "-"), MASTER, self.idx), "1:100-|1:5000+")
        self.assertEqual(F.match_junction(("1", 5600, "-"), ("1", 700, "+"), MASTER, self.idx), "1:100-|1:5000+")
        self.assertIsNone(F.match_junction(("1", 100, "+"), ("1", 6200, "-"), MASTER, self.idx))

    def test_strand_must_match(self):
        self.assertIsNone(F.match_junction(("1", 100, "-"), ("1", 5000, "+"), MASTER, self.idx))

    def test_cross_chromosome(self):
        self.assertEqual(F.match_junction(("1", 100, "-"), ("2", 70, "+"), MASTER, self.idx), "1:100+|2:70-")

    def test_three_cells_shared_and_graph_only_ids(self):
        # connections: DEL 1:100->1:5000 in cells a,b (a is off by 40 bp); a private junction 3:10->3:9000 in a and b
        ga = [(1, ("1", 100, "+"), ("1", 5000, "-"), "DEL-like", 1.0), (2, ("3", 10, "+"), ("3", 9000, "-"), "DEL-like", 1.0)]
        gb = [(1, ("1", 140, "+"), ("1", 5000, "-"), "DEL-like", 1.0), (2, ("3", 600, "+"), ("3", 9100, "-"), "DEL-like", 1.0)]
        gc = [(1, ("1", 100, "+"), ("1", 5000, "-"), "DEL-like", 0.0)]  # weight 0 ignored
        cells, um, titles, gm = F.match_cells({"a": ga, "b": gb, "c": gc}, MASTER)
        self.assertEqual(cells["a"]["1"], "1:100-|1:5000+")
        self.assertEqual(cells["b"]["1"], "1:100-|1:5000+")
        self.assertEqual(cells["c"], {})
        self.assertEqual(um, {"a": [2], "b": [2]})
        self.assertEqual(cells["a"]["2"], cells["b"]["2"])  # one graph-derived id shared by the two cells
        self.assertIn(cells["a"]["2"], gm)
        self.assertEqual(len(gm), 1)


class Linking(unittest.TestCase):
    def test_small_events_share_one_id(self):
        self.assertTrue(F.link({"a", "b"}, {"b", "c"}))
        self.assertFalse(F.link({"a"}, {"c"}))
        self.assertFalse(F.link(set(), {"c"}))

    def test_large_events_need_jaccard(self):
        big = {f"j{i}" for i in range(12)}
        F.BIG_RULE = "jaccard"
        self.assertFalse(F.link(big, {"j0", "j1"}))
        F.BIG_RULE = "containment"
        self.assertTrue(F.link(big, {"j0", "j1"}))
        self.assertTrue(F.link(big, {f"j{i}" for i in range(8)}))

    def test_components_and_twins(self):
        js = [{"a"}, {"a", "b"}, {"z"}, {"q"}, {"r"}]
        comps = F.families(js, twin_groups=[[3, 4]])
        self.assertEqual(sorted(map(sorted, comps)), [[0, 1], [2], [3, 4]])

    def test_footprint_union(self):
        self.assertEqual(F.merge_footprint([[("1", 1, 10), ("2", 5, 6)], [("1", 8, 20), ("1", 100, 200)]]),
                         [("1", 1, 20), ("1", 100, 200), ("2", 5, 6)])


class EventJunctions(unittest.TestCase):
    def test_calling_cell_plus_shared_carrier_ids(self):
        alts = {"a": [(1, ("1", 100, "+"), ("1", 900, "-"), "t", 1.0), (2, ("1", 100, "+"), ("5", 50, "-"), "t", 1.0)],
                "b": [(1, ("1", 100, "+"), ("1", 900, "-"), "t", 1.0)]}
        cm = {"a": {"1": "J1", "2": "J2"}, "b": {"1": "J1"}}
        ev = {"Genome_Location": "1:50-1000", "sample": "a", "cell_ids": "a,b"}
        ids, nun, shared = F.event_junctions(ev, alts, cm)
        self.assertEqual(ids, ["J1"])   # J2 leaves the footprint
        self.assertEqual(shared, ["J1"])

    def test_build_families_adds_twin_and_counts(self):
        ev = [dict(type="Complex SV", vartype="TIC", Variant="TIC (1,2)", cell_ids="a,b", n_cells=2, Genome_Location="1:1-9", twin_group="P:1"),
              dict(type="Complex SV", vartype="rDelDup", Variant="rDelDup (1,2)", cell_ids="a,b", n_cells=2, Genome_Location="1:1-9", twin_group="P:1"),
              dict(type="Complex SV", vartype="DM", Variant="DM (1)", cell_ids="c", n_cells=1, Genome_Location="3:1-9")]
        fams = F.build_families("P", ev, [0, 1, 2], {0: ["x"], 1: ["y"], 2: []}, {})
        self.assertEqual(len(fams), 2)
        self.assertEqual(fams[0]["n_cells"], 2)
        self.assertEqual(fams[0]["classifications"], {"TIC (1,2)": 2, "rDelDup (1,2)": 2})
        self.assertFalse(fams[1]["resolved"])


class Guard(unittest.TestCase):
    def test_tiny_event_does_not_bridge_and_needs_3_shared_cells(self):
        # big1 {a,b} and big2 {b2,c}; tiny (1 cell) holds a and c -> must NOT merge big1 and big2
        cells = {f"c{i}" for i in range(5)}
        js = [{"a", "x"}, {"c", "y"}, {"a", "c"}]
        carr = [set(cells), set(cells), {"c0"}]
        has = {c: {"a", "x", "c", "y"} for c in cells}
        comps, links = F.guarded_families(js, carr, has)
        fam_of = {i: k for k, cp in enumerate(comps) for i in cp}
        self.assertNotEqual(fam_of[0], fam_of[1])
        self.assertIn(fam_of[2], (fam_of[0], fam_of[1]))  # attached to one family only
        self.assertTrue(any(l["kind"] == "weak_attach" for l in links))

    def test_shared_id_needs_three_common_carrier_cells(self):
        js = [{"a"}, {"a"}]
        carr = [{"c0", "c1", "c2", "c3"}, {"c2", "c3", "c4", "c5"}]
        has = {c: {"a"} for c in ["c0", "c1", "c2", "c3", "c4", "c5"]}
        comps, _ = F.guarded_families(js, carr, has)  # only 2 common cells
        self.assertEqual(len(comps), 2)
        carr[1] = {"c1", "c2", "c3", "c5"}
        comps, links = F.guarded_families(js, carr, has)
        self.assertEqual(len(comps), 1)
        self.assertEqual(links[0]["n_cells"], 3)


if __name__ == "__main__":
    unittest.main()
