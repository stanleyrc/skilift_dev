"""Synthetic fixtures only (no patient data): python3 -m pytest test_jake_trees.py  (or python3 test_jake_trees.py)"""
import gzip
import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gos_sc_jake_trees as j  # noqa: E402

GOS = "((c_a:1,c_b:1):1,((c_c:1,c_d:1):1,(c_e:1,c_f:1):1):1);"
JAKE = "((x_A:1,x_B:1)0.9:1,((x_C:1,x_D:1)0.8:1,(x_E:1,x_F:1,x_G:1)0.7:2)0.6:1);"


def write_gz(path, text):
    with gzip.open(path, "wt") as f:
        f.write(text)


class JakeTrees(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.src = os.path.join(self.tmp, "src", "P1")
        self.live = os.path.join(self.tmp, "live")
        os.makedirs(self.src)
        os.makedirs(os.path.join(self.live, "P1", "trees"))
        open(os.path.join(self.live, "P1", "tree.nwk"), "w").write(GOS)
        json.dump([{"name": "gos", "label": "Pipeline (current)"}], open(os.path.join(self.live, "P1", "trees", "index.json"), "w"))
        open(os.path.join(self.tmp, "alias.tsv"), "w").write("alias\tgos_id\n" + "".join("x-%s\tc_%s\n" % (c, c.lower()) for c in "ABCDEF"))
        for suffix in ("mapped.TBE", "mapped.tumour.TBE"):
            open(os.path.join(self.src, "P1.%s.nwk" % suffix), "w").write(JAKE)
        cols = "\tx_A\tx_B\tx_C\tx_D\tx_E\tx_F\tx_G\n"
        write_gz(os.path.join(self.src, "P1.counts.alt.tsv.gz"), "id" + cols + "chr1:100_C>T\t1\t2\t0\t0\t0\t0\t5\nchr2:50_G>A\t0\t0\t3\t1\t0\t0\t0\n")
        write_gz(os.path.join(self.src, "P1.counts.depth.tsv.gz"), "id" + cols + "chr1:100_C>T\t3\t2\t4\t0\t6\t1\t9\nchr2:50_G>A\t2\t2\t3\t3\t0\t5\t0\n")
        hdr = "chrom\tpos\tref\talt\tclass\tnode\tclade_size\tfunc\tgene\texonic_func\taachange\tcell\tref_reads\talt_reads\tvaf\n"
        rows = ["chr2\t50\tG\tA\tsubclonal\tN3\t2\texonic\tGENE2\tnonsynonymous_SNV\tp.A1B\tx-C\t0\t3\t1.0\n",
                "chr1\t100\tC\tT\ttruncal\tN1\t6\tintronic\tGENE1\t.\t.\tx-A\t2\t1\t0.3\n",
                "chr1\t100\tC\tT\ttruncal\tN1\t6\tintronic\tGENE1\t.\t.\tx-B\t0\t2\t1.0\n"]
        write_gz(os.path.join(self.src, "P1.per_cell_long_table.tsv.gz"), hdr + "".join(rows))

    def test_newick_roundtrip_and_prune(self):
        root = j.parse_newick(JAKE)
        self.assertEqual(j.write_newick(root), JAKE)
        alias = j.read_alias(os.path.join(self.tmp, "alias.tsv"))
        new, n, matched, unmatched = j.rename_tree(root, alias)
        self.assertEqual((n, matched, unmatched), (7, ["c_a", "c_b", "c_c", "c_d", "c_e", "c_f"], ["x_G"]))
        self.assertIn("(c_e:1,c_f:1)0.7:2", j.write_newick(new))

    def test_unary_collapse_keeps_child_label_and_adds_length(self):
        root = j.parse_newick("((a:1,b:1)0.9:2,c:1,d:1);")
        new = j.prune(root, {"a", "c", "d"})
        self.assertEqual(j.write_newick(new), "(a:3.0,c:1,d:1);")

    def test_support_keys_match_gos_convention(self):
        # same key as treeSupport.js: canonical side excludes the reference leaf (shared[0])
        root = j.parse_newick("((c_a:1,c_b:1)0.9:1,((c_c:1,c_d:1)0.8:1,(c_e:1,c_f:1)0.7:1)0.6:1);")
        e = j.support_entry(root)
        self.assertEqual((e["n_leaves_shared"], e["ref_leaf"]), (6, "c_a"))
        keys = {s["key"]: s for s in e["splits"]}
        self.assertEqual(len(keys), 3)  # (a,b) is a split of size 4 on the other side, (c,d) and (e,f) size 2; the (c,d,e,f) clade is the complement of (a,b)
        self.assertIn(j.split_key(["c_c", "c_d"]), keys)
        self.assertIn(j.split_key(["c_c", "c_d", "c_e", "c_f"]), keys)
        self.assertEqual(keys[j.split_key(["c_c", "c_d"])]["tbe"], 0.8)
        self.assertTrue(all(s["fbp"] is None for s in e["splits"]))
        self.assertEqual(j.fnv1a32(""), "811c9dc5")

    def test_export_trees_and_matrix(self):
        out = os.path.join(self.tmp, "out")
        rep = j.run("P1", self.src, os.path.join(self.tmp, "alias.tsv"), self.live, out)
        idx = json.load(open(os.path.join(out, "P1", "trees", "index.json")))
        self.assertEqual([e["name"] for e in idx], ["gos", "jake_tbe_261026", "jake_tbe_tumour_261026"])
        self.assertEqual(rep["trees"]["jake_tbe_261026"]["unmatched"], ["x_G"])
        sup = json.load(open(os.path.join(out, "P1", "trees", "support.json")))
        self.assertIn("jake_tbe_261026", sup["trees"])
        m = json.load(open(os.path.join(out, "P1", j.SET_FILE)))
        self.assertEqual([v["id"] for v in m["variants"]], ["chr1_100_C_T", "chr2_50_G_A"])
        self.assertEqual(m["cells"]["c_a"], [[0, 2, 1, None], [1, 2, 0, None]])  # ref = depth - alt; zero-depth cells omitted
        self.assertEqual(m["cells"]["c_d"], [[1, 2, 1, None]])
        self.assertNotIn("x_G", m["cells"])
        sets = json.load(open(os.path.join(out, "P1", "snv_sets.json")))
        self.assertEqual([s["name"] for s in sets], ["full", j.SET_NAME])
        self.assertIsNone(sets[1]["n_variants"])  # keeps the full matrix the default set
        self.assertEqual(open(os.path.join(self.live, "P1", "tree.nwk")).read(), GOS)

    def test_refuses_symlinked_stage(self):
        out = os.path.join(self.tmp, "out")
        os.makedirs(os.path.join(out, "P1"))
        os.symlink(os.path.join(self.live, "P1", "trees"), os.path.join(out, "P1", "trees"))
        with self.assertRaises(RuntimeError):
            j.run("P1", self.src, os.path.join(self.tmp, "alias.tsv"), self.live, out)


if __name__ == "__main__":
    unittest.main()
