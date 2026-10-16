import unittest, gos_sc_tree_support as g

class T(unittest.TestCase):
    def test_parse_leaves(self):
        t = g.parse_newick('((a:1,b:2)90:0.5,(c,d)x,e);')
        self.assertEqual(sorted(g.leaves(t)), list('abcde'))

    def test_fbp_and_tbe(self):
        tree = g.parse_newick('(((a,b),c),(d,(e,f)));')
        reps = [g.parse_newick(x) for x in [
            '(((a,b),c),(d,(e,f)));',   # identical
            '(((a,b),c),(d,(e,f)));',
            '(((a,c),b),(d,(e,f)));',   # loses ab
            '(((a,b),d),(c,(e,f)));',   # loses abc|def, keeps ab, ef
        ]]
        r = g.support(tree, reps)
        self.assertEqual(r['n_leaves_shared'], 6)
        by = {s['size']: s for s in r['splits']}
        k_ab = g.split_key(['c', 'd', 'e', 'f'])
        s_ab = [s for s in r['splits'] if s['key'] == k_ab][0]
        self.assertAlmostEqual(s_ab['fbp'], 3 / 4)
        k_ef = g.split_key(['e', 'f'])
        self.assertAlmostEqual([s for s in r['splits'] if s['key'] == k_ef][0]['fbp'], 1.0)
        # canonical side excludes ref_leaf 'a': abc|def -> def
        k_def = g.split_key(['d', 'e', 'f'])
        s_def = [s for s in r['splits'] if s['key'] == k_def][0]
        self.assertAlmostEqual(s_def['fbp'], 3 / 4)
        for s in r['splits']:
            self.assertTrue(0 <= s['tbe'] <= 1)
            self.assertGreaterEqual(s['tbe'], 0.0)
        self.assertGreaterEqual(s_ab['tbe'], s_ab['fbp'] - 1e-9)

    def test_leaf_restriction(self):
        tree = g.parse_newick('(((a,b),c),((d,x),(e,f)));')
        reps = [g.parse_newick('(((a,b),c),(d,(e,f)));')]
        r = g.support(tree, reps)
        self.assertEqual(r['excluded_leaves'], ['x'])
        self.assertEqual(r['n_leaves_shared'], 6)
        self.assertTrue(all(s['fbp'] == 1.0 for s in r['splits']))

    def test_key_stable(self):
        self.assertEqual(g.fnv1a32('a'), 'e40c292c')
        self.assertEqual(g.split_key(['b', 'a']), g.split_key(['a', 'b']))

if __name__ == '__main__':
    unittest.main()
