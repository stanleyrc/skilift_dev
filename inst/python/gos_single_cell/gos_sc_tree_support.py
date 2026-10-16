#!/usr/bin/env python3
"""Bootstrap support (Felsenstein bootstrap proportion, FBP, and optionally
transfer bootstrap expectation, TBE) for every tree shown in the gOS Compare Trees tab.

Output: <data>/<P>/trees/support.json

Encoding (so the front end can map splits onto drawn nodes without shipping leaf lists):
  * all tip names are normalised first: lower case, '-' -> '_' (the front end's tipKey).
  * shared leaf set S = leaves of the tree  intersect  leaves of the bootstrap replicates
    (tree leaves missing from the replicates are listed in `excluded_leaves`).
  * ref_leaf = lexicographically smallest name in S.
  * every internal edge is one bipartition of S. Its canonical side is the one NOT
    containing ref_leaf. key = "<size>:<fnv1a32 hex>" where size = |canonical side| and the hash is
    FNV-1a 32-bit over the UTF-8 bytes of the canonical side's names, sorted (code-point order)
    and joined with "|". Trivial splits (side size < 2 or > |S|-2) are skipped.
  * the root's two child edges give the same bipartition (same key; same value).
FBP = fraction (0-1) of replicates containing the bipartition (replicates restricted to S).
TBE = 1 - mean over replicates of min transfer distance / (size - 1), size = min side size.
"""
import argparse, glob, json, os, sys, subprocess

def norm(name):
    """Tip-name normalisation shared with the front end (tipKey): lower case, '-' -> '_'."""
    return name.lower().replace('-', '_')

def parse_newick(s):
    """Return nested dict nodes {name, children}. Handles quoted-free labels, lengths, internal labels."""
    s = s.strip().rstrip(';')
    root = {'name': None, 'children': []}
    stack = [root]
    i, n = 0, len(s)
    cur = root
    while i < n:
        c = s[i]
        if c == '(':
            node = {'name': None, 'children': []}
            stack[-1]['children'].append(node)
            stack.append(node)
            i += 1
        elif c == ',':
            i += 1
        elif c == ')':
            cur = stack.pop()
            i += 1
            j = i
            while j < n and s[j] not in ',():;':
                j += 1
            i = j  # internal label ignored
            if i < n and s[i] == ':':
                j = i + 1
                while j < n and s[j] not in ',()':
                    j += 1
                i = j
        else:
            j = i
            while j < n and s[j] not in ',():;':
                j += 1
            name = norm(s[i:j])
            i = j
            if i < n and s[i] == ':':
                j = i + 1
                while j < n and s[j] not in ',()':
                    j += 1
                i = j
            stack[-1]['children'].append({'name': name, 'children': []})
    # root wrapper holds exactly the outer node
    return root['children'][0] if len(root['children']) == 1 and root['children'][0]['children'] else root

def leaves(node):
    out, st = [], [node]
    while st:
        x = st.pop()
        if x['children']:
            st.extend(x['children'])
        else:
            out.append(x['name'])
    return out

def clade_leafsets(node):
    """List of leaf-name lists for every internal non-root node (post-order, iterative)."""
    res = []
    def rec(x):
        if not x['children']:
            return [x['name']]
        acc = []
        for ch in x['children']:
            acc.extend(rec(ch))
        res.append(acc)
        return acc
    sys.setrecursionlimit(100000)
    rec(node)
    res.pop()  # last appended is the root (all leaves)
    return res

def fnv1a32(text):
    h = 0x811c9dc5
    for b in text.encode('utf-8'):
        h ^= b
        h = (h * 0x01000193) & 0xffffffff
    return '%08x' % h

def split_key(names):
    nm = sorted(names)
    return '%d:%s' % (len(nm), fnv1a32('|'.join(nm)))

def tree_splits(tree, shared):
    """dict canonical-mask -> sorted names of canonical side. Bit i = shared sorted index i; bit 0 = ref."""
    idx = {nm: i for i, nm in enumerate(shared)}
    full = (1 << len(shared)) - 1
    n = len(shared)
    out = {}
    for ls in clade_leafsets(tree):
        m = 0
        for nm in ls:
            if nm in idx:
                m |= 1 << idx[nm]
        if m & 1:
            m = full ^ m
        k = bin(m).count('1')
        if k < 2 or k > n - 2:
            continue
        out[m] = None
    return out

def mask_names(m, shared):
    return [shared[i] for i in range(len(shared)) if m >> i & 1]

def support(tree, reps, shared_extra=None, do_tbe=True):
    tl = set(leaves(tree))
    rl = set(leaves(reps[0]))
    for r in reps[1:]:
        rl &= set(leaves(r))
    shared = sorted(tl & rl)
    n = len(shared)
    res = {'n_tree_leaves': len(tl), 'n_leaves_shared': n,
           'excluded_leaves': sorted(tl - rl), 'ref_leaf': shared[0] if shared else None, 'splits': []}
    if n < 4:
        return res
    ts = tree_splits(tree, shared)
    rep_sets = [set(tree_splits(r, shared)) for r in reps]
    nrep = len(reps)
    rep_lists = [list(s) for s in rep_sets]
    pop = lambda x: bin(x).count('1')
    for m in sorted(ts):
        size = pop(m)
        fbp = sum(1 for s in rep_sets if m in s) / nrep
        d = {'key': split_key(mask_names(m, shared)), 'size': size, 'fbp': round(fbp, 4)}
        if do_tbe:
            tot = 0.0
            for lst in rep_lists:
                best = size - 1
                # transfer distance to the best matching replicate branch; compare both orientations
                for b in lst:
                    x = m ^ b
                    dist = min(pop(x), n - pop(x))
                    if dist < best:
                        best = dist
                        if best == 0:
                            break
                tot += best
            d['tbe'] = round(1 - tot / nrep / (size - 1), 4)
        res['splits'].append(d)
    return res

def read_trees(path):
    with open(path) as f:
        return [parse_newick(l) for l in f if l.strip()]

def find_candidates(analysis, patient):
    pat = os.path.join(analysis, patient, 'genotyping_pipeline', '**', '*.raxml.bootstraps')
    return sorted(glob.glob(pat, recursive=True))

def pick_run(cands, live):
    """Choose the bootstrap file whose sibling supportFBP tree best matches the live tree."""
    best = None
    for b in cands:
        base = b[:-len('.bootstraps')]
        sf = base + '.supportFBP'
        if not os.path.exists(sf):
            continue
        st = read_trees(sf)[0]
        shared = sorted(set(leaves(live)) & set(leaves(st)))
        if len(shared) < 4:
            continue
        a, c = set(tree_splits(live, shared)), set(tree_splits(st, shared))
        frac = len(a & c) / max(1, len(a))
        score = (len(shared), frac)
        if best is None or score > best[0]:
            best = (score, b, sf, frac, len(shared))
    return best

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--patient', required=True)
    ap.add_argument('--data', default='/gpfs/commons/projects/imielinski_web/sclarke/gOS_GBM_sc/data')
    ap.add_argument('--analysis', default='/gpfs/commons/groups/imielinski_lab/projects/single_cell_GBM/analysis')
    ap.add_argument('--bootstraps', help='override: explicit .raxml.bootstraps file')
    ap.add_argument('--out', help='override output path')
    ap.add_argument('--min-agreement', type=float, default=0.1, help='below this fraction of live-tree clades also in the replicates\' search tree, support is declared not applicable')
    ap.add_argument('--no-tbe', action='store_true')
    ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args()
    P = a.patient
    ddir = os.path.join(a.data, P)
    out = a.out or os.path.join(ddir, 'trees', 'support.json')
    cands = [a.bootstraps] if a.bootstraps else find_candidates(a.analysis, P)
    if not cands:
        doc = {'format': 'gos-tree-support-v1', 'available': False, 'reason': 'no bootstrap replicates'}
    else:
        live = read_trees(os.path.join(ddir, 'tree.nwk'))[0]
        pick = pick_run(cands, live)
        if pick is None:
            pick = (None, cands[0], None, None, None)
        _, bfile, sfile, frac, nsh = pick
        reps = read_trees(bfile)
        if frac is not None and frac < a.min_agreement:
            doc = {'format': 'gos-tree-support-v1', 'available': False,
                   'reason': 'bootstrap replicates come from a different genotyping run (only %.0f%% of the live tree clades occur in their search tree)' % (100 * frac),
                   'source': {'bootstraps': bfile, 'search_tree_with_support': sfile}}
            return finish(a, out, doc, P)
        trees = {'gos': os.path.join(ddir, 'tree.nwk')}
        idx_path = os.path.join(ddir, 'trees', 'index.json')
        if os.path.exists(idx_path):
            for e in json.load(open(idx_path)):
                if e['name'] != 'gos':
                    trees[e['name']] = os.path.join(ddir, 'trees', e['name'] + '.nwk')
        doc = {'format': 'gos-tree-support-v1', 'available': True,
               'encoding': 'split key = "<size>:<fnv1a32 of sorted canonical-side names joined by |>"; canonical side excludes ref_leaf; see gos_sc_tree_support.py',
               'source': {'bootstraps': bfile, 'search_tree_with_support': sfile},
               'n_replicates': len(reps),
               'live_tree_split_agreement_with_source_search_tree': None if frac is None else round(frac, 4),
               'trees': {}}
        for name, path in trees.items():
            if not os.path.exists(path):
                continue
            t = read_trees(path)[0]
            doc['trees'][name] = support(t, reps, do_tbe=not a.no_tbe)
    finish(a, out, doc, P)

def finish(a, out, doc, P):
    if a.dry_run:
        print(P, doc.get('available'), doc.get('n_replicates'), {k: (v['n_leaves_shared'], len(v['splits'])) for k, v in doc.get('trees', {}).items()})
        return
    os.makedirs(os.path.dirname(out), exist_ok=True)
    if os.path.exists(out):
        os.replace(out, out + '.bak_261016')
    with open(out, 'w') as f:
        json.dump(doc, f, separators=(',', ':'))
    print('wrote', out)

if __name__ == '__main__':
    main()
