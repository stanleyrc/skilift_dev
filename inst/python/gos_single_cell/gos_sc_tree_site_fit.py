#!/usr/bin/env python3
"""Clade-violating mutations per tree (vignette V1) -> <out>/<P>/tree_site_fit.json.

For every tree in <P>/ (tree.nwk = "current", every trees/*.nwk) and every SNV site of the CURRENT matrix <P>/snv_matrix.json that has at
least `--min-carriers` carrier cells (genotype 1) on the tree: the best-matching-clade F1 = max over clades (internal nodes and leaves) of
2 * TP / (clade size + carriers), as gOS cladeFit.js and the V1 figure; flag `violating` = F1 < 0.5, `consistent` = F1 >= 0.8.
Tips = tree tips that are cells of the matrix. Also the SBS96 channel counts (context of snv_matrix.json, pyrimidine-referenced) of violating,
intermediate and consistent sites per tree, and the % of sites fitting a clade, so old and new trees compare on the same matrix.
Layout: sites = columnar (id, category, context, gene, n_carriers per tree); trees[name] = {tips, f1[], best_clade[], violating[]}.

usage: gos_sc_tree_site_fit.py [--data DIR] [--out DIR] [--patients P ...] [--min-carriers 3]
"""
import argparse
import glob
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gos_sc_figlib as F  # noqa: E402

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
VIOL, OK = 0.5, 0.8


def sbs96(context):
    """'T[T>C]C' -> pyrimidine-referenced channel; None when malformed."""
    m = re.match(r"^([ACGT])\[([ACGT])>([ACGT])\]([ACGT])$", context or "")
    if not m:
        return None
    l, r, a, rr = m.groups()
    if r in "AG":
        l, r, a, rr = COMP[rr], COMP[r], COMP[a], COMP[l]
    return "%s[%s>%s]%s" % (l, r, a, rr)


def best_clade_f1(layout, tips, carrier_matrix):
    """layout leaves restricted to `tips` (list in the order of carrier_matrix columns); returns (f1, best size) per site.
    carrier_matrix: sites x tips bool. Clades are the internal nodes and leaves of the tree pruned to tips."""
    idx = {t: i for i, t in enumerate(tips)}
    leaf_rows = [idx.get(l, -1) for l in layout["leaves"]]
    nL = len(layout["leaves"])
    # clade membership over the full leaf order, then restricted to tips
    sets = {}
    for nd in layout["nodes"]:
        members = [leaf_rows[k] for k in range(nd["firstLeaf"], nd["lastLeaf"] + 1) if leaf_rows[k] >= 0]
        if members:
            sets[tuple(members)] = members  # identical clades after pruning collapse
    M = np.zeros((len(sets), len(tips)), dtype=np.float32)
    for r, members in enumerate(sets.values()):
        M[r, members] = 1
    sizes = M.sum(1)
    C = carrier_matrix.astype(np.float32)
    k = C.sum(1)
    f1 = np.zeros(len(C))
    best = np.zeros(len(C), dtype=np.int64)
    for a in range(0, len(C), 2000):
        tp = C[a:a + 2000] @ M.T
        s = 2 * tp / (sizes[None, :] + k[a:a + 2000, None])
        j = s.argmax(1)
        f1[a:a + 2000] = s[np.arange(len(j)), j]
        best[a:a + 2000] = sizes[j].astype(np.int64)
    return f1, best


def export_patient(data, outdir, patient, min_carriers):
    sp = os.path.join(data, patient, "snv_matrix.json")
    if not os.path.exists(sp):
        print(patient, "no snv_matrix.json: not available", flush=True)
        return None
    snv = F.jload(sp)
    variants = snv["variants"]
    cells = list((snv.get("cells") or {}).keys())
    cidx = {c: i for i, c in enumerate(cells)}
    carriers = {}  # site -> set of cell idx
    for c, entries in snv["cells"].items():
        for e in entries:
            if len(e) >= 4 and e[3] == 1:
                carriers.setdefault(e[0], set()).add(cidx[c])
    trees = {}
    for name, path in [("current", os.path.join(data, patient, "tree.nwk"))] + [
            (os.path.splitext(os.path.basename(p))[0], p) for p in sorted(glob.glob(os.path.join(data, patient, "trees", "*.nwk")))]:
        if os.path.exists(path):
            trees[name] = path
    if "current" not in trees:
        print(patient, "no tree.nwk: not available", flush=True)
        return None
    res = {}
    for name, path in trees.items():
        with open(path) as fh:
            layout = F.layout_tree(F.parse_newick(fh.read()))
        tips = [t for t in layout["leaves"] if t in cidx]
        tip_cols = {cidx[t]: j for j, t in enumerate(tips)}
        n = np.zeros(len(variants), dtype=np.int64)
        rows = []
        for s in range(len(variants)):
            cols = [tip_cols[c] for c in carriers.get(s, ()) if c in tip_cols]
            n[s] = len(cols)
            rows.append(cols)
        use = np.where(n >= min_carriers)[0]
        Cm = np.zeros((len(use), len(tips)), dtype=bool)
        for r, s in enumerate(use):
            Cm[r, rows[s]] = True
        f1, best = best_clade_f1(layout, tips, Cm) if len(use) else (np.zeros(0), np.zeros(0, dtype=np.int64))
        res[name] = {"tips": len(tips), "n": n, "use": use, "f1": f1, "best": best, "path": path}
    cur_use = res["current"]["use"]
    sites = {"id": [], "category": [], "context": [], "gene": [], "driver": []}
    for s in range(len(variants)):
        v = variants[s]
        sites["id"].append(v["id"])
        sites["category"].append(v.get("category"))
        sites["context"].append(sbs96(v.get("context")))
        sites["gene"].append(v.get("gene") or "")
        sites["driver"].append(bool(v.get("driver")))
    out = {"format": "gos-sc-tree-site-fit/1", "patient": patient, "min_carriers": min_carriers, "f1_violating_below": VIOL,
           "f1_consistent_from": OK, "n_sites_matrix": len(variants), "sites": sites, "trees": {}}
    for name, r in res.items():
        use, f1, best = r["use"], r["f1"], r["best"]
        f1_full = [None] * len(variants)
        best_full = [None] * len(variants)
        for s, a, b in zip(use, f1, best):
            f1_full[s] = round(float(a), 4)
            best_full[s] = int(b)
        viol = [None if x is None else bool(x < VIOL) for x in f1_full]
        spec = {"violating": {}, "intermediate": {}, "consistent": {}}
        for s, a in zip(use, f1):
            ch = sites["context"][s]
            if ch is None:
                continue
            cls = "violating" if a < VIOL else ("consistent" if a >= OK else "intermediate")
            spec[cls][ch] = spec[cls].get(ch, 0) + 1
        nu = max(1, len(use))
        out["trees"][name] = {
            "tips": r["tips"], "source": r["path"], "n_sites_scored": int(len(use)),
            "pct_consistent": round(100.0 * float((f1 >= OK).sum()) / nu, 2) if len(use) else None,
            "pct_violating": round(100.0 * float((f1 < VIOL).sum()) / nu, 2) if len(use) else None,
            "n_carriers": r["n"].tolist(), "f1": f1_full, "best_clade": best_full, "violating": viol, "sbs96": spec,
        }
    pdir = os.path.join(outdir, patient)
    os.makedirs(pdir, exist_ok=True)
    with open(os.path.join(pdir, "tree_site_fit.json"), "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    print(patient, {k: (v["n_sites_scored"], v["pct_consistent"]) for k, v in out["trees"].items()}, flush=True)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", default=os.path.join(F.DATA_ROOT, "data"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-carriers", type=int, default=3)
    ap.add_argument("--patients", nargs="*", default=["BWH69", "BWH70", "BWH75", "MGH285", "MGH301", "MGH302", "MGH303", "MGH304", "MGH305"])
    a = ap.parse_args()
    for p in a.patients:
        export_patient(a.data, a.out, p, a.min_carriers)


if __name__ == "__main__":
    main()
