#!/usr/bin/env python3
"""Phylogenetic autocorrelation (Moran's I on the patient tree) of per-cell variables, paper Fig 3E -> <out>/<P>/phylo_signal.json.

A port of gOS figures.js phyloSignal + phyloSignalPanel.js: Moran's I with inverse patristic-distance weights between tumour cells
(tree.nwk, branch lengths summed to the lowest common ancestor), z-scored against 199 permutations of the values over the leaves
drawn with the same PRNG as the browser (mulberry32, seed 7), permutation p = (n_ge + 1) / (199 + 1). Variables, as in the panel:
  controls  Seq. depth (qc_depth), Ploidy, SNVs / cell (snv_count) from datafiles.json
  walks     log1p(copies) of every walk in a driver-gene group of walks.json (walkGroups defaults: curated walks, >= 3 carrier cells)
Extra, not in the gOS panel (source = "extra"):
  Mean somatic VAF   per cell, mean alt / (ref + alt) over its called tree sites (snv_matrix.json, genotype 1)
and for every variable `p_perm10000`, a permutation p from 10,000 random permutations (numpy, seed 7) for a finer p.
n = tumour cells (not Normal clones) on the tree; patients with < 6 such cells, or no tree.nwk, are not available.

usage: gos_sc_phylo_signal.py [--data DIR] [--out DIR] [--patients P ...]
"""
import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gos_sc_figlib as F  # noqa: E402
import gos_sc_walk_segment_cn as W  # noqa: E402  (walk_groups)

N_PERM = 199
SEED = 7
N_PERM_FINE = 10000
CONTROLS = [("Seq. depth", "qc_depth"), ("Ploidy", "ploidy"), ("SNVs / cell", "snv_count")]
M32 = 0xFFFFFFFF


def _imul(a, b):
    return ((a & M32) * (b & M32)) & M32


class Mulberry32:
    def __init__(self, seed):
        self.a = seed & M32

    def __call__(self):
        self.a = (self.a + 0x6D2B79F5) & M32
        t = self.a
        t = _imul(t ^ (t >> 15), t | 1)
        t = (t ^ ((t + _imul(t ^ (t >> 7), t | 61)) & M32)) & M32
        return ((t ^ (t >> 14)) & M32) / 4294967296.0


def leaf_distances(layout):
    """patristic distances between leaves (float32 n x n), each pair set at its LCA (gOS leafDistances)."""
    nodes = layout["nodes"]
    n = len(layout["leaves"])
    depth = np.zeros(n)
    for nd in nodes:
        if nd["isLeaf"]:
            depth[nd["firstLeaf"]] = nd["x"]
    D = np.zeros((n, n), dtype=np.float32)
    for nd in nodes:
        if nd["isLeaf"] or len(nd["children"]) < 2:
            continue
        kids = [nodes[c] for c in nd["children"]]
        for a in range(len(kids)):
            for b in range(a + 1, len(kids)):
                ia = np.arange(kids[a]["firstLeaf"], kids[a]["lastLeaf"] + 1)
                ib = np.arange(kids[b]["firstLeaf"], kids[b]["lastLeaf"] + 1)
                d = (depth[ia][:, None] + depth[ib][None, :] - 2 * nd["x"]).astype(np.float32)
                D[np.ix_(ia, ib)] = d
                D[np.ix_(ib, ia)] = d.T
    return D


def morans_i(values, Wm, n):
    mean = values.sum() / n
    d = values - mean
    den = float((d * d).sum())
    nz = Wm != 0
    wsum = float(Wm[nz].sum())
    num = float(d @ Wm @ d)
    return (n / wsum) * (num / den) if den > 0 and wsum > 0 else float("nan")


def phylo_signal(layout, cell_ids, variables, n_perm=N_PERM, seed=SEED, n_fine=N_PERM_FINE):
    """variables: [(key, {cell_id: value})]; returns [{key, I, z, p, n, p_perm10000}] (gOS phyloSignal)."""
    if not layout or not layout["leaves"]:
        return []
    row_of = {cid: i for i, cid in enumerate(layout["leaves"])}
    rows = [row_of[c] for c in cell_ids if c in row_of]
    n = len(rows)
    if n < 6:
        return []
    full = leaf_distances(layout)
    sub = full[np.ix_(rows, rows)].astype(np.float64)
    with np.errstate(divide="ignore"):
        Wm = np.where(sub > 0, 1.0 / sub, 0.0).astype(np.float32).astype(np.float64)
    np.fill_diagonal(Wm, 0.0)
    ids = [layout["leaves"][r] for r in rows]
    out = []
    for key, getter in variables:
        vals = np.array([_f(getter.get(i)) for i in ids], dtype=np.float64)
        bad = ~np.isfinite(vals)
        if bad.any():
            fin = vals[~bad]
            vals[bad] = fin.mean() if len(fin) else 0.0
        I = morans_i(vals, Wm, n)
        if not math.isfinite(I):
            out.append({"key": key, "I": None, "z": None, "p": None, "p_perm10000": None, "n": n})
            continue
        rand = Mulberry32(seed)
        perm = vals.copy()
        Is = np.zeros(n_perm)
        P = np.zeros((n_perm, n))
        for k in range(n_perm):
            for i in range(n - 1, 0, -1):
                j = int(math.floor(rand() * (i + 1)))
                perm[i], perm[j] = perm[j], perm[i]
            P[k] = perm
        D = P - P.sum(1, keepdims=True) / n
        den = (D * D).sum(1)
        nz = Wm != 0
        wsum = Wm[nz].sum()
        Is = (n / wsum) * (((D @ Wm) * D).sum(1) / den)
        mean = Is.mean()
        sd = math.sqrt(max(0.0, (Is ** 2).sum() / n_perm - mean * mean))
        rng = np.random.default_rng(seed)
        ge = 0
        done = 0
        while done < n_fine:
            b = min(1000, n_fine - done)
            Pb = np.array([rng.permutation(vals) for _ in range(b)])
            Db = Pb - Pb.sum(1, keepdims=True) / n
            Ib = (n / wsum) * (((Db @ Wm) * Db).sum(1) / (Db * Db).sum(1))
            ge += int((Ib >= I).sum())
            done += b
        out.append({"key": key, "I": I, "z": (I - mean) / sd if sd > 0 else None, "p": (int((Is >= I).sum()) + 1) / (n_perm + 1),
                    "p_perm10000": (ge + 1) / (n_fine + 1), "n": n})
    return out


def _f(x):
    try:
        v = float(x)
        return v
    except (TypeError, ValueError):
        return float("nan")


def mean_vaf(snv):
    out = {}
    for cell, entries in (snv.get("cells") or {}).items():
        v = [e[2] / (e[1] + e[2]) for e in entries if len(e) >= 4 and e[3] == 1 and e[1] is not None and e[2] is not None and (e[1] + e[2]) > 0]
        if v:
            out[cell] = float(np.mean(v))
    return out


def export_patient(data, outdir, patient, datafiles):
    tp = os.path.join(data, patient, "tree.nwk")
    if not os.path.exists(tp):
        print(patient, "no tree.nwk: not available", flush=True)
        return None
    with open(tp) as fh:
        layout = F.layout_tree(F.parse_newick(fh.read()))
    cells = F.patient_cells(datafiles, patient)
    rec = {c["cell_id"]: c for c in cells}
    wp = os.path.join(data, patient, "walks.json")
    walks = F.jload(wp) if os.path.exists(wp) else {}
    cell_ids = walks.get("cells") or [c["cell_id"] for c in cells]
    tumour = [i for i in cell_ids if i in rec and not F.is_normal_clone(rec[i]["clone_id"])]
    variables, meta = [], []
    for key, field in CONTROLS:
        variables.append((key, {i: rec[i].get(field) for i in tumour}))
        meta.append({"key": key, "kind": "control", "source": "gOS", "field": field})
    for gkey, ws in W.walk_groups(walks.get("walks") or [], cell_ids).items():
        for w in ws:
            key = w.get("label") or w.get("name") or str(w.get("id"))
            variables.append((key, {i: math.log1p(W._n((w.get("cells") or {}).get(i))) for i in tumour}))
            meta.append({"key": key, "kind": "walk", "source": "gOS", "gene_set": gkey, "walk_id": str(w.get("id")), "curated": w.get("curated")})
    sp = os.path.join(data, patient, "snv_matrix.json")
    if os.path.exists(sp):
        variables.append(("Mean somatic VAF", mean_vaf(F.jload(sp))))
        meta.append({"key": "Mean somatic VAF", "kind": "control", "source": "extra", "field": "snv_matrix.json alt/(ref+alt), genotype 1"})
    res = phylo_signal(layout, tumour, variables)
    if not res:
        print(patient, "fewer than 6 tumour cells on the tree: not available", flush=True)
        return None
    rows = []
    for r, m in zip(res, meta):
        rows.append(dict(m, I=r["I"], z=r["z"], p=r["p"], p_perm10000=r["p_perm10000"], n=r["n"]))
    out = {"format": "gos-sc-phylo-signal/1", "patient": patient, "n_perm": N_PERM, "seed": SEED, "n_perm_fine": N_PERM_FINE,
           "weights": "1 / patristic distance (tree.nwk)", "n_tumour_cells_on_tree": res[0]["n"], "rows": rows}
    pdir = os.path.join(outdir, patient)
    os.makedirs(pdir, exist_ok=True)
    with open(os.path.join(pdir, "phylo_signal.json"), "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    print(patient, len(rows), "variables", flush=True)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", default=os.path.join(F.DATA_ROOT, "data"))
    ap.add_argument("--datafiles", default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--patients", nargs="*", default=["BWH69", "BWH70", "BWH75", "MGH285", "MGH301", "MGH302", "MGH303", "MGH304", "MGH305"])
    a = ap.parse_args()
    datafiles = F.jload(a.datafiles or os.path.join(os.path.dirname(a.data.rstrip("/")), "datafiles.json"))
    for p in a.patients:
        export_patient(a.data, a.out, p, datafiles)


if __name__ == "__main__":
    main()
