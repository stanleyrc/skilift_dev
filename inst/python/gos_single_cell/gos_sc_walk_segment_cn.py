#!/usr/bin/env python3
"""Per-cell copy number at the segments of each ecDNA walk + the segment-correlation triangles of paper Fig 5C -> <out>/<P>/walk_segment_cn.json.

A port of gOS triangleCorrelation.js / figures.js columnCorrelations / cohortFiguresPanel.js lociDomains + walkGroups, so the numbers are
what the browser computes live. For every region (each walk, and each driver-gene group of walks as in walkGroups with its defaults):
  windows   the walk's loci merged to one padded window per chromosome (pad = max(1.5 Mb, 25% of the span), clipped to the chromosome)
  positions K = 96 sample positions spread over the windows (>= 2 per window)
  position_cn   cells x positions: total CN of the cell's complex.json segment covering each position (null = none)
  node_cn       cells x walk nodes: length-weighted mean total CN over each walk node (null = none)
  carriers / others   tumour cells with / without the region's walk copies >= 1 (groups: summed copies); 'others' needs >= 8 cells
  corr_carriers, corr_others   positions x positions Pearson r (same window only; null elsewhere) of log1p(max(0, CN)) after
                               centring every cell on its own mean (removes ploidy / depth), as drawn in the triangle
Tumour cells = cells of walks.json `cells` (else datafiles.json) that are not Normal clones and have a complex.json.

usage: gos_sc_walk_segment_cn.py [--data DIR] [--out DIR] [--patients P ...]
"""
import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gos_sc_figlib as F  # noqa: E402

K = 96
MIN_GROUP = 8
PAD = 1.5e6
CORR_DIGITS = 3


def js_round(x):
    return math.floor(x + 0.5)


def walk_loci(walks):
    return [(n["chromosome"], float(n["start"]), float(n["end"])) for w in walks for n in (w.get("nodes") or [])]


def loci_domains(loci, pad=PAD):
    """gOS lociDomains: [(a, b, chr)] with a, b global; chromosome-local coordinates are kept as (la, lb)."""
    by = {}
    for chr_, s, e in loci:
        key = F.norm_chr(chr_)
        if key is None or not math.isfinite(s):
            continue
        cur = by.get(key)
        by[key] = [min(cur[0], s), max(cur[1], e)] if cur else [s, e]
    out = []
    for chr_, (s, e) in by.items():
        p = max(pad, 0.25 * (e - s))
        la = max(1, s - p)
        lb = min(F.LEN[chr_], e + p)
        a, b = F.OFFSET[chr_] + 1 + la, F.OFFSET[chr_] + 1 + lb
        if b > a:
            out.append((a, b, chr_, la, lb))
    out.sort(key=lambda d: d[0])
    return out


def positions_for(domains):
    total = sum(d[1] - d[0] for d in domains)
    pos = []
    for k, d in enumerate(domains):
        a, b = d[0], d[1]
        n = max(2, js_round(K * (b - a) / max(1, total)))
        for i in range(n):
            pos.append((a + ((i + 0.5) / n) * (b - a), k))
    return pos


def column_correlations(M):
    """gOS columnCorrelations on a rows x m matrix, NaN where a column has no variance."""
    rows, m = M.shape
    M = M.astype(np.float64)
    mean = M.sum(0) / max(1, rows)
    D = M - mean
    sd = np.sqrt((D ** 2).sum(0))
    C = np.full((m, m), np.nan)
    ok = sd > 0
    G = D.T @ D
    den = np.outer(sd, sd)
    sel = np.outer(ok, ok)
    C[sel] = G[sel] / den[sel]
    return C


def triangle_matrix(cn_rows, positions):
    """cells x positions of centred log1p(CN), as triangleCorrelation.js corr()"""
    g = np.array([p[0] for p in positions])
    M = np.zeros((len(cn_rows), len(positions)), dtype=np.float32)
    for i, row in enumerate(cn_rows):
        v = row.at(g)
        x = np.where(np.isfinite(v), np.log1p(np.maximum(0, v)), 0.0).astype(np.float32)
        M[i] = (x.astype(np.float64) - x.astype(np.float64).sum() / max(1, len(positions))).astype(np.float32)
    return M


def gene_set_key(w):
    dg = w.get("driver_genes")
    lst = [] if dg is None else (dg if isinstance(dg, list) else [dg])
    genes = sorted(set(str(x) for x in lst if str(x)))
    return "+".join(genes) if genes else "Other"


def walk_groups(walks, cell_ids, curated_only=True, include_other=False, min_cells=3, min_cn=1):
    """gOS walkGroups: gene-set groups of walks with summed per-cell copies."""
    any_curated = any(w.get("curated") is True for w in walks)
    keep = []
    for w in walks:
        if curated_only and any_curated and w.get("curated") is False:
            continue
        if not include_other and gene_set_key(w) == "Other":
            continue
        n = sum(1 for c in cell_ids if (_n((w.get("cells") or {}).get(c)) >= min_cn))
        if n >= min_cells:
            keep.append(w)
    by = {}
    for w in keep:
        by.setdefault(gene_set_key(w), []).append(w)
    return by


def _n(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else 0.0
    except (TypeError, ValueError):
        return 0.0


def region_record(label, kind, key, walks, tumour_ids, rows_by_id):
    domains = loci_domains(walk_loci(walks))
    if not domains:
        return None
    positions = positions_for(domains)
    cn = {c: sum(_n((w.get("cells") or {}).get(c)) for w in walks) for c in tumour_ids}
    carriers = [c for c in tumour_ids if cn[c] >= 1]
    others = [c for c in tumour_ids if cn[c] < 1]
    g = np.array([p[0] for p in positions])
    rec = {
        "region": key, "kind": kind, "label": label, "walk_ids": [str(w.get("id")) for w in walks],
        "windows": [{"chr": d[2], "start": int(round(d[3])), "end": int(round(d[4]))} for d in domains],
        "k": len(positions), "window_of_position": [p[1] for p in positions],
        "position_chr": [domains[p[1]][2] for p in positions],
        "position_bp": [int(round(p[0] - F.OFFSET[domains[p[1]][2]] - 1)) for p in positions],
        "cells": tumour_ids, "walk_copies": [cn[c] for c in tumour_ids], "carriers": carriers, "others": others,
        "min_group": MIN_GROUP,
    }
    pos_cn = []
    for c in tumour_ids:
        v = rows_by_id[c].at(g)
        pos_cn.append([None if not np.isfinite(x) else int(round(x)) for x in v])
    rec["position_cn"] = pos_cn
    nodes = [(w, n) for w in walks for n in (w.get("nodes") or [])]
    rec["nodes"] = [{"walk": str(w.get("id")), "chr": n["chromosome"], "start": n["start"], "end": n["end"]} for w, n in nodes]
    nc = []
    for c in tumour_ids:
        row = []
        for w, n in nodes:
            ch = F.norm_chr(n["chromosome"])
            if ch is None:
                row.append(None)
                continue
            a = F.OFFSET[ch] + 1 + min(float(n["start"]), float(n["end"]))
            b = F.OFFSET[ch] + 1 + max(float(n["start"]), float(n["end"]))
            m = rows_by_id[c].mean_over(a, b)
            row.append(None if not math.isfinite(m) else round(m, 3))
        nc.append(row)
    rec["node_cn"] = nc

    def corr(ids):
        if not ids:
            return None
        C = column_correlations(triangle_matrix([rows_by_id[c] for c in ids], positions))
        wnd = np.array([p[1] for p in positions])
        same = wnd[:, None] == wnd[None, :]
        return [[None if (not np.isfinite(C[i, j]) or not same[i, j]) else round(float(C[i, j]), CORR_DIGITS) for j in range(len(positions))] for i in range(len(positions))]

    rec["corr_carriers"] = corr(carriers)
    rec["corr_others"] = corr(others) if len(others) >= MIN_GROUP else None
    if rec["corr_others"] is None:
        rec["corr_others_note"] = "not available (fewer than %d other tumour cells)" % MIN_GROUP
    return rec


def export_patient(data, outdir, patient, datafiles):
    wp = os.path.join(data, patient, "walks.json")
    if not os.path.exists(wp):
        print(patient, "no walks.json: not available", flush=True)
        return None
    w = F.jload(wp)
    walks = w.get("walks") or []
    cells = F.patient_cells(datafiles, patient)
    rec_by = {c["cell_id"]: c for c in cells}
    cell_ids = w.get("cells") or [c["cell_id"] for c in cells]
    tumour, rows = [], {}
    for cid in cell_ids:
        if cid not in rec_by or F.is_normal_clone(rec_by[cid]["clone_id"]):
            continue
        cx = os.path.join(data, cid, "complex.json")
        if not os.path.exists(cx):
            continue
        rows[cid] = F.total_cn_row(F.jload(cx))
        tumour.append(cid)
    regions = []
    for wk in walks:
        r = region_record("%s (walk %s)" % (wk.get("label") or wk.get("name"), wk.get("id")), "walk", "w:%s" % wk.get("id"), [wk], tumour, rows)
        if r:
            r["curated"] = wk.get("curated")
            r["driver_genes"] = wk.get("driver_genes")
            regions.append(r)
    for key, ws in walk_groups(walks, cell_ids).items():
        r = region_record("ec%s walks" % key, "group", "g:%s" % key, ws, tumour, rows)
        if r:
            regions.append(r)
    out = {"format": "gos-sc-walk-segment-cn/1", "patient": patient, "genome": "hg38", "k": K, "n_tumour_cells": len(tumour),
           "source": "complex.json total CN per cell; port of gOS triangleCorrelation.js", "regions": regions}
    pdir = os.path.join(outdir, patient)
    os.makedirs(pdir, exist_ok=True)
    with open(os.path.join(pdir, "walk_segment_cn.json"), "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    print(patient, len(regions), "regions,", len(tumour), "tumour cells", flush=True)
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
