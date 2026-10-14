#!/usr/bin/env python3
"""Cross-patient associations of driver / copy-number events with cell state and RNA programs.

For every recurrent event (a gene + class such as "EGFR amp", "TP53 snv", "CDKN2A homdel", "DM" in at least
two patients, or any amp / homdel of a GBM gene) and every patient that has it, compare carriers with
non-carriers among that patient's tumour cells that have RNA (the same quantities the patient Drivers
drawer shows):

  - cell-state membership (MES / AC / NPC / OPC) and cycling (Phase != G1): log odds ratio, Wald SE, Fisher p
  - S and G2M scores, and the gene-set program scores of gOS (mean per-gene z-score over the patient's RNA
    cells, helpers/singleCell/geneSets.js): Hedges g, SE, Mann-Whitney p

Per-patient effects are combined by inverse-variance random effects (DerSimonian-Laird; fixed effect also
stored) and by Stouffer's z of the signed per-patient z (weights sqrt(n cells)). BH q-values over all
combined tests (events seen in >= 2 patients). Needs >= 3 carriers and >= 3 non-carriers with RNA
in a patient for it to count.

Events: data/<P>/filtered.events.json (carriers = cell_ids, upstream-only fusions dropped, tier <= 2 or
GBM gene or complex SV). Tumour cells: root clade of data/<P>/clades.json. RNA: data/<P>/rna/{cells.json,
genes.tsv,matrix.*}. Writes data/_cohort/event_associations.json.

usage: gos_sc_cohort_event_assoc.py [--data DIR] [--out DIR] [--patients P ...] [--min-cells N]
"""
import argparse
import datetime
import json
import math
import os
from collections import defaultdict

import numpy as np
from scipy import sparse, stats

ROOT = "/gpfs/commons/projects/imielinski_web/sclarke/gOS_GBM_sc/data"
PATIENTS = ["BWH69", "BWH70", "BWH75", "MGH285", "MGH301", "MGH302", "MGH303", "MGH304", "MGH305"]
GBM_GENES = set("EGFR PTEN TP53 CDKN2A CDKN2B CDK4 CDK6 MDM2 MDM4 RB1 NF1 PIK3CA PIK3R1 PDGFRA MET MYCN MYC TERT IDH1 IDH2 ATRX BRAF FGFR3 LZTR1 STAG2 SPTA1 GABRA6 QKI SOX2 OLIG2".split())
STATES = ["MES", "AC", "NPC", "OPC"]
# gene-set programs (copied from gOS src/helpers/singleCell/geneSets.js)
GENE_SETS = {
    "MES-like": "CHI3L1 ANXA2 ANXA1 CD44 VIM MT2A C1S NAMPT EFEMP1 C1R SOD2 IFITM3 TIMP1 SPP1 A2M S100A11 MT1X S100A10 FN1 LGALS1",
    "AC-like": "CST3 S100B SLC1A3 HEPN1 HOPX MT3 SPARCL1 MLC1 GFAP FABP7 BCAN PON2 METTL7B SPARC GATM RAMP1 PMP2 AQP4 DBI EDNRB",
    "OPC-like": "BCAN PLP1 GPR17 FIBIN LHFPL3 OLIG1 PSAT1 SCRG1 OMG APOD SIRT2 TNR THY1 PHYHIPL SOX2-OT NKAIN4 LPPR1 PTPRZ1 VCAN DBI",
    "NPC-like": "DLL3 DLL1 SOX4 TUBB3 HES6 TAGLN3 NEU4 MARCKSL1 CD24 STMN1 TCF12 BEX1 OLIG1 MAP2 FXYD6 PTPRS MLLT11 NPPA BCAN MEST",
    "G2/M genes": "MKI67 TOP2A CENPF CCNB1 CDK1 NUSAP1 UBE2C BIRC5 TPX2 AURKB CDC20 PTTG1 HMGB2 CCNB2 KIF2C",
    "S-phase genes": "PCNA MCM2 MCM3 MCM4 MCM5 MCM6 MCM7 RRM1 RRM2 TYMS FEN1 GINS2 HELLS UNG CDC6",
    "Hypoxia": "VEGFA CA9 ADM NDRG1 SLC2A1 PGK1 LDHA BNIP3 ANKRD37 P4HA1 ENO1 HILPDA EGLN3 PDK1 ALDOA",
    "Interferon": "ISG15 IFI6 IFIT1 IFIT3 IFI27 MX1 OAS1 STAT1 IRF7 BST2 IFITM1 XAF1 RSAD2 IFI44L ISG20",
}
GENE_SETS = {k: v.split() for k, v in GENE_SETS.items()}


def event_class(e):
    vt = str(e.get("vartype", "")).lower()
    ty = str(e.get("type", "")).lower()
    if ty == "complex sv":
        return "complex"
    if vt == "amp" or "amp" in ty:
        return "amp"
    if vt == "homdel":
        return "homdel"
    if "fusion" in vt or "fusion" in ty:
        return "fusion"
    if vt == "snv" or ty in ("missense", "trunc", "splice") or "frameshift" in ty or "nonsense" in ty:
        return "snv"
    return None


def event_id(e, cls):
    """(row id, label, gene) of an event; complex events are pooled by class (DM, TIC, ...)."""
    if cls == "complex":
        vt = str(e.get("vartype", "complex"))
        return f"{vt} (any)", f"{vt} (any locus)", ""
    if cls == "fusion":
        g = str(e.get("fusion_genes") or e.get("gene") or "")
        return f"{g} fusion", f"{g} fusion", g
    g = str(e.get("gene") or "")
    return f"{g} {cls}", f"{g} {cls}", g


def tier_ok(e, cls):
    gene = str(e.get("gene") or "")
    if cls == "complex":
        return True
    try:
        t = float(e.get("Tier"))
    except (TypeError, ValueError):
        t = 9
    return t <= 2 or gene in GBM_GENES


def load_patient(root, p, min_cells):
    d = os.path.join(root, p)
    try:
        clades = json.load(open(os.path.join(d, "clades.json")))["clades"]
        cells = json.load(open(os.path.join(d, "rna", "cells.json")))["cells"]
        events = json.load(open(os.path.join(d, "filtered.events.json")))
    except FileNotFoundError as err:
        print(f"[{p}] skipped: {err}")
        return None
    roots = [c for c in clades if c.get("parent") in (None, "")]
    tumour = set(max(roots, key=lambda c: c["n_cells"])["cells"]) if roots else set()
    rna = [(i, c) for i, c in enumerate(cells) if c.get("cell_id") and c["cell_id"] in tumour]
    if len(rna) < 2 * min_cells:
        print(f"[{p}] skipped: {len(rna)} tumour cells with RNA")
        return None
    idx = np.array([i for i, _ in rna])
    ids = [c["cell_id"] for _, c in rna]
    # program scores: mean per-gene z over the patient's tumour RNA cells
    genes = [g.strip() for g in open(os.path.join(d, "rna", "genes.tsv"))]
    gi = {g: k for k, g in enumerate(genes)}
    indptr = np.fromfile(os.path.join(d, "rna", "matrix.indptr.i32"), dtype=np.int32)
    indices = np.fromfile(os.path.join(d, "rna", "matrix.indices.i32"), dtype=np.int32)
    data = np.fromfile(os.path.join(d, "rna", "matrix.data.f32"), dtype=np.float32)
    mat = sparse.csr_matrix((data, indices, indptr), shape=(len(genes), len(cells)))
    scores = {}
    for name, gl in GENE_SETS.items():
        acc = np.zeros(len(idx))
        used = 0
        for g in gl:
            if g not in gi:
                continue
            v = np.asarray(mat[gi[g]].toarray()).ravel()[idx].astype(float)
            sd = v.std()
            if sd > 0:
                acc += (v - v.mean()) / sd
                used += 1
        if used >= 3:
            scores[name] = acc / used
    cont = dict(scores)
    for key, src in (("S score", "S_Score"), ("G2M score", "G2M_Score")):
        v = np.array([c.get(src, np.nan) if c.get(src) is not None else np.nan for _, c in rna], dtype=float)
        if np.isfinite(v).sum() >= 2 * min_cells:
            cont[key] = v
    binary = {f"{s} state": np.array([c.get("state") == s for _, c in rna]) for s in STATES}
    binary["cycling (not G1)"] = np.array([c.get("Phase") not in (None, "G1") for _, c in rna])
    # events: id -> carriers (tumour cells with RNA)
    pos = {c: k for k, c in enumerate(ids)}
    evs = {}
    for e in events:
        cls = event_class(e)
        if not cls or not tier_ok(e, cls):
            continue
        if str(e.get("upstream_only")).lower() == "true":
            continue
        eid, label, gene = event_id(e, cls)
        if not gene and cls != "complex":
            continue
        car = {pos[c] for c in str(e.get("cell_ids", "")).split(",") if c in pos}
        rec = evs.setdefault(eid, {"label": label, "gene": gene, "cls": cls, "carriers": set()})
        rec["carriers"] |= car
    return {"ids": ids, "cont": cont, "binary": binary, "events": evs}


def hedges(a, b):
    n1, n0 = len(a), len(b)
    sp = math.sqrt(((n1 - 1) * a.var(ddof=1) + (n0 - 1) * b.var(ddof=1)) / (n1 + n0 - 2))
    if not sp > 0:
        return None
    j = 1 - 3 / (4 * (n1 + n0) - 9)
    g = j * (a.mean() - b.mean()) / sp
    se = math.sqrt((n1 + n0) / (n1 * n0) + g * g / (2 * (n1 + n0)))
    return g, se


def one_patient(kind, carrier_mask, values):
    """effect, se, p, signed z for carriers vs non-carriers; None if undefined."""
    if kind == "binary":
        a = int(values[carrier_mask].sum()); b = int(carrier_mask.sum() - a)
        c = int(values[~carrier_mask].sum()); d = int((~carrier_mask).sum() - c)
        eff = math.log(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)))
        se = math.sqrt(1 / (a + .5) + 1 / (b + .5) + 1 / (c + .5) + 1 / (d + .5))
        p = stats.fisher_exact([[a, b], [c, d]])[1]
        return eff, se, p, [a, b, c, d]
    x, y = values[carrier_mask], values[~carrier_mask]
    ok1, ok0 = np.isfinite(x), np.isfinite(y)
    x, y = x[ok1], y[ok0]
    if len(x) < 2 or len(y) < 2:
        return None
    h = hedges(x, y)
    if h is None:
        return None
    p = stats.mannwhitneyu(x, y, alternative="two-sided").pvalue
    return h[0], h[1], p, [len(x), len(y)]


def signed_z(eff, p):
    z = stats.norm.isf(min(max(p, 1e-300), 1) / 2)
    return math.copysign(z, eff) if eff != 0 else 0.0


def meta(rows):
    eff = np.array([r["effect"] for r in rows]); se = np.array([r["se"] for r in rows])
    w = 1 / se ** 2
    fe = (w * eff).sum() / w.sum(); fe_se = math.sqrt(1 / w.sum())
    q = float((w * (eff - fe) ** 2).sum()); k = len(rows)
    tau2 = max(0.0, (q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 1 else 0.0
    wr = 1 / (se ** 2 + tau2)
    re = (wr * eff).sum() / wr.sum(); re_se = math.sqrt(1 / wr.sum())
    zre = re / re_se
    ws = np.array([math.sqrt(r["n"]) for r in rows]); zs = np.array([r["z"] for r in rows])
    zst = float((ws * zs).sum() / math.sqrt((ws ** 2).sum()))
    return {
        "k": k, "effect": re, "lo": re - 1.96 * re_se, "hi": re + 1.96 * re_se, "p": 2 * stats.norm.sf(abs(zre)),
        "fixed_effect": fe, "fixed_se": fe_se, "tau2": tau2, "q": q, "i2": max(0.0, (q - (k - 1)) / q) if q > 0 and k > 1 else 0.0,
        "stouffer_z": zst, "stouffer_p": 2 * stats.norm.sf(abs(zst)), "k_up": int((eff > 0).sum()),
    }


def bh(ps):
    n = len(ps)
    order = np.argsort(ps)
    q = np.empty(n)
    prev = 1.0
    for rank, i in enumerate(order[::-1]):
        r = n - rank
        prev = min(prev, ps[i] * n / r)
        q[i] = prev
    return q


def clean(x):
    if isinstance(x, dict):
        return {k: clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, (np.floating, float)):
        return float(x) if math.isfinite(x) else None
    if isinstance(x, np.integer):
        return int(x)
    return x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=ROOT)
    ap.add_argument("--out", default=None)
    ap.add_argument("--patients", nargs="*", default=PATIENTS)
    ap.add_argument("--min-cells", type=int, default=3)
    a = ap.parse_args()
    out_dir = a.out or os.path.join(a.data, "_cohort")
    os.makedirs(out_dir, exist_ok=True)
    pts = {}
    for p in a.patients:
        r = load_patient(a.data, p, a.min_cells)
        if r:
            pts[p] = r
            print(f"[{p}] {len(r['ids'])} tumour RNA cells, {len(r['events'])} candidate events")
    # recurrent: >= 2 patients with the event (>= min carriers with RNA), or amp / homdel of a GBM gene
    seen = defaultdict(list)
    for p, r in pts.items():
        for eid, ev in r["events"].items():
            if len(ev["carriers"]) >= a.min_cells and len(r["ids"]) - len(ev["carriers"]) >= a.min_cells:
                seen[eid].append(p)
    chosen = {}
    for eid, ps in seen.items():
        ev = next(pts[p]["events"][eid] for p in ps)
        if len(ps) >= 2 or (ev["cls"] in ("amp", "homdel") and ev["gene"] in GBM_GENES):
            chosen[eid] = ps
    results = []
    for eid, ps in sorted(chosen.items()):
        ev0 = next(pts[p]["events"][eid] for p in ps)
        outcomes = defaultdict(list)
        for p in ps:
            r = pts[p]
            mask = np.zeros(len(r["ids"]), dtype=bool)
            mask[list(r["events"][eid]["carriers"])] = True
            for name, v in list(r["binary"].items()):
                if v.sum() < 1:
                    continue
                res = one_patient("binary", mask, v)
                if res:
                    outcomes[("binary", name)].append((p, mask, res))
            for name, v in r["cont"].items():
                res = one_patient("cont", mask, np.asarray(v, dtype=float))
                if res:
                    outcomes[("cont", name)].append((p, mask, res))
        for (kind, name), lst in outcomes.items():
            per = [{"patient": p, "effect": res[0], "se": res[1], "p": res[2], "z": signed_z(res[0], res[2]),
                    "n": int(len(mask)), "n_carriers": int(mask.sum()), "table": res[3]} for p, mask, res in lst]
            row = {"event": eid, "label": ev0["label"], "gene": ev0["gene"], "cls": ev0["cls"], "outcome": name,
                   "measure": "log odds ratio" if kind == "binary" else "Hedges g", "n_patients_event": len(ps), "per_patient": per}
            row.update(meta(per) if len(per) >= 2 else {"k": len(per)})
            results.append(row)
    multi = [r for r in results if r["k"] >= 2]
    if multi:
        for key in ("p", "stouffer_p"):
            qs = bh(np.array([r[key] for r in multi]))
            for r, q in zip(multi, qs):
                r["q" if key == "p" else "stouffer_q"] = float(q)
    events = [{"event": eid, "label": next(pts[p]["events"][eid] for p in ps)["label"], "cls": next(pts[p]["events"][eid] for p in ps)["cls"],
               "gene": next(pts[p]["events"][eid] for p in ps)["gene"], "patients": ps,
               "carriers": {p: len(pts[p]["events"][eid]["carriers"]) for p in ps}} for eid, ps in sorted(chosen.items())]
    doc = {
        "format": "gos-sc-event-associations/1",
        "created": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "patients": {p: len(r["ids"]) for p, r in pts.items()},
        "params": {"min_cells": a.min_cells, "tests": "Fisher exact (binary), Mann-Whitney (continuous)", "meta": "DerSimonian-Laird random effects; Stouffer z weighted by sqrt(n)"},
        "events": events,
        "results": results,
        "caveats": [
            "Carriers vs non-carriers are compared within each patient, so patient-level confounding is removed, but carriers of one event are often also carriers of the others on the same clade.",
            "Cells with RNA and a DNA match only; small groups give wide intervals. Events in one patient only are listed without a combined estimate.",
        ],
    }
    path = os.path.join(out_dir, "event_associations.json")
    if os.path.exists(path):
        os.replace(path, path + ".bak_261014") if not os.path.exists(path + ".bak_261014") else None
    with open(path, "w") as f:
        json.dump(clean(doc), f, allow_nan=False)
    print(f"wrote {path}: {len(events)} events, {len(results)} tests ({len(multi)} combined)")


if __name__ == "__main__":
    main()
