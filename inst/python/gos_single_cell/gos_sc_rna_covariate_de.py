#!/usr/bin/env python3
"""Cross-patient RNA differential on the per-cell covariates (rna/covariates.json) -> _cohort/rna_covariate_de.json.

Three families of tests, each run WITHIN every patient that has enough tumour cells with RNA, then combined across patients:

  driver   carriers vs non-carriers of each recurrent driver (the same event in >= 2 patients) or of any amp / homdel of
           a GBM gene (single patient tests are listed but cannot be combined)
  cn       high vs low per-cell copy number of an oncogene (EGFR MDM2 CDK4 PDGFRA MET CDK6 MYCN), split at the patient's
           median ("CN: <gene>" covariate)
  clade    a subclonal clade vs its sister clade(s) (cells of the parent not in the clade): the four high-burden clades of
           the clade table (MGH303 n209, MGH304 n121, MGH305 n185, BWH70 n135), the next clades of each patient by relative
           burden, and a combined test of each patient's highest-burden eligible clade

Per patient and gene the statistic is the Wilcoxon rank-sum test on the log-normalised expression the gOS compare tab uses
(helpers/singleCell/rnaStats.js differentialExpression: log2 fold change of mean(expm1)+1, pct expressing, minimum pct 0.1
in either group), here computed STRATIFIED by plate x cycling (Phase != G1) so that plate / region and S-phase copy-number
inflation are not mistaken for a result: the rank test is run inside each stratum with both groups present and the
stratum z scores are combined (weights sqrt(nA*nB/(nA+nB))). If the strata would discard more than 30% of the smaller
group the stratification falls back to cycling only, then to none (recorded as `adjust` per patient).

Across patients a gene's signed z scores are combined by Stouffer (weights sqrt of the effective n) and by Fisher (as the
gOS cohort differential does), with BH q-values over the genes of a test and the number of patients whose fold change has
the same sign (`n_conc`). A gene is a ROBUST hit when q (Stouffer) < 0.05, |mean log2FC| >= 0.5, in >= 2 patients and
concordant in every patient tested. Enrichment per test: rank-sum of the combined z of the genes in a set against all other
tested genes, for MSigDB Hallmark, the Neftel MES / AC / OPC / NPC programs and the 3CA metaprograms.

usage: gos_sc_rna_covariate_de.py [--data DIR] [--out DIR] [--patients P ...] [--genesets DIR]
"""
import argparse
import datetime
import json
import math
import os
import re
import sys
from collections import defaultdict

import numpy as np
from scipy import sparse, stats

ROOT = "/gpfs/commons/projects/imielinski_web/sclarke/gOS_GBM_sc/data"
GENESETS = "/gpfs/commons/projects/imielinski_web/sclarke/gOS_dev_ase/genesets"
PATIENTS = ["BWH69", "BWH70", "BWH75", "MGH285", "MGH301", "MGH302", "MGH303", "MGH304", "MGH305"]
GBM_GENES = set("EGFR PTEN TP53 CDKN2A CDKN2B CDK4 CDK6 MDM2 MDM4 RB1 NF1 PIK3CA PIK3R1 PDGFRA MET MYCN MYC TERT IDH1 IDH2 ATRX BRAF FGFR3 LZTR1 STAG2 SPTA1 GABRA6 QKI SOX2 OLIG2".split())
CN_GENES = ["EGFR", "MDM2", "CDK4", "PDGFRA", "MET", "CDK6", "MYCN"]
NAMED_CLADES = ["MGH303_n209", "MGH304_n121", "MGH305_n185", "BWH70_n135"]
EXTRA_CLADES_PER_PATIENT = 2
MIN_GROUP = 3        # cells per group per patient
HARD_MIN = 15
DATA_ROOT = ROOT
MIN_CLADE = 5        # cells per side for a clade test
MIN_PCT = 0.1
TUMOUR_ONLY = True
QC_STRATA = True
LFC_ROBUST = 0.5
TOP_TABLE = 50
VOLCANO_N = 800
NEFTEL = {
    "MES-like (Neftel)": "CHI3L1 ANXA2 ANXA1 CD44 VIM MT2A C1S NAMPT EFEMP1 C1R SOD2 IFITM3 TIMP1 SPP1 A2M S100A11 MT1X S100A10 FN1 LGALS1",
    "AC-like (Neftel)": "CST3 S100B SLC1A3 HEPN1 HOPX MT3 SPARCL1 MLC1 GFAP FABP7 BCAN PON2 METTL7B SPARC GATM RAMP1 PMP2 AQP4 DBI EDNRB",
    "OPC-like (Neftel)": "BCAN PLP1 GPR17 FIBIN LHFPL3 OLIG1 PSAT1 SCRG1 OMG APOD SIRT2 TNR THY1 PHYHIPL SOX2-OT NKAIN4 LPPR1 PTPRZ1 VCAN DBI",
    "NPC-like (Neftel)": "DLL3 DLL1 SOX4 TUBB3 HES6 TAGLN3 NEU4 MARCKSL1 CD24 STMN1 TCF12 BEX1 OLIG1 MAP2 FXYD6 PTPRS MLLT11 NPPA BCAN MEST",
}


def log(*a):
    print(*a, flush=True)


# ----------------------------------------------------------------------------------------------------------------
# data
# ----------------------------------------------------------------------------------------------------------------
class Patient:
    pass


def load_patient(root, p):
    d = os.path.join(root, p)
    try:
        cov = json.load(open(os.path.join(d, "rna", "covariates.json")))
        cells = json.load(open(os.path.join(d, "rna", "cells.json")))["cells"]
    except FileNotFoundError as err:
        log(f"[{p}] skipped: {err}")
        return None
    try:
        clades = json.load(open(os.path.join(d, "clades.json")))["clades"]
    except FileNotFoundError:
        clades = []
    cvals = cov["cells"]
    keep = [i for i, c in enumerate(cells) if c.get("rna_id") in cvals]
    n_cov = len(keep)
    dropped = {}
    if TUMOUR_ONLY:
        # malignant tumour cells only, on both sides: Seurat Cell_Type == Malignant and a DNA clone that is not Normal
        def malignant(i):
            c = cells[i]
            clone = str(cvals[c["rna_id"]].get("DNA clone") or "")
            if c.get("Cell_Type") != "Malignant":
                return str(c.get("Cell_Type") or "unannotated")
            if clone.lower() == "normal":
                return "DNA clone Normal"
            return None
        why = {i: malignant(i) for i in keep}
        for r in why.values():
            if r:
                dropped[r] = dropped.get(r, 0) + 1
        keep = [i for i in keep if not why[i]]
    if len(keep) < 2 * MIN_GROUP:
        log(f"[{p}] skipped: {len(keep)} tumour cells with covariates")
        return None
    keep = np.array(keep)
    P = Patient()
    P.n_covariate_cells = n_cov
    P.keep = keep
    P.dropped = dropped
    lib = np.log10(np.array([float(cells[i].get("nCount_RNA") or 0) + 1 for i in keep]))
    mt = np.array([float(cells[i].get("percent_mt") if cells[i].get("percent_mt") is not None else np.nan) for i in keep])
    mt_hi = np.nan_to_num(mt, nan=0.0) > np.nanpercentile(mt, 75) if np.isfinite(mt).any() else np.zeros(len(keep), dtype=bool)
    P.qc = np.array([f"lib{int(l > np.median(lib))}mt{int(m)}" for l, m in zip(lib, mt_hi)])
    P.lib = lib
    P.mtv = mt
    P.name = p
    P.n = len(keep)
    P.rna_ids = [cells[i]["rna_id"] for i in keep]
    P.cell_ids = [cells[i].get("cell_id") for i in keep]
    P.cov = [cvals[r] for r in P.rna_ids]
    P.cycling = np.array([cells[i].get("Phase") not in (None, "G1") for i in keep])
    P.plate = np.array([str(c.get("Plate") or c.get("Region") or "NA") for c, i in zip(P.cov, keep)])
    P.field_meta = {f["name"]: f for f in cov["fields"]}
    # expression: genes x tumour cells, dense float32 on genes seen in >= 3 cells
    genes = [g.strip() for g in open(os.path.join(d, "rna", "genes.tsv"))]
    indptr = np.fromfile(os.path.join(d, "rna", "matrix.indptr.i32"), dtype=np.int32)
    indices = np.fromfile(os.path.join(d, "rna", "matrix.indices.i32"), dtype=np.int32)
    data = np.fromfile(os.path.join(d, "rna", "matrix.data.f32"), dtype=np.float32)
    mat = sparse.csr_matrix((data, indices, indptr), shape=(len(genes), len(cells))).tocsc()[:, keep].tocsr()
    seen = np.asarray((mat > 0).sum(axis=1)).ravel() >= 3
    P.genes = np.array(genes)[seen]
    P.X = mat[np.where(seen)[0]].toarray().astype(np.float32)
    P.E = np.expm1(P.X)
    P.nz = P.X > 0
    P.rank_cache = {}
    # clades -> masks over tumour cells (by DNA cell id)
    pos = {c: k for k, c in enumerate(P.cell_ids) if c}
    P.clades = {}
    for c in clades:
        m = np.zeros(P.n, dtype=bool)
        for cid in c.get("cells") or []:
            if cid in pos:
                m[pos[cid]] = True
        P.clades[c["id"]] = {"mask": m, "parent": c.get("parent"), "burden_rel": c.get("burden_rel"), "n_private_snv": c.get("n_private_snv"), "is_trunk": c.get("is_trunk"), "depth": c.get("depth")}
    log(f"[{p}] {P.n} tumour RNA cells, {len(P.genes)} genes, {len(P.clades)} clades")
    return P


def drivers_of(P):
    """event -> carrier mask (covariate 'Driver: X' = carrier)"""
    out = {}
    for name in P.field_meta:
        if not name.startswith("Driver: "):
            continue
        m = np.array([c.get(name) == "carrier" for c in P.cov])
        non = np.array([c.get(name) == "non-carrier" for c in P.cov])
        out[name[len("Driver: "):]] = (m, non)
    return out


def event_gene_class(ev):
    mt = re.match(r"^(\S+) (amp|homdel|snv|fusion)$", ev)
    return (mt.group(1), mt.group(2)) if mt else (None, "complex")


# ----------------------------------------------------------------------------------------------------------------
# stratified Wilcoxon
# ----------------------------------------------------------------------------------------------------------------
def stratum_z(P, idx, a_mask):
    """z (A - B) of the rank-sum per gene inside the cells idx; a_mask: bool over idx (True = group A)."""
    key = idx.tobytes()
    ranks = P.rank_cache.get(key)
    if ranks is None:
        sub = P.X[:, idx]
        ranks = stats.rankdata(sub, axis=1).astype(np.float32)
        zeros = (sub == 0).sum(axis=1).astype(np.float64)
        P.rank_cache[key] = (ranks, zeros)
        ranks = P.rank_cache[key]
    ranks, zeros = ranks
    nA = int(a_mask.sum())
    nB = len(idx) - nA
    n = nA + nB
    ra = ranks[:, a_mask].sum(axis=1, dtype=np.float64)
    u = ra - nA * (nA + 1) / 2.0
    mu = nA * nB / 2.0
    tie = zeros ** 3 - zeros
    var = nA * nB / 12.0 * ((n + 1) - tie / (n * (n - 1)))
    sd = np.sqrt(np.maximum(var, 1e-12))
    diff = u - mu
    z = np.where(diff == 0, 0.0, (diff - 0.5 * np.sign(diff)) / sd)
    z[var <= 1e-9] = 0.0
    return z


def strata_plan(P, a, b):
    """list of (idx, a_mask_in_idx, weight) and the adjustment used"""
    used = a | b
    chain = [("plate+cycling+qc", [P.plate, P.cycling, P.qc]), ("cycling+qc", [P.cycling, P.qc])] if QC_STRATA else []
    chain += [("plate+cycling", [P.plate, P.cycling]), ("cycling", [P.cycling]), ("none", [])]
    for label, keys in chain:
        if keys:
            lab = np.array(["|".join(str(k[i]) for k in keys) for i in range(P.n)])
        else:
            lab = np.zeros(P.n, dtype=int)
        plan = []
        for s in np.unique(lab[used]):
            sel = used & (lab == s)
            na, nb = int((sel & a).sum()), int((sel & b).sum())
            if na >= 2 and nb >= 2:
                idx = np.where(sel)[0]
                plan.append((idx, a[idx], math.sqrt(na * nb / (na + nb)), na, nb))
        ka = sum(x[3] for x in plan)
        kb = sum(x[4] for x in plan)
        if ka >= MIN_GROUP and kb >= MIN_GROUP and ka >= 0.7 * a.sum() and kb >= 0.7 * b.sum():
            return plan, label
        if label == "none":
            return (plan, label) if ka >= MIN_GROUP and kb >= MIN_GROUP else (None, label)
    return None, "none"


def de_patient(P, a, b):
    """per-gene arrays for one patient comparison (A vs B) or None"""
    plan, adjust = strata_plan(P, a, b)
    if not plan:
        return None
    zs = np.zeros(len(P.genes))
    wsum = 0.0
    w2 = 0.0
    for idx, am, w, _, _ in plan:
        zs += w * stratum_z(P, idx, am)
        wsum += w
        w2 += w * w
    z = zs / math.sqrt(w2)
    # restrict to cells actually used in the strata
    ua = np.zeros(P.n, dtype=bool)
    ub = np.zeros(P.n, dtype=bool)
    for idx, am, _, _, _ in plan:
        ua[idx[am]] = True
        ub[idx[~am]] = True
    nA, nB = int(ua.sum()), int(ub.sum())
    pa = P.nz[:, ua].mean(axis=1)
    pb = P.nz[:, ub].mean(axis=1)
    # log2FC as the gOS differential, computed per stratum and weighted like the z scores
    lfc = np.zeros(len(P.genes))
    for idx, am, w, _, _ in plan:
        ea = P.E[:, idx[am]].mean(axis=1)
        eb = P.E[:, idx[~am]].mean(axis=1)
        lfc += w * (np.log2(ea + 1) - np.log2(eb + 1))
    lfc /= wsum
    ok = np.maximum(pa, pb) >= MIN_PCT
    p = 2 * stats.norm.sf(np.abs(z))
    d_mt = float(np.nanmean(P.mtv[ua]) - np.nanmean(P.mtv[ub])) if np.isfinite(P.mtv).any() else None
    d_lib = float(P.lib[ua].mean() - P.lib[ub].mean())
    return {"ia": P.keep[ua].tolist(), "ib": P.keep[ub].tolist(), "d_mt": d_mt, "d_lib": d_lib, "z": z, "lfc": lfc, "p": p, "pa": pa, "pb": pb, "ok": ok, "nA": nA, "nB": nB, "adjust": adjust, "n_strata": len(plan), "w": math.sqrt(nA * nB / (nA + nB))}


# ----------------------------------------------------------------------------------------------------------------
# meta + enrichment
# ----------------------------------------------------------------------------------------------------------------
def bh(ps):
    ps = np.asarray(ps, dtype=float)
    n = len(ps)
    if not n:
        return ps
    order = np.argsort(ps)
    q = np.empty(n)
    prev = 1.0
    for rank, i in enumerate(order[::-1]):
        r = n - rank
        prev = min(prev, ps[i] * n / r)
        q[i] = prev
    return q


def combine(units):
    """units: [(patient label, de result, gene names)]; returns list of gene rows (unsorted) + universe"""
    by = defaultdict(list)
    for label, r, genes in units:
        ok = np.where(r["ok"])[0]
        for i in ok:
            by[genes[i]].append((label, r["z"][i], r["lfc"][i], r["p"][i], r["pa"][i], r["pb"][i], r["w"]))
    rows = []
    kmin = 2 if len(units) >= 2 else 1
    for g, lst in by.items():
        if len(lst) < kmin:
            continue
        w = np.array([x[6] for x in lst])
        z = np.array([x[1] for x in lst])
        zst = float((w * z).sum() / math.sqrt((w * w).sum()))
        lfc = float(np.mean([x[2] for x in lst]))
        ps = [max(x[3], 1e-300) for x in lst]
        fisher = float(stats.chi2.sf(-2 * sum(math.log(x) for x in ps), 2 * len(ps)))
        rows.append({
            "gene": g, "avg_log2FC": lfc, "pct_1": float(np.mean([x[4] for x in lst])), "pct_2": float(np.mean([x[5] for x in lst])),
            "z": zst, "p_val": float(2 * stats.norm.sf(abs(zst))), "fisher_p": fisher, "n_patients": len(lst),
            "n_concordant": int(sum(1 for x in lst if np.sign(x[1]) == np.sign(zst) and x[1] != 0)),
            "per": {x[0]: round(float(x[2]), 3) for x in lst},
        })
    if rows:
        q = bh([r["p_val"] for r in rows])
        qf = bh([r["fisher_p"] for r in rows])
        for r, a, b in zip(rows, q, qf):
            r["q_val"] = float(a)
            r["fisher_q"] = float(b)
    return rows


def is_robust(r):
    return r["n_patients"] >= 2 and r["n_concordant"] == r["n_patients"] and r["q_val"] < 0.05 and abs(r["avg_log2FC"]) >= LFC_ROBUST


def load_sets(gdir):
    sets = []
    def read(fn, collection):
        path = os.path.join(gdir, fn)
        if not os.path.exists(path):
            return
        for line in open(path):
            f = line.rstrip("\n").split("\t")
            if len(f) > 2 and f[0]:
                sets.append((collection, f[0], sorted(set(x for x in f[2:] if x))))
    read("hallmark.gmt", "Hallmark")
    read("gbm_3ca.gmt", "3CA metaprogram")
    for k, v in NEFTEL.items():
        sets.append(("Neftel", k, v.split()))
    return sets


def enrich(rows, sets):
    if len(rows) < 200:
        return {"up": [], "down": []}
    zmap = {r["gene"]: r["z"] for r in rows}
    genes = np.array(list(zmap))
    z = np.array([zmap[g] for g in genes])
    index = {g: i for i, g in enumerate(genes)}
    out = []
    for coll, name, gl in sets:
        ii = np.array([index[g] for g in gl if g in index])
        if len(ii) < 5 or len(ii) > 500:
            continue
        mask = np.zeros(len(genes), dtype=bool)
        mask[ii] = True
        inside, outside = z[mask], z[~mask]
        p = stats.mannwhitneyu(inside, outside, alternative="two-sided").pvalue
        lead = ii[np.argsort(-np.abs(z[ii]))][:6]
        out.append({"collection": coll, "term": name, "n": int(len(ii)), "mean_z": float(inside.mean()), "p": float(p), "leading": [str(genes[i]) for i in lead]})
    if not out:
        return {"up": [], "down": []}
    q = bh([o["p"] for o in out])
    for o, qq in zip(out, q):
        o["q"] = float(qq)
    up = sorted([o for o in out if o["mean_z"] > 0], key=lambda o: o["p"])[:8]
    down = sorted([o for o in out if o["mean_z"] < 0], key=lambda o: o["p"])[:8]
    return {"up": up, "down": down}


def clean(x):
    if isinstance(x, dict):
        return {k: clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, (np.floating, float)):
        return float(f"{x:.5g}") if math.isfinite(x) else None
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


# ----------------------------------------------------------------------------------------------------------------
# tests
# ----------------------------------------------------------------------------------------------------------------
def build_specs(pts):
    """test id -> {family, label, desc, per: {patient: (maskA, maskB, meta)}}"""
    specs = {}

    # (a) drivers
    ev_by_patient = {p: drivers_of(P) for p, P in pts.items()}
    events = defaultdict(dict)
    for p, d in ev_by_patient.items():
        for ev, (car, non) in d.items():
            if car.sum() >= MIN_GROUP and non.sum() >= MIN_GROUP:
                events[ev][p] = (car, non, {})
    for ev, per in sorted(events.items()):
        gene, cls = event_gene_class(ev)
        if len(per) >= 2 or (cls in ("amp", "homdel") and gene in GBM_GENES):
            specs[f"driver:{ev}"] = {"family": "driver", "label": f"{ev}: carriers vs non-carriers", "a": f"{ev} carriers", "b": "non-carriers",
                                     "gene": gene if cls in ("amp", "homdel", "snv") else None, "per": per, "units": "patients"}

    # (b) CN level
    for g in CN_GENES:
        per = {}
        for p, P in pts.items():
            name = f"CN: {g}"
            if name not in P.field_meta:
                continue
            v = np.array([c.get(name) if isinstance(c.get(name), (int, float)) else np.nan for c in P.cov], dtype=float)
            okc = np.isfinite(v)
            if okc.sum() < 2 * MIN_GROUP or len(np.unique(v[okc])) < 2:
                continue
            med = np.median(v[okc])
            hi, lo = okc & (v > med), okc & (v <= med)
            if hi.sum() == 0 or lo.sum() == 0:
                hi, lo = okc & (v >= med), okc & (v < med)
            if hi.sum() >= MIN_GROUP and lo.sum() >= MIN_GROUP:
                per[p] = (hi, lo, {"median_cn": float(med)})
        if per:
            specs[f"cn:{g}"] = {"family": "cn", "label": f"{g} copy number: high vs low (median split)", "a": f"{g} CN above median", "b": f"{g} CN at or below median",
                                "gene": g, "per": per, "units": "patients"}

    # (c) clades
    cohort_rows = {}
    try:
        for r in json.load(open(os.path.join(DATA_ROOT, "_cohort", "clades.json")))["rows"]:
            cohort_rows[(r["patient"], r["clade"])] = r
    except (FileNotFoundError, KeyError):
        pass
    top_unit = {}       # unrestricted: best candidate with >= 5 cells each side (first version, kept as a sensitivity row)
    top_hard = {}       # hardened: SNV-based burden > 0 and >= 15 RNA tumour cells each side
    info = {}
    for p, P in pts.items():
        cands = []
        allsub = []
        for cid, c in P.clades.items():
            if c["is_trunk"] or not c["parent"] or c["parent"] not in P.clades:
                continue
            m = c["mask"]
            sister = P.clades[c["parent"]]["mask"] & ~m
            allsub.append((cid, c.get("burden_rel") or 0, m, sister))
            if m.sum() >= MIN_CLADE and sister.sum() >= MIN_CLADE:
                cands.append((cid, c.get("burden_rel") or 0, m, sister))
        cands.sort(key=lambda x: -x[1])
        allsub.sort(key=lambda x: -x[1])
        chosen = [x for x in cands if x[0] in NAMED_CLADES]
        chosen += [x for x in cands if x[0] not in NAMED_CLADES][:EXTRA_CLADES_PER_PATIENT]
        if cands:
            top_unit[p] = cands[0]
        # hardened choice and the reason for each patient
        if not allsub or allsub[0][1] <= 0:
            info[p] = {"patient": p, "included": False, "reason": "no SNV-based clade burden (burden_rel is 0 or missing: no SNV matrix for this patient)"}
        else:
            ok = [x for x in allsub if x[1] > 0 and x[2].sum() >= HARD_MIN and x[3].sum() >= HARD_MIN]
            if ok:
                cid, br, m, sis = ok[0]
                top_hard[p] = ok[0]
                row = cohort_rows.get((p, cid), {})
                info[p] = {"patient": p, "included": True, "reason": "highest burden_rel among subclonal clades with >= %d RNA tumour cells on each side" % HARD_MIN, "clade": cid, "burden_rel": br,
                           "n_private_snv": row.get("n_private_snv"), "top_signature": row.get("top_signature"), "n_a": int(m.sum()), "n_b": int(sis.sum())}
            else:
                cid, br, m, sis = allsub[0]
                info[p] = {"patient": p, "included": False, "reason": "no subclonal clade has >= %d RNA tumour cells on each side (highest-burden clade %s: %d vs %d)" % (HARD_MIN, cid, int(m.sum()), int(sis.sum()))}
        for cid, br, m, s_ in chosen:
            specs[f"clade:{cid}"] = {"family": "clade", "label": f"{cid}: clade vs sister", "a": f"{cid}", "b": "sister clade(s)", "gene": None, "units": "single clade",
                                     "per": {p: (m, s_, {"burden_rel": br, "clade": cid, "n_private_snv": cohort_rows.get((p, cid), {}).get("n_private_snv"), "top_signature": cohort_rows.get((p, cid), {}).get("top_signature")})}, "named": cid in NAMED_CLADES}
    def pack(unit):
        return {p: (m, sis, {"clade": cid, "burden_rel": br, "n_private_snv": cohort_rows.get((p, cid), {}).get("n_private_snv"), "top_signature": cohort_rows.get((p, cid), {}).get("top_signature")})
                for p, (cid, br, m, sis) in unit.items()}
    if len(top_hard) >= 2:
        specs["clade:top"] = {"family": "clade", "label": "Highest-burden subclonal clade vs sister, combined over patients", "a": "highest-burden clade", "b": "sister clade(s)",
                              "gene": None, "per": pack(top_hard), "units": "patients", "clade_info": list(info.values())}
    if len(top_unit) >= 2:
        specs["clade:top_unrestricted"] = {"family": "clade", "label": "Highest-burden subclonal clade vs sister, unrestricted (sensitivity)", "a": "highest-burden clade", "b": "sister clade(s)",
                                           "gene": None, "per": pack(top_unit), "units": "patients", "sensitivity": True}
    return specs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=ROOT)
    ap.add_argument("--out", default=None)
    ap.add_argument("--patients", nargs="*", default=PATIENTS)
    ap.add_argument("--genesets", default=GENESETS)
    ap.add_argument("--all-cells", action="store_true", help="sensitivity: keep non-malignant cells (the first version)")
    ap.add_argument("--no-qc", action="store_true", help="sensitivity: no library-size / percent_mt strata")
    ap.add_argument("--max-tests", type=int, default=0, help="debug: stop after N tests")
    a = ap.parse_args()
    global TUMOUR_ONLY, QC_STRATA, DATA_ROOT
    DATA_ROOT = a.data
    TUMOUR_ONLY = not a.all_cells
    QC_STRATA = not a.no_qc
    out_dir = a.out or os.path.join(a.data, "_cohort")
    os.makedirs(out_dir, exist_ok=True)

    # pass 1: masks only (load one patient at a time would lose the cross-patient eligibility, so load all; ~100 MB each)
    pts = {}
    for p in a.patients:
        P = load_patient(a.data, p)
        if P:
            pts[p] = P
    specs = build_specs(pts)
    if a.max_tests:
        specs = dict(list(specs.items())[: a.max_tests])
    log(f"{len(specs)} tests: " + ", ".join(f"{f}={sum(1 for s in specs.values() if s['family'] == f)}" for f in ("driver", "cn", "clade")))

    # pass 2: per patient differential for every test
    res = defaultdict(dict)
    for p, P in pts.items():
        todo = [(tid, s["per"][p]) for tid, s in specs.items() if p in s["per"]]
        for k, (tid, (am, bm, meta)) in enumerate(todo):
            r = de_patient(P, am, bm)
            if r is not None:
                r["meta"] = meta
                res[tid][p] = r
        log(f"[{p}] {len(todo)} tests done")
        P.rank_cache.clear()

    sets = load_sets(a.genesets)
    log(f"{len(sets)} gene sets")
    tests = []
    cells_doc = {}
    for tid, s in specs.items():
        units = [(p, r, pts[p].genes) for p, r in res[tid].items()]
        if not units:
            continue
        rows = combine(units)
        if not rows:
            continue
        cells_doc[tid] = {p: {"a": r["ia"], "b": r["ib"]} for p, r in res[tid].items()}
        rows.sort(key=lambda r: (r["p_val"], -abs(r["avg_log2FC"])))
        robust = [r for r in rows if is_robust(r)]
        sig = [r for r in rows if r["q_val"] < 0.05]
        k = len(units)
        gene = s.get("gene")
        t = {
            "id": tid, "family": s["family"], "clade_info": s.get("clade_info"), "sensitivity": bool(s.get("sensitivity")), "label": s["label"], "group_a": s["a"], "group_b": s["b"], "gene": gene,
            "k_patients": k, "combinable": k >= 2, "units": s["units"],
            "patients": [{"patient": p, "n_a": r["nA"], "n_b": r["nB"], "adjust": r["adjust"], "n_strata": r["n_strata"], "d_mt": r.get("d_mt"), "d_lib": r.get("d_lib"), **{kk: vv for kk, vv in r["meta"].items()}} for p, r in res[tid].items()],
            "n_tested": len(rows), "n_sig_q05": len(sig), "n_robust": len(robust) if k >= 2 else 0,
            "n_up_robust": sum(1 for r in robust if r["avg_log2FC"] > 0) if k >= 2 else 0,
            "n_down_robust": sum(1 for r in robust if r["avg_log2FC"] < 0) if k >= 2 else 0,
            "top_up": [r["gene"] for r in rows if r["avg_log2FC"] > 0 and r["q_val"] < 0.05][:8],
            "top_down": [r["gene"] for r in rows if r["avg_log2FC"] < 0 and r["q_val"] < 0.05][:8],
            "robust_up": [r["gene"] for r in sorted(robust, key=lambda r: -r["z"]) if r["avg_log2FC"] > 0 and r["gene"] != gene][:8],
            "robust_down": [r["gene"] for r in sorted(robust, key=lambda r: r["z"]) if r["avg_log2FC"] < 0 and r["gene"] != gene][:8],
            "genes": rows[:TOP_TABLE],
            "volcano": {"genes": [r["gene"] for r in rows[:VOLCANO_N]], "lfc": [round(r["avg_log2FC"], 3) for r in rows[:VOLCANO_N]], "p": [float(f"{r['p_val']:.3g}") for r in rows[:VOLCANO_N]], "q": [float(f"{r['q_val']:.3g}") for r in rows[:VOLCANO_N]]},
            "enrichment": enrich(rows, sets),
        }
        if gene:
            hit = next((r for r in rows if r["gene"] == gene), None)
            if hit:
                t["self_gene"] = {"gene": gene, "avg_log2FC": hit["avg_log2FC"], "q_val": hit["q_val"], "n_patients": hit["n_patients"]}
        tests.append(t)
        log(f"{tid}: k={k} tested={len(rows)} sig={len(sig)} robust={t['n_robust']}")
    tests.sort(key=lambda t: (-t["n_robust"], -t["n_sig_q05"]))

    doc = {
        "format": "gos-sc-rna-covariate-de/1",
        "created": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "patients": {p: P.n for p, P in pts.items()},
        "tumour_filter": {p: {"cells_with_covariates": P.n_covariate_cells, "kept": P.n, "dropped": P.dropped} for p, P in pts.items()},
        "params": {"test": "stratified Wilcoxon rank-sum (strata = plate x cycling x library size / percent_mt), z combined by Stouffer; Fisher also given", "tumour_only": TUMOUR_ONLY, "qc_strata": QC_STRATA, "min_pct": MIN_PCT,
                   "min_group_cells": MIN_GROUP, "robust": f"q < 0.05, |log2FC| >= {LFC_ROBUST}, >= 2 patients, concordant in all", "cn_split": "median within patient",
                   "enrichment": "rank-sum of combined z, in-set vs other tested genes"},
        "tests": tests,
        "caveats": [
            "Cells are not independent samples of a patient; p values are cell-level and optimistic. Concordance across patients is the safer evidence.",
            "Driver carriers are often one clade or one plate; stratifying by plate x cycling removes only part of that, and the adjustment used is recorded per patient.",
            "Copy-number tests include the gene itself and its neighbours on the amplicon (cis dosage); they are a positive control, not a discovery.",
            "Clade tests have one patient each and cannot be combined; only the combined highest-burden test has more than one patient.",
            "Only malignant tumour cells with a DNA match and RNA are used (Cell_Type Malignant, DNA clone not Normal); non-malignant cells are dropped on both sides.",
            "Strata include library size (log nCount, above/below the patient median) and high percent_mt (top quartile) besides plate and cycling; the adjustment used is recorded per patient.",
        ],
    }
    path = os.path.join(out_dir, "rna_covariate_de.json")
    if os.path.exists(path) and not os.path.exists(path + ".bak_261014"):
        os.replace(path, path + ".bak_261014")
    with open(path, "w") as f:
        json.dump(clean(doc), f, allow_nan=False, separators=(",", ":"))
    cpath = os.path.join(out_dir, "rna_covariate_de_cells.json")
    if os.path.exists(cpath) and not os.path.exists(cpath + ".bak_261014"):
        os.replace(cpath, cpath + ".bak_261014")
    with open(cpath, "w") as f:
        json.dump({"format": "gos-sc-rna-covariate-de-cells/1", "note": "row indices into rna/cells.json of each patient: cells compared in each test (group a / group b)", "tests": cells_doc}, f, separators=(",", ":"))
    log(f"wrote {cpath}: {os.path.getsize(cpath) / 1e6:.1f} MB")
    log(f"wrote {path}: {len(tests)} tests, {os.path.getsize(path) / 1e6:.1f} MB")


if __name__ == "__main__":
    sys.exit(main())
