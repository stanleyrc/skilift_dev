#!/usr/bin/env python3
"""Cohort TMB associations -> analysis/tmb_assoc/out/tmb_associations.json (copied to _cohort/tmb_associations.json).

Relates per-cell TMB (DNA, <P>/tmb.json) to the RNA, cell-cycle, copy-number and driver state of the same cell
(DNA-RNA linked malignant cells, <P>/rna/cells.json + covariates.json), and clade TMB to clade-level signature
exposure and cell-cycle composition.

Per-cell TMB measures (alt calls at tree sites / callable Mb, divided by the cell's truncal-site sensitivity = dropout
corrected):  total, subclonal (sites shared by an internal clade), private (cell-private sites; a lower bound).
Within a patient total TMB mostly encodes clade membership, so every test is also run adjusting for the DNA clone.

Per patient test = partial Spearman (rank-rank, residualised on a design matrix):
   raw        no adjustment
   qc         cell QC: sensitivity, callable fraction, ADO, MAPD, log median het depth
   qc_clone   qc + DNA clone indicators (within-clone association)
   qc_sdna    qc + DNA replication score (rt_z), so an RNA S-phase effect is not a DNA replication artefact
Across patients: Stouffer z weighted by sqrt(n) (Fisher-z of the partial r), sign concordance count.
Categorical covariates (Phase, RNA state): Kruskal-Wallis on QC-residualised ranks per patient, Fisher across patients.
Genome-wide: the same partial Spearman for every gene (>= 10% cells expressing) vs each TMB measure.
Clade level: per patient Spearman across clades (clade TMB vs signature fraction / mean S_Score / state fractions), combined.

usage: gos_sc_tmb_assoc.py [--data DIR] [--out DIR] [--patients P ...]
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
GENESETS = "/gpfs/commons/projects/imielinski_web/sclarke/gOS_dev_ase/genesets"
PATIENTS = ["BWH69", "BWH70", "BWH75", "MGH285", "MGH301", "MGH302", "MGH303", "MGH304", "MGH305"]
MEASURES = ["total", "subclonal", "private", "subclonal_pc", "private_pc"]
GENE_MEASURES = ["total", "subclonal", "private"]
MEASURE_LABEL = {"total": "Total TMB", "subclonal": "Subclonal TMB", "private": "Private TMB (lower bound)",
                 "subclonal_pc": "Subclonal TMB per copy (/ ploidy x 2)", "private_pc": "Private TMB per copy (/ ploidy x 2)"}
N_PERM = 2000
MIN_CLADE_CELLS = 5
ROBUST_MEASURE = "subclonal"
ROBUST_ADJ = "qc_clone_region"
ROBUST_FOCUS = ["Gene CN burden (mean |log2 CN/ploidy|)", "Module: Hypoxia (Hallmark)", "Module: Glycolysis (Hallmark)", "Module: MTORC1 signalling (Hallmark)"]
ADJ = ["raw", "qc", "qc_region", "qc_clone", "qc_clone_region", "qc_sdna"]
MODULES = {"MTORC1 signalling (Hallmark)": "HALLMARK_MTORC1_SIGNALING", "Glycolysis (Hallmark)": "HALLMARK_GLYCOLYSIS", "Hypoxia (Hallmark)": "HALLMARK_HYPOXIA",
           "MYC targets v1 (Hallmark)": "HALLMARK_MYC_TARGETS_V1", "Oxidative phosphorylation (Hallmark)": "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
           "Stress (3CA)": "3CA_5_STRESS", "Stress in vitro (3CA)": "3CA_7_STRESS_IN_VITRO"}
MIN_CELLS = 20
MIN_GROUP = 3
MIN_PCT = 0.10
TOP_GENES = 40
GENE_VARIANTS = [("qc_clone_region", "qc + clone + region"), ("qc_clone_region_cycle", "qc + clone + region + S/G2M score")]
SIGS_MERGE = {"SBS40a": "SBS40", "SBS40b": "SBS40", "SBS40c": "SBS40"}
NEFTEL = {
    "MES-like (Neftel)": "CHI3L1 ANXA2 ANXA1 CD44 VIM MT2A C1S NAMPT EFEMP1 C1R SOD2 IFITM3 TIMP1 SPP1 A2M S100A11 MT1X S100A10 FN1 LGALS1",
    "AC-like (Neftel)": "CST3 S100B SLC1A3 HEPN1 HOPX MT3 SPARCL1 MLC1 GFAP FABP7 BCAN PON2 METTL7B SPARC GATM RAMP1 PMP2 AQP4 DBI EDNRB",
    "OPC-like (Neftel)": "BCAN PLP1 GPR17 FIBIN LHFPL3 OLIG1 PSAT1 SCRG1 OMG APOD SIRT2 TNR THY1 PHYHIPL SOX2-OT NKAIN4 LPPR1 PTPRZ1 VCAN DBI",
    "NPC-like (Neftel)": "DLL3 DLL1 SOX4 TUBB3 HES6 TAGLN3 NEU4 MARCKSL1 CD24 STMN1 TCF12 BEX1 OLIG1 MAP2 FXYD6 PTPRS MLLT11 NPPA BCAN MEST",
}
SIG_LABEL = {"_clock": "SBS1+SBS5 (clock-like)", "_apobec": "SBS2+SBS13 (APOBEC)"}


def log(*a):
    print(*a, flush=True)


# ---------------------------------------------------------------------------------------------------------------
# statistics (numpy/scipy only)
# ---------------------------------------------------------------------------------------------------------------
def rank(a):
    return stats.rankdata(a, axis=-1).astype(np.float64)


def design(cols, n):
    """orthonormal basis of [intercept, z-scored columns]"""
    parts = [np.ones((n, 1))]
    for c in cols:
        c = np.asarray(c, dtype=np.float64).reshape(n, -1)
        sd = c.std(axis=0)
        keep = sd > 0
        c = (c - c.mean(axis=0)) / np.where(keep, sd, 1)
        parts.append(c[:, keep])
    D = np.hstack(parts)
    u, s, _ = np.linalg.svd(D, full_matrices=False)
    return u[:, s > 1e-8 * s[0]]


def partial_corr(y, X, Q):
    """partial Spearman of y (n,) with each row of X (g, n) given basis Q (n, k) incl. intercept. returns r, z, p, df"""
    n = len(y)
    X = np.atleast_2d(X)
    yr = rank(y)
    Xr = rank(X)
    yr = yr - Q @ (Q.T @ yr)
    Xr = Xr - (Xr @ Q) @ Q.T
    ny = np.linalg.norm(yr)
    nx = np.linalg.norm(Xr, axis=1)
    with np.errstate(all="ignore"):
        r = (Xr @ yr) / (nx * ny)
    r = np.where((nx > 1e-9) & (ny > 1e-9), r, np.nan)
    k = Q.shape[1]
    df = n - k - 1
    if df < 3:
        return r, np.full_like(r, np.nan), np.full_like(r, np.nan), df
    rr = np.clip(r, -0.999999, 0.999999)
    z = np.arctanh(rr) * math.sqrt(max(n - k - 2, 1))
    p = 2 * stats.norm.sf(np.abs(z))
    return r, z, p, df


def bh(ps):
    ps = np.asarray(ps, dtype=float)
    out = np.full(len(ps), np.nan)
    idx = np.where(np.isfinite(ps))[0]
    if not len(idx):
        return out
    order = idx[np.argsort(ps[idx])]
    m = len(order)
    prev = 1.0
    for rk in range(m, 0, -1):
        i = order[rk - 1]
        prev = min(prev, ps[i] * m / rk)
        out[i] = prev
    return out


def stouffer(units):
    """units: [(patient, z, n, r)] -> combined dict (weights sqrt(n)) or None"""
    units = [u for u in units if u[1] is not None and np.isfinite(u[1])]
    if not units:
        return None
    w = np.array([math.sqrt(u[2]) for u in units])
    z = np.array([u[1] for u in units])
    zc = float((w * z).sum() / math.sqrt((w * w).sum()))
    mean_r = float(np.average([u[3] for u in units], weights=w))
    sgn = np.sign(zc)
    return {"z": zc, "p": float(2 * stats.norm.sf(abs(zc))), "mean_r": mean_r, "n_patients": len(units),
            "n_concordant": int(sum(1 for u in units if np.sign(u[1]) == sgn and u[1] != 0)), "n_cells": int(sum(u[2] for u in units))}


def fisher(ps):
    ps = [max(p, 1e-300) for p in ps if p is not None and np.isfinite(p)]
    if not ps:
        return None
    return float(stats.chi2.sf(-2 * sum(math.log(p) for p in ps), 2 * len(ps)))


def clean(x):
    if isinstance(x, dict):
        return {k: clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, np.ndarray):
        return clean(x.tolist())
    if isinstance(x, (np.floating, float)):
        return float(f"{x:.5g}") if math.isfinite(x) else None
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


# ---------------------------------------------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------------------------------------------
class P:
    pass


def jload(path):
    with open(path) as f:
        return json.load(f)


def num(v):
    try:
        v = float(v)
        return v if math.isfinite(v) else np.nan
    except (TypeError, ValueError):
        return np.nan


def load_patient(root, p):
    d = os.path.join(root, p)
    if not os.path.exists(os.path.join(d, "tmb.json")):
        return None, "not available (no tmb.json: no per-cell SNV tree / TMB export for this patient)"
    if not os.path.exists(os.path.join(d, "rna", "covariates.json")):
        return None, "not available (no RNA covariate export)"
    tmb = jload(os.path.join(d, "tmb.json"))
    cells = jload(os.path.join(d, "rna", "cells.json"))["cells"]
    cov = jload(os.path.join(d, "rna", "covariates.json"))
    cvals = cov["cells"]
    tcell = {c["cell_id"]: c for c in tmb["cells"]}
    qc, sph = {}, {}
    for fn, store in (("qc.json", qc), ("sphase.json", sph)):
        fp = os.path.join(d, "precompute", fn)
        if os.path.exists(fp):
            for c in jload(fp)["cells"]:
                store[c["cell_id"]] = c
    keep = []
    for i, c in enumerate(cells):
        cid, rid = c.get("cell_id"), c.get("rna_id")
        if not cid or cid not in tcell or rid not in cvals:
            continue
        if c.get("Cell_Type") != "Malignant":
            continue
        if str(cvals[rid].get("DNA clone") or "").lower() == "normal":
            continue
        keep.append(i)
    if len(keep) < MIN_CELLS:
        return None, f"not available ({len(keep)} DNA-RNA linked malignant cells with TMB; need {MIN_CELLS})"
    o = P()
    o.name = p
    o.noisy = bool(tmb.get("noisy"))
    o.keep = np.array(keep)
    o.n = len(keep)
    o.cid = [cells[i]["cell_id"] for i in keep]
    o.rid = [cells[i]["rna_id"] for i in keep]
    rc = [cells[i] for i in keep]
    cv = [cvals[r] for r in o.rid]
    tc = [tcell[c] for c in o.cid]
    sens = np.array([num(t.get("sensitivity")) for t in tc])
    o.meas = {
        "total": np.array([num(t.get("tmb_corrected")) for t in tc]),
        "subclonal": np.array([num(t.get("tmb_subclonal")) for t in tc]) / sens,
        "private": np.array([num(t.get("tmb_private")) for t in tc]) / sens,
    }
    q = [qc.get(c, {}) for c in o.cid]
    with np.errstate(all="ignore"):
        o.qcm = {
            "sensitivity": sens,
            "callable_frac": np.array([num(t.get("callable_frac")) for t in tc]),
            "ado": np.array([num(x.get("ado")) for x in q]),
            "mapd": np.array([num(x.get("mapd")) for x in q]),
            "log_het_depth": np.log(np.array([num(x.get("het_depth_med")) for x in q])),
        }
    s = [sph.get(c, {}) for c in o.cid]
    o.dna_s = np.array([num(x.get("rt_z")) for x in s])
    o.dna_s_prob = np.array([num(x.get("s_prob")) for x in s])
    o.clone = np.array([str(c.get("DNA clone") or "NA") for c in cv])
    o.region = np.array([str(c.get("Region") or "NA") for c in rc])
    o.state = np.array([str(c.get("state") or "NA") for c in rc])
    o.phase = np.array([str(c.get("Phase") or "NA") for c in rc])
    ncount = np.array([num(c.get("nCount_RNA")) for c in rc])
    o.num = {
        "S_Score": np.array([num(c.get("S_Score")) for c in rc]),
        "G2M_Score": np.array([num(c.get("G2M_Score")) for c in rc]),
        "DNA replication score (rt_z)": o.dna_s,
        "DNA S-phase probability": o.dna_s_prob,
        "Cycling (RNA Phase S/G2M)": np.array([(1.0 if c.get("Phase") != "G1" else 0.0) if c.get("Phase") else np.nan for c in rc]),
        "Ploidy": np.array([num(c.get("Ploidy")) for c in cv]),
        "log10 nCount_RNA (technical)": np.log10(ncount + 1),
        "percent_mt (technical)": np.array([num(c.get("percent_mt")) for c in rc]),
    }
    with np.errstate(all="ignore"):
        dna = o.num["Ploidy"] / 2.0
        o.meas["subclonal_pc"] = o.meas["subclonal"] / dna
        o.meas["private_pc"] = o.meas["private"] / dna
    # proxy for a MES-high area: fraction of MES cells among the RNA cells of the same sampling region
    o.num["Region MES fraction (proxy for MES-high area)"] = np.array([np.mean(o.state[o.region == r] == "MES") for r in o.region])
    # gene copy-number burden: mean |log2(CN / ploidy)| over the panel genes of the covariate export
    cn_names = [f["name"] for f in cov["fields"] if f["name"].startswith("CN: ")]
    if cn_names:
        pl = np.where(np.isfinite(o.num["Ploidy"]), o.num["Ploidy"], 2.0)
        B = np.array([[num(c.get(nm)) for nm in cn_names] for c in cv])
        with np.errstate(all="ignore"):
            L = np.abs(np.log2(np.clip(B, 0.1, None) / pl[:, None]))
            o.num["Gene CN burden (mean |log2 CN/ploidy|)"] = np.nanmean(L, axis=1)
    o.binary, o.bin_known = {}, {}
    ecn = [f["name"] for f in cov["fields"] if f["name"].startswith("ecDNA copies: ")]
    if ecn:
        E = np.array([[num(c.get(nm)) for nm in ecn] for c in cv])
        tot = np.nansum(E, axis=1)
        o.num["ecDNA copies (log1p, all amplicons)"] = np.log1p(tot)
        o.binary["ecDNA carrier (any)"] = tot > 0
        o.bin_known["ecDNA carrier (any)"] = np.ones(o.n, dtype=bool)
    for f in cov["fields"]:
        nm = f["name"]
        if nm.startswith("Driver: "):
            v = [c.get(nm) for c in cv]
            k = nm[len("Driver: "):]
            o.binary[k] = np.array([x == "carrier" for x in v])
            o.bin_known[k] = np.array([x in ("carrier", "non-carrier") for x in v])
    # expression (log-normalised, as the gOS compare tab)
    genes = [g.strip() for g in open(os.path.join(d, "rna", "genes.tsv"))]
    indptr = np.fromfile(os.path.join(d, "rna", "matrix.indptr.i32"), dtype=np.int32)
    indices = np.fromfile(os.path.join(d, "rna", "matrix.indices.i32"), dtype=np.int32)
    data = np.fromfile(os.path.join(d, "rna", "matrix.data.f32"), dtype=np.float32)
    mat = sparse.csr_matrix((data, indices, indptr), shape=(len(genes), len(cells))).tocsc()[:, o.keep].tocsr()
    seen = np.asarray((mat > 0).sum(axis=1)).ravel() >= max(10, MIN_PCT * o.n)
    o.genes = np.array(genes)[seen]
    o.X = mat[np.where(seen)[0]].toarray().astype(np.float32)
    gi = {g: i for i, g in enumerate(o.genes)}
    for label, setname in MODULES.items():
        ii = [gi[g] for g in MODSETS.get(setname, []) if g in gi]
        if len(ii) >= 10:
            Z = o.X[ii].astype(np.float64)
            sd = Z.std(axis=1)
            Z = (Z - Z.mean(axis=1, keepdims=True))[sd > 0] / sd[sd > 0, None]
            o.num["Module: " + label] = Z.mean(axis=0)
        uset = [g for g in MODSETS.get(setname, []) if g in gi and g not in SHARED]
        if len(uset) >= 10:
            jj = [gi[g] for g in uset]
            Z = o.X[jj].astype(np.float64)
            sd = Z.std(axis=1)
            Z = (Z - Z.mean(axis=1, keepdims=True))[sd > 0] / sd[sd > 0, None]
            o.num["Module (unique genes): " + label] = Z.mean(axis=0)
    o.tmb_clades = tmb.get("clades", [])
    o.summary = tmb.get("patient_summary", {})
    try:
        o.clade_cells = {c["id"]: set(c.get("cells") or []) for c in jload(os.path.join(d, "clades.json"))["clades"]}
    except FileNotFoundError:
        o.clade_cells = {}
    try:
        o.sigjoint = jload(os.path.join(d, "signatures_joint.json"))
    except FileNotFoundError:
        o.sigjoint = None
    return o, None


# ---------------------------------------------------------------------------------------------------------------
# per-cell tests
# ---------------------------------------------------------------------------------------------------------------
def qc_columns(o, mask):
    return [k for k, v in o.qcm.items() if np.isfinite(v[mask]).mean() >= 0.95 and np.nanstd(v[mask]) > 0]


def dummies(labels, mask):
    lab = labels[mask]
    u = list(np.unique(lab))
    if len(u) < 2:
        return None
    return np.column_stack([(lab == x).astype(float) for x in u[1:]])


def basis(o, adj, mask, extra=None):
    """(Q, ok) design basis for the cells in `mask` under adjustment `adj`; ok refines mask to cells with all inputs"""
    ok = mask.copy()
    qcn = qc_columns(o, mask) if adj != "raw" else []
    for k in qcn:
        ok &= np.isfinite(o.qcm[k])
    if adj == "qc_sdna":
        ok &= np.isfinite(o.dna_s)
    for e in extra or []:
        ok &= np.isfinite(e)
    cols = [o.qcm[k][ok] for k in qcn]
    if adj == "qc_sdna":
        cols.append(o.dna_s[ok])
    if "clone" in adj:
        cd = dummies(o.clone, ok)
        if cd is not None:
            cols.append(cd)
    if "region" in adj:
        cd = dummies(o.region, ok)
        if cd is not None:
            cols.append(cd)
    for e in extra or []:
        cols.append(e[ok])
    return design(cols, int(ok.sum())), ok


def test_numeric(o, y, x, adj):
    m = np.isfinite(y) & np.isfinite(x)
    if m.sum() < MIN_CELLS or np.nanstd(x[m]) == 0:
        return None
    Q, ok = basis(o, adj, m)
    n = int(ok.sum())
    if n < MIN_CELLS:
        return None
    r, z, p, df = partial_corr(y[ok], x[ok][None, :], Q)
    if not np.isfinite(r[0]):
        return None
    return {"r": float(r[0]), "z": float(z[0]), "p": float(p[0]), "n": n}


def test_categorical(o, y, lab, groups):
    m = np.isfinite(y) & np.isin(lab, groups)
    if m.sum() < MIN_CELLS:
        return None
    Q, ok = basis(o, "qc", m)
    yr = rank(y[ok])
    res = yr - Q @ (Q.T @ yr)
    L = lab[ok]
    use = [g for g in groups if (L == g).sum() >= MIN_GROUP]
    if len(use) < 2:
        return None
    H, p = stats.kruskal(*[res[L == g] for g in use])
    return {"H": float(H), "p": float(p), "n": int(ok.sum()),
            "groups": {g: {"n": int((L == g).sum()), "median_tmb": float(np.median(y[ok][L == g]))} for g in use}}


def pack_row(family, cname, meas, adj, per, pts):
    units = [(pn, t["z"], t["n"], t["r"]) for pn, t in per.items()]
    return {"family": family, "covariate": cname, "measure": meas, "adjust": adj, "per_patient": per,
            "combined": stouffer(units), "combined_no_noisy": stouffer([u for u in units if not pts[u[0]].noisy])}


def per_cell_numeric(pts):
    rows = []
    names = []
    for o in pts.values():
        names += [k for k in o.num if k not in names]
    for cname in names:
        for meas in MEASURES:
            for adj in ADJ:
                per = {}
                for pn, o in pts.items():
                    if cname in o.num:
                        t = test_numeric(o, o.meas[meas], o.num[cname], adj)
                        if t:
                            per[pn] = t
                if per:
                    rows.append(pack_row("numeric", cname, meas, adj, per, pts))
    bnames = []
    for o in pts.values():
        bnames += [k for k in o.binary if k not in bnames]
    for cname in bnames:
        for meas in MEASURES:
            for adj in ADJ:
                per = {}
                for pn, o in pts.items():
                    if cname not in o.binary:
                        continue
                    kn, car = o.bin_known[cname], o.binary[cname]
                    if (car & kn).sum() < MIN_GROUP or (~car & kn).sum() < MIN_GROUP:
                        continue
                    t = test_numeric(o, o.meas[meas], np.where(kn, car.astype(float), np.nan), adj)
                    if t:
                        t["n_carrier"], t["n_non"] = int((car & kn).sum()), int((~car & kn).sum())
                        per[pn] = t
                if per:
                    rows.append(pack_row("binary", cname, meas, adj, per, pts))
    # BH q across all covariates tested in each measure x adjustment family
    for key in ("combined", "combined_no_noisy"):
        fam = defaultdict(list)
        for r in rows:
            if r[key]:
                fam[(r["measure"], r["adjust"])].append(r)
        for lst in fam.values():
            for r, qq in zip(lst, bh([r[key]["p"] for r in lst])):
                r[key]["q"] = float(qq)
    return rows


def per_cell_categorical(pts):
    rows = []
    for fam, attr, groups in (("Phase", "phase", ["G1", "S", "G2M"]), ("RNA state", "state", ["AC", "MES", "NPC", "OPC"])):
        for meas in MEASURES:
            per = {}
            for pn, o in pts.items():
                t = test_categorical(o, o.meas[meas], getattr(o, attr), groups)
                if t:
                    per[pn] = t
            ps = [t["p"] for t in per.values()]
            rows.append({"family": "categorical", "covariate": fam, "measure": meas, "adjust": "qc", "per_patient": per,
                         "fisher_p": fisher(ps), "fisher_p_no_noisy": fisher([t["p"] for pn, t in per.items() if not pts[pn].noisy]),
                         "n_patients": len(ps), "n_sig_patients": int(sum(1 for x in ps if x < 0.05))})
    return rows


def sphase_confounding(pts):
    """does RNA S-phase / DNA replication track DNA callability? Spearman with each cell QC metric"""
    out = []
    for xname in ("S_Score", "G2M_Score", "DNA replication score (rt_z)"):
        for qn in ("sensitivity", "callable_frac", "ado", "mapd", "log_het_depth"):
            per = {}
            for pn, o in pts.items():
                x, y = o.num[xname], o.qcm[qn]
                m = np.isfinite(x) & np.isfinite(y)
                if m.sum() < MIN_CELLS:
                    continue
                r, pv = stats.spearmanr(x[m], y[m])
                z = math.atanh(max(min(r, 0.999999), -0.999999)) * math.sqrt(max(m.sum() - 3, 1))
                per[pn] = {"r": float(r), "p": float(pv), "n": int(m.sum()), "z": z}
            units = [(pn, t["z"], t["n"], t["r"]) for pn, t in per.items()]
            if per:
                out.append({"x": xname, "qc": qn, "per_patient": per, "combined": stouffer(units),
                            "combined_no_noisy": stouffer([u for u in units if not pts[u[0]].noisy])})
    return out


# ---------------------------------------------------------------------------------------------------------------
# genome-wide
# ---------------------------------------------------------------------------------------------------------------
MODSETS = {}
SHARED = set()
SHARED_BY = {}


def load_modsets():
    for fn in ("hallmark.gmt", "gbm_3ca.gmt"):
        path = os.path.join(GENESETS, fn)
        if os.path.exists(path):
            for line in open(path):
                f = line.rstrip("\n").split("\t")
                if len(f) > 2 and f[0] in MODULES.values():
                    MODSETS[f[0]] = sorted(set(x for x in f[2:] if x))
    cnt = defaultdict(list)
    for k, v in MODSETS.items():
        for g in v:
            cnt[g].append(k)
    for g, ks in cnt.items():
        if len(ks) >= 2:
            SHARED.add(g)
    for label, setname in MODULES.items():
        SHARED_BY[label] = sorted(g for g in MODSETS.get(setname, []) if g in SHARED)


def load_sets():
    sets = []
    for fn, coll in (("hallmark.gmt", "Hallmark"), ("gbm_3ca.gmt", "3CA metaprogram")):
        path = os.path.join(GENESETS, fn)
        if os.path.exists(path):
            for line in open(path):
                f = line.rstrip("\n").split("\t")
                if len(f) > 2 and f[0]:
                    sets.append((coll, f[0], sorted(set(x for x in f[2:] if x))))
    for k, v in NEFTEL.items():
        sets.append(("Neftel", k, v.split()))
    return sets


_SETS = None


def enrich(rows):
    global _SETS
    if len(rows) < 500:
        return {"up": [], "down": []}
    if _SETS is None:
        _SETS = load_sets()
    genes = np.array([r["gene"] for r in rows])
    z = np.array([r["z"] for r in rows])
    index = {g: i for i, g in enumerate(genes)}
    res = []
    for coll, name, gl in _SETS:
        ii = np.array([index[g] for g in gl if g in index])
        if len(ii) < 5 or len(ii) > 500:
            continue
        mask = np.zeros(len(genes), dtype=bool)
        mask[ii] = True
        pv = stats.mannwhitneyu(z[mask], z[~mask], alternative="two-sided").pvalue
        lead = ii[np.argsort(-np.abs(z[ii]))][:6]
        res.append({"collection": coll, "term": name, "n": int(len(ii)), "mean_z": float(z[mask].mean()), "p": float(pv),
                    "leading": [str(genes[i]) for i in lead]})
    if not res:
        return {"up": [], "down": []}
    q = bh([r["p"] for r in res])
    for r, qq in zip(res, q):
        r["q"] = float(qq)
    return {"up": sorted([r for r in res if r["mean_z"] > 0], key=lambda r: r["p"])[:8],
            "down": sorted([r for r in res if r["mean_z"] < 0], key=lambda r: r["p"])[:8]}


def genewide(pts):
    out = {}
    for meas in GENE_MEASURES:
        for adj, adjlab in GENE_VARIANTS:
            by = defaultdict(list)
            n_by = {}
            for pn, o in pts.items():
                y = o.meas[meas]
                extra = [o.num["S_Score"], o.num["G2M_Score"]] if adj.endswith("cycle") else None
                if extra is not None and min(np.isfinite(e).sum() for e in extra) < MIN_CELLS:
                    continue
                Q, ok = basis(o, "qc_clone_region", np.isfinite(y), extra)
                if ok.sum() < MIN_CELLS:
                    continue
                r, z, p, df = partial_corr(y[ok], o.X[:, ok].astype(np.float64), Q)
                n_by[pn] = int(ok.sum())
                for gi in np.where(np.isfinite(z))[0]:
                    by[o.genes[gi]].append((pn, float(z[gi]), int(ok.sum()), float(r[gi])))
            rows = []
            for g, lst in by.items():
                if len(lst) < 2:
                    continue
                c = stouffer(lst)
                nn = stouffer([u for u in lst if not pts[u[0]].noisy])
                rows.append({"gene": g, "z": c["z"], "p": c["p"], "mean_r": c["mean_r"], "n_patients": c["n_patients"],
                             "n_concordant": c["n_concordant"], "z_no_noisy": nn["z"] if nn else None,
                             "per": {u[0]: round(u[3], 3) for u in lst}})
            if rows:
                for r, qq in zip(rows, bh([r["p"] for r in rows])):
                    r["q"] = float(qq)
            rows.sort(key=lambda r: r["p"])
            sig = [r for r in rows if r["q"] < 0.05]
            out[f"{meas}|{adj}"] = {"measure": meas, "adjust": adj, "adjust_label": adjlab, "n_genes": len(rows), "n_q05": len(sig),
                                    "n_q05_concordant": int(sum(1 for r in sig if r["n_concordant"] == r["n_patients"])),
                                    "n_cells": n_by, "up": [r for r in rows if r["z"] > 0][:TOP_GENES],
                                    "down": [r for r in rows if r["z"] < 0][:TOP_GENES], "enrichment": enrich(rows)}
            log(f"genes {meas} {adj}: {len(rows)} tested, q<0.05: {len(sig)}")
    return out


# ---------------------------------------------------------------------------------------------------------------
# clade + patient level
# ---------------------------------------------------------------------------------------------------------------
def sname(s):
    return SIGS_MERGE.get(s, s)


def frac_dict(activities):
    f = defaultdict(float)
    for a in activities:
        f[sname(a["signature"])] += float(a.get("fraction") or 0)
    f["_apobec"] = f.get("SBS2", 0.0) + f.get("SBS13", 0.0)
    f["_clock"] = f.get("SBS1", 0.0) + f.get("SBS5", 0.0)
    return f


def spearman_unit(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 6 or np.std(x[m]) == 0 or np.std(y[m]) == 0:
        return None
    r, p = stats.spearmanr(x[m], y[m])
    z = math.atanh(max(min(r, 0.999999), -0.999999)) * math.sqrt(max(m.sum() - 3, 1))
    return {"r": float(r), "p": float(p), "n": int(m.sum()), "z": z}


def clade_level(pts):
    clade_cov = {}
    sigfr = {}
    for pn, o in pts.items():
        cidx = {c: i for i, c in enumerate(o.cid)}
        rows = []
        for t in o.tmb_clades:
            if t.get("is_trunk"):
                continue
            ii = np.array([cidx[c] for c in o.clade_cells.get(t["id"], set()) if c in cidx], dtype=int)
            r = {"id": t["id"], "tmb": num(t.get("tmb_corrected")), "n_rna": int(len(ii))}
            if len(ii) >= 3:
                r["S_Score"] = float(np.nanmean(o.num["S_Score"][ii]))
                r["G2M_Score"] = float(np.nanmean(o.num["G2M_Score"][ii]))
                r["Cycling fraction (RNA Phase)"] = float(np.nanmean(o.num["Cycling (RNA Phase S/G2M)"][ii]))
                for s in ("AC", "MES", "NPC", "OPC"):
                    r["Fraction " + s] = float(np.mean(o.state[ii] == s))
                if np.isfinite(o.dna_s[ii]).any():
                    r["DNA replication score (rt_z)"] = float(np.nanmean(o.dna_s[ii]))
            rows.append(r)
        clade_cov[pn] = rows
        sigfr[pn] = {}
        if o.sigjoint:
            sset = {sname(x["signature"]) for x in o.sigjoint.get("signature_set", [])}
            for st in o.sigjoint.get("sets", []):
                if st.get("group") == "clade" and str(st.get("name", "")).startswith("clade: "):
                    sigfr[pn][st["name"][len("clade: "):]] = (frac_dict(st.get("activities", [])), sset, st.get("n"))
    sig_rows, comp_rows = [], []
    for sg in ["SBS1", "SBS5", "SBS8", "SBS40", "SBS37", "SBS34", "SBS46", "_clock", "_apobec"]:
        per = {}
        for pn in pts:
            xs, ys = [], []
            for r in clade_cov[pn]:
                e = sigfr[pn].get(r["id"])
                if e is None or not np.isfinite(r["tmb"]):
                    continue
                f, sset, _ = e
                if not sg.startswith("_") and sg not in sset:
                    continue
                xs.append(f.get(sg, 0.0))
                ys.append(r["tmb"])
            if len(xs) >= 6:
                t = spearman_unit(np.array(ys), np.array(xs))
                if t:
                    t["points"] = [[round(a, 4), round(b, 4)] for a, b in zip(ys, xs)]
                    per[pn] = t
        if per:
            units = [(pn, t["z"], t["n"], t["r"]) for pn, t in per.items()]
            sig_rows.append({"signature": SIG_LABEL.get(sg, sg), "per_patient": per, "combined": stouffer(units),
                             "combined_no_noisy": stouffer([u for u in units if not pts[u[0]].noisy])})
    for cn in ("S_Score", "G2M_Score", "Cycling fraction (RNA Phase)", "DNA replication score (rt_z)", "Fraction AC", "Fraction MES", "Fraction NPC", "Fraction OPC"):
        per = {}
        for pn in pts:
            rows = [r for r in clade_cov[pn] if cn in r and np.isfinite(r["tmb"]) and np.isfinite(r[cn])]
            if len(rows) >= 6:
                t = spearman_unit(np.array([r["tmb"] for r in rows]), np.array([r[cn] for r in rows]))
                if t:
                    t["points"] = [[round(r["tmb"], 4), round(r[cn], 4)] for r in rows]
                    per[pn] = t
        if per:
            units = [(pn, t["z"], t["n"], t["r"]) for pn, t in per.items()]
            comp_rows.append({"covariate": cn, "per_patient": per, "combined": stouffer(units),
                              "combined_no_noisy": stouffer([u for u in units if not pts[u[0]].noisy])})
    pat = []
    for pn, o in pts.items():
        al = {st["name"]: st for st in (o.sigjoint or {}).get("sets", [])}.get("all")
        if not al:
            continue
        s = o.summary
        pat.append({"patient": pn, "noisy": o.noisy, "cell_tmb_median": (s.get("cell_tmb_corrected") or {}).get("median"),
                    "subclonal": (s.get("subclonal") or {}).get("tmb_corrected"), "private": (s.get("private") or {}).get("tmb_corrected"),
                    "fractions": dict(frac_dict(al["activities"]))})
    tests = []
    for sg in ["SBS1", "SBS5", "SBS8", "SBS40", "SBS37", "_clock", "_apobec"]:
        for mk in ("cell_tmb_median", "subclonal"):
            for excl in (False, True):
                rows = [p for p in pat if not (excl and p["noisy"]) and p[mk] is not None]
                x = np.array([p[mk] for p in rows], dtype=float)
                y = np.array([p["fractions"].get(sg, 0.0) for p in rows])
                if len(rows) >= 5 and np.std(y) > 0:
                    r, pv = stats.spearmanr(x, y)
                    tests.append({"signature": SIG_LABEL.get(sg, sg), "tmb": mk, "excluding_noisy": excl, "rho": float(r), "p": float(pv), "n": len(rows)})
    return {"signatures": sig_rows, "composition": comp_rows, "patient_level": pat, "patient_level_tests": tests,
            "n_clades": {pn: len(clade_cov[pn]) for pn in pts}}


def highest_burden_clade(pts, key="Module: MTORC1 signalling (Hallmark)"):
    """highest-TMB non-trunk clade vs its sister cells (parent clade cells not in it), module score residualised on region (rank), rank-sum"""
    out = {}
    for pn, o in pts.items():
        if key not in o.num:
            continue
        cidx = {c: i for i, c in enumerate(o.cid)}
        par = {}
        try:
            par = {c["id"]: c.get("parent") for c in jload(os.path.join(ROOT, pn, "clades.json"))["clades"]}
        except FileNotFoundError:
            pass
        best = None
        for t in sorted(o.tmb_clades, key=lambda t: -(num(t.get("tmb_corrected")) if np.isfinite(num(t.get("tmb_corrected"))) else -1)):
            if t.get("is_trunk"):
                continue
            a = np.zeros(o.n, dtype=bool)
            for c in o.clade_cells.get(t["id"], set()):
                if c in cidx:
                    a[cidx[c]] = True
            pc = o.clade_cells.get(par.get(t["id"]), set())
            b = np.zeros(o.n, dtype=bool)
            for c in pc:
                if c in cidx:
                    b[cidx[c]] = True
            b &= ~a
            if a.sum() >= 5 and b.sum() >= 5:
                best = (t, a, b)
                break
        if not best:
            out[pn] = {"status": "not available (no clade with >= 5 RNA cells and >= 5 sister cells)"}
            continue
        t, a, b = best
        m = a | b
        Q = design([dummies(o.region, m)] if dummies(o.region, m) is not None else [], int(m.sum()))
        x = o.num[key][m]
        xr = rank(x)
        res = xr - Q @ (Q.T @ xr)
        g = a[m]
        U, p = stats.mannwhitneyu(res[g], res[~g], alternative="two-sided")
        out[pn] = {"clade": t["id"], "tmb": num(t.get("tmb_corrected")), "n_clade": int(a.sum()), "n_sister": int(b.sum()),
                   "mean_diff": float(o.num[key][a].mean() - o.num[key][b].mean()), "p_region_adjusted": float(p),
                   "direction": "up" if o.num[key][a].mean() > o.num[key][b].mean() else "down"}
    return out


# ---------------------------------------------------------------------------------------------------------------
# robustness: within-clade permutation, leaf-clade test, verdict
# ---------------------------------------------------------------------------------------------------------------
def deepest_clade_groups(o):
    """label per linked cell: the deepest called clade (smallest set) with >= MIN_CLADE_CELLS linked cells containing it, else 'trunk'"""
    cidx = {c: i for i, c in enumerate(o.cid)}
    sets = []
    for cid, cells in o.clade_cells.items():
        ii = frozenset(cidx[c] for c in cells if c in cidx)
        if MIN_CLADE_CELLS <= len(ii) < o.n:
            sets.append((len(ii), cid, ii))
    sets.sort()
    lab = np.array(["trunk"] * o.n, dtype=object)
    for _, cid, ii in reversed(sets):  # large first so that smaller (deeper) clades overwrite
        lab[list(ii)] = cid
    return lab


def perm_matrix(groups, idx, B, rng):
    """B x len(idx) index matrix: positions of idx shuffled within the clade group"""
    g = groups[idx]
    out = np.tile(np.arange(len(idx)), (B, 1))
    for lab in np.unique(g):
        pos = np.where(g == lab)[0]
        if len(pos) < 2:
            continue
        for b in range(B):
            out[b, pos] = pos[rng.permutation(len(pos))]
    return out


def permutation_tests(pts, covs, measures, adj="qc_region", B=N_PERM, seed=1):
    """shuffle the covariate among cells of the same deepest clade; statistic = Stouffer z of the partial Spearman (adjustment adj)"""
    rng = np.random.default_rng(seed)
    groups = {pn: deepest_clade_groups(o) for pn, o in pts.items()}
    out = {}
    cache = {}
    for cname in covs:
        for meas in measures:
            obs, perm, per = [], [], {}
            for pn, o in pts.items():
                if cname not in o.num:
                    continue
                y, x = o.meas[meas], o.num[cname]
                m = np.isfinite(y) & np.isfinite(x)
                if m.sum() < MIN_CELLS or np.nanstd(x[m]) == 0:
                    continue
                Q, ok = basis(o, adj, m)
                n = int(ok.sum())
                if n < MIN_CELLS:
                    continue
                idx = np.where(ok)[0]
                r0, z0, p0, _ = partial_corr(y[ok], x[ok][None, :], Q)
                if not np.isfinite(z0[0]):
                    continue
                key = (pn, idx.tobytes())
                if key not in cache:
                    cache[key] = perm_matrix(groups[pn], idx, B, rng)
                pm = cache[key]
                Xp = x[ok][pm]
                r, z, p, _ = partial_corr(y[ok], Xp, Q)
                w = math.sqrt(n)
                obs.append((w, float(z0[0])))
                perm.append((w, np.nan_to_num(z)))
                per[pn] = {"z": float(z0[0]), "p_perm": float((1 + np.sum(np.abs(z) >= abs(z0[0]))) / (1 + B)), "n_groups": int(len(np.unique(groups[pn][idx])))}
            if not obs:
                continue
            W = math.sqrt(sum(w * w for w, _ in obs))
            zo = sum(w * z for w, z in obs) / W
            zp = sum(w * z for w, z in perm) / W
            out[f"{cname}|{meas}"] = {"covariate": cname, "measure": meas, "adjust": adj, "z_obs": float(zo), "p_perm": float((1 + np.sum(np.abs(zp) >= abs(zo))) / (1 + B)),
                                      "n_perm": B, "n_patients": len(obs), "per_patient": per}
        log("perm", cname)
    return out


def leaf_clades(o):
    """leaf-most called clades (>= MIN_CLADE_CELLS linked cells, not the whole sample, no eligible clade strictly inside): [(id, index array)]"""
    cidx = {c: i for i, c in enumerate(o.cid)}
    el = []
    for cid, cells in o.clade_cells.items():
        ii = frozenset(cidx[c] for c in cells if c in cidx)
        if MIN_CLADE_CELLS <= len(ii) < o.n:
            el.append((cid, ii))
    leaves = [(cid, ii) for cid, ii in el if not any(j < ii for _, j in el)]
    return [(cid, np.array(sorted(ii))) for cid, ii in leaves]


def leaf_clade_tests(pts):
    """one point per leaf clade: clade mean TMB measure vs clade mean covariate; per patient Spearman, Stouffer across patients"""
    names = []
    for o in pts.values():
        names += [k for k in o.num if k not in names]
    leaves = {pn: leaf_clades(o) for pn, o in pts.items()}
    rows = []
    for cname in names:
        for meas in MEASURES:
            per = {}
            for pn, o in pts.items():
                if cname not in o.num:
                    continue
                xs, ys = [], []
                for cid, ii in leaves[pn]:
                    with np.errstate(all="ignore"):
                        x = np.nanmean(o.num[cname][ii]) if np.isfinite(o.num[cname][ii]).any() else np.nan
                        y = np.nanmean(o.meas[meas][ii]) if np.isfinite(o.meas[meas][ii]).any() else np.nan
                    xs.append(x)
                    ys.append(y)
                t = spearman_unit(np.array(ys), np.array(xs))
                if t:
                    t["points"] = [[round(float(a), 4), round(float(b), 4)] for a, b in zip(ys, xs) if np.isfinite(a) and np.isfinite(b)]
                    per[pn] = t
            if per:
                units = [(pn, t["z"], t["n"], t["r"]) for pn, t in per.items()]
                rows.append({"covariate": cname, "measure": meas, "per_patient": per, "combined": stouffer(units),
                             "combined_no_noisy": stouffer([u for u in units if not pts[u[0]].noisy])})
    for key in ("combined", "combined_no_noisy"):
        fam = defaultdict(list)
        for r in rows:
            if r[key]:
                fam[r["measure"]].append(r)
        for lst in fam.values():
            for r, qq in zip(lst, bh([r[key]["p"] for r in lst])):
                r[key]["q"] = float(qq)
    return {"min_cells": MIN_CLADE_CELLS, "n_leaves": {pn: len(v) for pn, v in leaves.items()}, "rows": rows}


def find_numeric(res, cname, meas, adj):
    for r in res["numeric"]:
        if r["covariate"] == cname and r["measure"] == meas and r["adjust"] == adj:
            return r
    return None


def robustness(res):
    """numbers for the robustness table. Main readout = clade-level test (one point per leaf clade); the cell-level partial
    Spearman (QC + clone + region), the per-copy and overlap-free variants and the within-clade permutation are shown beside it.
    The verdict is derived in the UI helper (helpers/singleCell/tmbAssoc.js robustnessVerdict) from these numbers only."""
    meas, adj = ROBUST_MEASURE, ROBUST_ADJ
    leaf = {(r["covariate"], r["measure"]): r for r in res["clade_level_covariates"]["rows"]}
    main = [r for r in res["numeric"] if r["measure"] == meas and r["adjust"] == adj and r["combined"] and r["family"] == "numeric"]
    top_cell = [r["covariate"] for r in sorted(main, key=lambda r: r["combined"]["p"])[:5]]
    top_clade = [r["covariate"] for r in sorted([r for r in res["clade_level_covariates"]["rows"] if r["measure"] == meas and r["combined"]], key=lambda r: r["combined"]["p"])[:5]]
    covs = [c for c in dict.fromkeys(ROBUST_FOCUS + top_clade + top_cell) if not c.startswith("Module (unique")]

    def cell(c, m):
        r = find_numeric(res, c, m, adj)
        cc = r and r["combined"]
        return {"r": cc["mean_r"], "p": cc["p"], "q": cc.get("q"), "n_patients": cc["n_patients"], "n_concordant": cc["n_concordant"]} if cc else None

    def clade(c, m):
        r = leaf.get((c, m))
        cc = r and r["combined"]
        return {"r": cc["mean_r"], "p": cc["p"], "q": cc.get("q"), "n_patients": cc["n_patients"], "n_concordant": cc["n_concordant"], "n_points": cc["n_cells"]} if cc else None

    rows = []
    for c in covs:
        pm = res["permutation"].get(f"{c}|{meas}")
        u = ("Module (unique genes): " + c[len("Module: "):]) if c.startswith("Module: ") else None
        rows.append({"covariate": c, "measure": meas, "adjust": adj,
                     "clade": clade(c, meas), "cell": cell(c, meas),
                     "per_copy": {"clade": clade(c, meas + "_pc"), "cell": cell(c, meas + "_pc")},
                     "permutation": pm and {"p": pm["p_perm"], "z_obs": pm["z_obs"], "n_perm": pm["n_perm"], "adjust": pm["adjust"]},
                     "overlap_free_applicable": u is not None,
                     "overlap_free": {"clade": clade(u, meas), "cell": cell(u, meas)} if u else None})
    return rows


def cell_table(pts):
    out = {}
    for pn, o in pts.items():
        cols = {"cell": o.cid, "total": o.meas["total"], "subclonal": o.meas["subclonal"], "private": o.meas["private"],
                "state": list(o.state), "phase": list(o.phase), "clone": list(o.clone)}
        cols.update(o.num)
        for k, v in o.qcm.items():
            cols["QC: " + k] = v
        for k, v in o.binary.items():
            cols["bin: " + k] = np.where(o.bin_known[k], v.astype(float), np.nan)
        out[pn] = cols
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=ROOT)
    ap.add_argument("--out", default="/gpfs/commons/groups/imielinski_lab/projects/single_cell_GBM/analysis/tmb_assoc/out")
    ap.add_argument("--patients", nargs="*", default=PATIENTS)
    a = ap.parse_args()
    load_modsets()
    pts, na = {}, {}
    for p in a.patients:
        o, why = load_patient(a.data, p)
        if o is None:
            na[p] = why
            log(f"[{p}] {why}")
        else:
            pts[p] = o
            log(f"[{p}] {o.n} linked malignant cells, {len(o.genes)} genes, noisy={o.noisy}")
    res = {"format": "gos-sc-tmb-assoc/1", "created": datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
           "patients": {p: {"n_cells": o.n, "noisy": o.noisy, "clones": sorted(set(o.clone))} for p, o in pts.items()},
           "not_available": na, "measures": MEASURE_LABEL,
           "adjustments": {"raw": "none", "qc": "cell QC (sensitivity, callable fraction, ADO, MAPD, log het depth)",
                           "qc_region": "qc + sampling region", "qc_clone": "qc + DNA clone", "qc_clone_region": "qc + DNA clone + sampling region", "qc_sdna": "qc + DNA replication score (rt_z)"}}
    res["numeric"] = per_cell_numeric(pts)
    res["categorical"] = per_cell_categorical(pts)
    res["sphase_confounding"] = sphase_confounding(pts)
    res["clade_level"] = clade_level(pts)
    res["genes"] = genewide(pts)
    res["mtorc1_clade_check"] = highest_burden_clade(pts)
    res["module_shared_genes_removed"] = SHARED_BY
    res["spatial"] = {"region": "sampling region (Region) is an adjustment (qc_region, qc_clone_region) and 'Region MES fraction' is a covariate (proxy for a MES-high area)",
                      "slidetag": "not available (BWH75 has no tmb.json; no slide-tag positions for the TMB patients)"}
    pcovs = [k for k in dict.fromkeys(k for o in pts.values() for k in o.num)]
    res["permutation"] = permutation_tests(pts, pcovs, ["subclonal", "private", "total", "subclonal_pc", "private_pc"])
    res["clade_level_covariates"] = leaf_clade_tests(pts)
    res["robustness"] = robustness(res)
    res["cells"] = cell_table(pts)
    res["methods"] = ("Per-cell TMB = alt calls at CellPhy tree sites / callable Mb, divided by the cell's truncal-site sensitivity. Total, subclonal (sites "
                      "shared by an internal clade) and private (cell-private; lower bound) are reported separately. Per patient: partial Spearman (rank-rank, "
                      "residualised on the listed adjustment). Across patients: Stouffer z weighted by sqrt(n) with sign concordance; BH within a gene table. "
                      "Linked cells: DNA cell with an RNA cell annotated Malignant and not DNA clone Normal.")
    res["caveats"] = [
        "Total and subclonal TMB mostly encode clade membership: the qc_clone adjustment (within-clone) separates TMB from lineage.",
        "Private TMB is a lower bound (cell-private sites are not in the candidate set) and noisy; BWH69 SNV calls are noisy (shown with and without it).",
        "Replicating (S-phase) cells have altered coverage and callability: see the S-phase confounding table and the qc_sdna adjustment.",
        "Nested clades share cells; clade-level p values are descriptive.",
        "Patient-level signature associations use at most 6 patients and are descriptive only.",
        "Patients without a TMB export or without DNA-RNA linked cells are listed as not available; nothing is imputed."]
    os.makedirs(a.out, exist_ok=True)
    fp = os.path.join(a.out, "tmb_associations.json")
    with open(fp, "w") as f:
        json.dump(clean(res), f, separators=(",", ":"))
    log("wrote", fp, os.path.getsize(fp) // 1024, "KB")


if __name__ == "__main__":
    main()
