#!/usr/bin/env python3
"""Clade-level TMB and its associations -> _cohort/tmb_clade_assoc.json  (unit = tree clade / branch, not cell).

Inputs per patient (all under DATA/<P>/): tmb.json (clades + cells), clades.json (listed clades: >= 5 cells, >= 10 private
SNVs; cells, parent, branch_length, private_snvs, private_events, signatures), snv_matrix.json (variants: category, node,
alt_cells, map_confidence, cellphy_input, anchor), tree.nwk (tip branch lengths), rna/cells.json + rna/covariates.json
(per-cell state, region, phase, ploidy, CN, ecDNA copies; DNA cell ids).

Rows written (one per clade):
  * listed clades  (kind "clade"):  SNVs mapped to the branch = clades.json private_snvs, looked up in snv_matrix.json
  * single-cell branches (kind "singleton"): cell-private SNVs of snv_matrix (category private, anchor = the cell)
For each row the SNV count is exported for four variant sets so the front end can toggle them:
  all        every SNV mapped to the branch
  conf       map_confidence >= 2 (the next-best tree placement is >= 2 log-likelihood units worse)
  multi      SNV with alt reads in >= 2 cells (drops SNVs seen in a single cell)
  confmulti  both
Branch TMB = SNVs / callable Mb / discovery power (clades: tmb.json discovery_power_detect; single-cell branches: the
cell's truncal-site sensitivity, so private TMB is a LOWER BOUND).  Cumulative TMB (root-to-clade) and the per-branch-length
rate are computed in the front end from the same rows; the statistics here (default_stats) mirror
src/helpers/singleCell/cladeTmbAssoc.js and are the cross-check for it.

Statistics: per patient partial Spearman (rank-rank, residualised on rank(n_cells)), normal approximation z = atanh(r) *
sqrt(n - 3), permutation p-value when n < 20; across patients Stouffer z weighted by sqrt(n) with sign concordance.
Nested clades are not independent, so the default test set is TERMINAL clades only (no other listed clade inside).

usage: gos_sc_tmb_clade.py [--data DIR] [--out FILE] [--patients P ...]
"""
import argparse
import datetime
import json
import math
import os
import re
import shutil
from collections import Counter, defaultdict

import numpy as np
from scipy import stats

DATA = "/gpfs/commons/projects/imielinski_web/sclarke/gOS_GBM_sc/data"
PATIENTS = ["BWH69", "BWH70", "BWH75", "MGH285", "MGH301", "MGH302", "MGH303", "MGH304", "MGH305"]
CONF_MIN = 2.0
MIN_LINKED = 3
N_PERM = 4000
PERM_BELOW = 20
MIN_N = 5
GENES = ["EGFR", "MYCN", "CDK4", "MDM2"]
VARIANT_SETS = ["all", "conf", "multi", "confmulti"]
MEASURES = ["branch", "cumulative", "rate"]
CLOCK = ("SBS1", "SBS5")
APOBEC = ("SBS2", "SBS13")
SIG_MERGE = {"SBS40a": "SBS40", "SBS40b": "SBS40", "SBS40c": "SBS40"}

COV_META = {
    "n_cells": ("Clade size (cells)", "structure (technical)", "cells"),
    "depth": ("Clade depth (listed ancestors)", "structure (technical)", "levels"),
    "branch_length": ("CellPhy branch length", "structure (technical)", "subst./site"),
    "drv_branch": ("Driver events gained on the branch", "drivers", "events"),
    "drv_cum": ("Driver events on the root-to-clade path", "drivers", "events"),
    "ploidy": ("Mean ploidy", "instability", "copies"),
    "cn_burden": ("Gene CN burden (mean |log2 CN/ploidy|, panel genes)", "instability", "log2"),
    "ecdna_frac": ("Fraction of cells with ecDNA (any)", "ecDNA / CN", "fraction"),
    "frac_AC": ("RNA state AC fraction", "RNA state", "fraction"),
    "frac_MES": ("RNA state MES fraction", "RNA state", "fraction"),
    "frac_NPC": ("RNA state NPC fraction", "RNA state", "fraction"),
    "frac_OPC": ("RNA state OPC fraction", "RNA state", "fraction"),
    "cycling": ("Cycling fraction (RNA Phase S/G2M)", "RNA state", "fraction"),
    "s_score": ("Mean S score", "RNA state", "score"),
    "region_entropy": ("Sampling-region diversity (Shannon entropy)", "region", "bits"),
    "region_top": ("Fraction of cells in the clade's main region", "region", "fraction"),
    "sig_clock": ("Branch signature: clock-like (SBS1+SBS5)", "signatures", "fraction of branch SNVs"),
    "sig_apobec": ("Branch signature: APOBEC (SBS2+SBS13)", "signatures", "fraction of branch SNVs"),
}
for g in GENES:
    COV_META["ecdna_" + g] = (g + " ecDNA copies (mean per cell)", "ecDNA / CN", "copies")
    COV_META["cn_" + g] = (g + " copy number (mean per cell)", "ecDNA / CN", "copies")


def log(*a):
    print(datetime.datetime.now().strftime("%H:%M:%S"), *a, flush=True)


def jload(p):
    with open(p) as f:
        return json.load(f)


def num(v):
    try:
        x = float(v)
        return x if math.isfinite(x) else float("nan")
    except (TypeError, ValueError):
        return float("nan")


def nanmean(a):
    a = np.asarray(a, dtype=float)
    a = a[np.isfinite(a)]
    return float(a.mean()) if len(a) else None


def rnd(x, k=5):
    return None if x is None or not math.isfinite(x) else round(float(x), k)


def tip_lengths(path):
    """tip name -> terminal branch length from a newick string"""
    txt = open(path).read()
    return {m.group(1): float(m.group(2)) for m in re.finditer(r"[(,]([^(),:;]+):([0-9.eE+-]+)", txt)}


def merge_sig(acts):
    out = defaultdict(float)
    for a in acts:
        out[SIG_MERGE.get(a["signature"], a["signature"])] += float(a["activity"])
    return out


def cell_table(d):
    """DNA cell id -> per-cell covariates (from the RNA-linked covariate export)"""
    fp = os.path.join(d, "rna", "covariates.json")
    if not os.path.exists(fp):
        return {}, [], "not available (no RNA covariate export)"
    cov = jload(fp)
    cells = jload(os.path.join(d, "rna", "cells.json"))["cells"]
    rna = {c["cell_id"]: c for c in cells if c.get("cell_id")}
    fields = [f["name"] for f in cov["fields"]]
    cn_names = [n for n in fields if n.startswith("CN: ")]
    ec = [n for n in fields if n.startswith("ecDNA copies: ")]
    out = {}
    for cid, r in rna.items():
        v = cov["cells"].get(r.get("rna_id"))
        if v is None or r.get("Cell_Type") != "Malignant" or str(v.get("DNA clone") or "").lower() == "normal":
            continue
        pl = num(v.get("Ploidy"))
        row = {"state": r.get("state"), "region": r.get("Region"), "phase": r.get("Phase"), "s": num(r.get("S_Score")), "ploidy": pl,
               "dm": v.get("Driver: DM (any)")}
        if cn_names and math.isfinite(pl):
            B = np.array([num(v.get(n)) for n in cn_names])
            with np.errstate(all="ignore"):
                row["cnb"] = float(np.nanmean(np.abs(np.log2(np.clip(B, 0.1, None) / max(pl, 0.5)))))
        for g in GENES:
            row["cn_" + g] = num(v.get("CN: " + g)) if ("CN: " + g) in v else float("nan")
            names = [n for n in ec if g in n[len("ecDNA copies: "):].rsplit(" ", 1)[0].split(",")]
            row["ec_" + g] = float(np.nansum([num(v.get(n)) for n in names])) if names else float("nan")
            row["ec_has_" + g] = bool(names)
        out[cid] = row
    return out, cn_names, None


def clade_covariates(cells, ct, min_linked):
    """mean/fraction covariates over the member cells with a covariate-export record"""
    rows = [ct[c] for c in cells if c in ct]
    cv = {"n_linked": len(rows)}
    if len(rows) < min_linked:
        return cv
    st = [r["state"] for r in rows if r["state"] in ("AC", "MES", "NPC", "OPC")]
    for s in ("AC", "MES", "NPC", "OPC"):
        cv["frac_" + s] = rnd(sum(1 for x in st if x == s) / len(st)) if st else None
    ph = [r["phase"] for r in rows if r["phase"]]
    cv["cycling"] = rnd(sum(1 for p in ph if p != "G1") / len(ph)) if ph else None
    cv["s_score"] = rnd(nanmean([r["s"] for r in rows]))
    cv["ploidy"] = rnd(nanmean([r["ploidy"] for r in rows]))
    cv["cn_burden"] = rnd(nanmean([r.get("cnb", float("nan")) for r in rows]))
    dm = [r["dm"] for r in rows if r["dm"] in ("carrier", "non-carrier")]
    cv["ecdna_frac"] = rnd(sum(1 for x in dm if x == "carrier") / len(dm)) if dm else None
    rg = Counter(r["region"] for r in rows if r["region"])
    if rg:
        p = np.array(list(rg.values()), dtype=float) / sum(rg.values())
        cv["region_entropy"] = rnd(float(-(p * np.log2(p)).sum()))
        cv["region_top"] = rnd(float(p.max()))
        cv["region_main"] = rg.most_common(1)[0][0]
    for g in GENES:
        if any(r["ec_has_" + g] for r in rows):
            cv["ecdna_" + g] = rnd(nanmean([r["ec_" + g] for r in rows]))
        v = nanmean([r["cn_" + g] for r in rows])
        cv["cn_" + g] = rnd(v) if v is not None else None
    return cv


def snv_counts(ids, V):
    """counts of the four variant sets for a list of variant ids (cellphy_input SNVs only)"""
    c = dict.fromkeys(VARIANT_SETS, 0)
    for i in ids:
        v = V.get(i)
        if v is None or not v.get("cellphy_input"):
            continue
        conf = num(v.get("map_confidence")) >= CONF_MIN
        multi = (v.get("alt_cells") or 0) >= 2
        c["all"] += 1
        c["conf"] += conf
        c["multi"] += multi
        c["confmulti"] += conf and multi
    return c


def build_patient(p, root):
    d = os.path.join(root, p)
    if not os.path.exists(os.path.join(d, "tmb.json")):
        return None, None, "not available (no tmb.json: no per-cell SNV tree / TMB export for this patient)"
    if not os.path.exists(os.path.join(d, "clades.json")) or not os.path.exists(os.path.join(d, "snv_matrix.json")):
        return None, None, "not available (no clades.json / snv_matrix.json)"
    tmb = jload(os.path.join(d, "tmb.json"))
    if not tmb.get("clades"):
        return None, None, "not available (tmb.json has no clades)"
    cl = jload(os.path.join(d, "clades.json"))["clades"]
    mat = jload(os.path.join(d, "snv_matrix.json"))
    V = {v["id"]: v for v in mat["variants"]}
    tclade = {c["id"]: c for c in tmb["clades"]}
    tcell = {c["cell_id"]: c for c in tmb["cells"]}
    ct, _, _ = cell_table(d)
    tips = tip_lengths(os.path.join(d, "tree.nwk")) if os.path.exists(os.path.join(d, "tree.nwk")) else {}
    sigfile = {}
    rows, recon = [], {}
    cellsets = {c["id"]: frozenset(c["cells"]) for c in cl}
    ids = [c["id"] for c in cl]
    # terminal = no other listed clade strictly inside
    terminal = {i: not any(cellsets[j] < cellsets[i] for j in ids if j != i) for i in ids}
    parent_of = {c["id"]: c.get("parent") for c in cl}
    depth = {}

    def dep(i):
        if i not in depth:
            depth[i] = 0 if not parent_of.get(i) or parent_of[i] not in cellsets else dep(parent_of[i]) + 1
        return depth[i]

    ev_branch = {}
    for c in cl:
        ev_branch[c["id"]] = sum(1 for e in (c.get("private_events") or []) if e.get("driver"))

    def ev_cum(i):
        n, seen = 0, set()
        while i and i in cellsets and i not in seen:
            seen.add(i)
            n += ev_branch[i]
            i = parent_of.get(i)
        return n

    # signature fractions per clade
    allsigs = Counter()
    for c in cl:
        if c.get("signatures"):
            for k in merge_sig(c["signatures"]):
                allsigs[k] += 1
    n_ok_checks = 0
    mapped_listed = set()
    for c in cl:
        t = tclade.get(c["id"])
        if t is None:
            continue
        cnt = snv_counts(c["private_snvs"], V)
        mapped_listed.update(i for i in c["private_snvs"] if V.get(i, {}).get("cellphy_input"))
        if t["n_private_snv"] != len(c["private_snvs"]) or t["n_private_snv"] != c["n_private_snv"]:
            recon.setdefault("count_mismatch", []).append(c["id"])
        else:
            n_ok_checks += 1
        row = {"patient": p, "id": c["id"], "kind": "clade", "node": c.get("node"), "parent": c.get("parent"), "clone": c.get("clone_id"),
               "n_cells": c["n_cells"], "terminal": terminal[c["id"]], "is_trunk": bool(c.get("is_trunk")),
               "callable_mb": rnd(t["callable_mb"], 4), "power": rnd(t["discovery_power_detect"], 4), "snv": cnt,
               "tmb_corrected_exported": t["tmb_corrected"], "burden_per_1000_sites": t.get("burden_per_1000_sites"),
               "drivers_branch_labels": [e.get("gene") or e.get("label") for e in (c.get("private_events") or []) if e.get("driver")][:6]}
        cv = {"n_cells": c["n_cells"], "depth": dep(c["id"]), "branch_length": rnd(c.get("branch_length"), 6),
              "drv_branch": ev_branch[c["id"]], "drv_cum": ev_cum(c["id"])}
        cv.update(clade_covariates(c["cells"], ct, MIN_LINKED))
        if c.get("signatures"):
            sg = merge_sig(c["signatures"])
            tot = sum(sg.values())
            if tot > 0:
                cv["sig_clock"] = rnd(sum(sg.get(k, 0) for k in CLOCK) / tot)
                cv["sig_apobec"] = rnd(sum(sg.get(k, 0) for k in APOBEC) / tot)
                for k, a in sg.items():
                    cv["sig:" + k] = rnd(a / tot)
        row["cov"] = cv
        rows.append(row)
    # single-cell branches from the private SNVs of the matrix, anchored on the cell
    priv = defaultdict(list)
    for v in mat["variants"]:
        if v.get("category") == "private" and v.get("cellphy_input") and v.get("anchor"):
            priv[v["anchor"]].append(v["id"])
    sizes = {i: len(cellsets[i]) for i in ids}
    term_cells = set().union(*[cellsets[i] for i in ids if terminal[i]]) if ids else set()
    mapped_priv = 0
    for cell, tc in tcell.items():
        if not tc.get("callable_mb") or not tc.get("sensitivity"):
            continue
        vid = [i for i in priv.get(cell, []) if i not in mapped_listed]
        cnt = snv_counts(vid, V)
        mapped_priv += cnt["all"]
        inside = [i for i in ids if cell in cellsets[i]]
        par = min(inside, key=lambda i: sizes[i]) if inside else None
        row = {"patient": p, "id": cell, "kind": "singleton", "node": None, "parent": par, "clone": tc.get("clone_id"), "n_cells": 1,
               "terminal": cell not in term_cells, "is_trunk": False, "callable_mb": rnd(tc["callable_mb"], 4),
               "power": rnd(tc["sensitivity"], 4), "snv": cnt, "drivers_branch_labels": []}
        cv = {"n_cells": 1, "depth": (dep(par) + 1) if par else 0, "branch_length": rnd(tips.get(cell), 6)}
        cv.update(clade_covariates([cell], ct, 1))
        row["cov"] = cv
        rows.append(row)
    # reconciliation against the matrix
    cat = Counter((v["category"], bool(v.get("cellphy_input"))) for v in mat["variants"])
    total_cp = sum(n for (k, ci), n in cat.items() if ci)
    recon.update({"n_listed_clades": len(cl), "n_singleton_branches": sum(1 for r in rows if r["kind"] == "singleton"),
                  "clade_count_checks_ok": n_ok_checks,
                  "matrix_cellphy_snvs": total_cp,
                  "snvs_on_listed_clades": len(mapped_listed),
                  "snvs_on_singleton_branches": mapped_priv,
                  "snvs_on_other_branches": total_cp - len(mapped_listed) - mapped_priv,
                  "note": "listed clades = tree clades with >= 5 cells and >= 10 private SNVs (clades.json); SNVs on smaller clades (2-4 cells or < 10 SNVs) and outside-tumour calls are not rows (snvs_on_other_branches)"})
    return rows, {"n_cells": len(tcell), "noisy": bool(tmb.get("noisy")), "n_linked_cells": len(ct), "reconciliation": recon}, None


# ---------------------------------------------------------------------------------------------------------------
# statistics (mirrored in src/helpers/singleCell/cladeTmbAssoc.js)
# ---------------------------------------------------------------------------------------------------------------
def rank(a):
    return stats.rankdata(a).astype(float)


def resid(y, Q):
    return y - Q @ (Q.T @ y)


def partial_spearman(x, y, size, rng):
    """partial Spearman of x,y given rank(size). returns dict r, z, p, p_perm, n or None"""
    x, y, size = (np.asarray(v, dtype=float) for v in (x, y, size))
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(size)
    n = int(ok.sum())
    if n < MIN_N or len(set(x[ok])) < 3 or len(set(y[ok])) < 3:
        return None
    rx, ry = rank(x[ok]), rank(y[ok])
    rs = rank(size[ok])
    Q = np.ones((n, 1))
    if np.std(rs) > 0:
        Q = np.linalg.qr(np.column_stack([np.ones(n), rs - rs.mean()]))[0]
    ex, ey = resid(rx, Q), resid(ry, Q)
    nx, ny = np.linalg.norm(ex), np.linalg.norm(ey)
    if nx < 1e-9 or ny < 1e-9:
        return None
    r = float(ex @ ey / (nx * ny))
    k = Q.shape[1]
    if n - k - 1 < 3:
        return None
    z = math.atanh(max(min(r, 0.999999), -0.999999)) * math.sqrt(max(n - k - 2, 1))
    out = {"r": r, "z": z, "p": float(2 * stats.norm.sf(abs(z))), "n": n}
    if n < PERM_BELOW:
        B = N_PERM
        perm = np.array([rng.permutation(ey) for _ in range(B)])
        rp = (perm @ ex) / (nx * ny)
        out["p_perm"] = float((1 + np.sum(np.abs(rp) >= abs(r) - 1e-12)) / (B + 1))
    return out


def stouffer(units):
    units = [u for u in units if u and math.isfinite(u["z"])]
    if not units:
        return None
    w = np.array([math.sqrt(u["n"]) for u in units])
    z = np.array([u["z"] for u in units])
    zc = float((w * z).sum() / math.sqrt((w * w).sum()))
    sg = np.sign(zc)
    return {"z": zc, "p": float(2 * stats.norm.sf(abs(zc))), "n_patients": len(units), "n": int(sum(u["n"] for u in units)),
            "mean_r": float(np.average([u["r"] for u in units], weights=w)),
            "n_concordant": int(sum(1 for u in units if np.sign(u["z"]) == sg and u["z"] != 0))}


def tmb_value(row, vs, power=True):
    n = row["snv"][vs]
    if not row["callable_mb"] or not row["power"]:
        return float("nan")
    return n / row["callable_mb"] / row["power"]


def add_derived(rows):
    """cumulative TMB per variant set along the listed-parent chain, rate per branch length"""
    by = {(r["patient"], r["id"]): r for r in rows}
    for r in rows:
        r["_tmb"] = {vs: tmb_value(r, vs) for vs in VARIANT_SETS}
    for r in rows:
        cum = {vs: 0.0 for vs in VARIANT_SETS}
        node, seen = r, set()
        while node is not None and id(node) not in seen:
            seen.add(id(node))
            for vs in VARIANT_SETS:
                cum[vs] += node["_tmb"][vs]
            node = by.get((r["patient"], node["parent"])) if node["parent"] else None
        r["_cum"] = cum


def measure_value(r, vs, meas):
    if meas == "branch":
        return r["_tmb"][vs]
    if meas == "cumulative":
        return r["_cum"][vs]
    bl = r["cov"].get("branch_length")
    return r["_tmb"][vs] / bl if bl and bl > 0 else float("nan")


def covariate_keys(rows):
    keys = [k for k in COV_META if any(r["cov"].get(k) is not None for r in rows)]
    sig = sorted({k for r in rows for k in r["cov"] if k.startswith("sig:")})
    return keys, sig


def run_stats(rows, patients, vs, meas, terminal, min_cells, drop_single, seed=1):
    keys, sigs = covariate_keys(rows)
    res = []
    for ck in keys + sigs:
        per, units = {}, []
        for pn in patients:
            rr = [r for r in rows if r["patient"] == pn and r["n_cells"] >= min_cells and (r["kind"] != "singleton" or not drop_single)
                  and (r["terminal"] or not terminal) and not r["is_trunk"]]
            if ck == "n_cells":
                continue
            x = [r["cov"].get(ck, float("nan")) if r["cov"].get(ck) is not None else float("nan") for r in rr]
            y = [measure_value(r, vs, meas) for r in rr]
            s = [r["n_cells"] for r in rr]
            t = partial_spearman(x, y, s, np.random.default_rng(seed))
            if t:
                per[pn] = {k: rnd(v, 5) if isinstance(v, float) else v for k, v in t.items()}
                units.append(t)
        comb = stouffer(units)
        if per:
            res.append({"covariate": ck, "per_patient": per, "combined": {k: rnd(v, 5) if isinstance(v, float) else v for k, v in comb.items()} if comb else None})
    res.sort(key=lambda r: (r["combined"]["p"] if r["combined"] else 2))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=DATA)
    ap.add_argument("--out", default=os.path.join(DATA, "_cohort", "tmb_clade_assoc.json"))
    ap.add_argument("--patients", nargs="*", default=PATIENTS)
    a = ap.parse_args()
    allrows, pinfo, na = [], {}, {}
    for p in a.patients:
        log("patient", p)
        rows, info, why = build_patient(p, a.data)
        if why:
            na[p] = why
            log("  ", why)
            continue
        allrows += rows
        pinfo[p] = info
        log("  rows", len(rows), "clades", sum(1 for r in rows if r["kind"] == "clade"), info["reconciliation"])
    add_derived(allrows)
    pats = list(pinfo)
    default_stats = {}
    for vs in VARIANT_SETS:
        for meas in MEASURES:
            for term in (True, False):
                key = f"{vs}|{meas}|{'terminal' if term else 'all'}"
                default_stats[key] = run_stats(allrows, pats, vs, meas, term, 5, True)
        log("stats", vs)
    keys, sigs = covariate_keys(allrows)
    cov_meta = [{"key": k, "label": COV_META[k][0], "group": COV_META[k][1], "unit": COV_META[k][2]} for k in keys if k != "n_cells"]
    cov_meta += [{"key": s, "label": "Branch signature " + s[4:] + " (fraction of branch SNVs)", "group": "signatures", "unit": "fraction of branch SNVs"} for s in sigs]
    for r in allrows:
        r.pop("_tmb", None)
        r.pop("_cum", None)
    out = {
        "format": "gos-sc-tmb-clade/1",
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "patients": pinfo,
        "not_available": na,
        "variant_sets": {"all": "all SNVs mapped to the branch", "conf": "confidently mapped SNVs (map_confidence >= %g: next-best placement >= %g log-lik units worse)" % (CONF_MIN, CONF_MIN),
                         "multi": "SNVs with alt reads in >= 2 cells (SNVs seen in a single cell removed)", "confmulti": "confidently mapped and seen in >= 2 cells"},
        "measures": {"branch": "Branch TMB (SNVs mapped to the branch / callable Mb / discovery power)",
                     "cumulative": "Cumulative TMB (sum of branch TMB along the listed root-to-clade path)",
                     "rate": "Branch TMB per CellPhy branch length (SNVs / Mb per unit branch length)"},
        "covariates": cov_meta,
        "params": {"conf_min": CONF_MIN, "min_linked_cells": MIN_LINKED, "n_perm": N_PERM, "perm_below_n": PERM_BELOW, "min_n": MIN_N,
                   "default": {"variant_set": "all", "measure": "branch", "terminal_only": True, "min_cells": 5, "remove_singletons": True}},
        "method": "Unit = tree clade (branch). Branch TMB = SNVs mapped to the branch / callable Mb of the clade / discovery power "
                  "(tmb.json); single-cell branches use the cell's callable Mb and truncal-site sensitivity and are a lower bound (a site needs "
                  "alt reads in >= 2 cells to enter the candidate set). Cumulative TMB sums branch TMB over the listed ancestors. A per-cell "
                  "normalisation (SNVs / clade cells) is not meaningful because every cell of the clade carries every branch SNV; clade size is "
                  "a covariate instead, and the rate per CellPhy branch length is offered (it shares information with the SNV count). "
                  "Covariates are means / fractions over the clade's cells with a DNA-linked RNA / CN record (>= %d cells; single-cell "
                  "branches use their own cell). Per patient: partial Spearman of covariate vs TMB adjusted for rank(clade size); z = atanh(r) "
                  "* sqrt(n-3); permutation p (%d permutations) when n < %d. Across patients: Stouffer z weighted by sqrt(n) with sign "
                  "concordance. Nested clades share cells, so tests default to TERMINAL clades only (no listed clade inside; trunk clades "
                  "excluded); plots show all clades." % (MIN_LINKED, N_PERM, PERM_BELOW),
        "caveats": ["Clade-level n is small (6-12 terminal clades per patient); p-values are exploratory and not multiplicity-corrected.",
                    "Signature fractions are computed from the same branch SNVs as the TMB; they are not independent of it.",
                    "SNVs on unlisted small clades (2-4 cells or < 10 SNVs) and outside-tumour calls are not attributed to any row (see patients.*.reconciliation).",
                    "Single-cell branches: private TMB is a lower bound and sensitivity-corrected per cell."],
        "default_stats": default_stats,
        "clades": allrows,
    }
    if os.path.exists(a.out):
        bak = a.out + ".bak_" + datetime.date.today().strftime("%y%m%d")
        if not os.path.exists(bak):
            shutil.copy2(a.out, bak)
    with open(a.out, "w") as f:
        json.dump(out, f, separators=(",", ":"), allow_nan=False)
    log("wrote", a.out, os.path.getsize(a.out), "bytes;", len(allrows), "rows")


if __name__ == "__main__":
    main()
