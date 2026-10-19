#!/usr/bin/env python3
"""Clade DE against a widened comparison group (sister, next outgroups, all other tumour cells) for ONE patient clade.

Reuses gos_sc_rna_covariate_de.py (same stratified Wilcoxon, same strata/QC fallbacks, same min pct, same enrichment sets).
Comparators for clade C with parent chain C -> A1 -> A2 ...:
  sister       cells of A1 not in C
  level_k      cells of A(k) not in A(k-1)   (ring k, A0 = C)
  outgroup     first ring with >= MIN_OUT RNA tumour cells on its own
  cumulative   rings 1..k, smallest k with >= MIN_OUT cells
  all_other    every tumour RNA cell not in C (C's subclades are inside C and so excluded)
Also: Hallmark module score (mean z-scored log-expression of the set) test and a confounder panel.

usage: gos_sc_clade_outgroup_de.py --patient MGH302 --clade MGH302_n234 --out DIR --staging DIR
"""
import argparse, json, math, os, sys, importlib.util
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("covde", os.path.join(HERE, "gos_sc_rna_covariate_de.py"))
cd = importlib.util.module_from_spec(spec); spec.loader.exec_module(cd)
MIN_OUT = 30


def rings(P, cid):
    """ordered list of (name, mask) for ancestors' rings"""
    out, cur, prev = [], cid, P.clades[cid]["mask"]
    k = 0
    while True:
        par = P.clades[cur]["parent"]
        if not par or par not in P.clades:
            break
        k += 1
        pm = P.clades[par]["mask"]
        out.append((f"ring{k}:{par}", pm & ~prev, par))
        prev, cur = pm, par
    return out


def comparators(P, cid, tumour_all=True):
    m = P.clades[cid]["mask"]
    rs = rings(P, cid)
    comps = {}
    sizes = [{"level": n, "n_rna_tumour": int(r.sum())} for n, r, _ in rs]
    if rs:
        comps["sister"] = (rs[0][1], f"sister ({rs[0][2]} minus clade)")
    out_done = cum_done = False
    cum = np.zeros(P.n, dtype=bool)
    for n, r, par in rs:
        cum |= r
        if not out_done and r.sum() >= MIN_OUT and n != rs[0][0]:
            comps["outgroup"] = (r, f"next outgroup ring {n}"); out_done = True
        elif not out_done and r.sum() >= MIN_OUT and n == rs[0][0]:
            out_done = True  # the sister itself is already wide enough
        if not cum_done and cum.sum() >= MIN_OUT:
            comps["cumulative"] = (cum.copy(), f"sister + outgroup rings up to {n} (cumulative)"); cum_done = True
    comps["all_other"] = (~m, "all other tumour RNA cells")
    return comps, sizes


def module_score(P, genes):
    idx = [i for i, g in enumerate(P.genes) if g in genes and not g.startswith("MT-")]
    X = P.X[idx]
    z = (X - X.mean(axis=1, keepdims=True)) / np.maximum(X.std(axis=1, keepdims=True), 1e-6)
    return z.mean(axis=0), len(idx)


def ols_group(y, g, covs):
    """t-test of group coefficient in y ~ 1 + g + covs"""
    A = np.column_stack([np.ones(len(y)), g] + covs)
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    res = y - A @ beta
    dof = len(y) - A.shape[1]
    s2 = res @ res / dof
    cov = s2 * np.linalg.pinv(A.T @ A)
    t = beta[1] / math.sqrt(cov[1, 1])
    return float(beta[1]), float(t), float(2 * stats.t.sf(abs(t), dof))


def score_test(P, a, b, score):
    import copy
    S = copy.copy(P)
    S.X = score[None, :].astype(np.float32); S.E = S.X; S.nz = np.ones_like(S.X, dtype=bool)
    S.genes = np.array(["score"]); S.rank_cache = {}
    r = cd.de_patient(S, a, b)
    if r is None:
        return None
    ia, ib = np.array(r["ia"]), np.array(r["ib"])
    pos = {k: i for i, k in enumerate(P.keep)}
    ua = np.array([pos[i] for i in ia]); ub = np.array([pos[i] for i in ib])
    return S, r, ua, ub


def confounders(P, a, b, score, genes_ox):
    used = np.where(a | b)[0]
    g = a[used].astype(float)
    y = score[used]
    mt = P.mtv[used]; lib = P.lib[used]; cyc = P.cycling[used].astype(float)
    ok = np.isfinite(mt)
    res = {"n_cells": int(len(used)), "n_a": int(a.sum()), "n_b": int(b.sum())}
    res["mean_diff_raw"] = float(y[g == 1].mean() - y[g == 0].mean())
    res["corr_score_percent_mt"] = float(stats.spearmanr(y[ok], mt[ok])[0])
    res["corr_score_log10_nCount"] = float(stats.spearmanr(y, lib)[0])
    res["corr_score_cycling"] = float(stats.spearmanr(y, cyc)[0])
    res["delta_percent_mt_a_minus_b"] = float(np.nanmean(mt[g == 1]) - np.nanmean(mt[g == 0]))
    res["delta_log10_nCount"] = float(lib[g == 1].mean() - lib[g == 0].mean())
    res["frac_cycling_a"] = float(cyc[g == 1].mean()); res["frac_cycling_b"] = float(cyc[g == 0].mean())
    fits = {}
    sd = lambda v: (v - v[ok].mean()) / max(v[ok].std(), 1e-9)
    fits["unadjusted"] = ols_group(y, g, [])
    fits["+percent_mt"] = ols_group(y[ok], g[ok], [sd(mt)[ok]])
    fits["+log10_nCount"] = ols_group(y, g, [sd(lib)])
    fits["+cycling"] = ols_group(y, g, [cyc])
    fits["+all three"] = ols_group(y[ok], g[ok], [sd(mt)[ok], sd(lib)[ok], cyc[ok]])
    res["ols_group_effect(beta,t,p)"] = {k: [round(x, 5) for x in v] for k, v in fits.items()}
    # same stratified Wilcoxon with weaker/other strata
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patient", required=True); ap.add_argument("--clade", required=True)
    ap.add_argument("--data", default=cd.ROOT); ap.add_argument("--out", required=True); ap.add_argument("--staging", required=True)
    ap.add_argument("--genesets", default="/gpfs/commons/projects/imielinski_web/sclarke/gOS_dev/genesets")
    ap.add_argument("--set", default="HALLMARK_OXIDATIVE_PHOSPHORYLATION")
    a_ = ap.parse_args()
    cd.DATA_ROOT = a_.data
    os.makedirs(a_.out, exist_ok=True); os.makedirs(a_.staging, exist_ok=True)
    P = cd.load_patient(a_.data, a_.patient)
    cid = a_.clade
    m = P.clades[cid]["mask"]
    comps, sizes = comparators(P, cid)
    sets = cd.load_sets(a_.genesets)
    ox = next(s for s in sets if s[1] == a_.set)[2]
    score, nox = module_score(P, set(ox))
    summary = {"patient": a_.patient, "clade": cid, "n_clade_rna_tumour": int(m.sum()), "ring_sizes": sizes, "oxphos_genes_used": nox, "comparators": {}}
    tests, cells = [], {}
    for name, (cm, desc) in comps.items():
        cm = cm & ~m
        s = {"description": desc, "n_b": int(cm.sum())}
        if cm.sum() < cd.MIN_GROUP:
            s["skipped"] = "too few cells"; summary["comparators"][name] = s; continue
        r = cd.de_patient(P, m, cm)
        if r is None:
            s["skipped"] = "strata plan failed"; summary["comparators"][name] = s; continue
        rows = cd.combine([(a_.patient, r, P.genes)])
        rows.sort(key=lambda x: (x["p_val"], -abs(x["avg_log2FC"])))
        sig = [x for x in rows if x["q_val"] < 0.05]
        strong = [x for x in sig if abs(x["avg_log2FC"]) >= cd.LFC_ROBUST]
        enr = cd.enrich(rows, sets)
        s.update({"n_a": r["nA"], "n_b_used": r["nB"], "adjust": r["adjust"], "n_tested": len(rows), "n_sig_q05": len(sig),
                  "n_sig_q05_abslfc_ge_0.5": len(strong), "n_robust_filter(>=2 patients)": 0,
                  "top_genes": [(x["gene"], round(x["avg_log2FC"], 3), float(f"{x['q_val']:.3g}")) for x in rows[:10]],
                  "enrich_up_top5": [(o["term"], o["n"], round(o["mean_z"], 3), float(f"{o['p']:.3g}"), float(f"{o['q']:.3g}")) for o in enr["up"][:5]],
                  "enrich_down_top5": [(o["term"], o["n"], round(o["mean_z"], 3), float(f"{o['p']:.3g}"), float(f"{o['q']:.3g}")) for o in enr["down"][:5]]})
        oh = next((o for o in enr["up"] + enr["down"] if o["term"] == a_.set), None)
        s["oxphos_enrichment_in_top8"] = oh
        # full-enrichment OXPHOS p regardless of top8
        zmap = {x["gene"]: x["z"] for x in rows}
        ins = np.array([zmap[g] for g in ox if g in zmap]); outs = np.array([v for g, v in zmap.items() if g not in set(ox)])
        s["oxphos_gene_set_enrichment"] = {"n": int(len(ins)), "mean_z": float(ins.mean()), "p": float(stats.mannwhitneyu(ins, outs).pvalue)}
        # module score under several strata schemes
        sc = {}
        for lab, qc in (("strata_plate+cycling+qc(live)", True), ("strata_plate+cycling", False)):
            cd.QC_STRATA = qc
            out = score_test(P, m, cm, score)
            if out:
                S, rr, ua, ub = out
                sc[lab] = {"z": float(rr["z"][0]), "p": float(rr["p"][0]), "adjust": rr["adjust"], "mean_diff": float(score[ua].mean() - score[ub].mean())}
        cd.QC_STRATA = True
        s["oxphos_score_stratified"] = sc
        s["confounders"] = confounders(P, m, cm, score, ox)
        summary["comparators"][name] = s
        # tables
        with open(os.path.join(a_.out, f"genes_{name}.tsv"), "w") as f:
            f.write("gene\tavg_log2FC\tpct_clade\tpct_comp\tz\tp\tq\n")
            for x in rows:
                f.write(f"{x['gene']}\t{x['avg_log2FC']:.4f}\t{x['pct_1']:.3f}\t{x['pct_2']:.3f}\t{x['z']:.3f}\t{x['p_val']:.3g}\t{x['q_val']:.3g}\n")
        # staged entry in the live test format
        tid = f"clade:{cid}|vs:{name}"
        n_up = [x["gene"] for x in rows if x["avg_log2FC"] > 0 and x["q_val"] < 0.05][:8]
        n_dn = [x["gene"] for x in rows if x["avg_log2FC"] < 0 and x["q_val"] < 0.05][:8]
        tests.append(cd.clean({
            "id": tid, "family": "clade", "clade_info": None, "sensitivity": True, "label": f"{cid}: clade vs {desc}", "group_a": cid, "group_b": desc, "gene": None,
            "k_patients": 1, "combinable": False, "units": "single clade",
            "patients": [{"patient": a_.patient, "n_a": r["nA"], "n_b": r["nB"], "adjust": r["adjust"], "n_strata": r["n_strata"], "d_mt": r.get("d_mt"), "d_lib": r.get("d_lib"), "clade": cid}],
            "n_tested": len(rows), "n_sig_q05": len(sig), "n_robust": 0, "n_up_robust": 0, "n_down_robust": 0, "top_up": n_up, "top_down": n_dn, "robust_up": [], "robust_down": [],
            "genes": rows[:cd.TOP_TABLE],
            "volcano": {"genes": [x["gene"] for x in rows[:cd.VOLCANO_N]], "lfc": [round(x["avg_log2FC"], 3) for x in rows[:cd.VOLCANO_N]], "p": [float(f"{x['p_val']:.3g}") for x in rows[:cd.VOLCANO_N]], "q": [float(f"{x['q_val']:.3g}") for x in rows[:cd.VOLCANO_N]]},
            "enrichment": enr}))
        cells[tid] = {a_.patient: {"a": r["ia"], "b": r["ib"]}}
        cd.log(name, s["n_b"], s["n_sig_q05"], s["oxphos_gene_set_enrichment"], sc)
    json.dump(cd.clean(summary), open(os.path.join(a_.out, "summary.json"), "w"), indent=1, allow_nan=False)
    json.dump({"format": "gos-sc-rna-covariate-de-staged/1", "tests": tests}, open(os.path.join(a_.staging, "tests.json"), "w"), separators=(",", ":"), allow_nan=False)
    json.dump({"format": "gos-sc-rna-covariate-de-cells/1", "tests": cells}, open(os.path.join(a_.staging, "cells.json"), "w"), separators=(",", ":"))


if __name__ == "__main__":
    sys.exit(main())
