"""ecDNA copy number as a source of plasticity.

Per patient x oncogene (amplified in >=10% of tumour cells), using two per-cell
copy-number sources side by side:
  walk  : summed copies of the ecDNA walks (walks.json) that carry the gene
  total : the gene's total CN from rna/covariates.json ("CN: <gene>", per-cell
          JaBbA / CN data; works where ecDNA calling is poor, e.g. MGH303)
For each source:
  * distribution: CV, dispersion (VMR), fraction of cells at high CN
  * heritability on the DNA tree: nearest-relative Spearman of CN, permutation
    null within clone (does CN track the tree beyond the clone labels?)
  * CN vs state: per outcome (Neftel state one-vs-rest, cycling, program scores)
    OLS  outcome ~ clone + cycling(S/G2M) + CN_z  -> beta per SD of CN, p, BH q
  * residual test (the plasticity test): CN minus its clade mean (leave-one-out;
    clone if no clade) -> same regression with the residual
  * dosage -> expression of the gene (Spearman, and clone-adjusted)
Caveats stated in the output: S-phase inflates CN (cycling is a covariate;
outcome "cycling" itself is confounded by construction); only RNA-linked cells
enter state tests.

Writes <out>/_cohort/ecdna_plasticity.json and <out>/<P>/ecdna_plasticity.json.
usage: gos_sc_ecdna_plasticity.py OUT_DIR [PATIENT ...]
"""
import gzip
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import sparse, stats

DATA = os.environ.get("GOS_SC", "/gpfs/commons/projects/imielinski_web/sclarke/gOS_GBM_sc")
STORY = os.environ.get("STORY_OUT", "/gpfs/commons/groups/imielinski_lab/projects/single_cell_GBM/analysis/story")
sys.path.insert(0, os.path.join(STORY, "scripts"))
from heritability import leaf_distances, nearest, parse_newick  # noqa: E402

PATIENTS = ["BWH69", "BWH70", "BWH75", "MGH285", "MGH301", "MGH302", "MGH303", "MGH304", "MGH305"]
ONCO = ["EGFR", "MDM2", "CDK4", "MYCN", "PDGFRA", "MET", "GLI1", "DDIT3", "MYC", "MDM4", "CCND1", "CCND2", "KIT", "HMGA2", "CDK6", "MYB"]
AMP_CN = 4        # total CN considered amplified
HIGH_CN = 10      # "high" copy number (matches the 10-49 / 50+ levels in covariates)
MIN_FRAC = 0.10
NPERM = 1000
rng = np.random.default_rng(11)
CAVEATS = [
    "S-phase cells carry inflated copy number; cycling (S/G2M) is a covariate in every CN-state model and the 'cycling' outcome is confounded by construction.",
    "State tests use only cells with both DNA and RNA (n given per row); CN distribution and heritability use every tumour cell with DNA.",
    "Total CN (JaBbA / per-cell CN) is gene-level and includes ecDNA plus any chromosomal amplification; walk copies cover only curated circular walks.",
    "Residuals are CN minus the leave-one-out clade mean (clone when no clade): a CN change that is not explained by lineage.",
]


def clean(o):
    """NaN / inf -> None so the files are valid JSON for the browser."""
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def bh(p):
    p = np.asarray(p, float)
    if len(p) == 0:
        return p
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def ols_beta(y, x, covars):
    """beta, p of rank(x) in y ~ 1 + covars + rank(x) (standardised); beta in outcome SD per SD of the CN rank."""
    n = len(y)
    x = stats.rankdata(x)  # rank-based: one outlier cell cannot carry the effect
    sd = np.std(x)
    if n < 15 or sd == 0 or np.std(y) == 0:
        return None
    z = (x - x.mean()) / sd
    X = np.column_stack([np.ones(n), covars, z]) if covars.shape[1] else np.column_stack([np.ones(n), z])
    # drop collinear covariate columns
    keep = [0]
    for j in range(1, X.shape[1] - 1):
        if np.linalg.matrix_rank(X[:, keep + [j]]) > len(keep):
            keep.append(j)
    keep.append(X.shape[1] - 1)
    X = X[:, keep]
    if n - X.shape[1] < 5:
        return None
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ coef
    s2 = res @ res / (n - X.shape[1])
    cov = s2 * np.linalg.pinv(X.T @ X)
    se = np.sqrt(cov[-1, -1])
    if se == 0:
        return None
    t = coef[-1] / se
    p = 2 * stats.t.sf(abs(t), n - X.shape[1])
    ysd = np.std(y)
    return float(coef[-1] / ysd), float(p)   # beta in outcome-SD per CN-SD


def dummies(labels):
    s = pd.Series(labels).astype(str)
    d = pd.get_dummies(s, drop_first=True).values.astype(float)
    return d if d.shape[1] else np.zeros((len(s), 0))


def nn_rho(v, nn):
    return float(stats.spearmanr(v, v[nn])[0])


def heritability(v, nn, clones):
    obs = nn_rho(v, nn)
    null = []
    for _ in range(NPERM):
        w = v.copy()
        for g in np.unique(clones):
            idx = np.where(clones == g)[0]
            w[idx] = rng.permutation(v[idx])
        null.append(nn_rho(w, nn))
    null = np.array(null)
    return {"rho": obs, "null_mean": float(null.mean()), "p_within_clone": float((1 + (null >= obs).sum()) / (1 + NPERM))}


def residual_vs_group(v, groups):
    """v minus leave-one-out group mean."""
    out = np.full(len(v), np.nan)
    for g in np.unique(groups):
        idx = np.where(groups == g)[0]
        if len(idx) < 3:
            continue
        s = v[idx].sum()
        out[idx] = v[idx] - (s - v[idx]) / (len(idx) - 1)
    return out


def load_rna(p):
    d = f"{DATA}/data/{p}/rna"
    cells = pd.DataFrame(json.load(open(f"{d}/cells.json"))["cells"])
    genes = [g.strip() for g in open(f"{d}/genes.tsv") if g.strip()]
    X = sparse.csc_matrix((np.fromfile(f"{d}/matrix.data.f32", "<f4"), np.fromfile(f"{d}/matrix.indices.i32", "<i4"),
                           np.fromfile(f"{d}/matrix.indptr.i32", "<i4")), shape=(len(cells), len(genes)))
    return cells, {g: i for i, g in enumerate(genes)}, X


def process(p, scores, mods):
    d = f"{DATA}/data/{p}"
    cov = json.load(open(f"{d}/rna/covariates.json"))
    cc = cov["cells"]
    names = set(f["name"] for f in cov["fields"])
    parent, leaf, length = parse_newick(open(f"{d}/tree.nwk").read())
    rna, gi, X = load_rna(p)
    rna["row"] = np.arange(len(rna))
    rna = rna[rna.cell_id.notna()].set_index("cell_id")
    walks = None
    if os.path.exists(f"{d}/walks.json"):
        W = json.load(open(f"{d}/walks.json"))
        walks = [w for w in W["walks"] if w.get("ncells", 0) >= 5]
    # tumour DNA cells on the tree
    ids = [c for c in cc if c in leaf and str(cc[c].get("DNA clone", "")).lower() not in ("normal", "none", "")]
    if len(ids) < 20:
        return None
    clone = np.array([str(cc[c].get("DNA clone")) for c in ids])
    clade = np.array([str(cc[c].get("Clade")) if cc[c].get("Clade") not in (None, "") else str(cc[c].get("DNA clone")) for c in ids])
    E, P = leaf_distances(parent, length, [leaf[c] for c in ids])
    nn = nearest(E, P)
    # candidate genes
    genes = [g for g in ONCO if f"CN: {g}" in names or (walks and any(g in w["genes"] for w in walks))]
    sc = scores[scores.patient == p].pivot_table(index="rna_id", columns="module", values="score")
    outcome_names = ["AC", "MES", "NPC", "OPC", "cycling"] + [m for m in sc.columns if mods[m]["meta"] or p in mods[m]["patients"]]
    labels = {"AC": "AC state", "MES": "MES state", "NPC": "NPC state", "OPC": "OPC state", "cycling": "Cycling (S/G2M)"}
    labels.update({m: f"{m} {mods[m]['label']}" for m in sc.columns})
    in_rna = np.array([c in rna.index for c in ids])
    rows, cells_out = [], {}
    for g in genes:
        srcs = {}
        if f"CN: {g}" in names:
            srcs["total"] = np.array([cc[c].get(f"CN: {g}", np.nan) if cc[c].get(f"CN: {g}") is not None else np.nan for c in ids], float)
        if walks:
            ws = [w for w in walks if g in w["genes"]]
            if ws:
                srcs["walk"] = np.array([sum(w["cells"].get(c, 0) for w in ws) for c in ids], float)
        amp = (srcs["total"] >= AMP_CN).mean() if "total" in srcs else 0
        wamp = (srcs["walk"] > 0).mean() if "walk" in srcs else 0
        if max(np.nan_to_num(amp), wamp) < MIN_FRAC:
            continue
        row = {"patient": p, "gene": g, "n_cells": len(ids), "n_rna": int(in_rna.sum()), "sources": {}}
        for src, v in srcs.items():
            ok = ~np.isnan(v)
            if ok.sum() < 20 or np.std(v[ok]) == 0:
                continue
            vv = np.where(ok, v, np.nanmedian(v))
            mean = float(np.mean(v[ok]))
            s = {"mean": mean, "median": float(np.median(v[ok])), "max": float(np.max(v[ok])),
                 "cv": float(np.std(v[ok]) / mean) if mean > 0 else None,
                 "vmr": float(np.var(v[ok]) / mean) if mean > 0 else None,
                 "frac_amplified": float((v[ok] >= AMP_CN).mean()) if src == "total" else float((v[ok] > 0).mean()),
                 "frac_high": float((v[ok] >= HIGH_CN).mean())}
            s["heritability"] = heritability(vv, nn, clone)
            # within-clone CV: variability not explained by clone identity
            wc = [np.std(vv[clone == k]) / np.mean(vv[clone == k]) for k in np.unique(clone) if (clone == k).sum() >= 8 and np.mean(vv[clone == k]) > 0]
            s["cv_within_clone"] = float(np.mean(wc)) if wc else None
            res = residual_vs_group(vv, clade)
            s["resid_sd"] = float(np.nanstd(res))
            # state models on RNA-linked cells
            sel = np.where(in_rna & ~np.isnan(res))[0]
            rid = [rna.loc[ids[i], "rna_id"] for i in sel]
            st = rna.loc[[ids[i] for i in sel], "state"].values
            ph = rna.loc[[ids[i] for i in sel], "Phase"].values
            cyc = np.isin(ph, ["S", "G2M"]).astype(float)
            tests = []
            for o in outcome_names:
                if o in ("AC", "MES", "NPC", "OPC"):
                    y = (st == o).astype(float)
                    if pd.isna(st).any():
                        continue
                elif o == "cycling":
                    y = cyc
                else:
                    y = np.array([sc[o].get(r, np.nan) for r in rid], float)
                ok2 = ~np.isnan(y)
                if ok2.sum() < 20 or y[ok2].std() == 0:
                    continue
                cv_ = [dummies(clone[sel][ok2])]
                if o != "cycling":
                    cv_.append(cyc[ok2][:, None])
                cov_m = np.column_stack(cv_)
                a = ols_beta(y[ok2], vv[sel][ok2], cov_m)
                r_ = ols_beta(y[ok2], res[sel][ok2], cov_m)
                if a is None or r_ is None:
                    continue
                tests.append({"outcome": o, "label": labels[o], "n": int(ok2.sum()), "beta": a[0], "p": a[1], "beta_resid": r_[0], "p_resid": r_[1]})
            s["tests"] = tests
            # dosage -> expression
            s["expr"] = None
            if g in gi:
                ri = np.where(in_rna)[0]
                expr = np.asarray(X[rna.loc[[ids[i] for i in ri], "row"].values, gi[g]].todense()).ravel()
                if (expr > 0).sum() >= 10 and np.std(expr) > 0 and np.std(vv[ri]) > 0:
                    rho, pr = stats.spearmanr(vv[ri], expr)
                    # clone-adjusted: rank-residualise both on clone
                    ex_r = pd.Series(stats.rankdata(expr)).groupby(clone[ri]).transform(lambda x: x - x.mean()).values
                    cn_r = pd.Series(stats.rankdata(vv[ri])).groupby(clone[ri]).transform(lambda x: x - x.mean()).values
                    rho_c, p_c = stats.pearsonr(cn_r, ex_r) if np.std(cn_r) > 0 and np.std(ex_r) > 0 else (np.nan, np.nan)
                    s["expr"] = {"rho": float(rho), "p": float(pr), "rho_within_clone": None if np.isnan(rho_c) else float(rho_c),
                                 "p_within_clone": None if np.isnan(p_c) else float(p_c), "detected": float((expr > 0).mean()), "n": len(ri)}
            row["sources"][src] = s
        if row["sources"]:
            rows.append(row)
            cells_out[g] = {"ids": ids, "clone": clone.tolist(), "clade": clade.tolist(),
                            **{src: v.tolist() for src, v in srcs.items()},
                            **{src + "_resid": residual_vs_group(np.where(np.isnan(v), np.nanmedian(v), v), clade).tolist() for src, v in srcs.items()}}
    outcomes = {}
    for o in outcome_names:
        vals = []
        for c, lk in zip(ids, in_rna):
            if not lk:
                vals.append(None)
            elif o in ("AC", "MES", "NPC", "OPC"):
                st_ = rna.loc[c, "state"]
                vals.append(None if pd.isna(st_) else float(st_ == o))
            elif o == "cycling":
                ph_ = rna.loc[c, "Phase"]
                vals.append(None if pd.isna(ph_) else float(ph_ in ("S", "G2M")))
            else:
                x_ = sc[o].get(rna.loc[c, "rna_id"], np.nan)
                vals.append(None if pd.isna(x_) else round(float(x_), 4))
        outcomes[o] = {"label": labels[o], "values": vals}
    return rows, cells_out, outcomes


def main():
    out_dir = sys.argv[1]
    pats = sys.argv[2:] or PATIENTS
    scores = pd.read_csv(f"{STORY}/out/module_scores.tsv.gz", sep="\t")
    mods = {m["id"]: m for m in json.load(open(f"{STORY}/out/modules.json"))["modules"]}
    allrows = []
    detail = {}
    for p in pats:
        try:
            r = process(p, scores, mods)
        except FileNotFoundError as e:
            print("skip", p, e, file=sys.stderr)
            continue
        if not r:
            continue
        rows, cells_out, outcomes = r
        allrows += rows
        print("==", p, [x["gene"] for x in rows], file=sys.stderr)
        detail[p] = {"outcomes": outcomes, "genes": cells_out, "rows": rows}
    # BH across all CN-state and residual tests (cohort-wide, per source)
    for key_p, key_q in (("p", "q"), ("p_resid", "q_resid")):
        flat = [(r, s, t) for r in allrows for s in r["sources"].values() for t in s["tests"]]
        for (r, s, t), q in zip(flat, bh([t[key_p] for _, _, t in flat])):
            t[key_q] = float(q)
    for r in allrows:
        for s in r["sources"].values():
            a = min(s["tests"], key=lambda t: t["p"]) if s["tests"] else None
            b = min(s["tests"], key=lambda t: t["p_resid"]) if s["tests"] else None
            s["top_state"] = a and {k: a[k] for k in ("outcome", "label", "beta", "p", "q", "n")}
            s["top_resid"] = b and {k: b[k] for k in ("outcome", "label", "beta_resid", "p_resid", "q_resid", "n")}
    for p, dd in detail.items():
        os.makedirs(f"{out_dir}/{p}", exist_ok=True)
        json.dump(clean({"patient": p, **dd}), open(f"{out_dir}/{p}/ecdna_plasticity.json", "w"))
    os.makedirs(f"{out_dir}/_cohort", exist_ok=True)
    json.dump(clean({"format": "ecdna_plasticity_v1", "amp_cn": AMP_CN, "high_cn": HIGH_CN, "min_frac": MIN_FRAC, "caveats": CAVEATS, "rows": allrows}),
              open(f"{out_dir}/_cohort/ecdna_plasticity.json", "w"))
    print("rows", len(allrows), file=sys.stderr)


if __name__ == "__main__":
    main()
