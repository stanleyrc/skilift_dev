"""ecDNA status from single-cell copy-number dispersion (ecDNA vs chromosomal / integrated).

Principle: ecDNA has no centromere and segregates randomly at mitosis, so its cell-to-cell
copy number varies far beyond measurement noise; chromosomal / integrated amplicons (HSR,
tandem dup, BFB) are inherited stably and stay at the noise floor within a clade.

Per patient:
  * focal amplicons from cn_binned 100 kb integer CN: CN >= max(5, 2 x cell ploidy,
    2 x the cell's median CN on that arm) in >= 5 tumour cells; bins merged across <= 500 kb.
    Driver-gene loci (+-50 kb minimum) inside amplicons are measured too.
  * per tumour cell, per region: mean of gOS 10 kb coverage bins (coverage.arrow, CN units);
    technical variance from alternating 20 kb half-blocks, Var_tech = (mean_A - mean_B)^2 / 4.
  * multiplicative noise floor CV2_ref: median over random same-size autosomal windows
    (non-amplified) in the same cells of (Var_obs - Var_tech) / mean^2.
  * D = (Var_obs - mean Var_tech - CV2_ref m^2) / m per group (all tumour cells, every SNV-tree
    clade); D_sib from tree cherries; 95% bootstrap CIs over cells / cherries.
Call per amplicon x clade (pre-registered, see README):
  not enough (< 10 carriers or median CN < 5); unstable (ecDNA-like): D CI low > 1 and
  D_sib > 1 when >= 5 cherries; stable (chromosomal / integrated): D CI high < 0.5
  (one round of random segregation from a uniform parent gives D = 0.5); else indeterminate.
Patient-level verdict per amplicon: ecDNA-like, chromosomal / integrated, mixed (unstable and
stable in disjoint clades), nested, indeterminate, not enough. gGnome complex-SV labels of the
carrier cells and curated circular walks are reported beside the call, not used by it.

Primary analysis excludes S-phase cells (precompute/sphase.json); a with-S-phase verdict is
reported as sensitivity. Needs numpy, pandas-free; pyarrow for coverage.arrow.

Writes <out>/<P>/ecdna_status.json (with per-cell region means `values`, aligned to `cells`, and
`clade_cells` index lists for the per-clade dot plot) and <out>/_cohort/ecdna_status.json.
usage: gos_sc_ecdna_status.py DATA_DIR OUT_DIR [PATIENT ...]
"""
import json
import os
import sys

import numpy as np

FORMAT = "gos-sc-ecdna-status/1"
CHR = [str(i) for i in range(1, 23)] + ["X"]
LEN = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
       145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
       101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468,
       156040895]
CEN = dict(zip(CHR, [123.4, 93.9, 90.9, 50.0, 48.8, 59.8, 60.1, 45.2, 43.0, 39.8, 53.4, 35.5,
                     17.7, 17.2, 19.0, 36.8, 25.1, 18.5, 26.2, 28.1, 12.0, 15.0, 60.6]))
OFF = np.concatenate([[0], np.cumsum(LEN)])[:-1]
BS = 10000
NB = [int(np.ceil(x / BS)) for x in LEN]
FIRST = np.concatenate([[0], np.cumsum(NB)])[:-1]
NTOT = int(sum(NB))
CI = {c: i for i, c in enumerate(CHR)}
DRIVERS = set("EGFR PDGFRA KIT KDR MET CDK4 CDK6 MDM2 MDM4 MYCN MYC MYCL CCND1 CCND2 CCND3 SOX2 GLI1 "
              "DDIT3 HMGA2 FOXO1 SOX5 BRAF MAP3K1 AKT3 PIK3CA FGFR1 FGFR2 FGFR3 ERBB2 TERT CDK2 "
              "YEATS4 FRS2 CPM LRIG3 WIF1 KRAS RICTOR IGF1R AKT1 NOTCH2 ID2 E2F3".split())
MIN_CARRIERS, MIN_CN, B = 10, 5.0, 1000
WSIZES = (10, 20, 50, 100, 200, 500)
DM_TYPES = {"DM", "Cpxdm"}


def gbin(chrom, pos):
    return int(FIRST[CI[chrom]] + (int(pos) - 1) // BS)


# ---------------------------------------------------------------- pure functions (tested)
def cherries_from_newick(s):
    """Pairs of leaves joined directly at one internal node."""
    s = s.strip().rstrip(";")
    stack, out, i, n = [[]], [], 0, len(s)
    while i < n:
        ch = s[i]
        if ch == "(":
            stack.append([])
            i += 1
        elif ch == ")":
            kids = stack.pop()
            if len(kids) == 2 and all(k[0] == "leaf" for k in kids):
                out.append((kids[0][1], kids[1][1]))
            i += 1
            while i < n and s[i] not in ",()":
                i += 1
            stack[-1].append(("node", None))
        elif ch == ",":
            i += 1
        else:
            j = i
            while j < n and s[j] not in ",()":
                j += 1
            stack[-1].append(("leaf", s[i:j].split(":")[0].strip().strip("'")))
            i = j
    return out


def dispersion(x, vt, cv2):
    """D over the last axis: (Var - mean Var_tech - cv2 m^2) / m. Returns (D, m)."""
    m = x.mean(-1)
    return (x.var(-1, ddof=1) - vt.mean(-1) - cv2 * m ** 2) / m, m


def cv2_ref(W, WV, idx):
    """Median and 90th percentile over windows of (Var - Var_tech) / m^2 for cell rows idx."""
    Wx, Vx = W[idx], WV[idx]
    with np.errstate(all="ignore"):
        c = (np.nanvar(Wx, -2, ddof=1) - np.nanmean(Vx, -2)) / np.nanmean(Wx, -2) ** 2
    return np.nanmedian(c, -1), np.nanpercentile(c, 90, axis=-1)


def call_group(n_car, med, d_lo, d_hi, n_cherries, d_sib):
    if n_car < MIN_CARRIERS or not np.isfinite(med) or med < MIN_CN:
        return "not enough"
    if np.isfinite(d_lo) and d_lo > 1:
        if n_cherries >= 5:
            return "unstable" if (np.isfinite(d_sib) and d_sib > 1) else "indeterminate"
        return "unstable (sister untested)"
    if np.isfinite(d_hi) and d_hi < 0.5:
        return "stable"
    return "indeterminate"


def verdict(group_calls, clade_cells):
    """group_calls: {group_id: call}; clade_cells: {group_id: set(cells)} (ALL excluded from 'mixed')."""
    un = [g for g, c in group_calls.items() if c.startswith("unstable")]
    st = [g for g, c in group_calls.items() if c == "stable"]
    pairs = [(u, s) for u in un for s in st if u != "ALL" and s != "ALL"
             and not (clade_cells.get(u, set()) & clade_cells.get(s, set()))]
    if pairs:
        return "mixed", pairs
    if un and not st:
        return "ecDNA-like", []
    if st and not un:
        return "chromosomal / integrated", []
    if un and st:
        return "nested", []
    if any(c == "indeterminate" for c in group_calls.values()):
        return "indeterminate", []
    return "not enough", []


def detect_amplicons(cn, ploidy, tumour, chroms, min_cells=5, gap=5):
    """cn: cells x 100 kb bins (NaN missing); chroms: list of {chr, first_bin, n_bins}.
    Returns list of dicts (chrom, start, end, core_bins (100 kb indices within chromosome), peak_cn, n_amp_cells)."""
    armmed = np.full_like(cn, np.nan)
    for ch in chroms:
        c, f, n = str(ch["chr"]), ch["first_bin"], ch["n_bins"]
        if c not in CEN:
            continue
        cb = min(int(CEN[c] * 1e6) // 100000, n)
        for a, b in ((f, f + cb), (f + cb, f + n)):
            if b > a:
                with np.errstate(all="ignore"):
                    armmed[:, a:b] = np.nanmedian(cn[:, a:b], 1)[:, None]
    thr = np.fmax(np.maximum(5, 2 * ploidy)[:, None], 2 * armmed)
    with np.errstate(invalid="ignore"):
        amp = (cn >= thr) & tumour[:, None]
    cand = amp.sum(0) >= min_cells
    out = []
    for ch in chroms:
        c, f, n = str(ch["chr"]), ch["first_bin"], ch["n_bins"]
        if c not in CI:
            continue
        idx = np.where(cand[f:f + n])[0]
        if not len(idx):
            continue
        groups = [[idx[0]]]
        for i in idx[1:]:
            (groups[-1].append(i) if i - groups[-1][-1] <= gap else groups.append([i]))
        for g in groups:
            sub = cn[:, f + g[0]:f + g[-1] + 1]
            carr = amp[:, f + g[0]:f + g[-1] + 1].any(1)
            with np.errstate(all="ignore"):
                med = np.nanmedian(sub[carr], 0)
            pk = np.nanmax(med)
            core = [g[0] + j for j in range(len(med)) if med[j] >= 0.5 * pk and cand[f + g[0] + j]]
            out.append(dict(chrom=c, start=g[0] * 100000 + 1, end=(g[-1] + 1) * 100000, core_bins=core,
                            peak_cn=float(pk), n_amp_cells=int(carr.sum())))
    return out, amp


def group_stats(x, vt, W, WV, ci, gi, cherries_idx, xall, vall, rng, cache_key=None, cache=None):
    """Statistics for one group. x/vt: carrier values; ci: carrier row indices; gi: all group rows."""
    if cache is not None and cache_key in cache:
        c0, c90, cb = cache[cache_key]
    else:
        c0, c90 = cv2_ref(W, WV, gi)
        cb, _ = cv2_ref(W, WV, gi[rng.integers(0, len(gi), size=(200, len(gi)))])
        if cache is not None:
            cache[cache_key] = (c0, c90, cb)
    D, m = dispersion(x, vt, c0)
    bi = rng.integers(0, len(x), size=(B, len(x)))
    Db, _ = dispersion(x[bi], vt[bi], cb[rng.integers(0, len(cb), size=B)])
    d_lo, d_hi = np.nanpercentile(Db, [2.5, 97.5])
    ds = ds_lo = ds_hi = np.nan
    ch = cherries_idx
    if len(ch) >= 2:
        d2 = (xall[ch[:, 0]] - xall[ch[:, 1]]) ** 2 / 2
        vtc = (vall[ch[:, 0]] + vall[ch[:, 1]]) / 2
        mc = (xall[ch[:, 0]] + xall[ch[:, 1]]) / 2

        def f(ii):
            mm = mc[ii].mean(-1)
            return (d2[ii].mean(-1) - vtc[ii].mean(-1) - c0 * mm ** 2) / mm
        ds = float(f(np.arange(len(ch))))
        ds_lo, ds_hi = np.nanpercentile(f(rng.integers(0, len(ch), size=(B, len(ch)))), [2.5, 97.5])
    med = float(np.median(x))
    return dict(n_carriers=int(len(x)), median=med, mean=float(m), cv=float(x.std(ddof=1) / m),
                cv2_ref=float(c0), D=float(D), D_lo=float(d_lo), D_hi=float(d_hi),
                D_cv2ref90=float(dispersion(x, vt, c90)[0]), n_cherries=int(len(ch)),
                D_sib=ds, D_sib_lo=float(ds_lo), D_sib_hi=float(ds_hi),
                cv_bio=float(np.sqrt(max(D, 0) / m)), cv_floor_random_segregation=float(np.sqrt(0.5 / m)),
                call=call_group(len(x), med, d_lo, d_hi, len(ch), ds))


# ---------------------------------------------------------------- I/O
def read_coverage(fn):
    import pyarrow.ipc as ipc  # lazy: only the exporter run needs pyarrow
    t = ipc.open_file(fn).read_all()
    x = t.column("x").to_numpy().astype(np.float64)
    y = t.column("y").to_numpy().astype(np.float64)
    ch = np.searchsorted(OFF, x - 0.5, side="right") - 1
    gb = FIRST[np.minimum(ch, len(CHR) - 1)] + np.round((x - OFF[ch] - 1) / BS).astype(np.int64)
    ok = (ch < len(CHR)) & np.isfinite(y) & (gb >= 0) & (gb < NTOT)
    v = np.full(NTOT, np.nan)
    v[gb[ok]] = y[ok]
    return v


def overlaps(chrom, s, e, loc):
    for part in str(loc).split(","):
        try:
            c, rr = part.split(":")
            a, b = rr.split("-")
        except ValueError:
            continue
        if c == chrom and int(a) <= e and int(b) >= s:
            return True
    return False


def r4(v):
    if isinstance(v, float):
        return None if not np.isfinite(v) else float(f"{v:.4g}")
    return v


def run_patient(data, P, seed=20261019):
    rng = np.random.default_rng(seed)
    cnb = json.load(open(f"{data}/{P}/cn_binned.json"))
    cells, normal = cnb["cells"], np.array(cnb["normal"], bool)
    sph = {}
    fn = f"{data}/{P}/precompute/sphase.json"
    if os.path.exists(fn):
        sph = {c["cell_id"]: bool(c.get("s_call")) for c in json.load(open(fn))["cells"]}
    ploidy = []
    for c in cells:
        try:
            m = json.load(open(f"{data}/{c}/metadata.json"))
            m = m[0] if isinstance(m, list) else m
            ploidy.append(float(m.get("ploidy") or 2))
        except Exception:
            ploidy.append(2.0)
    ploidy = np.array(ploidy)
    lv = cnb["levels"]["100kb"]
    nb = int(lv["n_bins"])
    cn = np.fromfile(f"{data}/{P}/{lv['files']['cn']}", dtype="<u2").reshape(len(cells), nb).astype(float)
    cn[cn == 65535] = np.nan
    amps, amp100 = detect_amplicons(cn, ploidy, ~normal, lv["chromosomes"])
    genes = json.load(open(f"{data}/_cohort/gene_spans.json"))["genes"]

    regions = []
    for k, a in enumerate(amps):
        aid = f"amp{k + 1}_chr{a['chrom']}_{a['start'] // 1000000}Mb"
        gin = sorted({g["gene"] for g in genes if g["chrom"] == a["chrom"] and g["end"] >= a["start"]
                      and g["start"] <= a["end"] and g["gene"] in DRIVERS})
        a.update(id=aid, genes=gin, broad=(a["end"] - a["start"]) > 30e6)
        idx = np.unique(np.concatenate([np.arange(gbin(a["chrom"], b * 100000 + 1), gbin(a["chrom"], (b + 1) * 100000) + 1)
                                        for b in a["core_bins"]]))
        regions.append(dict(id=aid, kind="amplicon", amp=aid, chrom=a["chrom"], start=a["start"], end=a["end"], idx=idx,
                            label=",".join(gin)))
        for g in genes:
            if g["gene"] in gin and g["chrom"] == a["chrom"]:
                mid = (g["start"] + g["end"]) // 2
                s0, e0 = max(1, min(g["start"], mid - 50000)), max(g["end"], mid + 50000)
                regions.append(dict(id=f"{aid}:{g['gene']}", kind="gene", amp=aid, chrom=g["chrom"], start=g["start"],
                                    end=g["end"], idx=np.arange(gbin(g["chrom"], s0), gbin(g["chrom"], min(e0, LEN[CI[g['chrom']]])) + 1),
                                    label=g["gene"]))
    mask = np.zeros(NTOT, bool)
    for a in amps:
        mask[gbin(a["chrom"], max(1, a["start"] - 2000000)):gbin(a["chrom"], min(LEN[CI[a["chrom"]]], a["end"] + 2000000)) + 1] = True
    for ch in lv["chromosomes"]:
        c = str(ch["chr"])
        if c in CI:
            for j in np.where(amp100.sum(0)[ch["first_bin"]:ch["first_bin"] + ch["n_bins"]] >= 2)[0]:
                a0 = gbin(c, j * 100000 + 1)
                mask[a0:a0 + 10] = True
    ok = ~mask
    ok[FIRST[CI["X"]]:] = False
    wins = {}
    for size in WSIZES:
        lst, tries = [], 0
        while len(lst) < 150 and tries < 20000:
            tries += 1
            a0 = int(rng.integers(0, FIRST[CI["X"]] - size))
            chn = CHR[int(np.searchsorted(FIRST, a0, side="right") - 1)]
            if ok[a0:a0 + size].all() and gbin(chn, LEN[CI[chn]]) >= a0 + size:
                lst.append(np.arange(a0, a0 + size))
        wins[size] = lst

    blockA = ((np.arange(NTOT) // 2) % 2) == 0
    tum = [i for i in range(len(cells)) if not normal[i] and os.path.exists(f"{data}/{cells[i]}/coverage.arrow")]
    tcells = [cells[i] for i in tum]
    X = np.full((len(tum), len(regions)), np.nan)
    V = np.full_like(X, np.nan)
    WX = {s: np.full((len(tum), len(wins[s])), np.nan) for s in WSIZES}
    WV = {s: np.full_like(WX[s], np.nan) for s in WSIZES}

    def meas(v, idx):
        vv = v[idx]
        g = np.isfinite(vv)
        if g.sum() < 4:
            return np.nan, np.nan
        a, b = blockA[idx] & g, (~blockA[idx]) & g
        if not a.any() or not b.any():
            return vv[g].mean(), np.nan
        return vv[g].mean(), (vv[a].mean() - vv[b].mean()) ** 2 / 4

    for r, i in enumerate(tum):
        v = read_coverage(f"{data}/{cells[i]}/coverage.arrow")
        for k, reg in enumerate(regions):
            X[r, k], V[r, k] = meas(v, reg["idx"])
        for s in WSIZES:
            for k, w in enumerate(wins[s]):
                WX[s][r, k], WV[s][r, k] = meas(v, w)
    V = np.nan_to_num(V)
    pl = ploidy[tum]
    sp = np.array([sph.get(c, False) for c in tcells])
    clades = json.load(open(f"{data}/{P}/clades.json"))["clades"]
    pos = {c: i for i, c in enumerate(tcells)}
    groups = [("ALL", set(tcells))] + [(c["id"], set(c["cells"])) for c in clades]
    cl_cells = dict(groups)
    cher = [(pos[a], pos[b]) for a, b in cherries_from_newick(open(f"{data}/{P}/tree.nwk").read()) if a in pos and b in pos]
    ev = json.load(open(f"{data}/{P}/filtered.events.json"))
    walks = json.load(open(f"{data}/{P}/walks.json"))["walks"] if os.path.exists(f"{data}/{P}/walks.json") else []

    out = []
    cache = {}
    for k, reg in enumerate(regions):
        ws = min(WSIZES, key=lambda s: abs(np.log(s / max(len(reg["idx"]), 1))))
        lab = {}
        for e in ev:
            if e.get("type") == "Complex SV" and overlaps(reg["chrom"], reg["start"], reg["end"], e.get("Genome_Location")):
                for c in str(e.get("cell_ids") or "").split(","):
                    lab.setdefault(c, set()).add(e["vartype"])
        res = {}
        for mode, keep in (("primary", ~sp), ("with_sphase", np.ones(len(tcells), bool))):
            rows, calls = [], {}
            x_all, v_all = X[:, k], V[:, k]
            carrier = np.isfinite(x_all) & (x_all >= np.maximum(4, 1.5 * pl))
            for gid, gcells in groups:
                gi = np.array([pos[c] for c in gcells if c in pos and keep[pos[c]] and np.isfinite(x_all[pos[c]])], int)
                if len(gi) < 5:
                    continue
                ci = gi[carrier[gi]]
                row = dict(group=gid, n_group=int(len(gi)), loss_frac=float(1 - len(ci) / len(gi)))
                if len(ci) >= 3:
                    cs = set(ci.tolist())
                    chi = np.array([p for p in cher if p[0] in cs and p[1] in cs], int).reshape(-1, 2)
                    row.update(group_stats(x_all[ci], v_all[ci], WX[ws], WV[ws], ci, gi, chi, x_all, v_all, rng,
                                           cache_key=(mode, gid, ws), cache=cache))
                    cc = [tcells[i] for i in ci]
                    row["gg_dm_cells"] = sum(1 for c in cc if lab.get(c, set()) & DM_TYPES)
                    row["gg_other_cells"] = sum(1 for c in cc if lab.get(c) and not (lab[c] & DM_TYPES))
                else:
                    row.update(n_carriers=int(len(ci)), call="not enough")
                calls[gid] = row["call"]
                if mode == "primary" and (row["call"] != "not enough" or gid == "ALL"):
                    rows.append({kk: r4(vv) for kk, vv in row.items()})
            vd, pairs = verdict(calls, cl_cells)
            res[mode] = dict(verdict=vd, mixed_pairs=[list(p) for p in pairs[:10]], groups=rows)
        wl = [dict(label=w["label"], ncells=w["ncells"]) for w in walks if w.get("curated") and w.get("circular")
              and any(str(nd["chromosome"]) == reg["chrom"] and nd["start"] <= reg["end"] and nd["end"] >= reg["start"]
                      for nd in w["nodes"])]
        gg = {}
        for c, s in lab.items():
            for t in s:
                gg[t] = gg.get(t, 0) + 1
        out.append(dict(id=reg["id"], kind=reg["kind"], amplicon=reg["amp"], label=reg["label"], chrom=reg["chrom"],
                        start=int(reg["start"]), end=int(reg["end"]), n_bins_10kb=int(len(reg["idx"])),
                        verdict=res["primary"]["verdict"], verdict_with_sphase=res["with_sphase"]["verdict"],
                        mixed_pairs=res["primary"]["mixed_pairs"], groups=res["primary"]["groups"],
                        gGnome_cells_by_type=gg, circular_walks=wl,
                        values=[None if not np.isfinite(v) else round(float(v), 1) for v in X[:, k]]))
    used = {g["group"] for a in out for g in a["groups"]}
    return dict(format=FORMAT, patient=P, n_tumour_cells=len(tcells), n_sphase_excluded=int(sp.sum()),
                cells=tcells, sphase=[int(v) for v in sp],
                clade_cells={g: sorted(pos[c] for c in cl_cells[g] if c in pos) for g in sorted(used) if g != "ALL"},
                rule=RULE, params=dict(min_carriers=MIN_CARRIERS, min_cn=MIN_CN, bootstrap=B, unstable_D_lo=1,
                                       stable_D_hi=0.5, window_sizes_10kb=list(WSIZES), seed=seed),
                amplicons=out)


RULE = ("D = (Var_obs - Var_tech - CV2_ref m^2) / m over carrier cells of a clade (S-phase excluded). "
        "unstable (ecDNA-like): D 95% CI low > 1 and sister-pair D_sib > 1 (when >= 5 cherries); "
        "stable (chromosomal / integrated): D 95% CI high < 0.5 (below one round of random segregation); "
        "not enough: < 10 carriers or median CN < 5; else indeterminate. "
        "Amplicon verdict: mixed if unstable and stable in disjoint clades. gGnome DM/cpxdm and circular walks are "
        "reported beside the call, not used by it.")


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    data, out = argv[0], argv[1]
    pats = argv[2:] or sorted(d for d in os.listdir(data) if os.path.exists(f"{data}/{d}/cn_binned.json")
                              and os.path.exists(f"{data}/{d}/clades.json"))
    cohort = []
    for P in pats:
        res = run_patient(data, P)
        os.makedirs(f"{out}/{P}", exist_ok=True)
        with open(f"{out}/{P}/ecdna_status.json", "w") as f:
            json.dump(res, f, separators=(",", ":"))
        for a in res["amplicons"]:
            g = next((x for x in a["groups"] if x["group"] == "ALL"), {})
            cohort.append(dict(patient=P, id=a["id"], kind=a["kind"], label=a["label"], chrom=a["chrom"], start=a["start"], end=a["end"], verdict=a["verdict"],
                               verdict_with_sphase=a["verdict_with_sphase"], n_carriers=g.get("n_carriers"),
                               median=g.get("median"), D=g.get("D"), D_lo=g.get("D_lo"), D_hi=g.get("D_hi"),
                               stable_clades=[x["group"] for x in a["groups"] if x["call"] == "stable"],
                               unstable_clades=[x["group"] for x in a["groups"] if str(x["call"]).startswith("unstable")],
                               gGnome_cells_by_type=a["gGnome_cells_by_type"]))
        print(P, len(res["amplicons"]), flush=True)
    os.makedirs(f"{out}/_cohort", exist_ok=True)
    with open(f"{out}/_cohort/ecdna_status.json", "w") as f:
        json.dump(dict(format=FORMAT, rule=RULE, rows=cohort), f, separators=(",", ":"))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
