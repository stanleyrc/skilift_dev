#!/usr/bin/env python3
"""Cohort subclonal drivers and subclonal signature shifts -> <data>/_cohort/subclonal_drivers.json

Private analysis exporter for the gOS single-cell "Subclonal" cohort view (scview=subclonal).
Reads only the per-patient gOS files that already exist (no new genotyping):

  filtered.events.json         tiered driver SNVs with carrier cells (Tier 1-2, damaging tier-3 missense)
  snv_matrix.json              CellPhy mutation map: every tree site, its category (truncal / subclonal /
                               private), node, VEP consequence and SBS96 context
  precompute/region_calls.json per-cell alt / depth at the driver sites (dropout-aware clonality)
  tree.nwk                     carriers vs best clade (F1, same score as gOS cladeFit.js)
  clades.json                  called clades with their private (branch) SNVs
  signatures_joint.json        the patient's signature set and the tree-set fits (truncal / subclonal / private)

Sections written:
  hits            one row per driver SNV (gene, consequence, carriers, dropout-aware fraction, tree node,
                  clonal / subclonal call, clade F1, signature posterior)
  genes           length-aware recurrence of subclonal truncating (and non-synonymous) SNVs per gene
  pathways        patients x pathway with (subclonal / clonal) hits
  signatures      per patient truncal vs subclonal contrasts (copied from signatures_joint.json)
  clade_signatures per called clade: NNLS refit of the clade's branch SNVs, bootstrap CI over sites, test vs
                  parent branch set and vs the truncal set, BH q across the cohort

Recurrence method (stated in the JSON):
  Background per patient p: n_p = number of non-truncal (subclonal + private) truncating SNVs (stop_gained,
  splice donor / acceptor, start_lost) anywhere in the patient's CellPhy mutation map. Each is assumed to land in
  the coding territory in proportion to bp (GENCODE v24 CCDS CDS union + 4 bp per coding intron for the
  splice sites), so the expected number in gene g is mu_gp = n_p * L_g / L_total. Observed: distinct non-truncal
  truncating SNVs in gene g. Tests: hits ~ Poisson(sum_p mu_gp) and patients-with-a-hit ~ Poisson-binomial
  with P_p = 1 - exp(-mu_gp). BH over every gene in the territory (genes with no hit have p = 1). The
  non-synonymous version uses missense + truncating SNVs and CDS bp only. This does not model trinucleotide
  context, copy number, expression / replication timing or per-gene callability; huge genes (PRKDC, LRP1B)
  are length-corrected only.

Usage:
  python3 gos_sc_subclonal_drivers.py --data <gOS data dir> --gtf gencode.v24.annotation.gtf \
      --cosmic COSMIC_v3.4_SBS96_GRCh38.txt --work <analysis dir> [--patients P1,P2] [--n-boot 500]
"""
import argparse
import datetime
import gzip
import json
import math
import os
import re
import shutil
import sys
from collections import Counter, defaultdict

import numpy as np
from scipy.optimize import nnls
from gos_sc_noisy import is_noisy
from clonality import event_clonality, METHOD as CLONALITY_METHOD

TRUNC_RE = re.compile(r"stop_gained|splice_acceptor|splice_donor|start_lost|frameshift")
MIN_CLADE_SITES = 50
APOBEC = ["SBS2", "SBS13"]
CLONAL_FRACTION = 0.85
MIN_CARRIER_VAF = 0.15  # below this the carriers' alt reads are not het-like (CN-2 het ~ 0.5): dropout model not used

PATHWAYS = {
    "DNA repair: NHEJ": ["PRKDC", "XRCC4", "XRCC5", "XRCC6", "LIG4", "DCLRE1C", "NHEJ1", "PAXX", "POLM", "POLL"],
    "DNA repair: HR": ["BRCA1", "BRCA2", "PALB2", "RAD51", "RAD51B", "RAD51C", "RAD51D", "RAD50", "MRE11", "NBN", "ATM", "ATR", "CHEK1", "CHEK2", "BARD1", "BRIP1", "RAD54L", "BLM", "WRN", "RBBP8"],
    "DNA repair: MMR": ["MLH1", "MLH3", "MSH2", "MSH3", "MSH6", "PMS1", "PMS2", "EXO1", "POLE", "POLD1"],
    "DNA repair: FA": ["FANCA", "FANCB", "FANCC", "FANCD2", "FANCE", "FANCF", "FANCG", "FANCI", "FANCL", "FANCM", "SLX4", "ERCC4", "UBE2T"],
    "Chromatin remodelling": ["ARID1A", "ARID1B", "ARID2", "ARID4A", "ARID4B", "SMARCA2", "SMARCA4", "SMARCB1", "SMARCC1", "SMARCC2", "PBRM1", "ATRX", "DAXX", "CHD2", "CHD4", "CHD7", "CHD8", "KMT2A", "KMT2B", "KMT2C", "KMT2D", "SETD2", "NSD1", "KDM5C", "KDM6A", "EP300", "CREBBP", "ASXL1", "ASXL2", "BCOR", "SPEN", "H3-3A", "H3F3A", "DNMT3A", "TET2", "IDH1", "IDH2"],
    "RTK / RAS / PI3K": ["EGFR", "PDGFRA", "MET", "FGFR1", "FGFR2", "FGFR3", "ERBB2", "NF1", "KRAS", "NRAS", "HRAS", "BRAF", "PTPN11", "PIK3CA", "PIK3CB", "PIK3R1", "PIK3R2", "PTEN", "AKT1", "AKT3", "MTOR", "TSC1", "TSC2", "SH2B3", "LZTR1", "PTPN2"],
    "p53 / RB / cell cycle": ["TP53", "MDM2", "MDM4", "CDKN2A", "CDKN2B", "CDKN2C", "CDK4", "CDK6", "RB1", "CCND1", "CCND2", "CCNE1", "PPM1D", "ATM", "CDC73"],
}


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def load(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


# ---------------------------------------------------------------- gene territory

def gene_territory(gtf, cache):
    """gene -> (cds_bp, splice_bp) from CCDS-tagged protein-coding transcripts (union per gene)."""
    if os.path.exists(cache):
        out = {}
        with open(cache) as f:
            next(f)
            for line in f:
                g, c, s, chrom, start, end = line.rstrip("\n").split("\t")
                out[g] = {"cds": int(c), "splice": int(s), "chrom": chrom, "start": int(start), "end": int(end)}
        return out
    log("building gene territory from", gtf)
    cds = defaultdict(list)  # gene -> [(chrom,s,e)]
    tx_cds = defaultdict(list)  # (gene, tx) -> [(s,e)]
    op = gzip.open if gtf.endswith(".gz") else open
    with op(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t", 8)
            if p[2] != "CDS" or 'tag "CCDS"' not in p[8] or 'gene_type "protein_coding"' not in p[8]:
                continue
            g = re.search(r'gene_name "([^"]+)"', p[8]).group(1)
            tx = re.search(r'transcript_id "([^"]+)"', p[8]).group(1)
            s, e = int(p[3]), int(p[4])
            cds[g].append((p[0], s, e))
            tx_cds[(g, tx)].append((s, e))
    splice = defaultdict(set)
    for (g, tx), iv in tx_cds.items():
        iv.sort()
        for (s1, e1), (s2, e2) in zip(iv, iv[1:]):
            splice[g].update([e1 + 1, e1 + 2, s2 - 2, s2 - 1])
    out = {}
    for g, iv in cds.items():
        chrom = iv[0][0]
        iv = sorted((s, e) for c, s, e in iv if c == chrom)
        tot, cs, ce = 0, None, None
        for s, e in iv:
            if cs is None:
                cs, ce = s, e
            elif s <= ce + 1:
                ce = max(ce, e)
            else:
                tot += ce - cs + 1
                cs, ce = s, e
        tot += ce - cs + 1
        out[g] = {"cds": tot, "splice": len(splice[g]), "chrom": chrom, "start": iv[0][0], "end": max(e for s, e in iv)}
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    with open(cache, "w") as f:
        f.write("gene\tcds_bp\tsplice_bp\tchrom\tstart\tend\n")
        for g in sorted(out):
            o = out[g]
            f.write(f"{g}\t{o['cds']}\t{o['splice']}\t{o['chrom']}\t{o['start']}\t{o['end']}\n")
    return out


# ---------------------------------------------------------------- tree

def parse_newick_clades(text):
    """List of leaf-name sets, one per node (internal nodes and leaves)."""
    text = text.strip().rstrip(";")
    stack, clades, cur, i = [[]], [], None, 0
    leaves_all = []
    while i < len(text):
        ch = text[i]
        if ch == "(":
            stack.append([])
            i += 1
        elif ch in ",)":
            if ch == ")":
                kids = stack.pop()
                s = set().union(*kids) if kids else set()
                clades.append(s)
                stack[-1].append(s)
            i += 1
            if ch == ")":  # skip label / length after ')' (not after ',', which would swallow the next leaf name)
                m = re.match(r"[^(),;]*", text[i:])
                i += len(m.group(0))
        else:
            m = re.match(r"[^(),:;]+", text[i:])
            name = m.group(0)
            i += len(name)
            m2 = re.match(r"(:[^(),;]*)?", text[i:])
            i += len(m2.group(0))
            s = {name}
            leaves_all.append(name)
            clades.append(s)
            stack[-1].append(s)
    return clades, leaves_all


def best_clade_f1(carriers, clades):
    carriers = set(carriers)
    if not carriers:
        return None
    best = (-1, None)
    for s in clades:
        k = len(s & carriers)
        if not k:
            continue
        f1 = 2 * k / (len(s) + len(carriers))
        if f1 > best[0]:
            best = (f1, s)
    return best


# ---------------------------------------------------------------- stats

def poisson_sf(k, mu):
    """P(X >= k), X ~ Poisson(mu)."""
    if k <= 0:
        return 1.0
    if mu <= 0:
        return 0.0
    term = math.exp(-mu)
    cdf = term
    for i in range(1, k):
        term *= mu / i
        cdf += term
    return max(0.0, 1.0 - cdf) if cdf < 1 - 1e-12 else poisson_sf_tail(k, mu)


def poisson_sf_tail(k, mu):
    # direct tail sum for tiny p
    s, term = 0.0, math.exp(-mu + k * math.log(mu) - math.lgamma(k + 1))
    i = k
    while term > 1e-300 and i < k + 500:
        s += term
        i += 1
        term *= mu / i
    return s


def poisson_binomial_sf(k, ps):
    """P(X >= k), X = sum of independent Bernoulli(ps)."""
    dist = np.zeros(len(ps) + 1)
    dist[0] = 1.0
    for p in ps:
        dist[1:] = dist[1:] * (1 - p) + dist[:-1] * p
        dist[0] *= 1 - p
    return float(min(1.0, dist[k:].sum()))


def bh(ps, m=None):
    ps = np.asarray(ps, float)
    n = len(ps)
    m = m or n
    if not n:
        return []
    order = np.argsort(ps)
    q = np.empty(n)
    prev = 1.0
    for rank in range(n - 1, -1, -1):
        i = order[rank]
        prev = min(prev, ps[i] * m / (rank + 1))
        q[i] = prev
    return [float(min(1, x)) for x in q]


# ---------------------------------------------------------------- signatures

def read_cosmic(path):
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        rows = [l.rstrip("\n").split("\t") for l in f]
    channels = [r[0] for r in rows]
    M = np.array([[float(x) for x in r[1:]] for r in rows])
    return channels, {s: M[:, j] for j, s in enumerate(head[1:])}


def counts96(contexts, chan_index):
    v = np.zeros(len(chan_index))
    for c in contexts:
        k = chan_index.get(c)
        if k is not None:
            v[k] += 1
    return v


def fit_fractions(v, S):
    if v.sum() == 0:
        return np.zeros(S.shape[1])
    a, _ = nnls(S, v)
    t = a.sum()
    return a / t if t > 0 else a


def boot_fractions(idx_list, S, n_boot, rng):
    """idx_list: channel index per site; returns (point, boot matrix n_boot x k)."""
    idx = np.asarray(idx_list)
    n = len(idx)
    K = S.shape[0]
    point = fit_fractions(np.bincount(idx, minlength=K).astype(float), S)
    B = np.empty((n_boot, S.shape[1]))
    for b in range(n_boot):
        samp = idx[rng.integers(0, n, n)]
        B[b] = fit_fractions(np.bincount(samp, minlength=K).astype(float), S)
    return point, B


def boot_p(diff_boot, point_diff):
    """Two-sided bootstrap p for a difference of fractions (percentile, +1 smoothing)."""
    n = len(diff_boot)
    lo = (np.sum(diff_boot <= 0) + 1) / (n + 1)
    hi = (np.sum(diff_boot >= 0) + 1) / (n + 1)
    return float(min(1.0, 2 * min(lo, hi)))


# ---------------------------------------------------------------- per patient

def cons_class(x):
    t = (x.get("type") or "").lower()
    if t.startswith("trunc"):
        return "trunc"
    if t.startswith("splice"):
        return "splice"
    if t.startswith("missense"):
        return "missense"
    return t or "other"


def site_key(chrom, pos, ref, alt):
    c = f"{chrom}"
    c = c if c.startswith("chr") else f"chr{c}"
    return f"{c}_{pos}_{ref}_{alt}"


def keep_event(e):
    if e.get("vartype") != "SNV":
        return False
    try:
        tier = int(float(e.get("Tier")))
    except (TypeError, ValueError):
        return False
    if tier <= 2:
        return True
    return tier == 3 and cons_class(e) == "missense" and e.get("driver_class") == "likely"


def dropout_clonality(rc, sid, normals):
    """Fraction of tumour cells carrying the site, dropout-aware (gOS siteEvidence.js clonality())."""
    if not rc:
        return None
    rid = None
    for s in rc.get("sites", []):
        if s["id"] == sid:
            rid = s
            break
    if not rid:
        return None
    v = rc["dna"].get(sid)
    if not v:
        return None
    alt, dp = v.get("alt") or [], v.get("dp") or []
    ta = td = 0
    rows = []
    for k, cell in enumerate(rc["cells"]):
        if cell in normals:
            continue
        a, d = (alt[k] or 0), (dp[k] or 0)
        ta += a
        td += d
        rows.append((a, d))
    if not td or not ta:
        return {"fraction": None, "observed": 0, "expected": None, "vaf": None}
    vaf = ta / td
    covered = [r for r in rows if r[1] > 0]
    observed = sum(1 for a, d in covered if a > 0)
    # carrier fraction f and per-carrier VAF v with f * v = pooled VAF: expected carriers showing an alt read
    # g(f) = sum_cells f * (1 - (1 - v)^d) is increasing in f; solve g(f) = observed by bisection on [VAF, 1]
    def g(f):
        v = min(0.999, vaf / f)
        return sum(f * (1 - (1 - v) ** d) for a, d in covered)
    lo_f, hi_f = max(vaf, 1e-4), 1.0
    if observed >= g(hi_f):
        f = 1.0
    elif observed <= g(lo_f):
        f = lo_f
    else:
        for _ in range(60):
            mid = (lo_f + hi_f) / 2
            lo_f, hi_f = (mid, hi_f) if g(mid) < observed else (lo_f, mid)
        f = (lo_f + hi_f) / 2
    expected = g(1.0)
    na = sum(rc["dna"][sid].get("alt", [0])[k] or 0 for k, c in enumerate(rc["cells"]) if c in normals)
    nd = sum(rc["dna"][sid].get("dp", [0])[k] or 0 for k, c in enumerate(rc["cells"]) if c in normals)
    return {
        "fraction": round(f, 3),
        "per_carrier_vaf": round(min(1.0, vaf / f), 3),
        "fraction_if_all_at_pooled_vaf": round(min(1.0, observed / expected), 3) if expected else None,
        "observed": observed,
        "expected": round(expected, 1),
        "vaf": round(vaf, 3),
        "tumour_alt": ta,
        "tumour_dp": td,
        "normal_alt": na,
        "normal_dp": nd,
        "covered_tumour": len(covered),
    }


def patient_block(data, P, cosmic, chan_index, n_boot, rng, territory):
    pd_ = os.path.join(data, P)
    snv = load(os.path.join(pd_, "snv_matrix.json"))
    events = load(os.path.join(pd_, "filtered.events.json")) or []
    rc = load(os.path.join(pd_, "precompute", "region_calls.json"))
    clades_doc = load(os.path.join(pd_, "clades.json"))
    sj = load(os.path.join(pd_, "signatures_joint.json"))
    tree_txt = open(os.path.join(pd_, "tree.nwk")).read() if os.path.exists(os.path.join(pd_, "tree.nwk")) else None
    out = {"patient": P, "available": bool(snv), "reason": None if snv else "no SNV mutation map (snv_matrix.json) for this patient"}
    if not snv:
        return out, [], [], None
    clades, leaves = parse_newick_clades(tree_txt) if tree_txt else ([], [])
    normals = set()
    if rc:
        normals = {c for c, k in zip(rc["cells"], rc["clone"]) if re.match(r"(?i)^normal$", f"{k}")}
    tumour_leaves = [l for l in leaves if l not in normals]
    n_tumour = len(tumour_leaves) or None
    full = set(leaves)
    root_tumour = set(tumour_leaves)
    clade_by_set = {}
    if clades_doc:
        for c in clades_doc.get("clades", []):
            clade_by_set[frozenset(c.get("cells") or [])] = c
    variants = snv["variants"]
    vby = {v["id"]: v for v in variants}

    # signature activities per tree set (joint fit)
    sets = {s["name"]: s for s in (sj or {}).get("sets", [])}

    def posterior(context, category):
        if not sj or context not in chan_index:
            return None
        setname = category if category in ("truncal", "subclonal", "private") and sets.get(category) else "all"
        acts = [(a["signature"], float(a.get("activity") or 0)) for a in (sets.get(setname) or {}).get("activities", [])]
        acts = [(s, a) for s, a in acts if a > 0 and s in cosmic]
        if not acts:
            return None
        k = chan_index[context]
        w = np.array([a * cosmic[s][k] for s, a in acts])
        if w.sum() <= 0:
            return None
        w = w / w.sum()
        ranked = sorted(zip([s for s, _ in acts], w), key=lambda x: -x[1])
        return {"set": setname, "ranked": [{"signature": s, "prob": round(float(p), 3)} for s, p in ranked if p >= 0.02][:4]}

    hits = []
    for e in events:
        if not keep_event(e):
            continue
        chrom, pos = e.get("seqnames"), e.get("start")
        m = re.search(r"([ACGT])>([ACGT])", e.get("Variant_g") or "")
        if not m:
            continue
        ref, alt = m.group(1), m.group(2)
        sid = site_key(chrom, pos, ref, alt)
        rcid = f"{sid.split('_')[0]}:{pos}:{ref}>{alt}"
        v = vby.get(sid)
        carriers = [c for c in f"{e.get('cell_ids') or ''}".split(",") if c]
        f1 = best_clade_f1(carriers, clades) if clades else None
        clade = None
        if f1 and f1[1] is not None:
            c = clade_by_set.get(frozenset(f1[1]))
            clade = {"size": len(f1[1]), "clade_id": c["id"] if c else None, "is_root": f1[1] >= root_tumour}
        clon = dropout_clonality(rc, rcid, normals)
        raw_frac = float(e["cell_fraction"]) if e.get("cell_fraction") not in (None, "", "None") else None
        low_vaf = bool(clon and clon.get("per_carrier_vaf") is not None and clon["per_carrier_vaf"] < MIN_CARRIER_VAF)
        frac = clon["fraction"] if clon and clon.get("fraction") is not None and not low_vaf else raw_frac
        frac_source = "raw (per-carrier VAF too low for the dropout model)" if low_vaf else ("dropout-aware" if clon and clon.get("fraction") is not None else "raw")
        cat = v["category"] if v else None
        tree_sub = cat in ("subclonal", "private") if cat else None
        ec = event_clonality(e)   # read-evidence rows: exported dropout-aware call (clonality.py), not the fraction threshold
        if ec is not None:
            frac_source = CLONALITY_METHOD + " detection test (clonality.py); carrier_fraction is " + frac_source
            if ec == "indeterminate":
                status = "indeterminate"
            elif (ec == "subclonal") == bool(tree_sub) or tree_sub is None:
                status = ec
            else:
                status = "discordant"
        elif frac is None:
            status = "unknown"
        elif frac < CLONAL_FRACTION and (tree_sub is None or tree_sub):
            status = "subclonal"
        elif frac >= CLONAL_FRACTION and (tree_sub is None or not tree_sub):
            status = "clonal"
        else:
            status = "discordant"
        hits.append({
            "patient": P,
            "gene": e.get("gene"),
            "site": rcid,
            "snv_id": sid,
            "chrom": f"chr{chrom}" if not f"{chrom}".startswith("chr") else chrom,
            "pos": int(pos),
            "ref": ref,
            "alt": alt,
            "protein": e.get("Variant") or "",
            "class": cons_class(e),
            "tier": int(float(e.get("Tier"))),
            "driver_class": e.get("driver_class"),
            "oncokb_effect": e.get("effect"),
            "n_carriers": int(float(e.get("n_cells") or 0)),
            "n_normal_carriers": int(float(e.get("n_normal_cells") or 0)) if e.get("n_normal_cells") not in (None, "None", "") else None,
            "n_tumour_cells": n_tumour,
            "carrier_ids": [c for c in carriers if c not in normals],
            "carrier_fraction_raw": raw_frac,
            "clonality": clon,
            "carrier_fraction": frac,
            "carrier_fraction_source": frac_source,
            "low_carrier_vaf": low_vaf,
            "tree": {"category": cat, "node": v.get("node"), "clade_cells": v.get("clade_cells"), "alt_cells": v.get("alt_cells"), "cellphy_input": v.get("cellphy_input"), "clade_score": v.get("clade_score")} if v else None,
            "clade_f1": round(f1[0], 3) if f1 else None,
            "best_clade": clade,
            "status": status,
            "clonality_test": e.get("clonality_stats") if ec is not None else None,
            "context": v.get("context") if v else None,
            "signature": posterior(v.get("context"), cat) if v else None,
        })
    # multi-hit per gene per patient
    cnt = Counter(h["gene"] for h in hits)
    for h in hits:
        h["hits_in_gene_patient"] = cnt[h["gene"]]

    # background: non-truncal truncating / missense / synonymous SNVs of the mutation map
    bg = Counter()
    gene_obs = defaultdict(lambda: {"trunc": [], "nonsyn": []})
    for v in variants:
        cs = v.get("consequence") or ""
        if v.get("category") not in ("subclonal", "private"):
            continue
        g = v.get("gene") or ""
        if TRUNC_RE.search(cs):
            bg["trunc"] += 1
            gene_obs[g]["trunc"].append(v["id"])
            gene_obs[g]["nonsyn"].append(v["id"])
        elif "missense" in cs:
            bg["missense"] += 1
            gene_obs[g]["nonsyn"].append(v["id"])
        elif "synonymous" in cs:
            bg["synonymous"] += 1
    truncal_bg = Counter()
    for v in variants:
        cs = v.get("consequence") or ""
        if v.get("category") != "truncal":
            continue
        truncal_bg["trunc" if TRUNC_RE.search(cs) else "missense" if "missense" in cs else "synonymous" if "synonymous" in cs else "other"] += 1
    out.update({
        "n_tumour_cells": n_tumour,
        "n_normal_cells": len(normals),
        "n_sites": len(variants),
        "site_categories": dict(Counter(v.get("category") for v in variants)),
        "background": {"subclonal": dict(bg), "truncal": {k: truncal_bg[k] for k in ("trunc", "missense", "synonymous")}},
        "noisy": is_noisy(P),
    })

    # ---- clade signatures
    sig_block = None
    if sj and clades_doc:
        base = [s["signature"] for s in sj.get("signature_set", []) if s["signature"] in cosmic]
        sigs = base + [s for s in APOBEC if s not in base]
        S = np.column_stack([cosmic[s] for s in sigs])

        def idx_of(ids):
            return [chan_index[vby[i]["context"]] for i in ids if i in vby and vby[i].get("context") in chan_index]

        truncal_ids = [v["id"] for v in variants if v.get("category") == "truncal" and v.get("cellphy_input")]
        tr_idx = idx_of(truncal_ids)
        tr_point, tr_boot = boot_fractions(tr_idx, S, n_boot, rng) if len(tr_idx) >= MIN_CLADE_SITES else (None, None)
        cl = {c["id"]: c for c in clades_doc.get("clades", [])}
        fits = {}
        rows = []
        for c in clades_doc.get("clades", []):
            ids = c.get("private_snvs") or []
            ix = idx_of(ids)
            r = {"clade": c["id"], "parent": c.get("parent"), "n_cells": c.get("n_cells"), "depth": c.get("depth"), "is_trunk": bool(c.get("is_trunk")) or c.get("parent") is None, "clone_id": c.get("clone_id"), "n_sites": len(ix), "private_drivers": (c.get("private_drivers") or "").strip()}
            if len(ix) < MIN_CLADE_SITES:
                r["fitted"] = False
                r["reason"] = f"n sites {len(ix)} < {MIN_CLADE_SITES}"
            else:
                pt, B = boot_fractions(ix, S, n_boot, rng)
                fits[c["id"]] = (pt, B)
                r["fitted"] = True
                r["fractions"] = [{"signature": s, "frac": round(float(pt[j]), 4), "lo": round(float(np.quantile(B[:, j], 0.025)), 4), "hi": round(float(np.quantile(B[:, j], 0.975)), 4)} for j, s in enumerate(sigs)]
            rows.append(r)
        # comparisons: vs truncal set, vs nearest fitted ancestor branch set
        tests = []
        for r in rows:
            if not r.get("fitted") or r["is_trunk"]:
                continue
            pt, B = fits[r["clade"]]
            comps = []
            if tr_point is not None:
                comps.append(("truncal", tr_point, tr_boot))
            par = r["parent"]
            while par and par in cl and par not in fits:
                par = cl[par].get("parent")
            if par and par in fits and not (cl[par].get("is_trunk") or cl[par].get("parent") is None):
                comps.append((par, fits[par][0], fits[par][1]))
            r["vs"] = []
            for name, p2, B2 in comps:
                for j, s in enumerate(sigs):
                    d = float(pt[j] - p2[j])
                    db = B[:, j] - B2[:, j]
                    t = {"clade": r["clade"], "vs": name, "signature": s, "diff": round(d, 4), "lo": round(float(np.quantile(db, 0.025)), 4), "hi": round(float(np.quantile(db, 0.975)), 4), "p": boot_p(db, d), "frac": round(float(pt[j]), 4), "frac_ref": round(float(p2[j]), 4)}
                    r["vs"].append(t)
                    tests.append(t)
        sig_block = {
            "signatures": sigs,
            "added": [s for s in APOBEC if s not in base],
            "truncal": {"n_sites": len(tr_idx), "fractions": [{"signature": s, "frac": round(float(tr_point[j]), 4), "lo": round(float(np.quantile(tr_boot[:, j], 0.025)), 4), "hi": round(float(np.quantile(tr_boot[:, j], 0.975)), 4)} for j, s in enumerate(sigs)] if tr_point is not None else None},
            "clades": rows,
            "_tests": tests,
        }
    contrasts = []
    if sj:
        for c in sj.get("contrasts", []):
            if c.get("kind") == "tree" or c["contrast"].startswith("subclonal") or c["contrast"].startswith("private"):
                contrasts.append({"contrast": c["contrast"], "n_a": c.get("n_a"), "n_b": c.get("n_b"), "rows": c.get("rows")})
        out["signature_sets"] = [{"name": s["name"], "n": s.get("n"), "group": s.get("group"), "activities": [{k: a.get(k) for k in ("signature", "fraction", "lo", "hi")} for a in s.get("activities", [])]} for s in sj.get("sets", []) if s.get("group") == "tree"]
    out["signature_contrasts"] = contrasts
    return out, hits, gene_obs, sig_block


CARRIER_VAF_MIN_RATIO = 0.5     # carrier pooled VAF below this fraction of the patient's real-subclonal baseline ...
CARRIER_VAF_MAX_P = 1e-3         # ... and binomially unlikely under that baseline -> flagged
SINGLE_READ_MAX_FRAC = 0.3       # fewer than this share of carriers with >= 2 alt reads (and >= 5 carriers) -> flagged


def binom_cdf(k, n, p):
    """P(X <= k), X ~ Binomial(n, p), summed in log space (n can be thousands of reads)."""
    from math import exp, lgamma, log
    if k >= n or p <= 0:
        return 1.0
    if p >= 1:
        return 0.0
    terms = [lgamma(n + 1) - lgamma(i + 1) - lgamma(n - i + 1) + i * log(p) + (n - i) * log(1 - p) for i in range(0, k + 1)]
    m = max(terms)
    return min(1.0, exp(m) * sum(exp(t - m) for t in terms))


def carrier_vaf_flags(data, hits, min_dp=4):
    """Adds h["carrier_vaf_check"] to each hit from snv_matrix.json per-cell [variant, depth, alt, call]:
    carrier alt / depth, share of carriers with >= 2 alt reads, the patient's baseline (pooled carrier VAF of subclonal tree
    sites with clade_score >= 0.8 and >= 3 alt cells), a binomial P of the carrier reads under that baseline, and flags."""
    by_p = {}
    for h in hits:
        by_p.setdefault(h["patient"], []).append(h)
    for P, hs in by_p.items():
        m = load(os.path.join(data, P, "snv_matrix.json"))
        if not m:
            continue
        variants, cells = m["variants"], m["cells"]
        vid = {v["id"]: i for i, v in enumerate(variants)}
        good = {i for i, v in enumerate(variants) if v.get("category") == "subclonal" and (v.get("clade_score") or 0) >= 0.8 and (v.get("alt_cells") or 0) >= 3}
        want = {vid[h["snv_id"]] for h in hs if h.get("snv_id") in vid}
        b_alt = b_dp = 0
        per = {i: [] for i in want}
        for recs in cells.values():
            for r in recs:
                i, dp, alt = r[0], r[1] or 0, r[2] or 0
                if alt < 1 or dp < min_dp:
                    continue
                if i in good:
                    b_alt += alt
                    b_dp += dp
                if i in per:
                    per[i].append((alt, dp))
        base = b_alt / b_dp if b_dp else None
        for h in hs:
            i = vid.get(h.get("snv_id"))
            cs = per.get(i) or []
            if not cs or base is None:
                h["carrier_vaf_check"] = {"available": False, "reason": "site not in snv_matrix or no baseline"}
                continue
            alt = sum(c[0] for c in cs)
            dp = sum(c[1] for c in cs)
            multi = sum(c[0] >= 2 for c in cs) / len(cs)
            pv = binom_cdf(alt, dp, base)
            flags = []
            if alt / dp < CARRIER_VAF_MIN_RATIO * base and pv < CARRIER_VAF_MAX_P:
                flags.append(f"carrier VAF {alt / dp:.2f} vs {base:.2f} for well-placed subclonal sites in {P} (copy number 2 expects ~0.5): likely amplification / sequencing error")
            if len(cs) >= 5 and multi < SINGLE_READ_MAX_FRAC:
                flags.append(f"single-read carriers: {sum(c[0] >= 2 for c in cs)} of {len(cs)} carriers have >= 2 alt reads")
            h["carrier_vaf_check"] = {"available": True, "carriers": len(cs), "alt": alt, "depth": dp, "vaf": round(alt / dp, 4),
                                      "baseline_vaf": round(base, 4), "p_below_baseline": pv, "carriers_2plus_alt_frac": round(multi, 3), "flags": flags}
    return {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--cosmic", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--patients", default="BWH69,BWH70,BWH75,MGH285,MGH301,MGH302,MGH303,MGH304,MGH305")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--extra", default=None, help="JSON with extra per-hit notes {site: {title, text, loci: [...]}}")
    ap.add_argument("--read-qc", default=None, help="JSON from subclonal/hit_read_qc.py {site: {...}}")
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    terr = gene_territory(a.gtf, os.path.join(a.work, "gene_territory_gencode_v24_ccds.tsv"))
    channels, cosmic = read_cosmic(a.cosmic)
    chan_index = {c: k for k, c in enumerate(channels)}
    patients, hits, sigblocks = [], [], {}
    obs = {}
    for P in a.patients.split(","):
        log("patient", P)
        blk, h, go, sb = patient_block(a.data, P, cosmic, chan_index, a.n_boot, rng, terr)
        patients.append(blk)
        hits.extend(h)
        if blk["available"]:
            obs[P] = go
        if sb:
            sigblocks[P] = sb

    # ---- read-level QC (hit_read_qc.py): artefact flags per hit, used as caveats and for a sensitivity test
    rq = load(a.read_qc) if a.read_qc else {}
    flagged_ids = set()
    for h in hits:
        r = (rq or {}).get(f"{h['patient']}|{h['site']}")
        if not r:
            continue
        c = r.get("counts") or {}
        na, nr = c.get("alt_reads", 0), c.get("ref_reads", 0)
        f = lambda k, n: (c.get(k, 0) / n) if n else None
        flags = []
        if na >= 3 and (f("alt_exon_split", na) or 0) >= 0.5:
            flags.append("spliced reads: alt sits on exon-exon reads (RNA / cDNA carry-over or processed copy)")
        elif na >= 3 and (f("alt_sa", na) or 0) >= 0.5 and (f("ref_sa", nr) or 0) < 0.25:
            flags.append("alt on split reads only")
        if na >= 3 and (f("alt_end10", na) or 0) >= 0.8 and (f("ref_end10", nr) or 0) <= 0.3:
            flags.append("alt at read-block ends")
        if r.get("other_patients_with_allele"):
            flags.append("same allele in " + ", ".join(r["other_patients_with_allele"]))
        h["read_qc"] = {
            "carriers_scanned": r.get("carriers_scanned"), "alt_reads": na, "ref_reads": nr,
            "alt_sa_frac": f("alt_sa", na), "ref_sa_frac": f("ref_sa", nr), "alt_exon_split_frac": f("alt_exon_split", na),
            "alt_end10_frac": f("alt_end10", na), "ref_end10_frac": f("ref_end10", nr),
            "other_patients_with_allele": r.get("other_patients_with_allele"), "flags": flags,
        }
        if flags:
            flagged_ids.add(h["snv_id"])

    # ---- carrier VAF vs real mutations (2026-10-16): in one cell a real het SNV at CN 2 gives VAF ~0.5, and the patient's
    # well-placed subclonal tree sites reach ~0.6-0.7 in their carrier cells (PTA imbalance, carriers conditioned on alt >= 1).
    # Hits whose carriers sit far below that, mostly with a single alt read, look like amplification / sequencing errors.
    carrier_vaf_flags(a.data, hits)
    for h in hits:
        cv = h.get("carrier_vaf_check")
        if not cv or not cv.get("flags"):
            continue
        rqh = h.setdefault("read_qc", {"flags": []})
        rqh.setdefault("flags", [])
        for fl in cv["flags"]:
            if fl not in rqh["flags"]:
                rqh["flags"].append(fl)
        flagged_ids.add(h["snv_id"])

    # ---- recurrence
    L_trunc_tot = sum(t["cds"] + t["splice"] for t in terr.values())
    L_cds_tot = sum(t["cds"] for t in terr.values())
    avail = [p for p in patients if p["available"]]
    genes_out = []
    for variant_set, key in (("trunc", "trunc"), ("nonsyn", None)):
        rows = []
        for g, t in terr.items():
            L = t["cds"] + t["splice"] if variant_set == "trunc" else t["cds"]
            Ltot = L_trunc_tot if variant_set == "trunc" else L_cds_tot
            per = []
            for p in avail:
                n_p = p["background"]["subclonal"].get("trunc", 0) if variant_set == "trunc" else p["background"]["subclonal"].get("trunc", 0) + p["background"]["subclonal"].get("missense", 0)
                mu = n_p * L / Ltot
                ids = obs[p["patient"]].get(g, {}).get(variant_set, [])
                per.append((p["patient"], mu, len(ids), ids, p.get("noisy"), len([i for i in ids if i not in flagged_ids])))
            hits_n = sum(x[2] for x in per)
            pats_n = sum(1 for x in per if x[2] > 0)
            mu_tot = sum(x[1] for x in per)
            if hits_n == 0:
                rows.append({"gene": g, "p_hits": 1.0, "p_patients": 1.0, "p_patients_no_bwh69": 1.0, "p_patients_read_qc_pass": 1.0, "_zero": True})
                continue
            p_h = poisson_sf(hits_n, mu_tot)
            p_p = poisson_binomial_sf(pats_n, [1 - math.exp(-x[1]) for x in per])
            per_nb = [x for x in per if not x[4]]
            pats_nb = sum(1 for x in per_nb if x[2] > 0)
            p_nb = poisson_binomial_sf(pats_nb, [1 - math.exp(-x[1]) for x in per_nb]) if pats_nb else 1.0
            pats_qc = sum(1 for x in per if x[5] > 0)
            p_qc = poisson_binomial_sf(pats_qc, [1 - math.exp(-x[1]) for x in per]) if pats_qc else 1.0
            rows.append({
                "gene": g, "set": variant_set, "length_bp": L, "cds_bp": t["cds"], "splice_bp": t["splice"],
                "hits": hits_n, "patients": pats_n, "expected_hits": round(mu_tot, 4),
                "ratio": round(hits_n / mu_tot, 1) if mu_tot else None,
                "p_hits": p_h, "p_patients": p_p, "patients_no_bwh69": pats_nb, "p_patients_no_bwh69": p_nb,
                "patients_read_qc_pass": pats_qc, "p_patients_read_qc_pass": p_qc,
                "per_patient": [{"patient": x[0], "expected": round(x[1], 4), "hits": x[2], "sites": x[3], "read_qc_flagged": [i for i in x[3] if i in flagged_ids]} for x in per if x[2] > 0],
            })
        m = len(rows)
        for k in ("p_hits", "p_patients", "p_patients_no_bwh69", "p_patients_read_qc_pass"):
            qs = bh([r[k] for r in rows])
            for r, q in zip(rows, qs):
                r[k.replace("p_", "q_")] = q
        keep = sorted([r for r in rows if not r.get("_zero")], key=lambda r: (r["p_patients"], r["p_hits"]))
        genes_out.append({"set": variant_set, "n_genes_tested": m, "rows": keep})

    # driver-hit gene table (per gene: patients with a subclonal hit, independent hits)
    by_gene = defaultdict(list)
    for h in hits:
        by_gene[h["gene"]].append(h)
    rec_t = {r["gene"]: r for r in genes_out[0]["rows"]}
    rec_n = {r["gene"]: r for r in genes_out[1]["rows"]}
    driver_genes = []
    for g, hs in by_gene.items():
        sub = [h for h in hs if h["status"] in ("subclonal", "discordant")]
        driver_genes.append({
            "gene": g,
            "patients": sorted({h["patient"] for h in hs}),
            "patients_subclonal": sorted({h["patient"] for h in sub}),
            "n_hits": len({(h["patient"], h["site"]) for h in hs}),
            "n_hits_subclonal": len({(h["patient"], h["site"]) for h in sub}),
            "multi_hit_patients": sorted({h["patient"] for h in hs if h["hits_in_gene_patient"] > 1}),
            "classes": dict(Counter(h["class"] for h in hs)),
            "length_bp": (terr.get(g) or {}).get("cds", 0) + (terr.get(g) or {}).get("splice", 0) or None,
            "recurrence_trunc": {k: rec_t[g].get(k) for k in ("hits", "patients", "expected_hits", "p_patients", "q_patients", "p_hits", "q_hits", "patients_no_bwh69", "q_patients_no_bwh69", "patients_read_qc_pass", "q_patients_read_qc_pass")} if g in rec_t else None,
            "recurrence_nonsyn": {k: rec_n[g].get(k) for k in ("hits", "patients", "expected_hits", "p_patients", "q_patients", "p_hits", "q_hits")} if g in rec_n else None,
        })
    driver_genes.sort(key=lambda r: (-len(r["patients_subclonal"]), -r["n_hits_subclonal"], r["gene"]))

    # pathways
    pw = []
    for name, gl in PATHWAYS.items():
        gs = set(gl)
        cells = {}
        for p in patients:
            P = p["patient"]
            if not p["available"]:
                cells[P] = {"available": False}
                continue
            hh = [h for h in hits if h["patient"] == P and h["gene"] in gs]
            other = []
            for g in gs:
                for sid in obs[P].get(g, {}).get("trunc", []):
                    if not any(h["snv_id"] == sid for h in hh):
                        other.append({"gene": g, "snv_id": sid})
            cells[P] = {
                "available": True,
                "subclonal": sorted({f"{h['gene']} {h['protein'] or h['class']}" for h in hh if h["status"] in ("subclonal", "discordant")}),
                "clonal": sorted({f"{h['gene']} {h['protein'] or h['class']}" for h in hh if h["status"] == "clonal"}),
                "other_truncating": other,
            }
        pw.append({"pathway": name, "genes": gl, "patients": cells})

    # clade signature BH across the cohort
    all_tests = [t for sb in sigblocks.values() for t in sb["_tests"]]
    qs = bh([t["p"] for t in all_tests])
    for t, q in zip(all_tests, qs):
        t["q"] = q
    sig_shifts = []
    for P, sb in sigblocks.items():
        for t in sb.pop("_tests"):
            if t["q"] < 0.1 and abs(t["diff"]) >= 0.05:
                sig_shifts.append({"patient": P, **t})
    sig_shifts.sort(key=lambda t: (t["q"], -abs(t["diff"])))

    extra = load(a.extra) if a.extra else {}
    for h in hits:
        if h["site"] in (extra or {}):
            h["note"] = extra[h["site"]]

    doc = {
        "format": "gos-sc-subclonal-drivers/1",
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "generator": "srctools inst/gos_single_cell/gos_sc_subclonal_drivers.py (feature/subclonal-drivers)",
        "definitions": {
            "hit": "filtered.events.json SNV with Tier 1-2, or tier-3 missense with driver_class 'likely' (damaging)",
            "carrier_fraction": "dropout-aware fraction f of tumour cells carrying the site, from precompute/region_calls.json: f and the per-carrier VAF v satisfy f x v = pooled tumour VAF and sum over covered tumour cells of f x (1 - (1 - v)^depth) = observed carrier cells (bisection); raw carrier cells / tumour cells when the site has no region calls or the per-carrier VAF is < 0.15 (not het-like: the alt reads look like low-level error, so a dropout correction would be meaningless; flagged low_carrier_vaf)",
            "subclonal": f"carrier fraction < {CLONAL_FRACTION} and the CellPhy mutation map places the site on a subclonal or private node; 'clonal' = fraction >= {CLONAL_FRACTION} and a truncal node; 'discordant' = the two disagree. Read-evidence rows (hotspot rescue, pooled scan, called-SNV rows) use the exported dropout-aware call instead "
                         "(clonality.py: detection vs the patient's truncal SNVs at matched depth, Poisson-binomial alpha 0.01, plus clade-restriction permutation alpha 0.05; "
                         "'indeterminate' = < 10 covered cells, no reference, or low detection without clade restriction)",
            "clade_f1": "best F1 between the carrier cells and any clade of tree.nwk (normals included, same as gOS cladeFit.js)",
            "signature_posterior": "P(signature | SBS96 context) proportional to the exposure of the site's tree set (truncal / subclonal / private) in signatures_joint.json x the COSMIC v3.4 profile",
        },
        "recurrence_method": (
            "Background per patient = non-truncal (subclonal + private) truncating SNVs (stop_gained, splice donor/acceptor, start_lost) in the patient's CellPhy mutation map; "
            "expected hits in a gene = sum over patients of background x gene territory bp / total territory bp (GENCODE v24 CCDS CDS union + 4 bp per coding intron; "
            f"{L_trunc_tot / 1e6:.1f} Mb over {len(terr)} genes). p_hits: Poisson on distinct hits; p_patients: Poisson-binomial on patients with >= 1 hit; BH q over all genes in the territory. "
            "'nonsyn' adds missense SNVs to both background and hits. Not modelled: trinucleotide context, copy number, expression / replication timing, per-gene callability, "
            "and SNV noise (patients flagged noisy via NOISY_PATIENTS, default none, are dropped in p_patients_no_bwh69, a legacy key name equal to p_patients when none is flagged; read-level artefact flags from hit_read_qc.py: see p/q_patients_read_qc_pass, which drops flagged driver sites). Hits come from the mutation map, so they include non-OncoKB genes; driver SNVs placed on a truncal node are not counted."
        ),
        "dnds_note": "Gene-level dN/dS is not fitted: subclonal truncating counts are 3-31 per patient. Cohort trunc:missense:synonymous counts per patient are in patients[].background.",
        "clade_signature_method": (
            f"Each called clade's branch SNVs (clades.json private_snvs with an SBS96 context) refitted by NNLS on the patient's signature set from signatures_joint.json plus {', '.join(APOBEC)} "
            f"(so an APOBEC gain can show); bootstrap over sites (n = {a.n_boot}) for 95% CIs; clades with < {MIN_CLADE_SITES} sites are not fitted. Each clade is compared with the truncal set (truncal CellPhy-input sites) "
            "and with its nearest fitted non-trunk ancestor; two-sided bootstrap p for the fraction difference; BH q across every clade x signature x comparison in the cohort. "
            "Shifts listed: q < 0.1 and |difference| >= 5 points."
        ),
        "n_boot": a.n_boot,
        "patients": patients,
        "hits": hits,
        "driver_genes": driver_genes,
        "recurrence": genes_out,
        "pathways": pw,
        "clade_signatures": sigblocks,
        "signature_shifts": sig_shifts,
    }
    out = os.path.join(a.data, "_cohort", "subclonal_drivers.json")
    if os.path.exists(out):
        shutil.copy(out, out + ".bak_" + datetime.date.today().strftime("%y%m%d"))
    tmp = out + ".tmp"
    with open(tmp, "w") as f:
        json.dump(doc, f, separators=(",", ":"), default=lambda o: float(o) if isinstance(o, np.floating) else int(o) if isinstance(o, np.integer) else str(o))
    os.replace(tmp, out)
    log("wrote", out, len(hits), "hits;", len(flagged_ids), "read-QC flagged")
    top = genes_out[0]["rows"][:12]
    for r in top:
        log(f"  {r['gene']:10s} hits={r['hits']} pats={r['patients']} exp={r['expected_hits']:.3f} p={r['p_patients']:.2e} q={r['q_patients']:.2e} q_noBWH69={r['q_patients_no_bwh69']:.2e} pats_qc={r['patients_read_qc_pass']} q_qc={r['q_patients_read_qc_pass']:.2e}")
    for s in sig_shifts[:15]:
        log("  shift", s["patient"], s["clade"], "vs", s["vs"], s["signature"], s["diff"], f"q={s['q']:.3g}")


if __name__ == "__main__":
    main()
