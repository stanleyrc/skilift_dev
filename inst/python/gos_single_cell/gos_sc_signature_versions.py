#!/usr/bin/env python3
"""gos_sc_signature_versions: SBS96 signature fits of one patient under several SNV-set definitions ("versions"),
per-SNV signature posteriors, version-to-version assignment flows and artifact diagnostics -> data/<P>/signature_versions.json.

Every version uses the SAME catalogue-choice rule (choose_catalogue):
  greedy forward selection from COSMIC (+ forced columns) on the version's 96-channel catalogue -- each step adds the
  signature with the largest multinomial log-likelihood gain if its likelihood-ratio test passes alpha / #candidates,
  it takes >= min_frac of the exposure and it raises the reconstruction cosine by >= min_cosine_gain (the first, best signature is always taken); then only signatures selected in >= min_stability of n_boot_select bootstrap
  resamples of the catalogue are kept. Kept signatures are refitted by NNLS with n_boot mutation-resampling bootstraps
  (exposure fractions with percentile 95% CIs, reconstruction cosine with CI), as in signatures_joint.
Fit modes: "joint" = one catalogue for the whole version; "per_set" = the rule applied to each set separately (legacy).

Per-SNV posterior: P(signature k | SNV in channel c) = a_k S[c,k] / sum_j a_j S[c,j] with a = the version's (or set's)
NNLS exposures; the assignment is the argmax. Because the posterior depends only on the channel, it is exported per
channel (96 x K) and per tree site as a compact code array; singletons are aggregated by channel.

Inputs (all real; nothing is imputed -- missing inputs give "not available (reason)" entries):
  data/<P>/snv_matrix.sites.json   tree sites (category, clade_score, alt_cells, map_confidence, context, cellphy_input)
  data/<P>/snv_qc.json             site QC (basic pass flag), optional
  --vcf-sites TSV.gz               per-site read evidence from the VCF (sv_parse_vcf.py in srctools): singletons + tree sites
  --cosmic TSV                     COSMIC SBS96 reference (MutationType + one column per signature)
  --artifact NAME=TSV              extra artifact signature columns (e.g. the SCAN2 PTA artifact signature)
  --fasta                          reference FASTA (+ .fai) for singleton trinucleotide contexts
Synthetic-fixture tests: tests/test_signature_versions.py.
"""
import argparse
import csv
import gzip
import json
import math
import os
import sys
from datetime import datetime

import numpy as np
from scipy.optimize import lsq_linear
from scipy.optimize import nnls as _nnls


def nnls(A, b):
    """scipy nnls with an empty-matrix guard and a bounded least-squares fallback when it does not converge."""
    if A.shape[1] == 0:
        return np.zeros(0), float(np.linalg.norm(b))
    try:
        return _nnls(A, b, maxiter=100 * A.shape[1])
    except RuntimeError:
        r = lsq_linear(A, b, bounds=(0, np.inf))
        return r.x, float(np.linalg.norm(A @ r.x - b))

BASES = "ACGT"
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
SBS96 = [f"{l}[{r}>{a}]{rr}" for r, alts in (("C", "AGT"), ("T", "ACG")) for a in alts for l in BASES for rr in BASES]
SBS96 = sorted(SBS96, key=lambda s: (s[2:5], s[0], s[6]))  # COSMIC order: C>A block, 5' base, 3' base
CH_INDEX = {c: i for i, c in enumerate(SBS96)}
FLAT_SIGS = ("SBS3", "SBS5", "SBS40a", "SBS40b", "SBS40c", "SBS40")

# ---------------------------------------------------------------- reference / contexts

def read_signature_table(path):
    """(names, 96 x K matrix in SBS96 order, columns sum to 1)."""
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    names = rows[0][1:]
    m = {r[0]: [float(x) for x in r[1:]] for r in rows[1:]}
    missing = [c for c in SBS96 if c not in m]
    if missing:
        raise ValueError(f"{path}: missing SBS96 channels {missing[:3]}")
    S = np.array([m[c] for c in SBS96], dtype=float)
    S = S / S.sum(axis=0, keepdims=True)
    return names, S


class Fasta:
    """Minimal .fai-indexed FASTA reader (no pysam)."""

    def __init__(self, path):
        self.fh = open(path, "rb")
        self.idx = {}
        with open(path + ".fai") as fh:
            for ln in fh:
                name, length, off, lb, lbytes = ln.split("\t")[:5]
                self.idx[name] = (int(length), int(off), int(lb), int(lbytes))

    def fetch(self, chrom, start, end):  # 1-based inclusive
        length, off, lb, lbytes = self.idx[chrom]
        start, end = max(1, start), min(length, end)
        s0, e0 = start - 1, end
        bs = off + (s0 // lb) * lbytes + s0 % lb
        be = off + (e0 // lb) * lbytes + e0 % lb
        self.fh.seek(bs)
        return self.fh.read(be - bs).replace(b"\n", b"").replace(b"\r", b"").decode().upper()


def sbs96_context(tri, ref, alt):
    """COSMIC channel for a trinucleotide (5'-ref-3' on the + strand) and ref>alt, pyrimidine-normalised; None if invalid."""
    if len(tri) != 3 or tri[1] != ref or ref not in BASES or alt not in BASES or ref == alt or any(b not in BASES for b in tri):
        return None
    if ref in "GA":
        tri = "".join(COMP[b] for b in reversed(tri)); ref, alt = COMP[ref], COMP[alt]
    return f"{tri[0]}[{ref}>{alt}]{tri[2]}"


def context_of_id(fa, mut):
    chrom, pos, ref, alt = mut.rsplit("_", 3)
    pos = int(pos)
    return sbs96_context(fa.fetch(chrom, pos - 1, pos + 1), ref, alt)

# ---------------------------------------------------------------- fitting core

def counts_of(channels):
    v = np.zeros(96)
    for c in channels:
        i = CH_INDEX.get(c)
        if i is not None:
            v[i] += 1
    return v


def cosine(a, b):
    d = math.sqrt(float(np.dot(a, a)) * float(np.dot(b, b)))
    return float(np.dot(a, b)) / d if d > 0 else float("nan")


def nnls_fit(counts, S, idx):
    e = np.zeros(S.shape[1])
    if len(idx):
        e[list(idx)] = nnls(S[:, list(idx)], counts)[0]
    return e


def loglik(counts, S, e):
    rec = S @ e
    tot = rec.sum()
    if tot <= 0:
        return -math.inf
    p = np.maximum(rec / tot, 1e-12)
    return float(np.dot(counts, np.log(p)))


def forward_select(counts, S, forced=(), alpha=0.05, min_frac=0.05, max_add=6, min_cosine_gain=0.02):
    """Greedy forward selection (see module docstring). Returns list of added column indices in order."""
    from scipy.stats import chi2
    cur = list(forced)
    added = []
    for _ in range(max_add):
        cands = [j for j in range(S.shape[1]) if j not in cur]
        if not cands:
            break
        if cur:
            e0 = nnls_fit(counts, S, cur); ll0 = loglik(counts, S, e0); cos0 = cosine(counts, S @ e0)
        else:
            ll0 = float(np.dot(counts, np.log(np.full(96, 1 / 96)))); cos0 = cosine(counts, np.ones(96))
        best = None
        for j in cands:
            e1 = nnls_fit(counts, S, cur + [j])
            tot = e1.sum()
            if tot <= 0 or e1[j] / tot < min_frac:
                continue
            ll1 = loglik(counts, S, e1)
            if best is None or ll1 > best[1]:
                best = (j, ll1, cosine(counts, S @ e1))
        if best is None:
            break
        if cur:  # the first signature is always taken (every catalogue gets >= 1); later ones must pass the tests
            delta = 2 * (best[1] - ll0)
            p = chi2.sf(delta, 1) if math.isfinite(delta) else 1.0
            gain = best[2] - (cos0 if math.isfinite(cos0) else 0.0)
            if not math.isfinite(delta) or p >= alpha / len(cands) or gain < min_cosine_gain:
                break
        cur.append(best[0]); added.append(best[0])
    return added


def resample(counts, rng):
    n = int(counts.sum())
    return rng.multinomial(n, counts / n).astype(float) if n > 0 else counts.copy()


def choose_catalogue(counts, S, names, forced=(), n_boot_select=200, min_stability=0.5, min_mutations=50, seed=1, **kw):
    """The shared catalogue-choice rule. Returns dict(signatures, stability, status, reason)."""
    n = int(counts.sum())
    if n < min_mutations:
        return {"signatures": [], "stability": {}, "status": "not_available", "reason": f"too few SNVs ({n} < {min_mutations})"}
    forced = [names.index(f) for f in forced]
    sel = forward_select(counts, S, forced, **kw)
    rng = np.random.default_rng(seed)
    hits = np.zeros(S.shape[1])
    for _ in range(n_boot_select):
        for j in forward_select(resample(counts, rng), S, forced, **kw):
            hits[j] += 1
    stab = {names[j]: round(hits[j] / n_boot_select, 3) for j in sel}
    keep = [j for j in sel if hits[j] / n_boot_select >= min_stability]
    if not keep and not forced and sel:
        keep = [sel[0]]  # nothing stable: keep the first (largest-gain) signature, flagged below
    out = forced + keep
    if not out:
        return {"signatures": [], "stability": stab, "status": "not_available", "reason": "no signature passed the selection"}
    return {"signatures": [names[j] for j in out], "stability": stab, "status": "ok",
            "reason": "" if (keep or forced) and all(stab.get(names[j], 1) >= min_stability for j in keep) else
            "no bootstrap-stable signature; kept the first selected"}


def refit_bootstrap(counts, S, names, signatures, n_boot=1000, seed=2):
    idx = [names.index(s) for s in signatures]
    S1 = S[:, idx]
    e = nnls(S1, counts)[0]
    rng = np.random.default_rng(seed)
    bf = np.full((n_boot, len(idx)), np.nan)
    bc = np.full(n_boot, np.nan)
    for b in range(n_boot):
        ct = resample(counts, rng)
        eb = nnls(S1, ct)[0]
        if eb.sum() > 0:
            bf[b] = eb / eb.sum()
        bc[b] = cosine(ct, S1 @ eb)
    frac = e / e.sum() if e.sum() > 0 else np.full(len(idx), np.nan)
    acts = []
    for k, s in enumerate(signatures):
        lo, hi = np.nanpercentile(bf[:, k], [2.5, 97.5])
        acts.append({"signature": s, "activity": round(float(e[k]), 2), "fraction": round(float(frac[k]), 4),
                     "lo": round(float(lo), 4), "hi": round(float(hi), 4)})
    acts.sort(key=lambda a: -a["activity"])
    lo, hi = np.nanpercentile(bc, [2.5, 97.5])
    return {"activities": acts, "cosine": round(cosine(counts, S1 @ e), 4), "cosine_lo": round(float(lo), 4),
            "cosine_hi": round(float(hi), 4), "exposure": e, "idx": idx}


def channel_posterior(S, idx, exposure):
    """96 x K posterior of each kept signature per channel (rows sum to 1; zero rows -> uniform)."""
    w = S[:, idx] * exposure[None, :]
    tot = w.sum(axis=1, keepdims=True)
    post = np.where(tot > 0, w / np.where(tot > 0, tot, 1), 1.0 / max(1, len(idx)))
    return post


def entropy_bits(counts):
    p = counts / counts.sum() if counts.sum() > 0 else counts
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())

# ---------------------------------------------------------------- SNV groups and versions

GROUPS = {
    "truncal": "tree site mapped to the tumour MRCA",
    "shared3": "tree site on an internal tumour branch, >= 3 alt tumour cells, clade_score >= 0.8, no alt read in a normal cell",
    "subclonal_other": "other internal-branch (subclonal) tree sites",
    "private_tree": "tree site mapped to a single tip (category 'private'); still >= 2 alt-call cells, so NOT a singleton",
    "outside_tumor": "tree site mapped outside the tumour clade",
    "singleton": "callable singleton: alt reads in exactly one sample (a tree tumour cell), >= 3 alt reads, VAF >= 0.2, alt on both strands, >= 10 other tumour cells covered",
    "singleton_lowconf": "other singletons of tree tumour cells (>= 2 alt reads in one sample only, failing a callable criterion)",
    "normal_singleton": "singleton in a normal (outgroup) cell: empirical single-cell artifact (+ normal somatic) spectrum",
}

VERSIONS = [
    {"id": "v1_legacy", "label": "1 Legacy: separate fits", "mode": "per_set",
     "sets": {"truncal": ["truncal"], "subclonal": ["shared3", "subclonal_other"], "private": ["private_tree"]},
     "definition": "truncal, subclonal and private tree sites fitted separately (each its own catalogue choice), as signatures.json"},
    {"id": "v2_current", "label": "2 Current: all tree sites jointly", "mode": "joint",
     "groups": ["truncal", "shared3", "subclonal_other", "private_tree", "outside_tumor"],
     "definition": "every tree (CellPhy-input) site in one fit, as the signatures_joint 'all' set"},
    {"id": "v3_ge2cells", "label": "3 Joint, >= 2 tumour cells", "mode": "joint", "min_alt_tumor": 2,
     "groups": ["truncal", "shared3", "subclonal_other", "private_tree"],
     "definition": "tumour tree sites with alt reads in >= 2 tumour cells (VCF AD_BAM)"},
    {"id": "v4_tree_cells", "label": "4 Phylogeny cells, incl. singletons", "mode": "joint",
     "groups": ["truncal", "shared3", "subclonal_other", "private_tree", "singleton"],
     "definition": "all SNVs of the cells in the tree: tumour tree sites + callable singletons of tree tumour cells"},
    {"id": "v5_subclonal_singletons", "label": "5 Subclonal + singletons", "mode": "joint",
     "groups": ["shared3", "subclonal_other", "private_tree", "singleton"],
     "definition": "non-truncal tree sites (subclonal + private) plus callable singletons"},
    {"id": "v6_subclonal_no_singletons", "label": "6 Subclonal, no singletons", "mode": "joint",
     "groups": ["shared3", "subclonal_other", "private_tree"],
     "definition": "non-truncal tree sites (subclonal + private), no singletons"},
    {"id": "v7_ge3cells", "label": "7 >= 3 cells (paper comparator)", "mode": "joint",
     "groups": ["truncal", "shared3"],
     "definition": "truncal + shared3 (the Fig 6 site sets)"},
    {"id": "v8_qc_pass", "label": "8 Site-QC pass", "mode": "joint", "qc": True,
     "groups": ["truncal", "shared3", "subclonal_other", "private_tree"],
     "definition": "tumour tree sites passing site QC (snv_qc basic; map_confidence only where snv_qc.json is absent) with map_confidence >= 1 (0 = ambiguous mapping)"},
    {"id": "v9_artifact_cols", "label": "9 Artifact columns (PTA)", "mode": "joint", "artifact": True,
     "groups": ["truncal", "shared3", "subclonal_other", "private_tree", "singleton"],
     "definition": "version 4 sites, catalogue = COSMIC + artifact signatures forced in (absorbs PTA artifacts)"},
]


def is_callable_singleton(r):
    alt, dp = r["alt"], r["depth"]
    return (alt >= 3 and dp > 0 and alt / dp >= 0.2 and r["adf_alt"] >= 1 and r["adr_alt"] >= 1 and r["n_cov_tumor"] >= 11)


def read_vcf_sites(path):
    """Yields dict rows of sv_parse_vcf.py output with numeric fields converted."""
    num = ("n_alt_tumor", "n_alt_normal", "n_alt_other", "n_cov_tumor", "n_cov_normal", "alt", "depth", "adf_alt", "adr_alt", "adf_ref", "adr_ref")
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for k in num:
                r[k] = int(r[k]) if r[k] not in ("", None) else 0
            r["med_vaf"] = float(r["med_vaf"]) if r["med_vaf"] not in ("", None) else None
            yield r


def load_tree_sites(data_dir, P):
    sj = json.load(open(os.path.join(data_dir, P, "snv_matrix.sites.json")))
    cols = sj["columns"]

    def col(name):
        c = cols.get(name)
        if c is None:
            return [None] * sj["n_sites"]
        if c.get("type") == "dict":
            d = c["dict"]
            return [d[i] if i is not None and i >= 0 else None for i in c["codes"]]
        if c.get("type") == "bool":
            return [bool(x) if x is not None else None for x in c["values"]]
        return c["values"]

    names = ["id", "category", "clade_score", "alt_cells", "map_confidence", "context", "cellphy_input", "driver", "gene"]
    data = {n: col(n) for n in names}
    return [dict(zip(names, vals)) for vals in zip(*(data[n] for n in names))], sj.get("cells", [])


def assign_group(site, ev):
    cat = site["category"]
    if cat == "truncal":
        return "truncal"
    if cat == "outside_tumor" or cat is None or cat == "unmapped":
        return "outside_tumor"
    if cat == "private":
        return "private_tree"
    if cat == "subclonal":
        n_norm = ev["n_alt_normal"] if ev else 0
        if (site["alt_cells"] or 0) >= 3 and (site["clade_score"] or 0) >= 0.8 and n_norm == 0:
            return "shared3"
        return "subclonal_other"
    return "outside_tumor"

# ---------------------------------------------------------------- per-patient driver

def build_snvs(P, data_dir, vcf_sites_path, fasta_path):
    """List of SNV dicts: id, group, channel, is_tree, ev (VCF read evidence or None), qc_pass."""
    tree, cells = load_tree_sites(data_dir, P)
    ev = {}
    singles = []
    if vcf_sites_path and os.path.exists(vcf_sites_path):
        for r in read_vcf_sites(vcf_sites_path):
            if r["kind"] in ("tree", "tree_singleton"):
                ev[r["mutation"]] = r
            if r["kind"] == "singleton":
                singles.append(r)
    qc = None
    qf = os.path.join(data_dir, P, "snv_qc.json")
    if os.path.exists(qf):
        q = json.load(open(qf))
        qc = {k: bool(v.get("basic")) for k, v in q.get("sites", {}).items()}
    snvs = []
    for s in tree:
        if s["cellphy_input"] is False:
            continue
        e = ev.get(s["id"])
        qc_pass = (qc is None or qc.get(s["id"], False)) and (s["map_confidence"] or 0) >= 1
        snvs.append({"id": s["id"], "group": assign_group(s, e), "channel": s["context"], "is_tree": True, "ev": e,
                     "qc_pass": qc_pass, "alt_tumor": e["n_alt_tumor"] if e else None,
                     "driver": bool(s.get("driver")), "gene": s.get("gene") or ""})
    fa = Fasta(fasta_path) if (singles and fasta_path) else None
    for r in singles:
        if r["carrier_role"] == "tree_tumor":
            g = "singleton" if is_callable_singleton(r) else "singleton_lowconf"
        elif r["carrier_role"] == "normal":
            g = "normal_singleton"
        else:
            continue
        snvs.append({"id": r["mutation"], "group": g, "channel": context_of_id(fa, r["mutation"]), "is_tree": False,
                     "ev": r, "qc_pass": None, "alt_tumor": r["n_alt_tumor"], "carrier": r["carrier"]})
    return snvs, cells, {"tree_sites": len(tree), "vcf_evidence_tree_sites": len(ev), "snv_qc": qc is not None}


def version_members(v, snvs):
    gs = set(v.get("groups", [])) | {g for gl in v.get("sets", {}).values() for g in gl}
    out = []
    for i, s in enumerate(snvs):
        if s["group"] not in gs or s["channel"] not in CH_INDEX:
            continue
        if v.get("min_alt_tumor") and s["is_tree"] and (s["alt_tumor"] is None or s["alt_tumor"] < v["min_alt_tumor"]):
            continue
        if v.get("qc") and not s["qc_pass"]:
            continue
        out.append(i)
    return out


def fit_version(v, snvs, S, names, art_names, n_boot, n_boot_select, seed):
    members = version_members(v, snvs)
    sets = v["sets"] if v["mode"] == "per_set" else {"all": v["groups"]}
    forced = art_names if v.get("artifact") else []
    if v.get("artifact") and not art_names:
        return {"id": v["id"], "label": v["label"], "definition": v["definition"], "mode": v["mode"], "status": "not_available",
                "reason": "no artifact signature reference given", "n": 0, "n_singleton": 0, "groups": [], "forced": [], "sets": []}, {}
    res_sets, assign = [], {}
    for sname, gl in sets.items():
        idx = [i for i in members if snvs[i]["group"] in gl]
        ct = counts_of(snvs[i]["channel"] for i in idx)
        cat = choose_catalogue(ct, S, names, forced=forced, n_boot_select=n_boot_select, seed=seed)
        entry = {"name": sname, "groups": gl, "n": len(idx), "n_singleton": sum(snvs[i]["group"] == "singleton" for i in idx),
                 "catalogue": [int(x) for x in ct], "selection": cat}
        if cat["status"] == "ok":
            fit = refit_bootstrap(ct, S, names, cat["signatures"], n_boot=n_boot, seed=seed + 1)
            post = channel_posterior(S, fit["idx"], fit["exposure"])
            sig_order = [names[j] for j in fit["idx"]]
            entry.update({k: fit[k] for k in ("activities", "cosine", "cosine_lo", "cosine_hi")})
            entry["signatures"] = sig_order
            entry["channel_posterior"] = [[round(float(x), 4) for x in row] for row in post]
            entry["channel_argmax"] = [sig_order[int(np.argmax(row))] for row in post]
            for i in idx:
                assign[i] = entry["channel_argmax"][CH_INDEX[snvs[i]["channel"]]]
        else:
            entry.update({"status": "not_available", "reason": cat["reason"]})
        res_sets.append(entry)
    out = {"id": v["id"], "label": v["label"], "definition": v["definition"], "mode": v["mode"], "status": "ok",
           "n": len(members), "n_singleton": sum(snvs[i]["group"] == "singleton" for i in members),
           "groups": sorted({g for s in res_sets for g in s["groups"]}), "forced": forced, "sets": res_sets}
    if all(s.get("status") == "not_available" for s in res_sets):
        out.update({"status": "not_available", "reason": "; ".join(s["reason"] for s in res_sets)})
    return out, assign


def flows(assigns, snvs, vids, split):
    """Flow tables between consecutive and all version pairs: {from, to, rows:[[sig_a, sig_b, n]]} per split."""
    sel = {"all": lambda s: True, "singleton": lambda s: s["group"] == "singleton",
           "shared": lambda s: s["group"] in ("truncal", "shared3", "subclonal_other"),
           "private_tree": lambda s: s["group"] == "private_tree"}[split]
    out = []
    for a in range(len(vids)):
        for b in range(a + 1, len(vids)):
            A, B = assigns[vids[a]], assigns[vids[b]]
            tab = {}
            for i in set(A) | set(B):
                if not sel(snvs[i]):
                    continue
                k = (A.get(i, "(not in version)"), B.get(i, "(not in version)"))
                tab[k] = tab.get(k, 0) + 1
            out.append({"from": vids[a], "to": vids[b], "rows": sorted([[x, y, n] for (x, y), n in tab.items()], key=lambda r: -r[2])})
    return out


def diagnostics(snvs, S, names, art_S, art_names, cells):
    """Per-group spectrum metrics and read evidence; per-cell singleton burden."""
    ref = {n: S[:, names.index(n)] for n in FLAT_SIGS if n in names}
    for k, n in enumerate(art_names):
        ref[n] = art_S[:, k]
    groups = {}
    norm_ct = counts_of(s["channel"] for s in snvs if s["group"] == "normal_singleton")
    for g in GROUPS:
        mem = [s for s in snvs if s["group"] == g and s["channel"] in CH_INDEX]
        ct = counts_of(s["channel"] for s in mem)
        d = {"n": len(mem), "definition": GROUPS[g], "spectrum": [int(x) for x in ct]}
        if len(mem) >= 30:
            d["cosine"] = {n: round(cosine(ct, r), 4) for n, r in ref.items()}
            if norm_ct.sum() >= 30 and g != "normal_singleton":
                d["cosine"]["normal_singleton"] = round(cosine(ct, norm_ct), 4)
            d["cosine"]["flat_uniform"] = round(cosine(ct, np.ones(96)), 4)
            d["entropy_bits"] = round(entropy_bits(ct), 4)
            d["entropy_max_bits"] = round(math.log2(96), 4)
            cpg = sum(ct[CH_INDEX[c]] for c in SBS96 if c[2:5] == "C>T" and c[6] == "G")
            d["cpg_ct_frac"] = round(float(cpg / ct.sum()), 4)
            ev = [s["ev"] for s in mem if s["ev"]]
            if ev:
                vaf = [e["alt"] / e["depth"] for e in ev if e["depth"] > 0]
                alt = [e["alt"] for e in ev]
                sb = [min(e["adf_alt"], e["adr_alt"]) / max(1, e["alt"]) for e in ev]
                one_strand = [e["adf_alt"] == 0 or e["adr_alt"] == 0 for e in ev]
                q = lambda x: [round(float(v), 4) for v in np.percentile(x, [10, 25, 50, 75, 90])] if len(x) else None
                d["reads"] = {"n_with_evidence": len(ev), "vaf_q": q(vaf), "alt_reads_q": q(alt),
                              "minor_strand_frac_q": q(sb), "frac_one_strand": round(float(np.mean(one_strand)), 4),
                              "unit": "VCF AD_BAM/ADF_BAM/ADR_BAM summed over alt-carrying cells (singletons: the one carrier)",
                              "read_end_distance": "not available (needs per-read positions; not in the VCF)"}
        groups[g] = d
    burden = {}
    for s in snvs:
        if s["group"] in ("singleton", "singleton_lowconf", "normal_singleton") and s.get("carrier"):
            b = burden.setdefault(s["carrier"], {"singleton": 0, "singleton_lowconf": 0, "normal_singleton": 0})
            b[s["group"]] += 1
    return {"groups": groups, "cell_burden": burden}


def write_site_table(path, snvs, assigns, vids):
    """Exact SNVs of every version: id, group, channel, carrier cell (singletons), then the assigned signature per version."""
    with gzip.open(path, "wt") as fh:
        fh.write("\t".join(["mutation", "group", "channel", "is_tree", "carrier"] + vids) + "\n")
        for i, s in enumerate(snvs):
            if not any(i in assigns[v] for v in vids):
                continue
            fh.write("\t".join([s["id"], s["group"], s["channel"] or "", str(int(s["is_tree"])), s.get("carrier", "")] +
                                [assigns[v].get(i, "") for v in vids]) + "\n")


def run_patient(P, data_dir, vcf_sites, cosmic, artifacts, fasta, n_boot=1000, n_boot_select=200, seed=1, tree_codes=True, site_table=None):
    names, S = read_signature_table(cosmic)
    art_names, art_cols = [], []
    for spec in artifacts or []:
        nm, path = spec.split("=", 1)
        an, aS = read_signature_table(path)
        for k, a in enumerate(an):
            art_names.append(a if len(an) > 1 else nm); art_cols.append(aS[:, k])
    art_S = np.column_stack(art_cols) if art_cols else np.zeros((96, 0))
    S_all = np.column_stack([S, art_S]) if art_cols else S
    all_names = names + art_names
    snvs, cells, src = build_snvs(P, data_dir, vcf_sites, fasta)
    if not (vcf_sites and os.path.exists(vcf_sites)):
        src["singletons"] = "not available (no VCF evidence table)"
    versions, assigns = [], {}
    for v in VERSIONS:
        # artifact columns are only candidates where forced (version 9); other versions see COSMIC only
        if v.get("artifact"):
            res, asg = fit_version(v, snvs, S_all, all_names, art_names, n_boot, n_boot_select, seed)
        else:
            res, asg = fit_version(v, snvs, S, names, [], n_boot, n_boot_select, seed)
        versions.append(res); assigns[v["id"]] = asg
    vids = [v["id"] for v in VERSIONS]
    for res in versions:
        asg = assigns[res["id"]]
        res["n_singleton_cells"] = len({snvs[i].get("carrier") for i in asg if snvs[i]["group"] == "singleton"})
    if site_table:
        write_site_table(site_table, snvs, assigns, vids)
    tree_idx = [i for i, s in enumerate(snvs) if s["is_tree"]]
    sig_codes = sorted({x for a in assigns.values() for x in a.values()})
    code = {s: k for k, s in enumerate(sig_codes)}
    out = {
        "format": "gos-sc-signature-versions/1", "patient": P, "created": datetime.now().isoformat(timespec="seconds"),
        "cosmic": os.path.basename(cosmic), "artifact_signatures": art_names,
        "rule": {"selection": "forward selection from COSMIC (+ forced artifact columns), LRT alpha 0.05/#candidates, "
                              "min fraction 0.05, min cosine gain 0.02, first signature always taken, max 6; keep signatures selected in >= 50% of bootstrap re-selections",
                 "n_boot_select": n_boot_select, "n_boot": n_boot, "min_mutations": 50,
                 "posterior": "a_k S[c,k] / sum_j a_j S[c,j]; assignment = argmax"},
        "sources": src, "group_definitions": GROUPS,
        "group_counts": {g: sum(1 for s in snvs if s["group"] == g) for g in GROUPS},
        "versions": versions, "signature_codes": sig_codes,
        "flows": {sp: flows(assigns, snvs, vids, sp) for sp in ("all", "singleton", "shared", "private_tree")},
        "diagnostics": diagnostics(snvs, S, names, art_S, art_names, cells),
    }
    if tree_codes:
        out["tree_sites"] = {"id": [snvs[i]["id"] for i in tree_idx], "group": [snvs[i]["group"] for i in tree_idx],
                             "assign": {vid: [code.get(assigns[vid].get(i), -1) for i in tree_idx] for vid in vids}}
    return out


def cohort(per_patient):
    """Pooled flows and per-version signature-fraction summaries across patients."""
    pooled = {}
    for P, d in per_patient.items():
        for sp, fl in d["flows"].items():
            for f in fl:
                key = (sp, f["from"], f["to"])
                t = pooled.setdefault(key, {})
                for a, b, n in f["rows"]:
                    t[(a, b)] = t.get((a, b), 0) + n
    flows_out = {}
    for (sp, a, b), t in pooled.items():
        flows_out.setdefault(sp, []).append({"from": a, "to": b, "rows": sorted([[x, y, n] for (x, y), n in t.items()], key=lambda r: -r[2])})
    summary = []
    for P, d in per_patient.items():
        for v in d["versions"]:
            for s in v.get("sets", []):
                summary.append({"patient": P, "version": v["id"], "set": s["name"], "n": s["n"], "status": s.get("status", "ok"),
                                "cosine": s.get("cosine"), "signatures": [[a["signature"], a["fraction"], a["lo"], a["hi"]] for a in s.get("activities", [])]})
    return {"format": "gos-sc-signature-versions-cohort/1", "created": datetime.now().isoformat(timespec="seconds"),
            "patients": sorted(per_patient), "versions": [{"id": v["id"], "label": v["label"], "definition": v["definition"]} for v in VERSIONS],
            "flows": flows_out, "summary": summary}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--patients", required=True)
    ap.add_argument("--data", required=True, help="gOS data dir (inputs)")
    ap.add_argument("--out", required=True, help="output dir (writes <out>/<P>/signature_versions.json, <out>/_cohort/)")
    ap.add_argument("--vcf-sites", default="", help="pattern with {P}, e.g. work/{P}/vcf_sites.tsv.gz")
    ap.add_argument("--cosmic", required=True)
    ap.add_argument("--artifact", action="append", default=[])
    ap.add_argument("--fasta", default="")
    ap.add_argument("--site-table", default="", help="pattern with {P}: TSV.gz of the exact SNVs and assignments per version")
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--boot-select", type=int, default=200)
    a = ap.parse_args(argv)
    per = {}
    for P in a.patients.split(","):
        vs = a.vcf_sites.format(P=P) if a.vcf_sites else ""
        d = run_patient(P, a.data, vs, a.cosmic, a.artifact, a.fasta, a.boot, a.boot_select,
                        site_table=a.site_table.format(P=P) if a.site_table else None)
        os.makedirs(os.path.join(a.out, P), exist_ok=True)
        with open(os.path.join(a.out, P, "signature_versions.json"), "w") as fh:
            json.dump(d, fh, separators=(",", ":"))
        per[P] = d
        print(P, {v["id"]: (v["n"], v["status"]) for v in d["versions"]}, file=sys.stderr)
    if len(per) > 1:
        os.makedirs(os.path.join(a.out, "_cohort"), exist_ok=True)
        with open(os.path.join(a.out, "_cohort", "signature_versions.json"), "w") as fh:
            json.dump(cohort(per), fh, separators=(",", ":"))


if __name__ == "__main__":
    main()
