#!/usr/bin/env python3
"""Dropout-aware clonality for read-evidence driver rows (hotspot rescue, pooled driver scan, called-SNV rows).

Why: a read-evidence row's cell_fraction is carriers / covered tumour cells. In PTA single cells a clonal heterozygous variant is
seen in far fewer than 100% of covered cells (allelic dropout / imbalance, 1-2 reads per cell at GC-rich sites such as the TERT
promoter), so "cell_fraction < 0.85 => subclonal" calls clonal hotspots subclonal.

Rule (dropout_aware_clonality):
  expected detection per covered cell = detection rate (fraction of covered tumour cells with >= 1 alt read) of the patient's
      truncal SNVs (snv_matrix.json category "truncal") at the same depth bin (truncal_detection_by_depth)
  p_detection         one-sided Poisson-binomial P(X <= observed carriers | expected per-cell detection)
  p_clade_restriction permutation P of the best tree-clade F1 of the carriers (carriers re-drawn among the covered cells, weighted
                      by their expected detection, so deeper cells are not mistaken for a clade)
  status  indeterminate  n_covered < MIN_COVERED, or no truncal reference / per-cell depth and observed < RAW_CLONAL,
                         or detection significantly low but carriers not clade-restricted (locus dropout / artefact)
          subclonal      p_detection < ALPHA_DETECTION and p_clade_restriction < ALPHA_CLADE
          clonal         p_detection >= ALPHA_DETECTION (consistent with clonal), or (no reference) observed >= RAW_CLONAL
                         (dropout can only lower detection)
The raw cell_fraction is never changed; rows get clonality, clonality_method = "dropout-aware" and clonality_stats.
d3-free, numpy only. annotate_events() re-annotates a patient's filtered.events.json rows in place (see gos_sc_clonality.py).
"""
import json
import os
import re

import numpy as np

METHOD = "dropout-aware"
DEPTH_BINS = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 6), (7, 9), (10, 14), (15, 24), (25, None)]
ALPHA_DETECTION = 0.01
ALPHA_CLADE = 0.05
MIN_COVERED = 10
MIN_BIN_OBS = 20        # truncal-site observations needed for a depth bin's rate; sparser bins borrow the nearest bin
RAW_CLONAL = 0.85       # fallback when no truncal reference / per-cell depth: observed >= this is clonal
N_PERM = 2000
READ_EVIDENCE_SOURCES = ("site evidence (not in SNV calls)", "pooled driver scan", "called SNV (not a tree site)")


# ---------------------------------------------------------------- reference

def depth_bin(d):
    """Index of DEPTH_BINS holding depth d (>= 1), else None."""
    if d is None or d < 1:
        return None
    for i, (lo, hi) in enumerate(DEPTH_BINS):
        if d >= lo and (hi is None or d <= hi):
            return i
    return None


def bin_label(i):
    lo, hi = DEPTH_BINS[i]
    return f"{lo}+" if hi is None else (str(lo) if lo == hi else f"{lo}-{hi}")


def truncal_detection_by_depth(snv_matrix, tumour_cells=None):
    """Detection rate of truncal SNVs per depth bin from a compact snv_matrix.json
    ({variants: [{category}], cells: {cell: [[variant index, ref, alt, call], ...]}}).
    -> {bins, labels, n, detected, rate, n_sites, n_cells, source} or None when the patient has no truncal sites."""
    if not snv_matrix:
        return None
    truncal = {i for i, v in enumerate(snv_matrix.get("variants") or []) if v.get("category") == "truncal"}
    if not truncal:
        return None
    n = np.zeros(len(DEPTH_BINS), int)
    det = np.zeros(len(DEPTH_BINS), int)
    used = 0
    cells = snv_matrix.get("cells") or {}
    for cell, entries in cells.items():
        if tumour_cells is not None and cell not in tumour_cells:
            continue
        used += 1
        for x in entries:
            if x[0] not in truncal:
                continue
            ref, alt = int(x[1] or 0), int(x[2] or 0)
            b = depth_bin(ref + alt)
            if b is None:
                continue
            n[b] += 1
            det[b] += alt > 0
    if not n.sum():
        return None
    return {"bins": [list(b) for b in DEPTH_BINS], "labels": [bin_label(i) for i in range(len(DEPTH_BINS))], "n": n.tolist(), "detected": det.tolist(),
            "rate": [round(float(d / k), 4) if k >= MIN_BIN_OBS else None for d, k in zip(det, n)],
            "n_sites": len(truncal), "n_cells": used, "source": "snv_matrix.json truncal sites (fraction of covered tumour cells with >= 1 alt read)"}


def expected_rate(ref, d):
    """Truncal detection rate at depth d (nearest bin with data when d's bin is sparse)."""
    b = depth_bin(d)
    if b is None or not ref:
        return None
    rates = ref["rate"]
    if rates[b] is not None:
        return rates[b]
    for off in range(1, len(rates)):
        for j in (b - off, b + off):
            if 0 <= j < len(rates) and rates[j] is not None:
                return rates[j]
    return None


# ---------------------------------------------------------------- tests

def poisson_binomial_cdf(k, ps):
    """P(X <= k), X = sum of independent Bernoulli(ps)."""
    dist = np.zeros(len(ps) + 1)
    dist[0] = 1.0
    for p in ps:
        dist[1:] = dist[1:] * (1 - p) + dist[:-1] * p
        dist[0] *= 1 - p
    return float(min(1.0, dist[: k + 1].sum()))


def _clade_matrix(cells, clade_of_cell):
    """Boolean (clades x cells) membership over the covered cells; clades with >= 2 and < all covered cells."""
    ids = {}
    for c in cells:
        for cl in (clade_of_cell or {}).get(c) or ():
            ids.setdefault(cl, len(ids))
    if not ids:
        return np.zeros((0, len(cells)), bool), []
    m = np.zeros((len(ids), len(cells)), bool)
    for j, c in enumerate(cells):
        for cl in (clade_of_cell or {}).get(c) or ():
            m[ids[cl], j] = True
    size = m.sum(1)
    keep = (size >= 2) & (size < len(cells))
    names = [k for k, _ in sorted(ids.items(), key=lambda kv: kv[1])]
    return m[keep], [n for n, k in zip(names, keep) if k]


def _best_f1(m, x):
    """Best F1 of carrier vector(s) x (cells or cells x perms) against the clade rows of m."""
    k = x.sum(0)
    inter = m.astype(float) @ x.astype(float)
    f1 = 2 * inter / (m.sum(1)[:, None] + k) if x.ndim == 2 else 2 * inter / (m.sum(1) + k)
    return f1.max(0)


def clade_restriction(carrier, weights, m, n_perm=N_PERM, seed=0):
    """(p, best F1, best clade row) for carrier (bool per covered cell) against clades m, by weighted permutation."""
    k = int(carrier.sum())
    if m.shape[0] == 0 or k == 0 or k == len(carrier):
        return None, None, None
    obs = _best_f1(m, carrier)
    best = int(np.argmax(2 * (m.astype(float) @ carrier.astype(float)) / (m.sum(1) + k)))
    rng = np.random.default_rng(seed)
    w = np.asarray(weights, float)
    w = np.where(w > 0, w, 1e-6)
    w = w / w.sum()
    perm = np.zeros((len(carrier), n_perm), bool)
    for j in range(n_perm):
        perm[rng.choice(len(carrier), k, replace=False, p=w), j] = True
    null = _best_f1(m, perm)
    p = (1 + int((null >= obs - 1e-12).sum())) / (1 + n_perm)
    return p, float(obs), best


def dropout_aware_clonality(carriers_by_cell, depth_by_cell, clade_of_cell, truncal_detection_by_depth, n_perm=N_PERM, seed=0):
    """carriers_by_cell {cell: bool (>= 1 alt read)}, depth_by_cell {cell: reads} (None = per-cell depth unknown),
    clade_of_cell {cell: [clade ids the cell belongs to]} (nested tree clades), truncal_detection_by_depth (see above; None = no reference).
    -> {status, expected_detection, observed_detection, p_detection, p_clade_restriction, n_covered, n_carriers, ...}."""
    ref = truncal_detection_by_depth
    if depth_by_cell is None:
        cells = sorted(carriers_by_cell)
    else:
        cells = sorted(c for c, d in depth_by_cell.items() if d and d > 0)
    carrier = np.array([bool(carriers_by_cell.get(c)) for c in cells], bool)
    n_cov, n_car = len(cells), int(carrier.sum())
    out = {"status": None, "method": METHOD, "n_covered": n_cov, "n_carriers": n_car,
           "observed_detection": round(n_car / n_cov, 4) if n_cov else None, "expected_detection": None,
           "p_detection": None, "p_clade_restriction": None, "best_clade": None, "best_clade_f1": None,
           "alpha_detection": ALPHA_DETECTION, "alpha_clade": ALPHA_CLADE, "reason": None}
    if depth_by_cell is not None and n_cov:
        out["median_depth"] = float(np.median([depth_by_cell[c] for c in cells]))
    if n_cov < MIN_COVERED:
        out.update(status="indeterminate", reason=f"only {n_cov} covered tumour cells (< {MIN_COVERED})")
        return out
    ps = None
    if ref is not None and depth_by_cell is not None:
        ps = [expected_rate(ref, depth_by_cell[c]) for c in cells]
        if any(p is None for p in ps):
            ps = None
    if ps is None:
        why = "no truncal SNV reference" if ref is None else "per-cell depth not available"
        if out["observed_detection"] >= RAW_CLONAL:
            out.update(status="clonal", reason=f"{why}; observed detection {out['observed_detection']:.2f} >= {RAW_CLONAL} (dropout can only lower detection)")
        else:
            out.update(status="indeterminate", reason=f"{why}; observed detection {out['observed_detection']:.2f} is a lower bound under allelic dropout")
        return out
    ps = np.clip(np.array(ps, float), 1e-6, 1 - 1e-6)
    out["expected_detection"] = round(float(ps.mean()), 4)
    out["p_detection"] = poisson_binomial_cdf(n_car, ps)
    m, names = _clade_matrix(cells, clade_of_cell)
    pc, f1, best = clade_restriction(carrier, ps, m, n_perm, seed)
    out.update(p_clade_restriction=pc, best_clade_f1=round(f1, 3) if f1 is not None else None, best_clade=names[best] if best is not None else None,
               n_clades_tested=int(m.shape[0]))
    if ref.get("scope") == "cohort":
        out["reference_scope"] = "cohort (no truncal SNVs in this patient)"
    if out["p_detection"] >= ALPHA_DETECTION:
        out.update(status="clonal", reason=f"detection {out['observed_detection']:.2f} consistent with truncal SNVs at the same depth (expected {out['expected_detection']:.2f}, P = {out['p_detection']:.3g})")
    elif pc is not None and pc < ALPHA_CLADE:
        out.update(status="subclonal", reason=f"detection {out['observed_detection']:.2f} below truncal expectation {out['expected_detection']:.2f} (P = {out['p_detection']:.2g}) and carriers clade-restricted (best clade F1 {f1:.2f}, permutation P = {pc:.3g})")
    else:
        out.update(status="indeterminate", reason=f"detection {out['observed_detection']:.2f} below truncal expectation {out['expected_detection']:.2f} (P = {out['p_detection']:.2g}) "
                                                  f"but carriers not clade-restricted ({'no tree clades' if pc is None else f'permutation P = {pc:.3g}'}): locus-specific dropout or artefact")
    for k in ("p_detection", "p_clade_restriction"):
        if out[k] is not None:
            out[k] = float(f"{out[k]:.4g}")
    return out


# ---------------------------------------------------------------- patient inputs

def _load(path):
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return json.load(fh)


def parse_newick_clades(text):
    """Leaf-name sets of the internal nodes of a newick tree."""
    text = text.strip().rstrip(";")
    stack, clades, i = [[]], [], 0
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
            if ch == ")":       # internal-node label / branch length
                i += len(re.match(r"[^(),;]*", text[i:]).group(0))
        else:
            name = re.match(r"[^(),:;]+", text[i:]).group(0)
            i += len(name)
            i += len(re.match(r"(:[^(),;]*)?", text[i:]).group(0))
            stack[-1].append({name})
    return clades


def site_id_of(row):
    """chr<chrom>:<pos>:<ref>>alt> of a SNV / hotspot row (Variant_g "5:1295113-1295113 G>A"), else None."""
    m = re.match(r"^(?:chr)?([^:]+):(\d+)-(\d+) ([ACGT])>([ACGT])$", str(row.get("Variant_g") or ""))
    if not m or m.group(2) != m.group(3):
        return None
    return f"chr{m.group(1)}:{m.group(2)}:{m.group(4)}>{m.group(5)}"


def tumour_cells(data, patient):
    calls = _load(os.path.join(data, patient, "precompute", "region_calls.json")) or {}
    clone = calls.get("clone") or []
    return {c for i, c in enumerate(calls.get("cells") or []) if not (i < len(clone) and str(clone[i]).lower() == "normal")} or None


def patient_reference(data, patient, tumour=None):
    """truncal_detection_by_depth of data/<P>/snv_matrix.json (None when absent or without truncal sites)."""
    snv = _load(os.path.join(data, patient, "snv_matrix.json"))
    return truncal_detection_by_depth(snv, tumour) if snv else None


def cohort_reference(refs):
    """Pooled truncal detection by depth of several patients' references ({patient: ref}); detection at a given depth is a property of
    the PTA chemistry and is close across patients (per-patient rates are kept in `patients`)."""
    refs = {p: r for p, r in refs.items() if r}
    if not refs:
        return None
    n = np.sum([r["n"] for r in refs.values()], 0)
    det = np.sum([r["detected"] for r in refs.values()], 0)
    return {"bins": [list(b) for b in DEPTH_BINS], "labels": [bin_label(i) for i in range(len(DEPTH_BINS))], "n": n.tolist(), "detected": det.tolist(),
            "rate": [round(float(d / k), 4) if k >= MIN_BIN_OBS else None for d, k in zip(det, n)],
            "n_sites": int(sum(r["n_sites"] for r in refs.values())), "n_cells": int(sum(r["n_cells"] for r in refs.values())),
            "source": "truncal SNVs of " + ", ".join(sorted(refs)) + " pooled (snv_matrix.json)",
            "patients": {p: dict(zip(r["labels"], r["rate"])) for p, r in sorted(refs.items())}}


class PatientContext:
    """Per-patient inputs: region_calls (cell order, normals), site_evidence reads, truncal reference, tree clades, allsnv pileups."""

    def __init__(self, data, patient, allsnv_work=None, reference=False, fallback_reference=None):
        """reference: precomputed patient reference (False = compute from snv_matrix.json); fallback_reference: used when the
        patient has no truncal SNVs (cohort_reference of the other patients), marked scope = "cohort"."""
        d = os.path.join(data, patient)
        self.patient = patient
        self.calls = _load(os.path.join(d, "precompute", "region_calls.json"))
        self.ev = _load(os.path.join(d, "precompute", "site_evidence.json")) or {}
        cells = (self.calls or {}).get("cells") or []
        clone = (self.calls or {}).get("clone") or []
        self.cells = cells
        self.tumour_idx = [i for i in range(len(cells)) if not (i < len(clone) and str(clone[i]).lower() == "normal")]
        self.tumour = {cells[i] for i in self.tumour_idx}
        self.reference = patient_reference(data, patient, self.tumour or None) if reference is False else reference
        if self.reference is not None:
            self.reference = dict(self.reference, scope="patient")
        elif fallback_reference is not None:
            self.reference = dict(fallback_reference, scope="cohort")
        self.clade_of_cell = {}
        tp = os.path.join(d, "tree.nwk")
        if os.path.exists(tp):
            with open(tp) as fh:
                for k, s in enumerate(parse_newick_clades(fh.read())):
                    for c in s:
                        self.clade_of_cell.setdefault(c, []).append(f"node{k}")
        self.pileup = None
        if allsnv_work:
            self.pileup = _load(os.path.join(allsnv_work, f"pileup_own_{patient}.json"))

    def evidence(self, sid):
        """(carriers_by_cell, depth_by_cell, source) for a site from the best per-cell read source, else (None, None, None)."""
        cells, t = self.cells, self.tumour_idx
        ev = (self.ev.get("sites") or {}).get(sid)
        if ev:
            reads = ev.get("reads") or {}
            pooled = (((ev.get("presets") or {}).get("all") or {}).get("tumor") or {})
            if reads.get("rows") and not reads.get("capped") and ev.get("sample_every", 1) == 1:
                tset = set(t)
                dp, alt = {}, {}
                for r in reads["rows"]:
                    c, al = r[0], r[1]
                    if c is None or c not in tset or al not in ("A", "R"):
                        continue
                    dp[cells[c]] = dp.get(cells[c], 0) + 1
                    alt[cells[c]] = alt.get(cells[c], False) or al == "A"
                if sum(alt.values()) == pooled.get("cells_alt") and len(dp) == pooled.get("cells_cov"):
                    return alt, dp, "site_evidence reads (all reads)"
        v = ((self.calls or {}).get("dna") or {}).get(sid)
        if v:
            a, d = v.get("alt") or [], v.get("dp") or []
            keep = [i for i in t if i < len(d) and (d[i] or 0) > 0]
            dp = {cells[i]: int(d[i]) for i in keep}
            return {cells[i]: i < len(a) and (a[i] or 0) > 0 for i in keep}, dp, "region_calls (caller-filtered reads)"
        s = ((self.pileup or {}).get("sites") or {}).get(sid)
        if s and s.get("cells"):
            a, d = s["cells"].get("alt") or [], s["cells"].get("dp") or []
            dp = {cells[i]: int(d[i]) for i in t if i < len(d) and d[i] > 0}
            return {cells[i]: a[i] > 0 for i in t if i < len(d) and d[i] > 0}, dp, "BAM pileup (gos_sc_allsnv)"
        return None, None, None

    def assess(self, row, n_perm=N_PERM):
        sid = site_id_of(row)
        car, dp, src = self.evidence(sid) if sid else (None, None, None)
        if car is None:
            # per-cell depth not kept (pooled scan): carriers / covered counts only
            ids = [c for c in str(row.get("cell_ids") or "").split(",") if c]
            try:
                cov = int(float(row.get("covered_tumour_cells")))
            except (TypeError, ValueError):
                cov = None
            car = {c: True for c in ids}
            for j in range(max(0, (cov or len(ids)) - len(ids))):
                car[f"__covered_no_alt_{j}"] = False
            src = "row carriers / covered counts (per-cell depth not kept)"
            res = dropout_aware_clonality(car, None, self.clade_of_cell, self.reference, n_perm)
        else:
            res = dropout_aware_clonality(car, dp, self.clade_of_cell, self.reference, n_perm, seed=sum(map(ord, sid)))
        res["site"] = sid
        res["evidence_source"] = src
        if self.reference:
            res["reference"] = {"scope": self.reference.get("scope"), "source": self.reference["source"], "n_truncal_sites": self.reference["n_sites"],
                                "rate_by_depth": dict(zip(self.reference["labels"], self.reference["rate"]))}
        return res


def is_read_evidence(row):
    return row.get("source") in READ_EVIDENCE_SOURCES


def apply(row, res):
    """Write the result onto an event row (raw cell_fraction untouched)."""
    if row.get("clonality") is not None and row.get("clonality_method") != METHOD:
        row["clonality_fraction_rule"] = row["clonality"]       # earlier label (carriers / covered >= 0.85), kept for the diff
    row["clonality"] = res["status"]
    row["clonality_method"] = METHOD
    row["clonality_stats"] = res
    return row


def annotate_events(events, ctx, n_perm=N_PERM):
    """Annotate every read-evidence row of a filtered.events list in place; -> list of (row, result)."""
    out = []
    for e in events:
        if is_read_evidence(e) and site_id_of(e):
            res = ctx.assess(e, n_perm)
            apply(e, res)
            out.append((e, res))
    return out


def event_clonality(e):
    """'clonal' | 'subclonal' | 'indeterminate' from the exported dropout-aware field, else None (caller falls back)."""
    if e and e.get("clonality_method") == METHOD and e.get("clonality") in ("clonal", "subclonal", "indeterminate"):
        return e["clonality"]
    return None


_COHORT = {}


def cohort_fallback(data, exclude=()):
    """cohort_reference of every patient under data with truncal SNVs (cached per data dir)."""
    if data not in _COHORT:
        pats = sorted(p for p in os.listdir(data) if p[:3] in ("BWH", "MGH") and p[3:].isdigit())
        _COHORT[data] = {p: patient_reference(data, p, tumour_cells(data, p)) for p in pats}
    return cohort_reference({p: r for p, r in _COHORT[data].items() if p not in exclude})


def annotate_rows(data, patient, rows, allsnv_work=None, n_perm=N_PERM):
    """Dropout-aware clonality for freshly built read-evidence rows (used by the row-building exporters); rows changed in place.
    A failure to load the inputs leaves the rows unannotated (the exporter's own fields stay) and is reported, never invented."""
    try:
        ctx = PatientContext(data, patient, allsnv_work)
        if ctx.reference is None:
            ctx.reference = cohort_fallback(data)
            if ctx.reference is not None:
                ctx.reference = dict(ctx.reference, scope="cohort")
    except (OSError, ValueError) as err:
        print(f"clonality: {patient}: inputs not readable ({err}); rows left without dropout-aware clonality")
        return []
    return annotate_events(rows, ctx, n_perm)
