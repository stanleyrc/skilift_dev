#!/usr/bin/env python3
"""Per-mutation signature assignment from ONE patient-level fit, and group (clade / truncal / any cell set) signature fractions.

Method (the same for every patient):
  1. The patient's signature set and activities come from one fit of ALL the patient's SNVs
     (signatures_joint.json: set "all", cosmic_fit + bootstrap-stable additions). Nothing is added or forced per group.
  2. Every mutation is assigned to the signatures by its posterior given its SBS96 channel k and the patient activities a_s:
         P(s | k) = a_s * S[k, s] / sum_t a_t * S[k, t]            (S = COSMIC profile, columns sum to 1)
  3. The fraction of a group of mutations (a clade's private SNVs, the truncal set, ...) in signature s is the mean of the
     per-mutation posteriors over the group; 95% CI by bootstrap over the group's mutations.
  4. A signature that is not in the patient set cannot appear in a group. A real group-specific process therefore shows up as
     a poor fit: for every group we report the cosine between the observed 96-channel spectrum and (a) the patient-level
     mixture and (b) the best non-negative refit inside the patient set, each with a parametric-simulation p (n mutations drawn
     from the patient mixture, same statistic; flag_set_misfit when p < 0.01). Because the set never reproduces the pooled
     spectrum exactly, the group flag (flag_poor_fit) is clade-specific: the group's cosine to the patient's pooled observed
     spectrum against n draws from that spectrum (p_vs_pooled < 0.01). We do NOT add signatures to explain a flagged group.
  5. Group differences (vs truncal and vs the nearest fitted non-trunk ancestor): bootstrap over mutations, two-sided p,
     BH across the cohort (done by the caller).

Library functions are pure numpy and are unit-tested on synthetic data (tests/test_signature_assignment.py).

Usage (writes a staged copy; never touches live data):
  python3 gos_sc_signature_assignment.py --data <gOS data dir> --cosmic COSMIC_v3.4_SBS96_GRCh38.txt --out <stage dir> \
      [--patients P1,P2] [--n-boot 1000] [--n-sim 1000]
Outputs: <out>/<P>/signature_assignment.json (per-mutation posteriors, group table) and <out>/_cohort/signature_assignment.json
(clade_signatures / signature_shifts in the schema of subclonal_drivers.json).
"""
import argparse
import datetime
import json
import os
import sys

import numpy as np
from scipy.optimize import nnls

MIN_CLADE_SITES = 50
FIT_FLAG_P = 0.01


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def load(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def read_cosmic(path):
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        rows = [l.rstrip("\n").split("\t") for l in f]
    channels = [r[0] for r in rows]
    M = np.array([[float(x) for x in r[1:]] for r in rows])
    return channels, {s: M[:, j] for j, s in enumerate(head[1:])}


def bh(ps):
    m = len(ps)
    order = sorted(range(m), key=lambda i: ps[i])
    q = [1.0] * m
    prev = 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        prev = min(prev, ps[i] * m / (rank + 1))
        q[i] = prev
    return [float(min(1, x)) for x in q]


def boot_p(diff_boot):
    n = len(diff_boot)
    lo = (np.sum(diff_boot <= 0) + 1) / (n + 1)
    hi = (np.sum(diff_boot >= 0) + 1) / (n + 1)
    return float(min(1.0, 2 * min(lo, hi)))


# ------------------------------------------------------------------ core

def channel_posterior(S, act):
    """S: 96 x m profiles (columns sum to 1); act: length m activities (>= 0). Returns 96 x m matrix P[k, s] = P(s | channel k).
    Channels where every active signature has zero probability get a uniform-by-activity fallback (never NaN)."""
    act = np.asarray(act, float)
    W = S * act[None, :]
    tot = W.sum(axis=1, keepdims=True)
    fallback = act / act.sum() if act.sum() > 0 else np.full(len(act), 1.0 / len(act))
    return np.where(tot > 0, W / np.where(tot > 0, tot, 1.0), fallback[None, :])


def counts96(idx, K=96):
    return np.bincount(np.asarray(idx, dtype=int), minlength=K).astype(float)


def group_fractions(counts, Pch):
    """Mean per-mutation posterior of a group with channel counts."""
    n = counts.sum()
    return (counts @ Pch) / n if n > 0 else np.zeros(Pch.shape[1])


def bootstrap_fractions(idx, Pch, n_boot, rng):
    idx = np.asarray(idx, dtype=int)
    n = len(idx)
    K = Pch.shape[0]
    cnt = counts96(idx, K)
    point = group_fractions(cnt, Pch)
    B = np.empty((n_boot, Pch.shape[1]))
    for b in range(n_boot):
        B[b] = group_fractions(counts96(idx[rng.integers(0, n, n)], K), Pch)
    return point, B


def cosine(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(a @ b / (na * nb)) if na > 0 and nb > 0 else 0.0


def refit_cosine(counts, S):
    """Best non-negative combination of the patient-set profiles for the observed spectrum; returns (cosine, residual fraction)."""
    if counts.sum() == 0:
        return 0.0, 1.0
    a, _ = nnls(S, counts)
    rec = S @ a
    return cosine(counts, rec), float(np.linalg.norm(counts - rec) / np.linalg.norm(counts))


def fit_quality(counts, S, act, n_sim, rng, pooled=None):
    """Goodness of fit of a group's spectrum to the patient set. Returns dict with cosines and simulation p-values
    (probability that n mutations drawn from the patient mixture fit as badly or worse). The patient set never fits the pooled
    spectrum perfectly (p_fit_* is then small for every large group), so when `pooled` (the patient's 96-channel counts) is given
    a clade-specific test is added: p_vs_pooled = probability that n mutations drawn from the pooled spectrum are at least as
    far (cosine) from it as the group; a group departing from the patient's overall spectrum is where a missing process would show."""
    n = int(counts.sum())
    act = np.asarray(act, float)
    mix = S @ act
    mix = mix / mix.sum()
    cos_mix = cosine(counts, mix)
    cos_ref, resid = refit_cosine(counts, S)
    sim_mix = np.empty(n_sim)
    sim_ref = np.empty(n_sim)
    for i in range(n_sim):
        c = rng.multinomial(n, mix).astype(float)
        sim_mix[i] = cosine(c, mix)
        sim_ref[i] = refit_cosine(c, S)[0]
    extra = {}
    if pooled is not None and pooled.sum() > 0:
        pp = pooled / pooled.sum()
        cp = cosine(counts, pp)
        simp = np.array([cosine(rng.multinomial(n, pp).astype(float), pp) for _ in range(n_sim)])
        extra = {"cosine_vs_pooled": round(cp, 4), "p_vs_pooled": float((np.sum(simp <= cp) + 1) / (n_sim + 1))}
    return {
        **extra,
        "n": n,
        "cosine_patient_mixture": round(cos_mix, 4),
        "cosine_best_refit": round(cos_ref, 4),
        "refit_residual": round(resid, 4),
        "p_fit_mixture": float((np.sum(sim_mix <= cos_mix) + 1) / (n_sim + 1)),
        "p_fit_refit": float((np.sum(sim_ref <= cos_ref) + 1) / (n_sim + 1)),
        "expected_cosine_best_refit_median": round(float(np.median(sim_ref)), 4),
    }


def r4(x):
    return round(float(x), 4)


def fraction_rows(sigs, point, B):
    return [{"signature": s, "frac": r4(point[j]), "lo": r4(np.quantile(B[:, j], 0.025)), "hi": r4(np.quantile(B[:, j], 0.975))} for j, s in enumerate(sigs)]


# ------------------------------------------------------------------ patient

def patient_set(sj, cosmic):
    """(signatures, activities) of the patient-level fit ('all' tree set); signatures absent from COSMIC are dropped and reported."""
    sets = {s["name"]: s for s in sj.get("sets", [])}
    allset = sets.get("all")
    if not allset:
        return None, None, ["no 'all' set in signatures_joint.json"]
    acts = [(a["signature"], float(a.get("activity") or 0)) for a in allset.get("activities", [])]
    missing = [s for s, a in acts if a > 0 and s not in cosmic]
    acts = [(s, a) for s, a in acts if a > 0 and s in cosmic]
    return [s for s, _ in acts], np.array([a for _, a in acts]), missing


def patient_assignment(P, snv, clades_doc, sj, cosmic, chan_index, n_boot, n_sim, rng):
    """Returns (patient doc with per-mutation posteriors and group table, tests list)."""
    out = {"patient": P, "available": False}
    if not (snv and sj):
        out["reason"] = "no snv_matrix.json" if not snv else "no signatures_joint.json"
        return out, []
    sigs, act, missing = patient_set(sj, cosmic)
    if not sigs:
        out["reason"] = "; ".join(missing)
        return out, []
    S = np.column_stack([cosmic[s] for s in sigs])
    Pch = channel_posterior(S, act)
    variants = snv["variants"]
    vby = {v["id"]: v for v in variants}
    ids, chan = [], []
    for v in variants:
        k = chan_index.get(v.get("context"))
        if k is not None:
            ids.append(v["id"])
            chan.append(k)
    chan = np.array(chan, dtype=int)
    probs = Pch[chan]
    pos = {i: n for n, i in enumerate(ids)}

    def idx_of(sids):
        return [int(chan[pos[i]]) for i in sids if i in pos]

    groups = []   # (name, kind, ids)
    groups.append(("truncal", "set", [v["id"] for v in variants if v.get("category") == "truncal" and v.get("cellphy_input")]))
    groups.append(("subclonal", "set", [v["id"] for v in variants if v.get("category") == "subclonal"]))
    groups.append(("private", "set", [v["id"] for v in variants if v.get("category") == "private"]))
    clades = (clades_doc or {}).get("clades", [])
    for c in clades:
        groups.append((c["id"], "clade", c.get("private_snvs") or []))
    allidx = idx_of(ids)
    cnt_all = counts96(allidx)
    pat_point = group_fractions(cnt_all, Pch)

    rows, fits, setrows = [], {}, {}
    meta = {c["id"]: c for c in clades}
    for name, kind, sids in groups:
        ix = idx_of(sids)
        r = {"group": name, "kind": kind, "n_sites": len(ix)}
        if kind == "clade":
            c = meta[name]
            r.update({"parent": c.get("parent"), "n_cells": c.get("n_cells"), "depth": c.get("depth"), "is_trunk": bool(c.get("is_trunk")) or c.get("parent") is None, "clone_id": c.get("clone_id"), "private_drivers": (c.get("private_drivers") or "").strip()})
        if len(ix) < MIN_CLADE_SITES:
            r["fitted"] = False
            r["reason"] = f"n sites {len(ix)} < {MIN_CLADE_SITES}"
        else:
            pt, B = bootstrap_fractions(ix, Pch, n_boot, rng)
            fits[name] = (pt, B)
            r["fitted"] = True
            r["fractions"] = fraction_rows(sigs, pt, B)
            q = fit_quality(counts96(ix), S, act, n_sim, rng, pooled=cnt_all)
            q["flag_set_misfit"] = bool(q["p_fit_refit"] < FIT_FLAG_P)
            q["flag_poor_fit"] = bool(q["p_vs_pooled"] < FIT_FLAG_P)   # clade spectrum departs from the patient's overall spectrum
            r["fit"] = q
        rows.append(r)

    tests = []
    for r in rows:
        if r["kind"] != "clade" or not r.get("fitted") or r["is_trunk"]:
            continue
        pt, B = fits[r["group"]]
        comps = []
        if "truncal" in fits:
            comps.append(("truncal", fits["truncal"][0], fits["truncal"][1]))
        par = r["parent"]
        while par and par in meta and par not in fits:
            par = meta[par].get("parent")
        if par and par in fits and not (meta[par].get("is_trunk") or meta[par].get("parent") is None):
            comps.append((par, fits[par][0], fits[par][1]))
        r["vs"] = []
        for name, p2, B2 in comps:
            for j, s in enumerate(sigs):
                db = B[:, j] - B2[:, j]
                t = {"clade": r["group"], "vs": name, "signature": s, "diff": r4(pt[j] - p2[j]), "lo": r4(np.quantile(db, 0.025)), "hi": r4(np.quantile(db, 0.975)), "p": boot_p(db), "frac": r4(pt[j]), "frac_ref": r4(p2[j])}
                r["vs"].append(t)
                tests.append(t)

    pq = fit_quality(cnt_all, S, act, n_sim, rng)
    out.update({
        "available": True,
        "patient_fit": pq,
        "signatures": sigs,
        "activities": [{"signature": s, "activity": r4(a), "fraction": r4(pat_point[j])} for j, (s, a) in enumerate(zip(sigs, act))],
        "dropped_not_in_cosmic": missing,
        "method": "P(sig | SBS96 channel) = activity x COSMIC profile, normalised over the patient signature set; group fraction = mean over its mutations",
        "n_mutations": len(ids),
        "groups": rows,
        "mutations": {
            "ids": ids,
            "channel": [int(k) for k in chan],
            "probs": [[round(float(x), 3) for x in row] for row in probs],
            "top": [sigs[int(j)] for j in probs.argmax(axis=1)],
        },
    })
    return out, tests


def cohort_blocks(per_patient, tests_by_patient):
    """clade_signatures + signature_shifts in the schema of subclonal_drivers.json; BH q across the cohort."""
    all_tests = [t for ts in tests_by_patient.values() for t in ts]
    qs = bh([t["p"] for t in all_tests])
    for t, q in zip(all_tests, qs):
        t["q"] = q
    blocks = {}
    for P, d in per_patient.items():
        if not d.get("available"):
            continue
        tr = next((g for g in d["groups"] if g["group"] == "truncal"), None)
        blocks[P] = {
            "signatures": d["signatures"],
            "added": [],
            "truncal": {"n_sites": tr["n_sites"], "fractions": tr.get("fractions"), "fit": tr.get("fit")} if tr else None,
            "clades": [dict(clade=g["group"], **{k: v for k, v in g.items() if k not in ("group", "kind")}) for g in d["groups"] if g["kind"] == "clade"],
        }
    shifts = [{"patient": P, **t} for P, ts in tests_by_patient.items() for t in ts if t["q"] < 0.1 and abs(t["diff"]) >= 0.05]
    shifts.sort(key=lambda t: (t["q"], -abs(t["diff"])))
    return blocks, shifts


METHOD_TEXT = (
    "One signature set per patient, fitted once on all of the patient's SNVs (signatures_joint.json, set 'all'); no signature is added or forced for any clade. "
    "Each mutation is assigned by its posterior P(signature | SBS96 channel) = activity x COSMIC v3.4 profile, normalised over the patient set; a clade's fraction is the mean posterior over its private SNVs "
    f"(clades with < {MIN_CLADE_SITES} sites are not fitted), 95% CI by bootstrap over mutations. A signature outside the patient set cannot appear in a clade; "
    "a clade whose spectrum departs from the patient's overall spectrum more than chance allows (cosine, simulation p < 0.01) is flagged instead of being given extra signatures; the cosine of every group to the patient set is reported. "
    "Each clade is compared with the truncal set and with its nearest fitted non-trunk ancestor: two-sided bootstrap p over mutations, BH q across the cohort; shifts listed: q < 0.1 and |difference| >= 5 points. "
    "Because the activities are shared, a clade differs from the truncal set only through its channel composition."
)


def run(data, cosmic_path, patients, n_boot, n_sim, seed):
    rng = np.random.default_rng(seed)
    channels, cosmic = read_cosmic(cosmic_path)
    chan_index = {c: k for k, c in enumerate(channels)}
    per, tests = {}, {}
    for P in patients:
        log("patient", P)
        pd_ = os.path.join(data, P)
        d, t = patient_assignment(P, load(os.path.join(pd_, "snv_matrix.json")), load(os.path.join(pd_, "clades.json")), load(os.path.join(pd_, "signatures_joint.json")), cosmic, chan_index, n_boot, n_sim, rng)
        per[P] = d
        tests[P] = t
    return per, tests


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--cosmic", required=True)
    ap.add_argument("--out", required=True, help="stage directory (created); live data is never written")
    ap.add_argument("--patients", default="BWH69,BWH70,BWH75,MGH285,MGH301,MGH302,MGH303,MGH304,MGH305")
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--n-sim", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    per, tests = run(a.data, a.cosmic, a.patients.split(","), a.n_boot, a.n_sim, a.seed)
    blocks, shifts = cohort_blocks(per, tests)
    for P, d in per.items():
        os.makedirs(os.path.join(a.out, P), exist_ok=True)
        with open(os.path.join(a.out, P, "signature_assignment.json"), "w") as f:
            json.dump({"format": "gos-sc-signature-assignment/1", "created": datetime.datetime.now().isoformat(timespec="seconds"), **d}, f, separators=(",", ":"))
    os.makedirs(os.path.join(a.out, "_cohort"), exist_ok=True)
    with open(os.path.join(a.out, "_cohort", "signature_assignment.json"), "w") as f:
        json.dump({"format": "gos-sc-signature-assignment-cohort/1", "created": datetime.datetime.now().isoformat(timespec="seconds"), "method": METHOD_TEXT, "n_boot": a.n_boot, "n_sim": a.n_sim,
                   "clade_signatures": blocks, "signature_shifts": shifts,
                   "unavailable": {P: d.get("reason") for P, d in per.items() if not d.get("available")}}, f, separators=(",", ":"))
    log("wrote", a.out)
    for s in shifts[:20]:
        log("  shift", s["patient"], s["clade"], "vs", s["vs"], s["signature"], s["diff"], f"q={s['q']:.3g}")


if __name__ == "__main__":
    main()
