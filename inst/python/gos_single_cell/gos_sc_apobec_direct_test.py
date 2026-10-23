#!/usr/bin/env python3
"""Direct test of an APOBEC (SBS13 = C>G at TCW; SBS2 = C>T at TCW) signal in clade-private SNVs -> <out>/_cohort/apobec_direct_test.json

No signature fitting: counts the SBS96 channels that define the signature, compares them with the truncal set and with the other
clades' private SNVs, then checks whether the supporting reads look like real somatic heterozygous SNVs or like a single-cell
amplification (PTA) / batch artefact.

For each (patient, clade) with >= MIN_SITES private SNVs:
  counts       TCW C>G (and TCW C>T) among the clade's private SNVs vs the truncal set vs other clades' private SNVs (Fisher, one-sided),
               and vs the number expected from the patient-level signature mixture (gos_sc_signature_assignment posteriors)
  normals      alt reads / cells with alt reads in the patient's normal cells at the clade's TCW C>G sites vs its other private sites
  read support pooled VAF and alt reads per supporting cell of the TCW C>G sites vs the clade's other private sites (Mann-Whitney)
  plates       supporting cells (call == 1) of the TCW C>G sites per plate vs the clade's cells per plate; chi-square homogeneity of
               the share of alt observations that are TCW C>G across plates; the plate holding the most supporting cells
Verdict (stated rules): 'supported' only if the TCW C>G excess over truncal AND over other clades is significant (larger one-sided Fisher p, BH q < 0.05 over every clade tested),
the normals do not carry more alt reads at these sites than at the clade's other sites, the sites' VAF is not lower than the other
private sites (Mann-Whitney p >= 0.05 or higher), and supporting cells are not concentrated in one plate beyond that plate's
share of the clade (binomial p >= 0.01); otherwise 'not supported' with the failing criteria listed.

Usage: python3 gos_sc_apobec_direct_test.py --data <gOS data dir> --out <stage dir> --patients BWH70,MGH305
"""
import argparse
import datetime
import json
import os
import re
import sys

import numpy as np
from scipy import stats

MIN_SITES = 50
TCW = ("A", "T")


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def load(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


CTX = re.compile(r"^([ACGT])\[([ACGT])>([ACGT])\]([ACGT])$")


def classify(ctx):
    """'cg_tcw' (SBS13 channel), 'ct_tcw' (SBS2 channel) or None."""
    m = CTX.match(ctx or "")
    if not m:
        return None
    five, ref, alt, three = m.groups()
    if ref == "C" and five == "T" and three in TCW:
        if alt == "G":
            return "cg_tcw"
        if alt == "T":
            return "ct_tcw"
    return None


def plate_of(cell, patient):
    """(run, plate) label from a cell name such as BWH70_MR_1_pl1_10a or pool-Q-BWH70-MR-B2-CD45pos-pl2-4c."""
    m = re.match(r"^(.*?)[_-](pl\d+)[_-]([A-Za-z0-9]+)$", cell)
    if not m:
        return "unparsed"
    run = re.sub(r"^(pool[_-]Q[_-])?", "", m.group(1))
    run = re.sub(rf"^{patient}[_-]?", "", run, flags=re.I)
    return f"{run}|{m.group(2)}" if run else m.group(2)


def bh(ps):
    m = len(ps)
    order = sorted(range(m), key=lambda i: ps[i])
    q, prev = [1.0] * m, 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        prev = min(prev, ps[i] * m / (rank + 1))
        q[i] = prev
    return [float(min(1, x)) for x in q]


def finalize_verdicts(res, alpha=0.05):
    """Enrichment criterion: the larger of the two one-sided p (vs truncal, vs other clades) per clade, BH across every clade tested."""
    rs = [r for r in res if r.get("available")]
    ps = []
    for r in rs:
        c = r["cg_tcw"]
        a, b = c["truncal"]["p_greater"], c["other_clades"]["p_greater"]
        ps.append(1.0 if a is None or b is None else max(a, b))
    for r, p, q in zip(rs, ps, bh(ps)):
        r["cg_tcw"]["p_enrichment"] = p
        r["cg_tcw"]["q_enrichment"] = q
        v = r["verdict"]
        v["failed"] = [f for f in v["failed"] if not f.startswith("no significant TCW C>G")]
        v["passed"] = []
        c = r["cg_tcw"]
        if q < alpha:
            v["passed"].append(f"TCW C>G enriched: {c['k']}/{r['n_private_snvs']} = {c['frac']:.3f} vs truncal {c['truncal']['frac']:.3f} (p={c['truncal']['p_greater']:.2g}) and other clades {c['other_clades']['frac']:.3f} (p={c['other_clades']['p_greater']:.2g}); q={q:.2g} over {len(rs)} clades")
        else:
            v["failed"].insert(0, f"no significant TCW C>G enrichment: {c['k']}/{r['n_private_snvs']} = {c['frac']:.3f} vs truncal {c['truncal']['frac']} and other clades {c['other_clades']['frac']}; q={q:.2g} over {len(rs)} clades")
        v["status"] = "supported" if not v["failed"] else "not supported"
    return res


def fisher_greater(k1, n1, k2, n2):
    if n1 == 0 or n2 == 0:
        return None
    return float(stats.fisher_exact([[k1, n1 - k1], [k2, n2 - k2]], alternative="greater")[1])


def site_reads(cells_doc, site_pos, cell_set):
    """For each site position, (alt_reads, depth, n_alt_cells, n_call_cells, n_cells_covered) pooled over cell_set."""
    acc = {p: [0, 0, 0, 0, 0] for p in site_pos}
    per_cell_alt = {p: [] for p in site_pos}
    for c in cell_set:
        for e in cells_doc.get(c, []):
            p = e[0]
            if p not in acc:
                continue
            ref, alt, call = e[1] or 0, e[2] or 0, e[3]
            a = acc[p]
            a[0] += alt
            a[1] += ref + alt
            a[4] += 1 if (ref + alt) > 0 else 0
            if alt > 0:
                a[2] += 1
                per_cell_alt[p].append(alt)
            if call == 1:
                a[3] += 1
    return acc, per_cell_alt


def test_clade(P, clade, variants, vpos, cells_doc, normals, all_cells, truncal_ids, other_private_ids, expected_frac):
    ids = [i for i in (clade.get("private_snvs") or []) if i in vpos]
    cls = {i: classify(variants[vpos[i]].get("context")) for i in ids}
    cg = [i for i in ids if cls[i] == "cg_tcw"]
    ct = [i for i in ids if cls[i] == "ct_tcw"]
    other = [i for i in ids if cls[i] is None]
    tr_cls = [classify(variants[vpos[i]].get("context")) for i in truncal_ids if i in vpos]
    ot_cls = [classify(variants[vpos[i]].get("context")) for i in other_private_ids if i in vpos]
    n = len(ids)
    out = {"patient": P, "clade": clade["id"], "n_cells": clade.get("n_cells"), "n_private_snvs": n}
    for key, lab in (("cg_tcw", "SBS13 channel: TCW C>G"), ("ct_tcw", "SBS2 channel: TCW C>T")):
        k = sum(1 for i in ids if cls[i] == key)
        ktr, ntr = sum(1 for c in tr_cls if c == key), len(tr_cls)
        kot, not_ = sum(1 for c in ot_cls if c == key), len(ot_cls)
        out[key] = {"label": lab, "k": k, "n": n, "frac": round(k / n, 4) if n else None,
                    "truncal": {"k": ktr, "n": ntr, "frac": round(ktr / ntr, 4) if ntr else None, "p_greater": fisher_greater(k, n, ktr, ntr)},
                    "other_clades": {"k": kot, "n": not_, "frac": round(kot / not_, 4) if not_ else None, "p_greater": fisher_greater(k, n, kot, not_)}}
    out["cg_tcw"]["expected_if_patient_wide_rate"] = round(expected_frac * n, 2) if expected_frac is not None else None
    out["cg_tcw"]["site_ids"] = cg

    clade_cells = [c for c in (clade.get("cells") or []) if c in cells_doc]
    norm_cells = [c for c in normals if c in cells_doc]

    def pooled(siteids, cellset):
        pos = [vpos[i] for i in siteids]
        acc, pca = site_reads(cells_doc, pos, cellset)
        alt = sum(a[0] for a in acc.values())
        dep = sum(a[1] for a in acc.values())
        altcells = sum(a[2] for a in acc.values())
        callcells = sum(a[3] for a in acc.values())
        return acc, pca, {"alt_reads": alt, "depth": dep, "vaf": round(alt / dep, 4) if dep else None, "alt_cell_site_obs": altcells, "call_cell_site_obs": callcells}

    acc_cg, pca_cg, s_cg = pooled(cg, clade_cells) if cg else (None, None, None)
    acc_ot, pca_ot, s_ot = pooled(other, clade_cells) if other else (None, None, None)
    out["clade_reads"] = {"cg_tcw": s_cg, "other_private": s_ot}
    # per-site VAF (pooled over the clade's cells) comparison
    if acc_cg and acc_ot:
        v_cg = [a[0] / a[1] for a in acc_cg.values() if a[1] > 0]
        v_ot = [a[0] / a[1] for a in acc_ot.values() if a[1] > 0]
        if len(v_cg) >= 3 and len(v_ot) >= 3:
            u = stats.mannwhitneyu(v_cg, v_ot, alternative="less")
            out["vaf_per_site"] = {"cg_tcw_median": round(float(np.median(v_cg)), 4), "other_median": round(float(np.median(v_ot)), 4), "p_cg_lower": float(u.pvalue)}
        # carriers: reads per alt cell, fraction of alt cells with >= 2 alt reads
        def multi(pca):
            xs = [x for v in pca.values() for x in v]
            return {"n_alt_cell_obs": len(xs), "median_alt_reads": float(np.median(xs)) if xs else None, "frac_ge2": round(sum(1 for x in xs if x >= 2) / len(xs), 3) if xs else None}
        out["alt_reads_per_cell"] = {"cg_tcw": multi(pca_cg), "other_private": multi(pca_ot)}
    # normals
    if norm_cells and cg:
        _, _, n_cg = pooled(cg, norm_cells)
        _, _, n_ot = pooled(other, norm_cells) if other else (None, None, None)
        rate = lambda s: (s["alt_reads"] / s["depth"]) if s and s["depth"] else None
        out["normals"] = {"n_normal_cells": len(norm_cells), "cg_tcw": n_cg, "other_private": n_ot,
                          "alt_read_rate_cg": rate(n_cg), "alt_read_rate_other": rate(n_ot)}
        if n_cg and n_ot and n_cg["depth"] and n_ot["depth"]:
            out["normals"]["p_cg_higher"] = float(stats.fisher_exact([[n_cg["alt_reads"], n_cg["depth"] - n_cg["alt_reads"]], [n_ot["alt_reads"], n_ot["depth"] - n_ot["alt_reads"]]], alternative="greater")[1])
    else:
        out["normals"] = {"n_normal_cells": len(norm_cells), "note": "not available (no normal cells with reads)" if not norm_cells else "no TCW C>G sites"}
    # plates
    if cg and acc_cg is not None:
        plates_all = {}
        for c in clade_cells:
            plates_all.setdefault(plate_of(c, P), []).append(c)
        sup_by_plate, alt_cg_by_plate, alt_ot_by_plate = {}, {}, {}
        for pl, cs in plates_all.items():
            a1, _, s1 = pooled(cg, cs)
            sup_by_plate[pl] = s1["call_cell_site_obs"]
            alt_cg_by_plate[pl] = s1["alt_cell_site_obs"]
            if other:
                _, _, s2 = pooled(other, cs)
                alt_ot_by_plate[pl] = s2["alt_cell_site_obs"]
        tot_sup = sum(sup_by_plate.values())
        rows = []
        for pl, cs in plates_all.items():
            rows.append({"plate": pl, "clade_cells": len(cs), "supporting_cell_site_obs": sup_by_plate[pl], "alt_cell_site_obs_cg": alt_cg_by_plate[pl], "alt_cell_site_obs_other": alt_ot_by_plate.get(pl)})
        rows.sort(key=lambda r: -r["supporting_cell_site_obs"])
        out["plates"] = {"rows": rows}
        if tot_sup > 0 and len(plates_all) > 1:
            top = rows[0]
            share_cells = top["clade_cells"] / len(clade_cells)
            out["plates"]["top_plate"] = {"plate": top["plate"], "share_of_supporting": round(top["supporting_cell_site_obs"] / tot_sup, 3), "share_of_clade_cells": round(share_cells, 3),
                                          "p_binom_excess": float(stats.binomtest(top["supporting_cell_site_obs"], tot_sup, share_cells, alternative="greater").pvalue)}
            tab = [[r["alt_cell_site_obs_cg"], r["alt_cell_site_obs_other"]] for r in rows if r["alt_cell_site_obs_other"] is not None]
            tab = np.array(tab, float)
            tab = tab[tab.sum(axis=1) > 0]
            if tab.shape[0] > 1 and (tab.sum(axis=0) > 0).all():
                chi = stats.chi2_contingency(tab)
                out["plates"]["homogeneity_chi2_p"] = float(chi[1])
    # verdict
    fails, passes = [], []
    c13 = out["cg_tcw"]
    pt, po = c13["truncal"]["p_greater"], c13["other_clades"]["p_greater"]
    if pt is not None and po is not None and pt < 0.05 and po < 0.05:
        passes.append(f"TCW C>G enriched ({c13['k']}/{n} = {c13['frac']:.3f}) vs truncal ({c13['truncal']['frac']:.3f}, p={pt:.3g}) and vs other clades ({c13['other_clades']['frac']:.3f}, p={po:.3g})")
    else:
        fails.append(f"no significant TCW C>G enrichment: {c13['k']}/{n} = {c13['frac']:.3f} vs truncal {c13['truncal']['frac']} (p={pt if pt is None else round(pt, 3)}) and other clades {c13['other_clades']['frac']} (p={po if po is None else round(po, 3)})")
    nm = out["normals"]
    if nm.get("p_cg_higher") is not None and nm["p_cg_higher"] < 0.05 and (nm["alt_read_rate_cg"] or 0) > (nm["alt_read_rate_other"] or 0):
        fails.append(f"normal cells carry more alt reads at the TCW C>G sites ({nm['cg_tcw']['alt_reads']}/{nm['cg_tcw']['depth']}) than at the clade's other sites ({nm['other_private']['alt_reads']}/{nm['other_private']['depth']}), p={nm['p_cg_higher']:.3g}")
    vp = out.get("vaf_per_site")
    if vp and vp["p_cg_lower"] < 0.05:
        fails.append(f"TCW C>G sites have lower VAF than the clade's other sites (median {vp['cg_tcw_median']} vs {vp['other_median']}, p={vp['p_cg_lower']:.3g})")
    tp = (out.get("plates") or {}).get("top_plate")
    if tp and tp["p_binom_excess"] < 0.01:
        fails.append(f"supporting cells concentrated on plate {tp['plate']} ({tp['share_of_supporting']:.0%} of support vs {tp['share_of_clade_cells']:.0%} of the clade's cells, p={tp['p_binom_excess']:.2g})")
    ph = (out.get("plates") or {}).get("homogeneity_chi2_p")
    if ph is not None and ph < 0.01:
        fails.append(f"share of alt observations that are TCW C>G differs between plates (chi-square p={ph:.2g})")
    out["verdict"] = {"status": "supported" if not fails else "not supported", "failed": fails, "passed": passes}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--patients", default="BWH70,MGH305")
    ap.add_argument("--clades", default=None, help="comma list of clade ids to restrict to (default: every clade with >= %d private SNVs)" % MIN_SITES)
    a = ap.parse_args()
    only = set(a.clades.split(",")) if a.clades else None
    res = []
    for P in a.patients.split(","):
        pd_ = os.path.join(a.data, P)
        snv, cl = load(os.path.join(pd_, "snv_matrix.json")), load(os.path.join(pd_, "clades.json"))
        rc = load(os.path.join(pd_, "precompute", "region_calls.json"))
        if not (snv and cl):
            res.append({"patient": P, "available": False, "reason": "no snv_matrix.json / clades.json"})
            continue
        variants = snv["variants"]
        vpos = {v["id"]: n for n, v in enumerate(variants)}
        normals = [c for c, k in zip(rc["cells"], rc["clone"]) if re.match(r"(?i)^normal$", f"{k}")] if rc else []
        truncal = [v["id"] for v in variants if v.get("category") == "truncal" and v.get("cellphy_input")]
        private_by_clade = {c["id"]: set(c.get("private_snvs") or []) for c in cl["clades"]}
        allc = [classify(v.get("context")) for v in variants if v.get("context")]
        exp_frac = sum(1 for c in allc if c == "cg_tcw") / len(allc) if allc else None   # null: no clade-specific enrichment
        for c in cl["clades"]:
            if only and c["id"] not in only and c["id"].split("_")[-1] not in only:
                continue
            if len(c.get("private_snvs") or []) < MIN_SITES or c.get("is_trunk") or c.get("parent") is None:
                continue
            others = set().union(*[s for k, s in private_by_clade.items() if k != c["id"]]) if len(private_by_clade) > 1 else set()
            log(P, c["id"])
            r = test_clade(P, c, variants, vpos, snv["cells"], normals, list(snv["cells"]), truncal, sorted(others), exp_frac)
            r["available"] = True
            res.append(r)
    finalize_verdicts(res)
    os.makedirs(os.path.join(a.out, "_cohort"), exist_ok=True)
    doc = {"format": "gos-sc-apobec-direct-test/1", "created": datetime.datetime.now().isoformat(timespec="seconds"), "min_sites": MIN_SITES,
           "definitions": {"SBS13_channel": "TCW C>G (T[C>G]A, T[C>G]T)", "SBS2_channel": "TCW C>T (T[C>T]A, T[C>T]T)", "supporting_cell": "call == 1 in snv_matrix.json cells",
                           "expected_note": "expected_if_patient_wide_rate is the patient-wide TCW C>G fraction over all of the patient's mutations x clade size (a null of no clade-specific enrichment)"},
           "results": res}
    with open(os.path.join(a.out, "_cohort", "apobec_direct_test.json"), "w") as f:
        json.dump(doc, f, separators=(",", ":"), default=lambda o: float(o) if isinstance(o, np.floating) else int(o) if isinstance(o, np.integer) else str(o))
    for r in res:
        if r.get("available"):
            log(r["patient"], r["clade"], r["verdict"]["status"], "; ".join(r["verdict"]["failed"])[:300])


if __name__ == "__main__":
    main()
