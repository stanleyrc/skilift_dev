#!/usr/bin/env python3
"""Recovery of clinical-assay alterations (Oncopanel / Genexus) in the single-cell data, real counts only.

Inputs : data/_cohort/clinical_assays.json (gos_sc_clinical.py), clinical_sites.json (gos_sc_clinical_sites.py),
         data/<P>/{filtered.events.json,snv_matrix.json,tmb.json,<cell>/complex.json,<cell>/metadata.json},
         pooled driver scan (analysis/driver_sensitivity_261016/candidates_evaluated.tsv, coverage_gene.tsv),
         targeted pooled pileup of the clinical sites (analysis/clinical_recovery_261017/<P>/s*/{obs,depth}.tsv.gz,
         ds_pileup.py of the driver-sensitivity branch), data/_cohort/hotspot_rescue.json (TERT promoter).
Output : data/_cohort/clinical_recovery.json
   { format, params, patients: { P: { tumour_cells, clinical_status, snv:[...], cnv:[...], sv:[...], tert_promoter:[...],
                                       tmb:{...}, missed:[...] } } }
Status vocabulary: recovered | partial | not_recovered | not_covered | not_assessable (reason always given).
Rules (all fractions are real counts; thresholds are stated in params):
 SNV/indel recovered = >=2 tumour cells carry alt reads, 0 alt reads in normal cells (or normals not covered);
   partial = exactly 1 carrier cell; not_recovered = covered (pooled tumour depth >= 20 in >= 3 cells) but no tumour alt read;
   not_covered = pooled tumour depth < 20 or < 3 cells covered.
 CNV: per tumour cell gene CN from the genome graph (complex.json) vs the cell's own baseline (length-weighted median CN of the
   graph). gain: CN > baseline, amp: CN >= baseline + 4, loss: CN < baseline, deep deletion: CN == 0; event value = median over its genes.
   recovered >= 80% of tumour cells, partial 10-80%, not_recovered < 10%.
 SV: a graph ALT junction whose two breakends lie within SV_TOL of the two reported breakpoints; small intragenic events are
   below the graph's resolution (not_assessable); fusions: DNA fusion call with both gene names, else RNA fusion call.
usage: gos_sc_clinical_recovery.py [--data DIR]
"""
import argparse
import csv
import glob
import gzip
import json
import os
import re
import statistics
import datetime
from collections import defaultdict

DATA = "/gpfs/commons/projects/imielinski_web/sclarke/gOS_GBM_sc/data"
PROJ = "/gpfs/commons/groups/imielinski_lab/projects/single_cell_GBM"
SCAN = PROJ + "/analysis/driver_sensitivity_261016"
PILE = PROJ + "/analysis/clinical_recovery_261017"
GTF = "/gpfs/commons/groups/landau_lab/jzinno/ref/refdata-gex-GRCh38-2024-A/genes/genes.gtf"
MIN_DP, MIN_COV, AMP_DELTA, SV_TOL = 20, 3, 4, 250000
MIN_CARRIERS, MIN_VAF = 3, 0.02
FORMAT = "gos-sc-clinical-recovery/1"


# ---------- pure helpers (unit tested) ----------
def snv_status(carriers, covered_cells, dp, normal_alt=0, normal_dp=0, vaf=None):
    """-> (status, reason) from tumour carrier cells / covered cells / pooled tumour depth / pooled VAF."""
    if dp is None or covered_cells is None or dp < MIN_DP or covered_cells < MIN_COV:
        return "not_covered", "pooled tumour depth %s in %s cells (needs >= %d reads in >= %d cells)" % (dp, covered_cells, MIN_DP, MIN_COV)
    if carriers >= MIN_CARRIERS and vaf is not None and vaf >= MIN_VAF:
        return "recovered", "%d of %d covered tumour cells carry it (pooled VAF %.1f%%)" % (carriers, covered_cells, 100 * vaf)
    if carriers >= 1:
        return "partial", "weak: %d carrier cell(s) of %d covered, pooled VAF %.2f%% (below %d cells / %.0f%% VAF, near the error rate)" % (
            carriers, covered_cells, 100 * (vaf or 0), MIN_CARRIERS, 100 * MIN_VAF)
    return "not_recovered", "covered (%d reads, %d cells) with no tumour alt read" % (dp, covered_cells)


def germline_like(normal_alt, normal_dp):
    """alt allele also in normal cells at a heterozygous-like fraction (>= 20% of >= 20 normal reads)"""
    return bool(normal_dp and normal_dp >= 20 and normal_alt / normal_dp >= 0.2)


def cell_state(cn, base):
    """gain/amp/loss/deep for a cell: CN vs the cell's baseline."""
    if cn is None or base is None:
        return None
    s = set()
    if cn > base:
        s.add("gain")
    if cn >= base + AMP_DELTA:
        s.add("amp")
    if cn < base:
        s.add("loss")
    if cn == 0:
        s.add("deep")
    return s


def clinical_class(cnv_type):
    t = (cnv_type or "").lower()
    if "amplification" in t or t == "amp":
        return "amp"
    if "two copy" in t or "homozygous" in t:
        return "deep"
    if "gain" in t:
        return "gain"
    if "loss" in t or "deletion" in t:
        return "loss"
    return None


def cnv_status(frac):
    if frac is None:
        return "not_covered"
    return "recovered" if frac >= 0.8 else "partial" if frac >= 0.1 else "not_recovered"


def length_weighted_median(pairs):
    """[(cn, length)] -> weighted median cn"""
    pairs = sorted(pairs)
    tot = sum(l for _, l in pairs)
    acc = 0
    for cn, l in pairs:
        acc += l
        if acc >= tot / 2:
            return cn
    return None


def junction_match(junctions, bp1, bp2, tol=SV_TOL):
    """junctions: [(chrom1,pos1,chrom2,pos2)]; bp: (chrom,pos). -> min max-distance over junctions matching both breakpoints (either order)."""
    best = None
    for c1, p1, c2, p2 in junctions:
        for (a, b) in (((c1, p1), (c2, p2)), ((c2, p2), (c1, p1))):
            if a[0] == bp1[0] and b[0] == bp2[0]:
                d = max(abs(a[1] - bp1[1]), abs(b[1] - bp2[1]))
                if d <= tol and (best is None or d < best):
                    best = d
    return best


def gene_tokens(s):
    return set(re.split(r"::|-|\s|,", s or "")) - {""}


# ---------- data loading ----------
def load_genes(genes):
    pat = re.compile(r'gene_name "([^"]+)"')
    out = {}
    with open(GTF) as fh:
        for line in fh:
            if "\tgene\t" not in line:
                continue
            m = pat.search(line)
            if m and m.group(1) in genes and "gene_type \"protein_coding\"" in line:
                f = line.split("\t")
                out[m.group(1)] = (f[0].replace("chr", ""), int(f[3]), int(f[4]))
    return out


def tumour_cells(data, patient):
    cells = []
    for d in sorted(glob.glob(os.path.join(data, patient + "_*"))):
        mf = os.path.join(d, "metadata.json")
        if not os.path.exists(mf) or not os.path.exists(os.path.join(d, "complex.json")):
            continue
        m = json.load(open(mf))
        m = m[0] if isinstance(m, list) else m
        cells.append((os.path.basename(d), str(m.get("clone_id", "")).lower().startswith("normal")))
    return cells


def cell_graph(path):
    g = json.load(open(path))
    iv = {i["iid"]: i for i in g["intervals"]}
    segs = defaultdict(list)
    for i in iv.values():
        segs[str(i["chromosome"])].append((i["startPoint"], i["endPoint"], i["y"]))
    allp = [(i["y"], i["endPoint"] - i["startPoint"] + 1) for i in iv.values() if str(i["chromosome"]) not in ("X", "Y", "MT") and i["y"] is not None]
    base = length_weighted_median(allp)
    alt = []
    for c in g["connections"]:
        if c.get("type") != "ALT":
            continue
        ends = []
        for k in ("source", "sink"):
            n = c[k]
            i = iv.get(abs(n))
            if not i:
                break
            ends.append((str(i["chromosome"]), i["startPoint"] if n < 0 else i["endPoint"]))
        if len(ends) == 2:
            alt.append((ends[0][0], ends[0][1], ends[1][0], ends[1][1]))
    return segs, base, alt


def cn_at(segs, chrom, start, end):
    """max CN over the graph intervals overlapping the gene [start, end] (not the interval at the midpoint)"""
    ys = [y for s, e, y in segs.get(chrom, []) if s <= end and e >= start and y is not None]
    return max(ys) if ys else None


# ---------- per-source SNV evidence ----------
def read_pileup(patient):
    obs, depth = [], defaultdict(lambda: [0, 0, 0, 0])
    for d in glob.glob(os.path.join(PILE, patient, "s*")):
        if os.path.exists(os.path.join(d, "depth.tsv.gz")):
            with gzip.open(os.path.join(d, "depth.tsv.gz"), "rt") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    x = depth[(r["chrom"], int(r["pos"]))]
                    for k, f in enumerate(("dp_t", "cells_t", "dp_n", "cells_n")):
                        x[k] += int(r[f])
            with gzip.open(os.path.join(d, "obs.tsv.gz"), "rt") as fh:
                obs.extend(r for r in csv.DictReader(fh, delimiter="\t"))
    return obs, depth


def pileup_eval(site, obs, depth):
    """site: resolved clinical_sites row. -> evidence dict from the pooled pileup (real counts) or None."""
    if site.get("kind") == "SNV":
        opts = [{"pos": site["start"], "ref": site["ref"], "alt": site["alt"]}] + [
            {"pos": c["pos"], "ref": c["ref"], "alt": c["alt"]} for c in site.get("alt_candidates", [])]
    elif site.get("kind") == "SNV_PROTEIN":
        opts = site["options"]
    else:
        opts = None
    best = None
    if opts is not None:
        for o in opts:
            key = ("chr" + site["chrom"].replace("chr", ""), o["pos"])
            dp = depth.get((site["chrom"], o["pos"]), [0, 0, 0, 0])
            rows = [r for r in obs if r["chrom"] == site["chrom"] and int(r["pos"]) == o["pos"] and r["type"] == "SNV" and r["alt"] == o["alt"]]
            t = [r for r in rows if r["normal"] in ("False", "0")]
            n = [r for r in rows if r["normal"] in ("True", "1")]
            ev = {"pos": o["pos"], "ref": o["ref"], "alt": o["alt"], "tumour_alt": len(t), "tumour_dp": dp[0], "covered_cells": dp[1],
                  "carriers": len({r["cell"] for r in t}), "carrier_ids": sorted({r["cell"] for r in t}), "normal_alt": len(n), "normal_dp": dp[2],
                  "pooled_vaf": round(len(t) / dp[0], 4) if dp[0] else None}
            if best is None or ev["tumour_alt"] > best["tumour_alt"]:
                best = ev
    elif site.get("kind") == "DEL":  # same deleted length; the anchor may shift by repeat left-alignment (<= 8 bp)
        ln = site["end"] - site["start"] + 1
        rows = [r for r in obs if r["chrom"] == site["chrom"] and r["type"] == "DEL" and len(r["ref"]) - len(r["alt"]) == ln and
                abs(int(r["pos"]) - (site["start"] - 1)) <= 8]
        by = defaultdict(list)
        for r in rows:
            by[int(r["pos"])].append(r)
        if by:
            anchor = max(by, key=lambda k: len(by[k]))
            rows = by[anchor]
        else:
            anchor = site["start"] - 1
        t = [r for r in rows if r["normal"] in ("False", "0")]
        n = [r for r in rows if r["normal"] in ("True", "1")]
        dp = depth.get((site["chrom"], anchor), [0, 0, 0, 0])  # depth at the anchor base (the deleted reads keep it)
        best = {"pos": anchor, "tumour_alt": len(t), "tumour_dp": dp[0], "covered_cells": dp[1], "carriers": len({r["cell"] for r in t}),
                "carrier_ids": sorted({r["cell"] for r in t}),
                "normal_alt": len(n), "normal_dp": dp[2], "pooled_vaf": round(min(1.0, len(t) / dp[0]), 4) if dp[0] else None,
                "note": "VAF = deletion reads / base-level depth at the anchor, capped at 1; deletion of the same length (anchor %d, within 8 bp of the reported position)" % anchor}
    else:
        return None
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=DATA)
    a = ap.parse_args()
    data = a.data
    cl = json.load(open(os.path.join(data, "_cohort", "clinical_assays.json")))
    sites = json.load(open(os.path.join(data, "_cohort", "clinical_sites.json")))["patients"]
    hot = json.load(open(os.path.join(data, "_cohort", "hotspot_rescue.json")))["patients"]
    old_f = os.path.join(data, "_cohort", "clinical_sites.json.bak_261017")
    old_sites = json.load(open(old_f))["patients"] if os.path.exists(old_f) else {}
    scan = list(csv.DictReader(open(SCAN + "/candidates_evaluated.tsv"), delimiter="\t"))
    scan_cov = defaultdict(dict)
    for r in csv.DictReader(open(SCAN + "/coverage_gene.tsv"), delimiter="\t"):
        scan_cov[r["patient"]][r["gene"]] = r
    cnv_genes = {g for r in cl["patients"].values() for x in r["assays"] for c in x["cnvs"] for g in c["genes"]}
    cnv_genes |= {g for r in cl["patients"].values() for x in r["assays"] for s in x["svs"] for g in
                  [s.get("left", {}).get("feature", "").split(" ")[0], s.get("right", {}).get("feature", "").split(" ")[0]] if g}
    coords = load_genes(cnv_genes)
    panel_genes = {g for r in cl["patients"].values() for x in r["assays"] if x["assay"] == "oncopanel" for c in x["cnvs"] for g in c["genes"]}
    out = {}
    for p, rep in cl["patients"].items():
        res = {"assays": [{"assay": x["assay"], "status": x["status"]} for x in rep["assays"]]}
        cells = tumour_cells(data, p)
        tum = [c for c, nrm in cells if not nrm]
        res["tumour_cells"] = len(tum)
        res["normal_cells"] = len(cells) - len(tum)
        ok = [x for x in rep["assays"] if x["status"] == "ok"]
        res["clinical_status"] = "failed" if any(x["status"] == "failed" for x in rep["assays"]) else \
            "nothing detected" if any(x["status"] == "nothing detected" for x in rep["assays"]) and not ok else "ok" if ok else "none"
        events = json.load(open(os.path.join(data, p, "filtered.events.json")))
        smf = os.path.join(data, p, "snv_matrix.json")
        smv = {v["id"]: v for v in json.load(open(smf))["variants"]} if os.path.exists(smf) else None
        obs, depth = read_pileup(p)
        have_pileup = bool(depth)
        # ---- SNV / indel
        res["snv"] = []
        srows = sites.get(p, [])
        k = 0
        for x in rep["assays"]:
            for v in x["variants"]:
                site = srows[k]; k += 1
                row = {"assay": x["assay"], "gene": v["gene"], "cdna": v["cdna"], "protein": v["protein"], "tier": v["tier"],
                       "clinical_vaf": v["vaf"], "clinical_reads": v["reads"], "clinical_raw": v["raw"], "evidence": []}
                if not site.get("ok"):
                    row["status"], row["reason"] = "not_assessable", "no genomic position (%s)" % site.get("reason")
                    row["outcome"] = "unresolvable"
                    if v["gene"] == "TERT":
                        te = [e for e in events if e["gene"] == "TERT" and e.get("vartype") == "Hotspot"]
                        if te:
                            row["note"] = "single cells: %s in %s covered tumour cells (%s)" % (te[0]["Variant"].split(" ")[0], te[0].get("cells"), te[0].get("source") or "called")
                    if site.get("cdna"):
                        row["reason"] += "; the report's c./p. does not reproduce on any GENCODE protein-coding transcript (RefSeq not available), so no site can be pileup-checked"
                    if v["gene"] == "TERT":
                        row["reason"] += "; TERT promoter hotspots are assessed under tert_promoter"
                    res["snv"].append(row)
                    continue
                row["site"] = {k2: site[k2] for k2 in ("chrom", "start", "end", "kind", "transcript", "ref_check")}
                # filtered events
                ev = [e for e in events if e["gene"] == v["gene"] and (
                    (v["protein"] and e.get("Variant") == v["protein"]) or
                    (e.get("Genome_Location", "").startswith(site["chrom"].replace("chr", "") + ":") and e["type"] != "SCNA" and
                     site["start"] <= int(e["Genome_Location"].split(":")[1].split("-")[0]) <= site["end"] and e["vartype"] in ("SNV", "INDEL", "DEL", "INS")))]
                for e in ev[:1]:
                    row["evidence"].append({"source": "filtered events" + ("" if not e.get("source") else " (%s)" % e["source"]),
                                            "carriers": e.get("n_cells"), "covered": e.get("covered_tumour_cells") or e.get("tumour_cells"),
                                            "vaf": e.get("VAF") if e.get("VAF") is not None else (round(e["alt"] / (e["alt"] + e["ref"]), 4) if e.get("alt") is not None and (e["alt"] + e["ref"]) else None),
                                            "alt_reads": e.get("alt"), "ref_reads": e.get("ref"), "tier": e.get("Tier")})
                if smv is not None and site["kind"] == "SNV":
                    vid = "chr%s_%d_%s_%s" % (site["chrom"].replace("chr", ""), site["start"], site["ref"], site["alt"])
                    sv = smv.get(vid)
                    if sv:
                        row["evidence"].append({"source": "SNV matrix (CellPhy tree site)", "carriers": sv["alt_cells"], "category": sv["category"],
                                                "clade_f1": sv.get("clade_score")})
                sc = [r for r in scan if r["patient"] == p and r["gene"] == v["gene"] and (
                    (v["protein"] and r["protein"] == v["protein"]) or (r["chrom"] == site["chrom"] and int(r["pos"]) == site["start"]))]
                for r in sc[:1]:
                    row["evidence"].append({"source": "pooled driver scan", "carriers": int(r["carriers"]), "covered": int(r["covered"]),
                                            "vaf": float(r["scan_vaf_t"]) if r["scan_vaf_t"] else None, "clade_f1": float(r["clade_f1"]) if r["clade_f1"] else None})
                pe = pileup_eval(site, obs, depth) if have_pileup else None
                if have_pileup and pe is None:
                    row["status"], row["reason"] = "not_assessable", "complex indel (%s): allele not matched by the pooled pileup" % site["kind"]
                    row["source"] = "not assessed"
                    res["snv"].append(row)
                    continue
                if pe:
                    row["evidence"].append({"source": "targeted pooled pileup (all cells)", "carriers": pe["carriers"], "covered": pe["covered_cells"],
                                            "vaf": pe["pooled_vaf"], "alt_reads": pe["tumour_alt"], "depth": pe["tumour_dp"],
                                            "normal_alt": pe["normal_alt"], "note": pe.get("note")})
                # decision: prefer the pileup (uniform, all cells), else the other evidence
                if pe:
                    row["status"], row["reason"] = snv_status(pe["carriers"], pe["covered_cells"], pe["tumour_dp"], pe["normal_alt"], pe["normal_dp"], pe["pooled_vaf"])
                    row["germline_like"] = germline_like(pe["normal_alt"], pe["normal_dp"])
                    row["normal_alt"], row["normal_dp"] = pe["normal_alt"], pe["normal_dp"]
                    row["sc_vaf"], row["sc_carriers"], row["sc_covered"] = pe["pooled_vaf"], pe["carriers"], pe["covered_cells"]
                elif ev:
                    e = row["evidence"][0]
                    row["status"], row["reason"] = "recovered", "present in the filtered events (%s cells)" % e["carriers"]
                    row["sc_vaf"], row["sc_carriers"], row["sc_covered"] = e["vaf"], e["carriers"], e["covered"]
                else:
                    row["status"], row["reason"] = "not_covered", "no single-cell pileup for this patient (not available)"
                row["site"]["protein_verified"] = site.get("protein_verified")
                row["site"]["resolved_pos"] = (pe or {}).get("pos", site["start"])
                row["site"]["alt_positions_tested"] = site.get("alt_positions", [])
                row["source"] = ((ev[0].get("source") or "called SNV") if ev else "SNV matrix tree site" if any(e["source"].startswith("SNV matrix") for e in row["evidence"]) else
                                 "pooled driver scan" if sc else "targeted pileup" if pe and pe["carriers"] else "not found")
                for e in row["evidence"]:
                    if e.get("clade_f1") is not None:
                        row["clade_f1"] = e["clade_f1"]
                        break
                res["snv"].append(row)
        for row in res["snv"]:
            if "outcome" in row:
                continue
            old = (old_sites.get(p) or [])
            prev = next((o for o in old if o.get("gene") == row["gene"] and o.get("cdna") == row["cdna"] and o.get("protein") == row["protein"]), None)
            moved = bool(prev and prev.get("ok") and row.get("site") and prev["start"] != row["site"]["start"])
            row["outcome"] = (("recovered after correction (position moved from %s to %s)" % (prev["start"], row["site"]["start"]) if moved else "recovered")
                              if row["status"] == "recovered" else
                              "weak" if row["status"] == "partial" else
                              "truly absent (resolved site verified, covered, no tumour alt reads)" if row["status"] == "not_recovered" and (row["site"].get("protein_verified")) else
                              "absent at the resolved site (no protein change to verify the transcript; base matches)" if row["status"] == "not_recovered" else
                              row["status"])
            if moved:
                row["corrected_from"] = prev["start"]
            if row["clinical_vaf"] is not None and row["clinical_vaf"] <= 0.10 and row["status"] != "recovered":
                row["note"] = "possible clonal haematopoiesis (CD45- sorted cells would not carry it)"
        # ---- CNV
        res["cnv"] = []
        tcells = []
        for c in tum:
            segs, base, alt = cell_graph(os.path.join(data, c, "complex.json"))
            tcells.append((c, segs, base, alt))
        for x in rep["assays"]:
            items = [(c, False) for c in x["cnvs"]] + [({"location": None, "scope": "gene", "type": "Amplification" if n["what"].startswith("amp") else n["what"],
                                                          "genes": [n["gene"]], "copies": None, "raw": n["raw"]}, True) for n in x["negatives"]]
            for c, neg in items:
                cls = clinical_class(c["type"])
                gs = [g for g in c["genes"] if g in coords]
                row = {"assay": x["assay"], "location": c["location"], "scope": c["scope"], "type": c["type"], "genes": c["genes"],
                       "clinical_copies": c.get("copies"), "clinical_raw": c["raw"], "negative_finding": neg}
                if not tcells or not gs or cls is None:
                    row.update(status="not_assessable", reason="no tumour genome graphs" if not tcells else "gene coordinates / class unresolved")
                    res["cnv"].append(row)
                    continue
                sel = gs if len(gs) <= 25 else gs[:: max(1, len(gs) // 25)]
                vals, hits, carrier_vals = [], 0, []
                for cid, segs, base, alt in tcells:
                    cns = [cn_at(segs, *coords[g]) for g in sel]
                    cns = [v for v in cns if v is not None]
                    if not cns:
                        continue
                    cn = statistics.median(cns)
                    vals.append(cn)
                    st = cell_state(cn, base)
                    if st and cls in st:
                        hits += 1
                        carrier_vals.append(cn)
                n = len(vals)
                frac = hits / n if n else None
                row.update(n_cells=n, cells_with_class=hits, frac_cells=round(frac, 4) if frac is not None else None,
                           mean_cn=round(statistics.mean(vals), 2) if vals else None, median_cn=statistics.median(vals) if vals else None,
                           max_cn=max(vals) if vals else None, genes_used=sel, rule=cls,
                           carrier_median_cn=statistics.median(carrier_vals) if carrier_vals else None,
                           carrier_q25_cn=sorted(carrier_vals)[len(carrier_vals) // 4] if carrier_vals else None,
                           carrier_q75_cn=sorted(carrier_vals)[(3 * len(carrier_vals)) // 4] if carrier_vals else None,
                           cn_method="per-cell max graph CN over the gene interval; event value = median over its genes")
                if cls == "deep":
                    hd = [e for e in events if e["type"] == "SCNA" and e["vartype"] == "HOMDEL" and e["gene"] in c["genes"]]
                    row["homdel_events"] = [{"gene": e["gene"], "n_cells": e["n_cells"], "broad_cn0": e.get("broad_cn0"),
                                             "cn0_segment_mb_median": e.get("cn0_segment_mb_median")} for e in hd]
                    row["broad_cn0"] = bool(hd) and all(e.get("broad_cn0") for e in hd)
                if neg:
                    row["status"] = "concordant_negative" if (frac is not None and frac < 0.1) else "discordant_negative"
                    row["reason"] = "report states no %s; %d of %d tumour cells show it" % (c["genes"][0], hits, n)
                else:
                    row["status"] = cnv_status(frac)
                    if row.get("broad_cn0") and row["status"] == "recovered":
                        row["status"] = "partial"
                    row["reason"] = "%d of %d tumour cells with %s (median CN %s, carrier cells %s%s)%s" % (
                        hits, n, cls, row["median_cn"], row["carrier_median_cn"], ", clinical estimate %d copies" % c["copies"] if c.get("copies") else "",
                        "; the CN-0 calls are broad (chromosome/arm-scale), not evidence of a focal homozygous deletion" if row.get("broad_cn0") else "")
                res["cnv"].append(row)
        # ---- SV
        res["sv"] = []
        juncs = defaultdict(int)
        for cid, segs, base, alt in tcells:
            for j in set(alt):
                juncs[j] += 1
        jl = list(juncs)
        for x in rep["assays"]:
            for s in x["svs"]:
                row = {"assay": x["assay"], "kind": s["kind"], "tier": s.get("tier"), "clinical_raw": s["raw"], "left": s.get("left"), "right": s.get("right")}
                lg, rg = (s["left"] or {}).get("feature", ""), (s["right"] or {}).get("feature", "")
                if "chrom" in (s["left"] or {}):
                    c1, p1 = s["left"]["chrom"].replace("chr", ""), s["left"]["pos"]
                    c2, p2 = s["right"]["chrom"].replace("chr", ""), s["right"]["pos"]
                    if c1 == c2 and abs(p1 - p2) < 20000:
                        row.update(status="not_assessable", reason="intragenic event of %d bp, below the genome-graph resolution" % abs(p1 - p2))
                    else:
                        d = junction_match(jl, (c1, p1), (c2, p2))
                        if d is None:
                            row.update(status="not_recovered", reason="no graph ALT junction within %d kb of both breakpoints in any tumour cell" % (SV_TOL // 1000))
                        else:
                            m = sum(v for (a1, b1, a2, b2), v in juncs.items() if junction_match([(a1, b1, a2, b2)], (c1, p1), (c2, p2)) is not None)
                            row.update(status="recovered" if m >= 0.1 * len(tcells) else "partial", n_cells=m, nearest_bp_distance=d,
                                       reason="graph junction within %d bp of both breakpoints in %d of %d tumour cells (cell-junction pairs counted per cell)" % (d, min(m, len(tcells)), len(tcells)))
                else:
                    g1, g2 = lg.split(" ")[0], rg.split(" ")[0]
                    fus = [e for e in events if e["type"] == "Fusion" and {g1, g2} <= gene_tokens(e["fusion_genes"])]
                    one = [e for e in events if e["type"] in ("Fusion", "Complex SV") and ({g1, g2} & gene_tokens(e["fusion_genes"] or e["gene"]))]
                    d = None
                    if g1 in coords and g2 in coords and coords[g1][0] == coords[g2][0]:
                        d = junction_match(jl, (coords[g1][0], (coords[g1][1] + coords[g1][2]) // 2), (coords[g2][0], (coords[g2][1] + coords[g2][2]) // 2), SV_TOL)
                    rf = os.path.join(data, p, "rna", "fusions.json")
                    rna_hit = None
                    if os.path.exists(rf):
                        for f in json.load(open(rf)).get("fusions", []):
                            if {g1, g2} <= set(re.findall(r'"([A-Za-z0-9.\-]+)"', json.dumps(f))):
                                rna_hit = f
                                break
                    row["rna_checked"] = os.path.exists(rf)
                    if fus:
                        row.update(status="recovered", n_cells=fus[0]["n_cells"], reason="DNA fusion call %s in %s cells" % (fus[0]["fusion_genes"], fus[0]["n_cells"]))
                    elif d is not None:
                        m = sum(v for (a1, b1, a2, b2), v in juncs.items() if junction_match([(a1, b1, a2, b2)], (coords[g1][0], (coords[g1][1] + coords[g1][2]) // 2), (coords[g2][0], (coords[g2][1] + coords[g2][2]) // 2)) is not None)
                        row.update(status="partial", n_cells=m, nearest_bp_distance=d, reason="no fusion call, but a graph junction joins the %s and %s loci (within %d bp) in %d tumour cells" % (g1, g2, d, min(m, len(tcells))))
                    elif rna_hit:
                        row.update(status="partial", reason="RNA fusion call joins %s and %s (no DNA fusion call)" % (g1, g2))
                    elif one:
                        row.update(status="partial", reason="only one partner gene appears in DNA fusion calls (%s)" % one[0]["fusion_genes"])
                    else:
                        row.update(status="not_recovered", reason="no DNA fusion call, graph junction%s joining %s and %s" % (" or RNA fusion call" if row["rna_checked"] else " (RNA not available)", g1, g2))
                res["sv"].append(row)
        # ---- TERT promoter
        res["tert_promoter"] = []
        hp = {s["name"]: s for s in (hot.get(p, {}).get("sites") or [])}
        for x in rep["assays"]:
            for nm, t in x["tert_promoter"].items():
                s = hp.get("TERT " + nm)
                ev = next((e for e in events if e["gene"] == "TERT" and e.get("vartype") == "Hotspot" and str(e.get("Variant", "")).startswith(nm)), None)
                row = {"assay": x["assay"], "site": nm, "clinical": t["status"], "clinical_raw": t["raw"]}
                if s and s.get("level") is not None:
                    sc_present = s["level"] in ("supported", "present", "likely")
                    row.update(source="hotspot_rescue.json", sc_level=s["level"], carriers=s.get("carriers"), covered=s.get("covered"), alt=s.get("alt"), dp=s.get("dp"),
                               status=("recovered" if sc_present else "not_recovered") if t["status"] == "present" else
                               ("concordant_negative" if s["level"] == "absent" else "discordant_negative"),
                               reason=s.get("verdict_summary"))
                elif ev is not None and t["status"] == "present":
                    he = (ev.get("hotspot_evidence") or {}).get("tumor") or {}
                    row.update(source="filtered events (%s)" % (ev.get("source") or "called"), carriers=ev.get("n_cells"), covered=ev.get("covered_tumour_cells"),
                               alt=he.get("alt", ev.get("alt")), dp=he.get("dp"), status="recovered",
                               reason="%s carrier cells of %s covered tumour cells, %s alt reads (%s)" % (ev.get("n_cells"), ev.get("covered_tumour_cells"), he.get("alt", ev.get("alt")), ev.get("source") or "called"))
                elif ev is not None:
                    row.update(source="filtered events", status="discordant_negative", reason="report says absent but single cells have %s carrier cells" % ev.get("n_cells"))
                else:
                    row.update(status="not_assessable", reason=(s or {}).get("reason") or "no TERT promoter assessment in hotspot_rescue.json or the filtered events")
                if row["status"] == "concordant_negative":
                    row["label"] = "%s absent in both (clinical and single-cell); not a TERT-negative call" % nm
                res["tert_promoter"].append(row)
        # ---- TMB
        tf = os.path.join(data, p, "tmb.json")
        tmb = {"clinical": next((x["tmb_per_mb"] for x in rep["assays"] if x.get("tmb_per_mb") is not None), None)}
        if os.path.exists(tf):
            ps = json.load(open(tf))["patient_summary"]
            tmb.update(sc_truncal=ps["truncal"]["tmb_corrected"], sc_subclonal=ps["subclonal"]["tmb_corrected"], sc_private=ps["private"]["tmb_corrected"],
                       sc_cell_median=ps["cell_tmb_corrected"]["median"], sc_unit="mutations per Mb of callable diploid territory (tmb.json)")
        else:
            tmb["note"] = "not available (no tmb.json)"
        res["tmb"] = tmb
        # ---- inverse: single-cell drivers (tier 1-2, >=30% of tumour cells) not reported clinically
        reported = {v["gene"] for x in ok for v in x["variants"]} | {g for x in ok for c in x["cnvs"] for g in c["genes"]} | \
                   {t for x in ok for s in x["svs"] for ft in (s.get("left", {}).get("feature", ""), s.get("right", {}).get("feature", "")) for t in [ft.split(" ")[0]]}
        miss, seen = [], set()
        if ok and tum:
            for e in events:
                try:
                    tier = int(e.get("Tier"))
                except (TypeError, ValueError):
                    continue
                if tier > 2 or not e["gene"] or (e["n_cells"] or 0) < 0.3 * len(tum) or e["type"] == "Complex SV":
                    continue
                gs = gene_tokens(e["gene"]) if e["type"] == "Fusion" else {e["gene"]}
                key = (e["gene"], e["type"], e.get("Variant"))
                if key in seen or gs & reported:
                    continue
                seen.add(key)
                miss.append({"gene": e["gene"], "type": e["type"], "variant": e.get("Variant"), "tier": tier, "n_cells": e["n_cells"],
                             "cell_fraction": round(e["n_cells"] / len(tum), 3), "location": e.get("Genome_Location"),
                             "panel": "gene appears on the Oncopanel CNV gene lists" if e["gene"] in panel_genes else
                             "not in any Oncopanel report gene list (panel membership unknown)"
                             if any(x["assay"] == "oncopanel" for x in ok) else "Genexus panel gene list not in the report"})
        res["missed"] = miss
        out[p] = res
    dest = os.path.join(data, "_cohort", "clinical_recovery.json")
    if os.path.exists(dest) and not os.path.exists(dest + ".bak_261017"):
        os.replace(dest, dest + ".bak_261017")
    json.dump({"format": FORMAT, "created": datetime.datetime.now().isoformat(timespec="seconds"),
               "params": {"min_dp": MIN_DP, "min_covered_cells": MIN_COV, "amp_delta": AMP_DELTA, "sv_tolerance_bp": SV_TOL, "snv_recovered_min_carriers": MIN_CARRIERS, "snv_recovered_min_pooled_vaf": MIN_VAF,
                          "sources": ["clinical_assays.json", "clinical_sites.json", "filtered.events.json", "snv_matrix.json", "driver_sensitivity_261016 pooled scan",
                                      "clinical_recovery_261017 targeted pileup", "per-cell complex.json", "hotspot_rescue.json", "tmb.json"]},
               "patients": out}, open(dest, "w"), indent=1)
    c = defaultdict(lambda: defaultdict(int))
    for p, r in out.items():
        for kind in ("snv", "cnv", "sv", "tert_promoter"):
            for x in r[kind]:
                c[kind][x["status"]] += 1
    print({k: dict(v) for k, v in c.items()})


if __name__ == "__main__":
    main()
