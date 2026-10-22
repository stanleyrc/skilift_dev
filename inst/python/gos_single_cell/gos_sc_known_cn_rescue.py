#!/usr/bin/env python3
"""Tier K ("known biology") rescue of driver SNVs in copy-number-altered regions, with dropout- and CN-aware clonality.

Why: the tree-building site filter only takes sites that are confidently genotyped in CN-2 cells (srctools select_cellphy_sites needs
allelic_total_cn == 2; the site-QC pool needs diploid_fraction >= 0.80), so clonal drivers on a gained or lost arm (e.g. PIK3CA E545A on a
clonal 3q gain) never reach the SNV matrix. Known drivers are checked with priors, so they are rescued here from read-level evidence; they
feed the driver / oncoprint / Known-biology views and are never added to the tree-building site set.

Two parts per patient:
  1. Re-annotate every read-evidence driver row already in filtered.events.json (hotspot rescue, pooled driver scan, called SNV not a
     tree site, earlier runs of this script) with the CN-aware clonality (clonality.cn_aware_clonality) and the site's CN state.
  2. Sweep the patient's genotyping VCF (per-cell AD_BAM recounts; --vcf-manifest) over the known-driver gene spans (--known, made from the
     project driver list) and add the known-driver sites that pass read QC and are not events yet (source KNOWN_CN_SOURCE).
Known-driver rule for a swept site (biallelic SNV; ANNOVAR fields of the VCF): TERT promoter C228T / C250T; a configured hotspot residue
(gbm_hotspots, e.g. PIK3CA E545); stop-gain / stop-loss / start-loss / splice-site in a driver gene; or a protein-altering change with
>= COSMIC_MIN COSMIC (cosmic70) occurrences. Anything else in a CN-altered region is Tier N (novel) and is not touched here.
Read QC: pooled normal cells >= MIN_NORMAL_DP reads with 0 alt reads; gnomAD (gnomad40_genome_AF) <= MAX_POP_AF; >= MIN_CARRIERS
tumour cells with an alt read and >= MIN_ALT alt reads; alt reads on both strands when >= MIN_ALT. Sequoia's branch_filter = "germline" is
reported but not used (it labels every trunk variant germline).
A swept site is added only when its CN-aware clonality is clonal, or subclonal with best tree-clade F1 >= MIN_SUBCLONAL_F1; indeterminate
sites and splice-site changes (RNA carry-over is not checked from VCF counts) are reported only.
Rows get: evidence = "known driver (rescued)", known_basis, site_cn {modal_cn, cn2_fraction, counts, cn_altered}, pipeline_status
(tree site / site-filter stage), clonality (CN-aware; previous call kept in clonality_dropout_only) and clonality_stats.cn_aware = True.

usage: gos_sc_known_cn_rescue.py --data <gOS data> --out <staging data dir> --known known_drivers.json --vcf-manifest manifest.tsv
       [--patients A,B] [--tabix tabix] [--tag bak_<yymmdd>] [--dry]
  --out: rewritten files go there (data/<P>/filtered.events.json, <P>/known_cn_rescue.json, _cohort/known_cn_rescue.json); --out == --data
         writes in place after a <file>.<tag> backup. The manifest is a TSV with columns patient, vcf and optionally site_candidates,
         pool_sites (site-filter tables, for pipeline_status) and sample_map (TSV vcf_sample, pair; VCF sample names
         differ from the gOS cell ids for most patients).
"""
import argparse
import copy
import csv
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chrx_flags as cf  # noqa: E402
import clonality as cl  # noqa: E402

SOURCE = cl.KNOWN_CN_SOURCE
EVIDENCE = "known driver (rescued)"
MIN_NORMAL_DP = 5
MAX_POP_AF = 1e-3
MIN_CARRIERS = 3
MIN_ALT = 5
MIN_MINOR_STRAND = 0.1
COSMIC_MIN = 3
CN_ALTERED_CN2 = 0.80       # same cut as the site-QC pool (diploid_fraction)
SPAN_PAD = 2000
MIN_SUBCLONAL_F1 = 0.6     # a subclonal rescue must follow a tree clade (best-clade F1), not a few scattered cells
ADD_STATUSES = ("clonal", "subclonal")   # an indeterminate CN-aware call (low detection, carriers not clade-restricted) is artefact-like: reported, not added
TRUNC = ("stopgain", "stoploss", "startloss")
PATIENTS = ["BWH69", "BWH70", "BWH75", "MGH285", "MGH301", "MGH302", "MGH303", "MGH304", "MGH305"]


def jload(p):
    with open(p) as fh:
        return json.load(fh)


def jdump(p, obj):
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p + ".tmp", "w") as fh:
        json.dump(obj, fh, allow_nan=False)
    os.replace(p + ".tmp", p)


# ---------------------------------------------------------------- known-driver rule

def norm_known(k):
    hot = {g: ([v] if isinstance(v, str) else list(v)) for g, v in (k.get("hotspot_residues") or {}).items()}
    return {"genes": set(k.get("genes") or []) | set(hot), "hot": hot, "tert": {(t["chrom"], int(t["pos"]), t["ref"], t["alt"]): t["name"] for t in k.get("tert_promoter") or []}}


def parse_info(s):
    out = {}
    for kv in s.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
        else:
            out[kv] = True
    return out


def cosmic_count(v):
    """Total OCCURENCE of an ANNOVAR cosmic70 value 'ID=COSM1;OCCURENCE=2(lung),1(skin)' (escaped as \\x3d / \\x3b in INFO)."""
    if not v or v == ".":
        return 0
    v = v.replace("\\x3d", "=").replace("\\x3b", ";")
    m = re.search(r"OCCURENCE=([^;]+)", v)
    return sum(int(x) for x in re.findall(r"(\d+)\(", m.group(1))) if m else 0


def proteins(aa):
    """[(gene, transcript, protein 'p.E545A')] from ANNOVAR AAChange.refGene."""
    out = []
    for t in (aa or "").replace("\\x3b", ";").split(","):
        f = t.split(":")
        if len(f) >= 5 and f[4].startswith("p."):
            out.append((f[0], f[1], f[4]))
    return out


def known_basis(chrom, pos, ref, alt, info, known):
    """(basis, tier, type, protein, gene) when the site is a known driver by the rule, else None."""
    t = known["tert"].get((chrom, pos, ref, alt))
    if t:
        return f"TERT promoter hotspot ({t})", 1, "Promoter", t.split(" ", 1)[1] + " (promoter)", "TERT"
    genes = [g for g in str(info.get("Gene.refGene", "")).replace("\\x3b", ";").split(";") if g in known["genes"]]
    if not genes:
        return None
    gene = genes[0]
    func, efunc = str(info.get("Func.refGene", "")), str(info.get("ExonicFunc.refGene", ""))
    prots = [p for g, _, p in proteins(info.get("AAChange.refGene")) if g == gene]
    prot = prots[0] if prots else None
    for p in prots:
        m = re.match(r"p\.([A-Z])(\d+)", p)
        if m and f"{m.group(1)}{m.group(2)}" in known["hot"].get(gene, []):
            return f"configured hotspot residue {gene} {m.group(1)}{m.group(2)} (GBM driver list)", 1, "Missense", p, gene
    if "splicing" in func.split(";"):
        return f"splice-site change in a driver gene ({gene})", 2, "Splice", prot or "splice", gene
    if efunc in TRUNC:
        return f"{efunc} in a driver gene ({gene})", 2, "Trunc", prot, gene
    if efunc.startswith("nonsynonymous") and cosmic_count(info.get("cosmic70")) >= COSMIC_MIN:
        return f"COSMIC-recurrent protein change ({cosmic_count(info.get('cosmic70'))} COSMIC v70 occurrences) in a driver gene ({gene})", 2, "Missense", prot, gene
    return None


# ---------------------------------------------------------------- VCF sweep

def gene_regions(data, known):
    spans = jload(os.path.join(data, "_cohort", "gene_spans.json"))["genes"]
    out = []
    for g in spans:
        if g["gene"] in known["genes"]:
            out.append((f"chr{g['chrom']}", max(1, int(g["start"]) - SPAN_PAD), int(g["end"]) + SPAN_PAD, g["gene"]))
    return sorted(set(out))


def sample_counts(fmt, sample):
    """(ref, alt, alt_f, alt_r) from AD_BAM / ADF_BAM / ADR_BAM (fallback AD) of one sample column."""
    keys = fmt.split(":")
    vals = sample.split(":")
    d = dict(zip(keys, vals))

    def pair(k):
        v = d.get(k)
        if not v or v == ".":
            return None
        x = [int(y) if y not in (".", "") else 0 for y in v.split(",")]
        return x if len(x) >= 2 else None
    ad = pair("AD_BAM") or pair("AD")
    if ad is None:
        return 0, 0, 0, 0
    f, r = pair("ADF_BAM"), pair("ADR_BAM")
    return ad[0], ad[1], (f[1] if f else 0), (r[1] if r else 0)


def read_sample_map(path):
    """{vcf sample: gOS cell} from a TSV with columns vcf_sample, pair (None when no map is given)."""
    if not path:
        return None
    with open(path, newline="") as fh:
        return {r["vcf_sample"]: r["pair"] for r in csv.DictReader(fh, delimiter="\t")}


def sweep(vcf, regions, known, cell_names, tabix="tabix", sample_map=None):
    """Known-driver candidate records of the VCF in the regions -> list of site dicts with per-cell counts (cells by gOS name).
    sample_map {vcf sample: cell}; without one the VCF name with '-' -> '_' is used."""
    hdr = subprocess.run([tabix, "-H", vcf], check=True, capture_output=True, text=True).stdout.splitlines()
    samples = hdr[-1].split("\t")[9:]
    norm = {s: (sample_map.get(s) if sample_map is not None else s.replace("-", "_")) for s in samples}
    idx = [(j, norm[s]) for j, s in enumerate(samples) if norm[s] in cell_names]
    seen, out = set(), []
    # one tabix call per chunk of regions (the region list is short)
    for k in range(0, len(regions), 50):
        chunk = regions[k:k + 50]
        txt = subprocess.run([tabix, vcf] + [f"{c}:{s}-{e}" for c, s, e, _ in chunk], check=True, capture_output=True, text=True).stdout
        for ln in txt.splitlines():
            f = ln.split("\t")
            chrom, pos, ref, alt = f[0], int(f[1]), f[3], f[4]
            if len(ref) != 1 or len(alt) != 1 or (chrom, pos, alt) in seen:
                continue
            info = parse_info(f[7])
            func = str(info.get("Func.refGene", ""))
            if not (known["tert"].get((chrom, pos, ref, alt)) or any(x in func.split(";") for x in ("exonic", "splicing"))):
                continue
            kb = known_basis(chrom, pos, ref, alt, info, known)
            if kb is None:
                continue
            seen.add((chrom, pos, alt))
            counts = {}
            for j, c in idx:
                counts[c] = sample_counts(f[8], f[9 + j])
            out.append(dict(chrom=chrom, pos=pos, ref=ref, alt=alt, basis=kb[0], tier=kb[1], type=kb[2], protein=kb[3], gene=kb[4],
                            info={x: info.get(x) for x in ("Func.refGene", "ExonicFunc.refGene", "AAChange.refGene", "cosmic70", "avsnp150",
                                                            "gnomad40_genome_AF", "branch_filter", "Germline_qval")},
                            counts=counts))
    return out


def pop_af(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def read_qc(site, tumour, normal):
    t = [site["counts"][c] for c in tumour if c in site["counts"]]
    n = [site["counts"][c] for c in normal if c in site["counts"]]
    tr, ta = sum(x[0] for x in t), sum(x[1] for x in t)
    taf, tar = sum(x[2] for x in t), sum(x[3] for x in t)
    nr, na = sum(x[0] for x in n), sum(x[1] for x in n)
    carriers = [c for c in tumour if c in site["counts"] and site["counts"][c][1] > 0]
    covered = [c for c in tumour if c in site["counts"] and sum(site["counts"][c][:2]) > 0]
    af = pop_af(site["info"].get("gnomad40_genome_AF"))
    checks = []

    def chk(key, ok, text):
        checks.append({"key": key, "status": "pass" if ok else "fail", "text": text})
    chk("normal", nr + na >= MIN_NORMAL_DP and na == 0, f"normal cells: {na} alt / {nr + na} reads (need >= {MIN_NORMAL_DP} reads, 0 alt)")
    chk("population", af is None or af <= MAX_POP_AF, "gnomAD AF " + ("not annotated" if af is None else f"{af:g}") + f" (germline if > {MAX_POP_AF:g})")
    chk("carriers", len(carriers) >= MIN_CARRIERS and ta >= MIN_ALT, f"{len(carriers)} tumour carrier cells, {ta} alt reads (need >= {MIN_CARRIERS} cells, >= {MIN_ALT} reads)")
    minor = min(taf, tar) / ta if ta else 0
    chk("strand", ta >= MIN_ALT and (taf + tar == 0 or minor >= MIN_MINOR_STRAND), f"alt reads forward {taf} / reverse {tar}" + (" (no strand counts)" if taf + tar == 0 else ""))
    return {"pass": all(c["status"] == "pass" for c in checks), "checks": checks,
            "tumor": {"alt": ta, "ref": tr, "alt_f": taf, "alt_r": tar, "cells_alt": len(carriers), "cells_cov": len(covered)},
            "normal": {"alt": na, "ref": nr, "cells_cov": sum(1 for x in n if x[0] + x[1] > 0)},
            "carriers": carriers, "covered": covered}


# ---------------------------------------------------------------- pipeline status

def grep_ids(path, ids):
    """Subset of ids (chr_pos_ref_alt) present in the first column of a gzipped TSV."""
    if not path or not os.path.exists(path) or not ids:
        return None
    pat = "\n".join(ids)
    with tempfile.NamedTemporaryFile("w", suffix=".ids", delete=False) as fh:
        fh.write(pat + "\n")
    try:
        r = subprocess.run(f"gzip -dc {shlex_quote(path)} | cut -f1 | grep -F -x -f {shlex_quote(fh.name)}", shell=True, capture_output=True, text=True)
    finally:
        os.remove(fh.name)
    return set(r.stdout.split()) & set(ids)


def shlex_quote(s):
    import shlex
    return shlex.quote(s)


# ---------------------------------------------------------------- rows

def site_key(chrom, pos):
    return (str(chrom).replace("chr", ""), str(int(pos)))


def event_key(e):
    return (str(e.get("seqnames", "")).replace("chr", ""), str(e.get("start")))


def new_row(patient, s, qc, tumour_n):
    chrom = s["chrom"].replace("chr", "")
    t = qc["tumor"]
    covered = len(qc["covered"])
    ref_only = [c for c in qc["covered"] if c not in set(qc["carriers"])]
    return {
        "gene": s["gene"], "fusion_genes": "", "vartype": "Hotspot" if s["tier"] == 1 else "SNV", "type": s["type"], "Variant": s["protein"] or f"{s['ref']}>{s['alt']}",
        "gene_summary": f"{s['gene']} known driver recovered from the genotyping VCF in a copy-number-altered region (not a tree site).",
        "Variant_g": f"{chrom}:{s['pos']}-{s['pos']} {s['ref']}>{s['alt']}", "Tier": s["tier"], "tier_basis": "rule: " + s["basis"],
        "variant_summary": s["basis"], "effect": None, "effect_description": s["basis"], "is_multi_hit_per_gene": False,
        "Genome_Location": f"{chrom}:{s['pos']}-{s['pos']}", "seqnames": chrom, "start": str(s["pos"]), "end": str(s["pos"]),
        "sample": qc["carriers"][0], "role": None, "fusion_cn": None, "fusion_gene_coords": None, "upstream_only": None, "upstream_partner": None,
        "therapeutics": None, "prognoses": None, "diagnoses": None, "id": patient,
        "alt": t["alt"], "ref": t["ref"], "estimated_altered_copies": None, "segment_cn": None,
        "n_cells": len(qc["carriers"]), "n_normal_cells": 0, "cell_ids": ",".join(qc["carriers"]),
        "VAF": round(t["alt"] / max(1, t["alt"] + t["ref"]), 4), "cells": f"{len(qc['carriers'])}/{covered}",
        "cell_fraction": round(len(qc["carriers"]) / covered, 3) if covered else None,
        "source": SOURCE, "rescued": True, "covered_ref_cell_ids": ",".join(ref_only), "tumour_cells": tumour_n, "covered_tumour_cells": covered,
        "coverage_fraction": round(covered / tumour_n, 3) if tumour_n else None, "low_coverage": (covered / tumour_n if tumour_n else 0) < 0.5,
        "hotspot_evidence": {"preset": "genotyping VCF AD_BAM (per-cell CRAM recount)", "cell_source": "genotyping VCF AD_BAM",
                             "tumor": t, "normal": qc["normal"],
                             "verdict": {"level": "supported", "label": "read QC pass", "summary": "; ".join(c["text"] for c in qc["checks"]),
                                         "reasons": [], "checks": qc["checks"]}},
        "annovar": s["info"],
    }


def mark(row, res, state, status, basis=None):
    prev = row.get("clonality")
    if prev is not None and not (row.get("clonality_stats") or {}).get("cn_aware"):
        row["clonality_dropout_only"] = prev
    row["clonality"] = res["status"]
    row["clonality_method"] = cl.METHOD
    row["clonality_stats"] = res
    row["evidence"] = EVIDENCE
    row["site_cn"] = state
    row["cn_altered"] = state.get("cn_altered")
    if basis:
        row["known_basis"] = basis
    if status is not None:
        row["pipeline_status"] = status
    return row


def run_patient(data, out, patient, known, man, tabix, tag, dry, log=print, chrx_flagged=None, sex=None):
    pdir = os.path.join(data, patient)
    events = jload(os.path.join(pdir, "filtered.events.json"))
    kept = [e for e in events if e.get("source") != SOURCE]
    ctx = cl.PatientContext(data, patient)
    if ctx.reference is None:
        ctx.reference = cl.cohort_fallback(data)
        if ctx.reference is not None:
            ctx.reference = dict(ctx.reference, scope="cohort")
    try:
        cn = cl.CellCN(data, patient, chrx_flagged=chrx_flagged)
    except (OSError, KeyError, ValueError) as err:
        cn, cn_err = None, str(err)
    else:
        cn_err = None
    tumour = sorted(ctx.tumour)
    tset = set(tumour)
    normal = [c for i, c in enumerate(ctx.cells) if i not in set(ctx.tumour_idx)]
    snv = cl._load(os.path.join(pdir, "snv_matrix.json")) or {}
    tree_ids = {v["id"] for v in snv.get("variants") or []}
    report = {"patient": patient, "cn_source": cn.source if cn else f"not available ({cn_err})", "reference": (ctx.reference or {}).get("source"),
              "existing": [], "swept": [], "added": [], "notes": []}

    def state_of(chrom, pos):
        if cn is None:
            return {"modal_cn": None, "cn2_fraction": None, "counts": {}, "cn_altered": None}
        st = cn.state(chrom, pos, tset)
        if st["cn2_fraction"] is None:
            st["cn_altered"] = None
        elif cf.is_chrx(chrom) and chrx_flagged:
            # chrX with artefact cells excluded: expected CN follows the patient's sex (male 1, female 2); unknown sex -> not assessed
            exp = cf.expected_chrx_cn(sex)
            if exp is None:
                st["cn_altered"] = None
            else:
                frac = st["counts"].get(str(exp), 0) / st["n_cells"]
                st["cn_altered"] = frac < CN_ALTERED_CN2 or st["modal_cn"] != exp
                st["expected_cn"] = exp
        else:
            st["cn_altered"] = st["cn2_fraction"] < CN_ALTERED_CN2 or st["modal_cn"] != 2
        if st["cn_altered"] is None and st.get("chrx_not_assessed"):
            st["cn_status"] = cf.NOT_ASSESSED
        return st

    # 1. existing read-evidence driver rows
    for e in kept:
        sid = cl.site_id_of(e)
        if not (cl.is_read_evidence(e) and sid and e.get("gene") in known["genes"]):
            continue
        chrom, pos = sid.split(":")[0], int(sid.split(":")[1])
        res = ctx.assess_cn(e, cn)
        st = state_of(chrom, pos)
        before = e.get("clonality") if not (e.get("clonality_stats") or {}).get("cn_aware") else e.get("clonality_dropout_only")
        mark(e, res, st, None)
        report["existing"].append({"gene": e["gene"], "variant": e.get("Variant"), "site": sid, "source": e.get("source"), "cn_altered": st["cn_altered"], "cn_status": st.get("cn_status"),
                                   "modal_cn": st["modal_cn"], "cn2_fraction": st["cn2_fraction"], "clonality_before": before, "clonality_after": res["status"],
                                   "allele_loss": res.get("allele_loss"), "reason": res.get("reason"), "evidence_source": res.get("evidence_source")})

    # 2. VCF sweep
    vcf = (man or {}).get("vcf")
    if not vcf or not os.path.exists(vcf):
        report["notes"].append(f"VCF sweep not run: no genotyping VCF for {patient}" + (f" ({vcf} missing)" if vcf else ""))
        sites = []
    else:
        smap = read_sample_map((man or {}).get("sample_map"))
        sites = sweep(vcf, gene_regions(data, known), known, set(ctx.cells), tabix, smap)
        n_matched = len({smap[x] for x in smap if smap[x] in set(ctx.cells)}) if smap is not None else None
        report["vcf_cells_matched"] = n_matched
        if smap is not None and n_matched < 0.5 * len(ctx.cells):
            report["notes"].append(f"only {n_matched} of {len(ctx.cells)} cells matched to VCF samples")
    ids = [f"{s['chrom']}_{s['pos']}_{s['ref']}_{s['alt']}" for s in sites]
    in_sf1 = grep_ids((man or {}).get("site_candidates"), ids)
    in_pool = grep_ids((man or {}).get("pool_sites"), ids)
    by_pos = {}
    for e in kept:
        if e.get("vartype") in ("SNV", "Hotspot"):
            by_pos.setdefault(event_key(e), []).append(e)
    tumour_n = len(tumour)
    for s, mid in zip(sites, ids):
        qc = read_qc(s, tumour, normal)
        st = state_of(s["chrom"], s["pos"])
        status = {"tree_site": mid in tree_ids,
                  "site_filter_candidate": None if in_sf1 is None else mid in in_sf1,
                  "site_qc_pool": None if in_pool is None else mid in in_pool,
                  "sequoia_branch_filter": s["info"].get("branch_filter")}
        ex = by_pos.get(site_key(s["chrom"], s["pos"])) or []
        rec = {"gene": s["gene"], "variant": s["protein"], "site": f"{s['chrom']}:{s['pos']}:{s['ref']}>{s['alt']}", "basis": s["basis"], "tier": s["tier"],
               "qc_pass": qc["pass"], "qc": {c["key"]: c["text"] for c in qc["checks"] if c["status"] == "fail"}, "tumor": qc["tumor"], "normal": qc["normal"],
               "cn_altered": st["cn_altered"], "modal_cn": st["modal_cn"], "cn2_fraction": st["cn2_fraction"], "pipeline_status": status,
               "existing_sources": sorted({str(e.get("source")) for e in ex}) or None}
        report["swept"].append(rec)
        if ex:
            for e in ex:
                if e.get("evidence") == EVIDENCE:
                    e["pipeline_status"] = status
                    e.setdefault("known_basis", s["basis"])
            rec["action"] = "already an event"
            continue
        if not qc["pass"]:
            rec["action"] = "read QC fail"
            continue
        if st["cn_altered"] is None and st.get("cn_status"):
            rec["action"] = st["cn_status"] + " (left for review)"
            continue
        if not st["cn_altered"]:
            rec["action"] = "not CN-altered (diploid site missed for another reason; left for review)"
            continue
        row = new_row(patient, s, qc, tumour_n)
        alt = {c: s["counts"][c][1] for c in qc["covered"]}
        dp = {c: s["counts"][c][0] + s["counts"][c][1] for c in qc["covered"]}
        res = ctx.assess_cn(row, cn, counts=(alt, dp, "genotyping VCF AD_BAM (per-cell CRAM recount)"))
        mark(row, res, st, status, s["basis"])
        why = None
        if res["status"] not in ADD_STATUSES:
            why = "not added: carriers neither consistent with clonal nor clade-restricted (artefact-like)"
        elif res["status"] == "subclonal" and (res.get("best_clade_f1") or 0) < MIN_SUBCLONAL_F1:
            why = f"not added: subclonal but best tree-clade F1 {res.get('best_clade_f1')} < {MIN_SUBCLONAL_F1}"
        elif s["type"] == "Splice":
            why = "review: splice-site change (RNA / cDNA carry-over on spliced reads not checked from the VCF counts)"
        if why:
            rec.update(action=why, clonality=res["status"], clonality_reason=res.get("reason"), best_clade_f1=res.get("best_clade_f1"))
            continue
        kept.append(row)
        rec.update(action="added", clonality=res["status"], clonality_reason=res.get("reason"), allele_loss=res.get("allele_loss"))
        report["added"].append(rec)
    report["created"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    if not dry:
        odir = os.path.join(out, patient)
        path = os.path.join(odir, "filtered.events.json")
        if os.path.abspath(out) == os.path.abspath(data) and os.path.exists(path) and not os.path.exists(f"{path}.{tag}"):
            shutil.copy2(path, f"{path}.{tag}")
        if os.path.islink(path):
            os.remove(path)
        jdump(path, kept)
        jdump(os.path.join(odir, "known_cn_rescue.json"), dict(format="gos-sc-known-cn-rescue/1", **report))
    log(f"{patient}: existing re-annotated {len(report['existing'])}, swept {len(sites)}, added {len(report['added'])}")
    return report


def read_manifest(path):
    out = {}
    if not path:
        return out
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["patient"]] = {k: (v or None) for k, v in r.items()}
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--known", required=True)
    ap.add_argument("--vcf-manifest")
    ap.add_argument("--patients", default=",".join(PATIENTS))
    ap.add_argument("--tabix", default=os.environ.get("TABIX", "tabix"))
    ap.add_argument("--tag", default="bak_known_cn")
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--chrx-flags", help="cells_chrx_flags.tsv: chrX CN of flagged cells is 'not assessed' (chrx_flags.py)")
    ap.add_argument("--sex-table", help="patient_sex_evidence.tsv: expected chrX CN (male 1, female 2) for the unflagged cells")
    a = ap.parse_args(argv)
    flags = cf.load_flagged(a.chrx_flags) if a.chrx_flags else {}
    sexes = cf.load_sex(a.sex_table) if a.sex_table else {}
    known = norm_known(jload(a.known))
    man = read_manifest(a.vcf_manifest)
    reps = {}
    for p in [x for x in a.patients.split(",") if x]:
        reps[p] = run_patient(a.data, a.out, p, known, man.get(p), a.tabix, a.tag, a.dry, chrx_flagged=flags.get(p), sex=sexes.get(p))
    if not a.dry:
        summ = {p: {"existing": r["existing"], "added": r["added"], "n_swept": len(r["swept"]), "notes": r["notes"], "cn_source": r["cn_source"]} for p, r in reps.items()}
        jdump(os.path.join(a.out, "_cohort", "known_cn_rescue.json"), {"format": "gos-sc-known-cn-rescue-cohort/1", "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
                                                                       "evidence": EVIDENCE, "source": SOURCE, "patients": summ})
    return 0


if __name__ == "__main__":
    sys.exit(main())
