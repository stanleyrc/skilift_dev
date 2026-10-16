#!/usr/bin/env python3
"""gos_sc_allsnv: driver events from ALL called SNVs in driver genes, not only the tree sites.

Pipeline (each stage is a MODE; see README_allsnv.md in this folder or the header of each script):
  1. gos_sc_allsnv_select.R      called SNVs (cellphy_site_selection_all.tsv.gz) in OncoKB-cancer-gene + GBM region-set gene bodies -> WORK/<P>.candidates.txt
  2. gos_sc_hotspot_annotate.R   SnpEff + OncoKB (skilift sc_annotate_snv_drivers / sc_cell_snv_oncotable) with WORK=<this WORK> -> WORK/annotation.tsv
  3. MODE=prep       this script: annotated protein-altering / splice sites not already events -> WORK/pileup_sites_<P>.tsv
  4. gos_sc_allsnv_pileup.py MODE=own     per-cell read pileups of those sites in the patient's cells (BAMs)
  5. MODE=shortlist  provisional verdicts (gOS siteEvidence.js rankFindings through hotspot_verdicts.js) -> WORK/shortlist.tsv (sites that pass)
  6. gos_sc_allsnv_pileup.py MODE=cohort  the shortlisted sites in every patient's cells (cross-patient background for the cohort check)
  7. MODE=rows       final verdicts (with the cohort background) + RNA / cDNA carry-over check -> event rows -> filtered.events.json (backup first)

A row (source "called SNV (not a tree site)") carries: cell_ids (tumour cells with >= 1 alt read: MAPQ >= 20, BQ >= 13, no duplicate),
covered_ref_cell_ids, cells / cell_fraction (carriers / covered tumour cells), clonality (clonal >= 0.85, else subclonal; "not assessed" when
< 50% of tumour cells are covered), the pooled tumour / normal alt and ref reads per strand and the verdict with its checks (hotspot_evidence), and
the OncoKB tier / effect. Sites at an existing event position (tree sites, hotspot rescue) are skipped.
"""
import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gos_sc_hotspot_rescue as hr  # noqa: E402

SOURCE = "called SNV (not a tree site)"
DATA = os.environ.get("DATA_DIR", "/gpfs/commons/projects/imielinski_web/sclarke/gOS_GBM_sc/data")
WORK = os.environ.get("WORK", "/gpfs/commons/groups/imielinski_lab/projects/single_cell_GBM/gos_datasets/gOS_GBM_sc/work/allsnv")
GOS_SRC = os.environ.get("GOS_SRC", "/gpfs/commons/groups/imielinski_lab/home/sclarke/git/gOS_dev_wt_integration")
PATIENTS = [p for p in os.environ.get("PATIENTS", "BWH69,BWH70,MGH302,MGH303,MGH304,MGH305").split(",") if p]
ALL_PATIENTS = ["BWH69", "BWH70", "BWH75", "MGH285", "MGH301", "MGH302", "MGH303", "MGH304", "MGH305"]     # cohort background: every patient's cells
MIN_CARRIERS = 3
CLUSTER_BP = 3          # passing SNVs this close with (nearly) the same carrier cells are one MNV / misalignment, not independent SNVs
CLUSTER_JACCARD = 0.7
PASS_LEVELS = ("supported", "weak")      # levels worth a second look; accept() decides
CLONAL = 0.85
LOW_COVERAGE = 0.5
# RNA / cDNA carry-over (or processed pseudogene): alt reads sit on spliced / split / clipped alignments much more often than reference reads do
RNA_MIN_SPLICED = 2
RNA_ALT_FRACTION = 0.15
RNA_FOLD = 3.0
INCLUDE_TIER3_MISSENSE = os.environ.get("INCLUDE_TIER3_MISSENSE", "0") == "1"


def jload(p):
    with open(p) as fh:
        return json.load(fh)


def rna_like(sp):
    """(flag, text) from the spliced-read counts of tumour caller-filtered reads."""
    a, asp, r, rsp = sp["alt"], sp["alt_spl"], sp["ref"], sp["ref_spl"]
    if a == 0:
        return False, None
    fa, fr = asp / a, (rsp / r if r else 0.0)
    flag = asp >= RNA_MIN_SPLICED and fa >= RNA_ALT_FRACTION and fa >= RNA_FOLD * fr
    return flag, f"{asp} of {a} alt reads ({fa:.0%}) are spliced / split / clipped alignments vs {rsp} of {r} ref reads ({fr:.0%})"


def accept(x, tier, o):
    """Verdict rule for a site that is not a tree site: 'supported'; 'weak' only when the single warning is the strand test on a deep site (as hotspot
    rescue's discovered sites) or the site is an OncoKB tier 1-2 variant (a known oncogenic variant is a prior; the weak verdict stays on the row)."""
    if x["level"] == "supported":
        return True
    if x["level"] != "weak":
        return False
    if hr.strand_only_warning(x, {"presets": o["presets"]}):
        return True
    return tier is not None and tier <= 2


def read_candidates(patient, events):
    """Annotated protein-altering / splice sites of the patient that are not events yet: id -> annotation + joint-call stats."""
    ann = {}
    with open(os.path.join(WORK, "annotation.tsv"), newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["patient"] == patient:
                r = {k: (None if v in ("NA", "") else v) for k, v in r.items()}
                if r.get("tier") is not None:
                    r["tier"] = int(float(r["tier"]))
                ann[r["mutation"]] = r
    stats = {}
    with gzip.open(os.path.join(WORK, patient + ".sites.tsv.gz"), "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["mutation"] in ann:
                stats[r["mutation"]] = r
    with open(os.path.join(WORK, "driver_genes.txt")) as fh:
        drivers = {ln.strip() for ln in fh if ln.strip()}
    ann = {m: a for m, a in ann.items() if a.get("gene") in drivers}      # SnpEff may name a neighbouring (non-driver) gene
    taken = {(str(e.get("seqnames", "")).replace("chr", ""), str(e.get("start"))) for e in events if e.get("source") != SOURCE and e.get("vartype") in ("SNV", "Hotspot")}
    out, skipped = {}, 0
    for m, a in ann.items():
        chrom, pos, ref, alt = m.split("_")
        if (chrom.replace("chr", ""), pos) in taken:
            skipped += 1
            continue
        out[m] = dict(a, chrom=chrom, pos=int(pos), ref=ref, alt=alt, stats=stats.get(m, {}))
    return out, skipped


def site_id(c):
    return f"{c['chrom']}:{c['pos']}:{c['ref']}>{c['alt']}"


def build_docs(patient, cands, own, cohort):
    """region_calls-like and site_evidence-like documents (same layout as precompute/) for hotspot_verdicts.js."""
    rc = jload(os.path.join(DATA, patient, "precompute", "region_calls.json"))
    n = len(rc["cells"])
    normal = {i for i, c in enumerate(rc.get("clone") or []) if str(c).lower() == "normal"}
    sites, dna, ev = [], {}, {}
    for m, c in cands.items():
        sid = site_id(c)
        o = own["sites"].get(sid)
        if o is None:
            continue
        alt, dp = o["cells"]["alt"], o["cells"]["dp"]
        t = o["presets"]["caller"]["tumor"]
        nm = o["presets"]["caller"]["normal"]
        sites.append(dict(id=sid, chrom=c["chrom"], pos=c["pos"], ref=c["ref"], alt=c["alt"], name=f"{c['gene']} {c.get('protein') or ''}".strip(), gene=c["gene"],
                          source=SOURCE, hotspot=False, filter="PASS", prevalence=None, alt_tumor=t["alt_f"] + t["alt_r"], dp_tumor=t["alt_f"] + t["alt_r"] + t["ref_f"] + t["ref_r"],
                          cells_alt=t["cells_alt"], cells_cov=t["cells_cov"], alt_normal=nm["alt_f"] + nm["alt_r"], dp_normal=nm["alt_f"] + nm["alt_r"] + nm["ref_f"] + nm["ref_r"],
                          n_mut=sum(1 for i in range(n) if alt[i] > 0 and i not in normal), error=0.001))
        dna[sid] = dict(alt=alt, dp=dp)
        coh = [dict(v, patient=q) for q, d in (cohort or {}).items() if q != patient for v in [d["sites"].get(sid)] if v]
        ev[sid] = dict(sample_every=1, gene=c["gene"], n_reads=o["n_reads"], presets=o["presets"], cohort=coh)
    calls = dict(format="gos-sc-region-calls/1", patient=patient, cells=rc["cells"], clone=rc.get("clone"), sites=sites, dna=dna)
    evdoc = dict(format="gos-sc-site-evidence/1", patient=patient, sites=ev)
    return calls, evdoc


def verdicts(patient, cands, own, cohort, tag):
    calls, evdoc = build_docs(patient, cands, own, cohort)
    d = os.path.join(WORK, "verdict_" + tag, patient, "precompute")
    os.makedirs(d, exist_ok=True)
    json.dump(calls, open(os.path.join(d, "region_calls.json"), "w"), separators=(",", ":"))
    json.dump(evdoc, open(os.path.join(d, "site_evidence.json"), "w"), separators=(",", ":"))
    out = subprocess.run(["node", os.path.join(HERE, "hotspot_verdicts.js"), GOS_SRC, os.path.dirname(d)], check=True, capture_output=True, text=True).stdout
    return json.loads(out), calls, evdoc


def mode_prep():
    for p in PATIENTS:
        events = jload(os.path.join(DATA, p, "filtered.events.json"))
        cands, skipped = read_candidates(p, events)
        with open(os.path.join(WORK, f"pileup_sites_{p}.tsv"), "w") as fh:
            fh.write("id\tchrom\tpos\tref\talt\n")
            for c in cands.values():
                fh.write(f"{site_id(c)}\t{c['chrom']}\t{c['pos']}\t{c['ref']}\t{c['alt']}\n")
        print(p, len(cands), "annotated protein-altering / splice sites to pile up;", skipped, "already events", flush=True)


def mode_shortlist():
    ids = {}
    for p in PATIENTS:
        events = jload(os.path.join(DATA, p, "filtered.events.json"))
        cands, _ = read_candidates(p, events)
        own = jload(os.path.join(WORK, f"pileup_own_{p}.json"))
        v, _, _ = verdicts(p, cands, own, None, "prov")
        keep = 0
        for sid, x in v.items():
            c = next(c for c in cands.values() if site_id(c) == sid)
            if x["carrierCells"] >= MIN_CARRIERS and accept(x, c["tier"], own["sites"][sid]):
                ids[sid] = (c["chrom"], c["pos"], c["ref"], c["alt"])
                keep += 1
        print(p, keep, "of", len(v), "pass the provisional verdict", flush=True)
    with open(os.path.join(WORK, "shortlist.tsv"), "w") as fh:
        fh.write("id\tchrom\tpos\tref\talt\n")
        for sid, (ch, ps, r, a) in sorted(ids.items()):
            fh.write(f"{sid}\t{ch}\t{ps}\t{r}\t{a}\n")
    print(len(ids), "sites in the cross-patient shortlist")


def make_row(patient, c, calls, o, v, rna):
    sid = site_id(c)
    cells = calls["cells"]
    clone = calls.get("clone") or []
    tumour = [i for i in range(len(cells)) if not (i < len(clone) and hr.is_normal(clone[i]))]
    alt, dp = o["cells"]["alt"], o["cells"]["dp"]
    carriers = [cells[i] for i in tumour if alt[i] > 0]
    ref_only = [cells[i] for i in tumour if alt[i] == 0 and dp[i] > 0]
    covered = len(carriers) + len(ref_only)
    cov_frac = covered / len(tumour) if tumour else 0
    frac = len(carriers) / covered if covered else 0
    pr = o["presets"]["caller"]
    t, nm = pr["tumor"], pr["normal"]
    a_t, r_t = t["alt_f"] + t["alt_r"], t["ref_f"] + t["ref_r"]
    clonality = "not assessed" if cov_frac < LOW_COVERAGE else "clonal" if frac >= CLONAL else "subclonal"
    chrom = str(c["chrom"]).replace("chr", "")
    pos = int(c["pos"])
    return {
        "gene": c["gene"], "fusion_genes": "", "vartype": "SNV", "type": c["type"], "Variant": c["protein"] or f"{c['ref']}>{c['alt']}",
        "gene_summary": c.get("gene_summary"), "Variant_g": f"{chrom}:{pos}-{pos} {c['ref']}>{c['alt']}", "Tier": c["tier"], "variant_summary": c.get("variant_summary"),
        "effect": c.get("effect"), "effect_description": c.get("effect_description") or c.get("variant_summary"), "is_multi_hit_per_gene": False,
        "Genome_Location": f"{chrom}:{pos}-{pos}", "seqnames": chrom, "start": str(pos), "end": str(pos), "sample": carriers[0], "role": c.get("role"),
        "fusion_cn": None, "fusion_gene_coords": None, "upstream_only": None, "upstream_partner": None,
        "therapeutics": c.get("therapeutics"), "prognoses": c.get("prognoses"), "diagnoses": c.get("diagnoses"), "id": patient,
        "alt": a_t, "ref": r_t, "estimated_altered_copies": None, "segment_cn": None, "n_cells": len(carriers), "n_normal_cells": nm["cells_alt"],
        "cell_ids": ",".join(carriers), "VAF": round(a_t / max(1, a_t + r_t), 4), "cells": f"{len(carriers)}/{covered}", "cell_fraction": round(frac, 3),
        "source": SOURCE, "rescued": True, "covered_ref_cell_ids": ",".join(ref_only), "tumour_cells": len(tumour), "covered_tumour_cells": covered,
        "coverage_fraction": round(cov_frac, 3), "low_coverage": cov_frac < LOW_COVERAGE, "clonality": clonality,
        "hotspot_evidence": {
            "preset": "caller", "cell_source": "BAM pileup (MAPQ >= 20, BQ >= 13, no duplicate / secondary / QC-fail reads)",
            "tumor": {k: t[k] for k in ("alt_f", "alt_r", "ref_f", "ref_r", "cells_alt", "cells_cov")},
            "normal": {k: nm[k] for k in ("alt_f", "alt_r", "ref_f", "ref_r", "cells_alt", "cells_cov")},
            "splice": o["splice"], "rna_like": bool(rna[0]),
            "joint_call": {k: c["stats"].get(k) for k in ("tumor_alt_read_cells", "tumor_ref_cells", "tumor_cells_observed", "tumor_cells_diploid", "normal_alt_read_cells", "normal_truth_class")},
            "verdict": {"level": v["level"], "label": v["label"], "summary": v["summary"],
                        "reasons": [x["text"] for x in v["checks"] if x["status"] in ("warn", "fail")],
                        "checks": [{"key": x["key"], "status": x["status"], "text": x["text"]} for x in v["checks"]]},
        },
    }


def mode_rows():
    cohort = {}
    for p in ALL_PATIENTS:
        f = os.path.join(WORK, f"pileup_cohort_{p}.json")
        if os.path.exists(f):
            cohort[p] = jload(f)
    tag = os.environ.get("BACKUP_TAG", "bak_261016")
    dry = os.environ.get("DRY_RUN") == "1"
    outdir = os.environ.get("OUT_DIR", os.path.join(WORK, "staged"))
    summary = {}
    for p in PATIENTS:
        path = os.path.join(DATA, p, "filtered.events.json")
        events = jload(path)
        kept = [e for e in events if e.get("source") != SOURCE]
        cands, skipped = read_candidates(p, kept)
        own = jload(os.path.join(WORK, f"pileup_own_{p}.json"))
        sf = os.path.join(WORK, f"pileup_split_{p}.json")      # split-read counts (re-pileup of the sites that pass the provisional verdict)
        if os.path.exists(sf):
            for sid, x in jload(sf)["sites"].items():
                if sid in own["sites"]:
                    own["sites"][sid]["splice"] = x["splice"]
        v, calls, _ = verdicts(p, cands, own, cohort, "final")
        rows, report, passing = [], [], []
        for m, c in cands.items():
            sid = site_id(c)
            x = v.get(sid)
            o = own["sites"].get(sid)
            rec = dict(site=sid, gene=c["gene"], variant=c["protein"], type=c["type"], tier=c["tier"], called=False)
            report.append(rec)
            if x is None or o is None:
                rec["reason"] = "no reads"
                continue
            rna = rna_like(o["splice"])
            rec.update(level=x["level"], carriers=x["carrierCells"], covered_tumour=x["coveredTumour"], alt=o["presets"]["caller"]["tumor"]["alt_f"] + o["presets"]["caller"]["tumor"]["alt_r"],
                       normal_alt=o["presets"]["caller"]["normal"]["alt_f"] + o["presets"]["caller"]["normal"]["alt_r"], rna_like=rna[0], summary=x["summary"])
            if x["carrierCells"] < MIN_CARRIERS:
                rec["reason"] = f"only {x['carrierCells']} carrier cells"
            elif not accept(x, c["tier"], o):
                rec["reason"] = f"verdict {x['level']}: {x['summary']}"
            elif rna[0]:
                rec["reason"] = "RNA / cDNA carry-over-like: " + rna[1]
            elif c["protein"] and ("ext*" in c["protein"] or c["protein"].endswith("ext*?")):
                rec["reason"] = "stop-loss (extension), not a truncating variant"
            elif c["tier"] is not None and c["tier"] >= 3 and c["type"] == "Missense" and not INCLUDE_TIER3_MISSENSE:
                rec["reason"] = "passes QC but tier 3 missense (not added; INCLUDE_TIER3_MISSENSE=1 adds them)"
                rec["qc_pass"] = True
            else:
                passing.append((c, x, o, rna, rec))
        # clustered substitutions: neighbouring SNVs seen in the same cells are one MNV or a misalignment (the amino-acid change of each alone is wrong)
        carriers = {id(t[4]): {i for i, a in enumerate(t[2]["cells"]["alt"]) if a > 0} for t in passing}
        bad = set()
        for i, (c, _, _, _, rec) in enumerate(passing):
            for c2, _, _, _, rec2 in passing[i + 1:]:
                if c2["chrom"] == c["chrom"] and abs(c2["pos"] - c["pos"]) <= CLUSTER_BP:
                    a, b = carriers[id(rec)], carriers[id(rec2)]
                    if a and b and len(a & b) / len(a | b) >= CLUSTER_JACCARD:
                        bad.update((id(rec), id(rec2)))
        for c, x, o, rna, rec in passing:
            if id(rec) in bad:
                rec["reason"] = f"clustered substitutions (another called SNV within {CLUSTER_BP} bp in the same cells: MNV or misalignment)"
                continue
            row = make_row(p, c, calls, o, x, rna)
            rows.append(row)
            rec.update(called=True, clonality=row["clonality"], cell_fraction=row["cell_fraction"], low_coverage=row["low_coverage"])
        summary[p] = dict(candidates=len(cands), already_events=skipped, added=len(rows), sites=report)
        print(p, "candidates", len(cands), "added", len(rows), flush=True)
        if not dry:
            os.makedirs(os.path.join(outdir, p), exist_ok=True)
            with open(os.path.join(outdir, p, "filtered.events.json"), "w") as fh:
                json.dump(kept + rows, fh)
    doc = dict(format="gos-sc-driver-all-snvs/1", created=time.strftime("%Y-%m-%dT%H:%M:%S"), source=SOURCE, patients=summary)
    with open(os.path.join(WORK, "driver_all_snvs.json"), "w") as fh:
        json.dump(doc, fh)


if __name__ == "__main__":
    {"prep": mode_prep, "shortlist": mode_shortlist, "rows": mode_rows}[os.environ.get("MODE", "rows")]()
