#!/usr/bin/env python3
"""Hotspot rescue: configured hotspot drivers (TERT promoter C228T / C250T, IDH1 R132, H3-3A K27M / G34, BRAF V600E)
that the joint SNV calling misses (GC-rich, low depth, non-coding) are added to the patient-level filtered.events.json
from the read-level site evidence (data/<P>/precompute/site_evidence.json + region_calls.json).

Usage (after the `events` and gos_sc_snv_evidence.R steps, which rewrite filtered.events.json):
  PATIENTS=MGH285,MGH302 python3 gos_sc_hotspot_rescue.py          # default: every patient with site_evidence.json
  DRY_RUN=1 ...                                                    # print the rescued list, write nothing
  DATA_DIR=.../gOS_GBM_sc/data  GOS_SRC=<gOS_dev worktree>  BACKUP_TAG=bak_261015c

A configured hotspot site (slices/regions.json `hotspots`, with the configured alt allele) is called when the front end's own
verdict code (siteEvidence.js rankFindings, run by hotspot_verdicts.js) says "supported" or "weak" and at least MIN_CARRIERS
tumour cells show an alt read. The row is labelled source = "site evidence (not in SNV calls)" and carries the read evidence:
  cell_ids          tumour cells with >= 1 alt read (carriers)
  covered_ref_cell_ids  tumour cells with reads but no alt read ("no alt read"; allele dropout is common, so not "wild type")
  cells / cell_fraction carriers / covered tumour cells (a lower bound for the fraction: dropout), coverage = covered / tumour cells
  clonality / clonality_method / clonality_stats  dropout-aware call (clonality.py: detection vs the patient's truncal SNVs at matched depth,
                    clade restriction); never derive clonality from cell_fraction alone
  hotspot_evidence  pooled tumour / normal alt and ref reads per strand (preset "all"), cells, verdict and its checks
Cells outside `covered` have no coverage. Re-running replaces earlier rescued rows; the previous file is kept as .<BACKUP_TAG>.
"""
import csv
import json
import os
import shutil
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import clonality as cl  # noqa: E402

SOURCE = "site evidence (not in SNV calls)"
MIN_CARRIERS = 3
PASS_LEVELS = ("supported", "weak")      # configured hotspots (prior: a known driver position)
MIN_MINOR_STRAND = 0.35                  # see strand_only_warning
DISCOVERED_LEVELS = ("supported",)       # other sites of driver genes: only the clean ones
LOW_COVERAGE = 0.5          # covered / tumour cells below this is flagged low coverage

# OncoKB-style annotation of the configured hotspots (name -> fields). Tier 1: TERT promoter (brief / OncoKB), IDH1 R132,
# BRAF V600E, H3-3A K27M; H3-3A G34 tier 2.
ANNOT = {
    "TERT C228T": dict(type="Promoter", variant="C228T (promoter, -124 C>T)", tier=1, role="Oncogene", effect="Gain-of-function",
                       summary="TERT promoter hotspot mutation C228T (chr5:1295113 G>A on the + strand) creates an ETS/TCF binding motif and upregulates TERT; "
                               "a defining alteration of IDH-wildtype glioblastoma.",
                       therapeutics=None),
    "TERT C250T": dict(type="Promoter", variant="C250T (promoter, -146 C>T)", tier=1, role="Oncogene", effect="Gain-of-function",
                       summary="TERT promoter hotspot mutation C250T (chr5:1295135 G>A on the + strand) creates an ETS/TCF binding motif and upregulates TERT.",
                       therapeutics=None),
    "IDH1 R132H": dict(type="Missense", variant="p.R132H", tier=1, role="Oncogene", effect="Gain-of-function", summary="IDH1 R132H hotspot (c.395G>A)."),
    "IDH1 R132C": dict(type="Missense", variant="p.R132C", tier=1, role="Oncogene", effect="Gain-of-function", summary="IDH1 R132C hotspot (c.394C>T)."),
    "H3-3A K27M": dict(type="Missense", variant="p.K27M", tier=1, role="Oncogene", effect="Likely Gain-of-function", summary="H3-3A (H3F3A) K27M hotspot (c.83A>T)."),
    "H3-3A G34R": dict(type="Missense", variant="p.G34R", tier=2, role="Oncogene", effect="Likely Gain-of-function", summary="H3-3A (H3F3A) G34R hotspot (c.103G>A)."),
    "H3-3A G34V": dict(type="Missense", variant="p.G34V", tier=2, role="Oncogene", effect="Likely Gain-of-function", summary="H3-3A (H3F3A) G34V hotspot (c.104G>T)."),
    "BRAF V600E": dict(type="Missense", variant="p.V600E", tier=1, role="Oncogene", effect="Gain-of-function", summary="BRAF V600E hotspot (c.1799T>A)."),
}


def is_normal(c):
    return str(c).lower() == "normal"


def site_cells(calls, ev_site, site_id):
    """(carrier ids, covered-without-alt ids, source) of tumour cells for one site.
    Prefers the per-read rows of site_evidence.json (the same reads as the pooled counts); falls back to the
    caller-filtered per-cell counts of region_calls.json when the rows are capped or do not reproduce cells_alt / cells_cov."""
    cells, clone = calls["cells"], calls.get("clone") or []
    tumour = [i for i, _ in enumerate(cells) if not (i < len(clone) and is_normal(clone[i]))]
    tset = set(tumour)
    reads = (ev_site or {}).get("reads") or {}
    alt, cov = set(), set()
    pooled = (((ev_site or {}).get("presets") or {}).get("all") or {}).get("tumor") or {}
    ok = bool(reads) and not reads.get("capped") and (ev_site or {}).get("sample_every", 1) == 1
    if ok:
        for r in reads["rows"]:
            c, al = r[0], r[1]
            if c is None or c not in tset:
                continue
            if al == "A":
                alt.add(c)
                cov.add(c)
            elif al == "R":
                cov.add(c)
        ok = len(alt) == pooled.get("cells_alt") and len(cov) == pooled.get("cells_cov")
    source = "site_evidence reads (all reads)"
    if not ok:
        v = (calls.get("dna") or {}).get(site_id) or {}
        a, d = v.get("alt") or [], v.get("dp") or []
        alt = {i for i in tumour if i < len(a) and (a[i] or 0) > 0}
        cov = {i for i in tumour if (i < len(d) and (d[i] or 0) > 0) or i in alt}
        source = "region_calls (caller-filtered reads)"
    carriers = [cells[i] for i in sorted(alt)]
    ref_only = [cells[i] for i in sorted(cov - alt)]
    return carriers, ref_only, len(tumour), source


def mutation_id(h):
    return f"{h['chrom']}_{int(h['pos'])}_{h['ref']}_{h['alt']}"


def strand_only_warning(v, ev):
    """'weak' only because the strand-imbalance test fires: at tens of thousands of reads (amplified EGFR) a 46 / 54 split is
    "significant". Accepted when that is the single warning, nothing failed and the minor strand carries >= MIN_MINOR_STRAND of the alt reads."""
    if v["level"] != "weak":
        return False
    st = {c["key"]: c["status"] for c in v["checks"]}
    if st.get("strand") != "warn" or any(x in ("fail", "warn") for k, x in st.items() if k != "strand") or st.get("signal") != "pass" or st.get("cells") != "pass":
        return False
    t = ev["presets"]["all"]["tumor"]
    n = t["alt_f"] + t["alt_r"]
    return n > 0 and min(t["alt_f"], t["alt_r"]) / n >= MIN_MINOR_STRAND


def discover(calls, ev_doc, verdicts, genes, hotspot_ids=()):
    """Genotyped sites in driver genes (exon regions of the region set) that are not configured hotspots and whose site-evidence
    verdict is 'supported': candidates for annotation (SnpEff + OncoKB) and rescue."""
    out = []
    for s_ in calls.get("sites", []):
        sid = s_["id"]
        ev = (ev_doc.get("sites") or {}).get(sid)
        v = verdicts.get(sid)
        gene = (ev or {}).get("gene")
        if not ev or not v or sid in hotspot_ids or gene not in genes or len(s_["ref"]) != 1 or len(s_["alt"]) != 1:
            continue
        if v["level"] in DISCOVERED_LEVELS or strand_only_warning(v, ev):
            out.append(dict(chrom=s_["chrom"], pos=int(s_["pos"]), ref=s_["ref"], alt=s_["alt"], gene=gene, name=f"{gene} {s_['chrom']}:{s_['pos']} {s_['ref']}>{s_['alt']}", discovered=True))
    return out


def build_rows(patient, calls, ev_doc, regions, verdicts, existing_events=(), discovered=(), annotations=None):
    """Rescued event rows for one patient plus a report of every configured hotspot and discovered site (called or not, with the reason).
    annotations: { mutation id: SnpEff + OncoKB annotation } for the discovered sites (sc_annotate_snv_drivers via gos_sc_hotspot_annotate.R)."""
    annotations = annotations or {}
    rows, report = [], []
    taken = {(str(e.get("seqnames", "")).replace("chr", ""), str(e.get("start"))) for e in existing_events if e.get("source") != SOURCE and e.get("vartype") in ("SNV", "Hotspot")}
    configured = [h for h in regions.get("hotspots", []) if h.get("kind", "hotspot") == "hotspot" and h.get("source", "config") == "config"]
    # (SNVs carried over from the filtered events are already events)
    for h in list(configured) + list(discovered):
        name = h["name"]
        sid = f"{h['chrom']}:{int(h['pos'])}:{h['ref']}>{h['alt']}"
        rec = dict(patient=patient, name=name, site=sid, called=False)
        report.append(rec)
        ev = (ev_doc.get("sites") or {}).get(sid)
        verdict = verdicts.get(sid)
        if ev is None or verdict is None:
            rec["reason"] = "configured allele not genotyped (the caller reported another alt allele at this position)"
            continue
        carriers, ref_only, n_tumour, cell_src = site_cells(calls, ev, sid)
        tum = ev["presets"]["all"]["tumor"]
        nrm = ev["presets"]["all"]["normal"]
        rec.update(level=verdict["level"], gene=h["gene"], carriers=len(carriers), covered=len(carriers) + len(ref_only), alt=tum["alt_f"] + tum["alt_r"], dp=tum["alt_f"] + tum["alt_r"] + tum["ref_f"] + tum["ref_r"],
                   tumour_cells=n_tumour, tumor={k: tum[k] for k in ("alt_f", "alt_r", "ref_f", "ref_r", "cells_alt", "cells_cov")},
                   normal={k: nrm[k] for k in ("alt_f", "alt_r", "ref_f", "ref_r", "cells_alt", "cells_cov")}, verdict_summary=verdict["summary"])
        if verdict["level"] not in (DISCOVERED_LEVELS if h.get("discovered") else PASS_LEVELS) and not (h.get("discovered") and strand_only_warning(verdict, ev)):
            rec["reason"] = f"verdict {verdict['level']}"
            continue
        if len(carriers) < MIN_CARRIERS:
            rec["reason"] = f"only {len(carriers)} carrier cells"
            continue
        if (str(h["chrom"]).replace("chr", ""), str(int(h["pos"]))) in taken:
            rec["reason"] = "already in the SNV calls"
            continue
        if h.get("discovered"):
            an = annotations.get(mutation_id(h))
            if not an:
                rec["reason"] = "no protein-altering / splice annotation (not annotated or silent, intronic, UTR)"
                continue
            a = dict(type=an["type"], variant=an["protein"] or f"{h['ref']}>{h['alt']}", tier=an["tier"], role=an.get("role"), effect=an.get("effect"),
                     summary=an.get("variant_summary"), gene_summary=an.get("gene_summary"), effect_description=an.get("effect_description"),
                     therapeutics=an.get("therapeutics"), prognoses=an.get("prognoses"), diagnoses=an.get("diagnoses"), classification=an.get("classification"))
        else:
            a = ANNOT.get(name, dict(type="Missense", variant=f"{h['ref']}>{h['alt']}", tier=2, role=None, effect="Unknown", summary=name))
        covered = len(carriers) + len(ref_only)
        cov_frac = covered / n_tumour if n_tumour else 0
        chrom = str(h["chrom"]).replace("chr", "")
        row = {
            "gene": h["gene"], "fusion_genes": "", "vartype": "SNV" if h.get("discovered") else "Hotspot", "type": a["type"], "Variant": a["variant"],
            "gene_summary": a.get("gene_summary") or f"{h['gene']} hotspot recovered from read-level site evidence.", "Variant_g": f"{chrom}:{int(h['pos'])}-{int(h['pos'])} {h['ref']}>{h['alt']}",
            "Tier": a["tier"], "variant_summary": a["summary"], "effect": a["effect"], "effect_description": a.get("effect_description") or a["summary"],
            "is_multi_hit_per_gene": False, "Genome_Location": f"{chrom}:{int(h['pos'])}-{int(h['pos'])}", "seqnames": chrom,
            "start": str(int(h["pos"])), "end": str(int(h["pos"])), "sample": carriers[0], "role": a["role"],
            "fusion_cn": None, "fusion_gene_coords": None, "upstream_only": None, "upstream_partner": None,
            "therapeutics": a.get("therapeutics"), "prognoses": a.get("prognoses"), "diagnoses": a.get("diagnoses"), "id": patient,
            "alt": tum["alt_f"] + tum["alt_r"], "ref": tum["ref_f"] + tum["ref_r"], "estimated_altered_copies": None, "segment_cn": None,
            "n_cells": len(carriers), "n_normal_cells": nrm["cells_alt"], "cell_ids": ",".join(carriers),
            "VAF": round((tum["alt_f"] + tum["alt_r"]) / max(1, tum["alt_f"] + tum["alt_r"] + tum["ref_f"] + tum["ref_r"]), 4),
            "cells": f"{len(carriers)}/{covered}", "cell_fraction": round(len(carriers) / covered, 3),
            # rescue marking and read evidence
            "source": SOURCE, "rescued": True,
            "covered_ref_cell_ids": ",".join(ref_only),
            "tumour_cells": n_tumour, "covered_tumour_cells": covered, "coverage_fraction": round(cov_frac, 3), "low_coverage": cov_frac < LOW_COVERAGE,
            "hotspot_evidence": {
                "preset": "all", "cell_source": cell_src,
                "tumor": {k: tum[k] for k in ("alt_f", "alt_r", "ref_f", "ref_r", "cells_alt", "cells_cov")},
                "normal": {k: nrm[k] for k in ("alt_f", "alt_r", "ref_f", "ref_r", "cells_alt", "cells_cov")},
                "verdict": {"level": verdict["level"], "label": verdict["label"], "summary": verdict["summary"],
                            "reasons": [c["text"] for c in verdict["checks"] if c["status"] in ("warn", "fail")],
                            "checks": [{"key": c["key"], "status": c["status"], "text": c["text"]} for c in verdict["checks"]]},
            },
        }
        rows.append(row)
        rec.update(called=True, low_coverage=row["low_coverage"], variant=a["variant"], tier=a["tier"], discovered=bool(h.get("discovered")))
    return rows, report


def load(path):
    with open(path) as fh:
        return json.load(fh)


WORK = os.environ.get("WORK", "/gpfs/commons/groups/imielinski_lab/projects/single_cell_GBM/gos_datasets/gOS_GBM_sc/work/hotspot_rescue")
# region set stays in srctools (private pipeline config); override the srctools location with SRCTOOLS_GOS_SC
REGION_SET = os.path.join(os.environ.get("SRCTOOLS_GOS_SC", "/gpfs/commons/home/sclarke/git/srctools/inst/gos_single_cell"), "precompute", "config", "regions", "gbm.tsv")
ANNOTATED_COLS = ("patient", "mutation", "gene", "classification", "protein", "type", "tier", "role", "effect", "effect_description", "gene_summary", "variant_summary", "therapeutics", "prognoses", "diagnoses")


def driver_genes():
    """Genes of the region set (gene and hotspot rows of precompute/config/regions/gbm.tsv)."""
    out = set()
    with open(REGION_SET) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if f[0] in ("gene", "hotspot") and len(f) > 2 and f[2]:
                out.add(f[2])
    return out


def load_annotations():
    """{ patient: { mutation: row } } from WORK/annotation.tsv (gos_sc_hotspot_annotate.R); empty before that step has run."""
    path = os.path.join(WORK, "annotation.tsv")
    out = {}
    if not os.path.exists(path):
        return out
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r = {k: (None if v in ("NA", "") else v) for k, v in r.items()}
            if r.get("tier") is not None:
                try:
                    r["tier"] = int(float(r["tier"]))
                except ValueError:
                    r["tier"] = None
            out.setdefault(r["patient"], {})[r["mutation"]] = r
    return out


def run_patient(data, patient, gos_src, tag, dry, mode="call", annotations=None):
    pre = os.path.join(data, patient, "precompute")
    if not (os.path.exists(os.path.join(pre, "site_evidence.json")) and os.path.exists(os.path.join(pre, "region_calls.json"))):
        return None, [dict(patient=patient, name="(all)", called=False, reason="no precompute/site_evidence.json (site evidence step not run)")]
    calls, ev_doc = load(os.path.join(pre, "region_calls.json")), load(os.path.join(pre, "site_evidence.json"))
    regions = load(os.path.join(data, patient, "slices", "regions.json"))
    verdicts = json.loads(subprocess.run(["node", os.path.join(os.path.dirname(os.path.abspath(__file__)), "hotspot_verdicts.js"), gos_src, os.path.join(data, patient)],
                                         check=True, capture_output=True, text=True).stdout)
    path = os.path.join(data, patient, "filtered.events.json")
    events = load(path)
    kept = [e for e in events if e.get("source") != SOURCE]
    hot_ids = {f"{h['chrom']}:{int(h['pos'])}:{h['ref']}>{h['alt']}" for h in regions.get("hotspots", [])}
    disc = discover(calls, ev_doc, verdicts, driver_genes(), hot_ids)
    if mode == "select":     # candidates for gos_sc_hotspot_annotate.R
        os.makedirs(WORK, exist_ok=True)
        with open(os.path.join(WORK, f"{patient}.candidates.txt"), "w") as fh:
            fh.write("\n".join(mutation_id(h) for h in disc) + "\n")
        return [], [dict(patient=patient, name=h["name"], called=False, reason="candidate for annotation") for h in disc]
    rows, report = build_rows(patient, calls, ev_doc, regions, verdicts, kept, discovered=disc, annotations=(annotations or {}).get(patient, {}))
    # dropout-aware clonality (clonality.py): cell_fraction is a lower bound under allelic dropout, not a clonality call
    for r, res in cl.annotate_rows(data, patient, rows):
        rec = next((x for x in report if x.get("called") and x.get("gene") == r["gene"] and x.get("variant") == r["Variant"]), None)
        if rec is not None:
            rec.update(clonality=res["status"], clonality_method=cl.METHOD, clonality_reason=res["reason"])
    if not dry:
        bak = f"{path}.{tag}"
        if not os.path.exists(bak):
            shutil.copy2(path, bak)
        tmp = path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(kept + rows, fh)
        os.replace(tmp, path)
    return rows, report


def write_summary(data, report):
    """data/_cohort/hotspot_rescue.json: for every patient the configured hotspots with their read evidence, called or not
    (the cohort headline counts patients with a rescued call among those assessed). Patients not in this run keep their entry."""
    path = os.path.join(data, "_cohort", "hotspot_rescue.json")
    doc = load(path) if os.path.exists(path) else {"patients": {}}
    by = {}
    for r in report:
        by.setdefault(r["patient"], []).append({k: v for k, v in r.items() if k != "patient"})
    for p, sites in by.items():
        assessed = not (len(sites) == 1 and sites[0]["name"] == "(all)")
        doc["patients"][p] = {"assessed": assessed, "reason": None if assessed else sites[0].get("reason"), "sites": sites if assessed else []}
    doc.update(format="gos-sc-hotspot-rescue/1", created=time.strftime("%Y-%m-%dT%H:%M:%S"), source=SOURCE)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if os.path.exists(path) and not os.path.exists(path + ".bak_261015c"):
        shutil.copy2(path, path + ".bak_261015c")
    with open(path + ".tmp", "w") as fh:
        json.dump(doc, fh)
    os.replace(path + ".tmp", path)


def main():
    data = os.environ.get("DATA_DIR", "/gpfs/commons/projects/imielinski_web/sclarke/gOS_GBM_sc/data")
    gos_src = os.environ.get("GOS_SRC", "/gpfs/commons/groups/imielinski_lab/home/sclarke/git/gOS_dev_wt_hotspotrescue")
    tag = os.environ.get("BACKUP_TAG", "bak_hotspot_rescue")
    dry = os.environ.get("DRY_RUN") == "1"
    pats = [p for p in os.environ.get("PATIENTS", "").split(",") if p] or sorted(
        p for p in os.listdir(data) if p[:3] in ("BWH", "MGH") and p[3:].isdigit())
    allrep = []
    mode = os.environ.get("MODE", "call")
    ann = load_annotations()
    for p in pats:
        rows, report = run_patient(data, p, gos_src, tag, dry, mode, ann)
        allrep += report
        for r in report:
            print(p, r["name"], "CALLED" if r.get("called") else "no", {k: v for k, v in r.items() if k not in ("patient", "name", "called")}, flush=True)
    if not dry and mode == "call":
        write_summary(data, allrep)
    print(json.dumps([r for r in allrep if r.get("called")]))


if __name__ == "__main__":
    sys.exit(main())
