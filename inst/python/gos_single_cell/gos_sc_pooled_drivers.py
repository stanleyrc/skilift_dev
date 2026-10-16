#!/usr/bin/env python3
"""Import high-confidence pooled-scan driver calls (srctools feature/driver-sensitivity) into data/<P>/filtered.events.json.

  python3 gos_sc_pooled_drivers.py --data <gOS data dir> --candidates new_driver_candidates.tsv [--patients A,B] [--dry-run]
         [--scan-dir analysis/driver_sensitivity_261016] [--out-json dryrun.json] [--tag bak_<yymmdd>_pooled]

Rows with verdict == "high-confidence" are appended with the field set of the hotspot-rescue rows. Carrier cell ids are recovered from the
scan's per-shard observation tables (<scan-dir>/<P>/s*/obs.tsv.gz: one row per alt read, columns cell, normal, chrom, pos, type, ref, alt),
i.e. tumour cells with >= 1 alt read at the site (carriers) and normal cells with >= 1 alt read (n_normal_cells).
OncoKB fields come from <scan-dir>/annot/annotation.tsv (patient + chr_pos_ref_alt); where absent they stay null (never invented).
Exception: the TERT promoter hotspot uses the same curated C228T text as the existing hotspot-rescue rows.
Idempotent: events with the same gene + Genome_Location + alt allele (not from this script) are skipped; re-running replaces
earlier source == "pooled driver scan" rows. The previous file is kept as filtered.events.json.<tag> (first write only).
type uses the vocabulary the front end groups on (Missense / Trunc / Splice / Promoter); the finer class is in `consequence`.
"""
import argparse
import csv
import glob
import gzip
import json
import os
import shutil
import sys
import time

SOURCE = "pooled driver scan"
SCAN_SHA = "02065e8"
DEFAULT_SCAN = "/gpfs/commons/groups/imielinski_lab/projects/single_cell_GBM/analysis/driver_sensitivity_261016"
TERT_C228T = dict(type="Promoter", variant="C228T (promoter, -124 C>T)", tier=1, role="Oncogene", effect="Gain-of-function",
                  summary="TERT promoter hotspot mutation C228T (chr5:1295113 G>A on the + strand) creates an ETS/TCF binding motif and upregulates TERT; "
                          "a defining alteration of IDH-wildtype glioblastoma.",
                  gene_summary="TERT hotspot recovered from read-level site evidence.")
AA3 = dict(Ala="A", Arg="R", Asn="N", Asp="D", Cys="C", Gln="Q", Glu="E", Gly="G", His="H", Ile="I", Leu="L", Lys="K", Met="M", Phe="F",
           Pro="P", Ser="S", Thr="T", Trp="W", Tyr="Y", Val="V", Ter="*")


def nz(v):
    return None if v is None or str(v).strip() in ("", "NA") else v


def short_protein(p):
    """p.Arg209fs -> p.R209fs"""
    if not p or not p.startswith("p."):
        return p
    s, out, i = p[2:], "p.", 0
    while i < len(s):
        if s[i:i + 3] in AA3:
            out += AA3[s[i:i + 3]]
            i += 3
        else:
            out += s[i]
            i += 1
    return out


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_annotations(scan):
    out = {}
    path = os.path.join(scan, "annot", "annotation.tsv")
    if os.path.exists(path):
        for r in read_tsv(path):
            r = {k: nz(v) for k, v in r.items()}
            if r.get("tier") is not None:
                r["tier"] = int(float(r["tier"]))
            out[(r["patient"], r["mutation"])] = r
    return out


def snpeff(scan, patient, chrom, pos, ref, alt):
    """First SnpEff ANN entry (HGVS.c, HGVS.p, consequence) of the scan's annotated sites VCF, or None."""
    p = os.path.join(scan, "annot", patient, f"{patient}.sites.ann.vcf")
    if not os.path.exists(p):
        return None
    key = f"{chrom.replace('chr', '')}\t{pos}\t"
    with open(p) as fh:
        for ln in fh:
            if ln.startswith(key) and ln.split("\t")[3] == ref and ln.split("\t")[4] == alt and "ANN=" in ln:
                f = ln.split("ANN=")[1].split("\t")[0].split(",")[0].split("|")
                return dict(consequence=f[1], hgvs_c=f[9] or None, hgvs_p=f[10] or None)
    return None


def carriers_from_obs(scan, patient, chrom, pos, vtype, ref, alt):
    """(tumour carrier ids, normal carrier ids, alt reads per strand {tumor: {alt_f, alt_r}, normal: {...}}) from the per-shard obs tables."""
    tum, nrm = set(), set()
    strand = {"tumor": {"alt_f": 0, "alt_r": 0}, "normal": {"alt_f": 0, "alt_r": 0}}
    for f in sorted(glob.glob(os.path.join(scan, patient, "s*", "obs.tsv.gz"))):
        with gzip.open(f, "rt") as fh:
            next(fh)
            for ln in fh:
                x = ln.rstrip("\n").split("\t")
                if x[3] == str(pos) and x[2] == chrom and x[4] == vtype and x[5] == ref and x[6] == alt:
                    (nrm if x[1] == "True" else tum).add(x[0])
                    if len(x) > 7 and x[7] in "+-":
                        strand["normal" if x[1] == "True" else "tumor"]["alt_f" if x[7] == "+" else "alt_r"] += 1
    return sorted(tum), sorted(nrm), strand


def patient_cells(data, scan, patient):
    cells = set()
    p = os.path.join(data, patient, "precompute", "region_calls.json")
    if os.path.exists(p):
        with open(p) as fh:
            cells |= set(json.load(fh).get("cells", []))
    for f in glob.glob(os.path.join(scan, patient, "s*", "cells.tsv")):
        cells |= {r["cell"] for r in read_tsv(f)}
    return cells


# tumour suppressors for the truncating-mutation rule when OncoKB has no entry for the exact variant
RULE_TSGS = {"TP53", "NF1", "PTEN", "RB1", "STAG2", "ATRX", "CDKN2A", "PIK3R1", "BCOR", "SETD2", "ARID1A", "KMT2D"}


def classify(c, snp):
    cons = c["consequence"]
    if c["gene"] == "TERT" and c["pos"] == "1295113":
        return "Promoter", "Promoter"
    if "splice" in cons:
        return "Splice", "Splice_Site"
    if "frameshift" in cons:
        return "Trunc", "Frame_Shift_Del" if c["type"] == "DEL" else "Frame_Shift_Ins"
    if c["protein"] and c["protein"].endswith("*"):
        return "Trunc", "Nonsense"
    return "Missense", "Missense"


def n_tumour_cells(scan, patient):
    """Tumour (non-normal) cells in the scan, from the per-shard cells.tsv tables; None if absent."""
    cells = {}
    for f in glob.glob(os.path.join(scan, patient, "s*", "cells.tsv")):
        for r in read_tsv(f):
            cells[r["cell"]] = r.get("normal")
    return sum(v not in ("1", "True", "true") for v in cells.values()) if cells else None


def build_row(c, carriers, normal, ann=None, snp=None, sha=SCAN_SHA, strand=None, n_tumour=None):
    """One event row from a candidates.tsv row `c` (dict of str), carrier / normal cell ids, annotation row, SnpEff entry."""
    chrom, pos, ref, alt = c["chrom"].replace("chr", ""), int(c["pos"]), c["ref"], c["alt"]
    typ, cons = classify(c, snp)
    start, end = pos, pos + len(ref) - 1
    if c["gene"] == "TERT" and typ == "Promoter":
        a = TERT_C228T
        ann = dict(tier=a["tier"], role=a["role"], effect=a["effect"], variant_summary=a["summary"], effect_description=a["summary"], gene_summary=a["gene_summary"])
        variant = a["variant"]
    else:
        ann = ann or {}
        if typ == "Splice":
            variant = (snp or {}).get("hgvs_c") or None
        else:
            variant = nz(c.get("protein")) or short_protein((snp or {}).get("hgvs_p")) or (snp or {}).get("hgvs_c")
        variant = variant or f"{ref}>{alt}"
    covered, nc = int(c["covered"]), len(carriers)
    t_alt, t_dp = int(c["tumour_alt"]), int(c["tumour_dp"])
    tier = ann.get("tier")
    rule = None
    if tier is None and typ == "Trunc" and c["gene"] in RULE_TSGS:
        # no OncoKB entry for this exact variant: OncoKB's rule for truncating mutations in a tumour suppressor, labelled as such
        tier, rule = 2, f"Not in the OncoKB lookup; {c['gene']} is a tumour suppressor and this is a truncating mutation (OncoKB rule: likely loss-of-function)."
        ann = dict(ann, role=ann.get("role") or "Tumor Suppressor", effect="Likely Loss-of-function", variant_summary=rule)
    return {
        "gene": c["gene"], "fusion_genes": "", "vartype": "Hotspot" if typ == "Promoter" else ("INDEL" if c["type"] in ("DEL", "INS") else "SNV"), "type": typ, "consequence": cons,
        "Variant": variant, "gene_summary": ann.get("gene_summary"), "Variant_g": f"{chrom}:{start}-{end} {ref}>{alt}", "Tier": tier,
        "variant_summary": ann.get("variant_summary"), "effect": ann.get("effect"), "effect_description": ann.get("effect_description"),
        "is_multi_hit_per_gene": False, "Genome_Location": f"{chrom}:{start}-{end}", "seqnames": chrom, "start": str(start), "end": str(end),
        "sample": carriers[0] if carriers else None, "role": ann.get("role"), "fusion_cn": None, "fusion_gene_coords": None, "upstream_only": None,
        "upstream_partner": None, "therapeutics": ann.get("therapeutics"), "prognoses": ann.get("prognoses"), "diagnoses": ann.get("diagnoses"),
        "id": c["patient"], "alt": t_alt, "ref": t_dp - t_alt, "estimated_altered_copies": None, "segment_cn": None,
        "n_cells": nc, "n_normal_cells": len(normal), "cell_ids": ",".join(carriers),
        "VAF": round(t_alt / t_dp, 4) if t_dp else None, "cells": f"{nc}/{covered}", "cell_fraction": round(nc / covered, 3) if covered else None,
        "source": SOURCE, "pooled_scan": True, "rescued": True, "annotation_rule": rule,
        # fields the gOS rescued-event views read (as gos_sc_hotspot_rescue writes them); per-cell depth is pooled, so there are no covered_ref_cell_ids
        "covered_tumour_cells": str(covered), "tumour_cells": str(n_tumour) if n_tumour else None,
        "coverage_fraction": str(round(covered / n_tumour, 3)) if n_tumour else None, "low_coverage": "False",
        "hotspot_evidence": {
            "preset": "pooled", "cell_source": "pooled driver scan (pooled depth; per-cell depth not kept)",
            "tumor": dict((strand or {}).get("tumor") or {}, alt=t_alt, dp=t_dp),
            "normal": dict((strand or {}).get("normal") or {}, alt=int(c["normal_alt"]), dp=int(c["normal_dp"])),
            "verdict": {"level": "supported", "label": "pooled scan: high-confidence", "summary": c.get("why"), "reasons": []},
        },
        "evidence": {
            "scan_sha": sha, "covered": covered, "carriers": int(c["carriers"]),
            "carrier_vaf_median": float(c["carrier_vaf_median"]) if nz(c.get("carrier_vaf_median")) else None,
            "carriers_2plus_frac": float(c["carriers_2plus_frac"]) if nz(c.get("carriers_2plus_frac")) else None,
            "clade_best": nz(c.get("clade_best")), "clade_f1": float(c["clade_f1"]) if nz(c.get("clade_f1")) else None,
            "normal_alt": int(c["normal_alt"]), "normal_dp": int(c["normal_dp"]),
            "rna_alt": int(c["rna_alt"]), "rna_dp": int(c["rna_dp"]), "rna_cells": int(c["rna_alt_cells"]),
            "other_patients_alt_cells": int(c["other_patients_alt_cells"]), "why": c.get("why"),
        },
    }


def alt_of(e):
    vg = e.get("Variant_g") or ""
    return vg.split(" ")[-1].split(">")[-1] if ">" in vg else None


def merge(events, new_rows):
    """(kept events + rows to add, list of (row, reason) skipped). Our earlier rows are dropped first."""
    kept = [e for e in events if e.get("source") != SOURCE]
    out, skipped = list(kept), []
    for r in new_rows:
        if any(e.get("gene") == r["gene"] and e.get("Genome_Location") == r["Genome_Location"] and alt_of(e) == alt_of(r) for e in out):
            skipped.append((r, "event with same gene + Genome_Location + alt exists"))
        else:
            out.append(r)
    return out, skipped


def run(data, candidates, patients=None, dry=False, scan=DEFAULT_SCAN, tag=None, sha=SCAN_SHA, log=print):
    tag = tag or f"bak_{time.strftime('%y%m%d')}_pooled"
    rows = [r for r in read_tsv(candidates) if r["verdict"] == "high-confidence" and (not patients or r["patient"] in patients)]
    anns = load_annotations(scan)
    result, bypat = [], {}
    for c in rows:
        mid = f"{c['chrom']}_{c['pos']}_{c['ref']}_{c['alt']}"
        tum, nrm, strand = carriers_from_obs(scan, c["patient"], c["chrom"], c["pos"], c["type"], c["ref"], c["alt"])
        row = build_row(c, tum, nrm, anns.get((c["patient"], mid)), snpeff(scan, c["patient"], c["chrom"], c["pos"], c["ref"], c["alt"]), sha,
                        strand=strand, n_tumour=n_tumour_cells(scan, c["patient"]))
        cells = patient_cells(data, scan, c["patient"])
        row["_check"] = dict(carriers_expected=int(c["carriers"]), carriers_found=len(tum), cell_ids_in_patient=all(x in cells for x in tum), n_unknown=sum(x not in cells for x in tum))
        bypat.setdefault(c["patient"], []).append(row)
    for p, new in bypat.items():
        path = os.path.join(data, p, "filtered.events.json")
        with open(path) as fh:
            events = json.load(fh)
        clean = [{k: v for k, v in r.items() if k != "_check"} for r in new]
        merged, skipped = merge(events, clean)
        for r in new:
            reason = next((s for x, s in skipped if x["Genome_Location"] == r["Genome_Location"] and x["gene"] == r["gene"]), None)
            result.append(dict(row=r, skipped=reason))
            log(p, r["gene"], r["Variant"], "SKIP: " + reason if reason else "ADD", "n_cells", r["n_cells"], r["_check"])
        if not dry:
            bak = f"{path}.{tag}"
            if not os.path.exists(bak):
                shutil.copy2(path, bak)
            with open(path + ".tmp", "w") as fh:
                json.dump(merged, fh)
            os.replace(path + ".tmp", path)
    return result


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--patients", default="")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--scan-dir", default=DEFAULT_SCAN)
    ap.add_argument("--out-json", default=None, help="write the would-be / written rows here")
    ap.add_argument("--tag", default=None)
    ap.add_argument("--sha", default=SCAN_SHA)
    a = ap.parse_args(argv)
    res = run(a.data, a.candidates, set(p for p in a.patients.split(",") if p), a.dry_run, a.scan_dir, a.tag, a.sha)
    if a.out_json:
        with open(a.out_json, "w") as fh:
            json.dump(res, fh, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
