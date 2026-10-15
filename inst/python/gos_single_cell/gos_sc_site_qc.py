#!/usr/bin/env python3
"""Export per-site QC scores to <outdir>/<P>/snv_qc.json, keyed by snv_matrix variant id.
usage: gos_sc_site_qc.py PATIENT SCORES_TSV_GZ SNV_MATRIX_JSON OUT_JSON
"""
import csv, gzip, json, math, sys
patient, scores, matrix, out = sys.argv[1:5]
FIELDS = ["sum_logodds", "n_alt", "n_ref", "n_cov", "frac_informative", "pooled_alt_vaf", "coh_q",
          "clonal_llr", "mean_llr_ref", "n_ok_m1", "n_ok_0", "n_ok_p1", "basic"]
ids = {v["id"] for v in json.load(open(matrix))["variants"]}
def num(s):
    if s in ("", "NA", "NaN", "Inf", "-Inf"): return None
    try:
        x = float(s)
    except ValueError:
        return None
    return round(x, 4) if math.isfinite(x) else None
sites = {}
with gzip.open(scores, "rt") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        m = r["mutation"]
        if m not in ids: continue
        rec = {}
        for f in FIELDS:
            v = r.get(f)
            if v is None: continue
            rec[f] = (v == "TRUE") if f == "basic" else num(v)
        sites[m] = rec
json.dump({"schemaVersion": 1, "patient": patient, "fields": FIELDS, "sites": sites}, open(out, "w"), separators=(",", ":"))
print(f"{patient}: {len(sites)} of {len(ids)} matrix sites have QC")
