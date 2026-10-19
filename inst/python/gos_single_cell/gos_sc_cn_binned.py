#!/usr/bin/env python3
"""Per-cell binned copy number (total, major, minor, BAF) for paper-figure CN / BAF heatmaps -> <out>/<P>/cn_binned.json.

Source = the per-cell files gOS draws its heatmaps from: <cell>/complex.json (total CN segments, gOS cnRowFromGenome) and
<cell>/allelic.json (major/minor CN per segment, gOS allelicRowFromAllelic). A bin takes the value of the segment that covers the
bin centre (as gOS binAt); a bin whose centre lies in no segment is missing. BAF = minor / (major + minor) of the same bin.
Two levels (hg38, chromosomes 1-22, X, Y; bins restart at each chromosome start, the last bin of a chromosome is shorter):
  1 Mb   whole genome, inline in cn_binned.json (arrays per cell: null = missing)
  100 kb whole genome, columnar binary next to it: cn_binned_100kb.{cn,major,minor}.u16 (little-endian uint16, row-major
         cells x bins in the order of cn_binned.json cells, 65535 = missing). Region subsets are slices of the bin axis.
All cells of the patient in datafiles.json that have a complex.json are written (the `normal` flag marks Normal clones).

usage: gos_sc_cn_binned.py [--data DIR] [--out DIR] [--patients P ...]
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gos_sc_figlib as F  # noqa: E402

MISSING = 65535
LEVELS = {"1Mb": 1_000_000, "100kb": 100_000}


def bin_axis(size):
    """per chromosome bins: (chrom, start, end, global centre) arrays plus chromosome index ranges"""
    chroms, starts, ends, centres, ranges = [], [], [], [], []
    k = 0
    for c in F.CHROMS:
        L = F.LEN[c]
        n = (L + size - 1) // size
        s = np.arange(n, dtype=np.int64) * size + 1
        e = np.minimum(s + size - 1, L)
        chroms += [c] * n
        starts.append(s)
        ends.append(e)
        centres.append(F.OFFSET[c] + 1 + (s + e) / 2.0)
        ranges.append({"chr": c, "first_bin": k, "n_bins": int(n)})
        k += n
    return chroms, np.concatenate(starts), np.concatenate(ends), np.concatenate(centres), ranges


def _u16(v):
    out = np.full(v.shape, MISSING, dtype=np.uint16)
    ok = np.isfinite(v)
    out[ok] = np.clip(np.rint(v[ok]), 0, MISSING - 1).astype(np.uint16)
    return out


def cell_rows(data, cell_id):
    cx = os.path.join(data, cell_id, "complex.json")
    if not os.path.exists(cx):
        return None
    tot = F.total_cn_row(F.jload(cx))
    ap = os.path.join(data, cell_id, "allelic.json")
    maj = mnr = None
    if os.path.exists(ap):
        maj, mnr = F.allelic_rows(F.jload(ap))
    return tot, maj, mnr


def export_patient(data, outdir, patient, datafiles):
    cells = F.patient_cells(datafiles, patient)
    ids, meta, rows = [], [], []
    for c in cells:
        r = cell_rows(data, c["cell_id"])
        if r is None:
            continue
        ids.append(c["cell_id"])
        meta.append(c)
        rows.append(r)
    if not ids:
        print(patient, "no cells with complex.json", flush=True)
        return None
    out = {
        "format": "gos-sc-cn-binned/1", "patient": patient, "genome": "hg38", "n_cells": len(ids), "cells": ids,
        "clone_id": [m["clone_id"] for m in meta], "normal": [F.is_normal_clone(m["clone_id"]) for m in meta],
        "value_source": "bin-centre value of complex.json (total) and allelic.json (major, minor) segments, as gOS binAt",
        "baf": "minor / (major + minor)", "has_allelic": [r[1] is not None for r in rows], "levels": {},
    }
    pdir = os.path.join(outdir, patient)
    os.makedirs(pdir, exist_ok=True)
    for name, size in LEVELS.items():
        chroms, starts, ends, centres, ranges = bin_axis(size)
        nb = len(centres)
        mats = {k: np.full((len(ids), nb), np.nan) for k in ("cn", "major", "minor")}
        for i, (tot, maj, mnr) in enumerate(rows):
            mats["cn"][i] = tot.at(centres)
            if maj is not None:
                mats["major"][i] = maj.at(centres)
                mats["minor"][i] = mnr.at(centres)
        lv = {"bin_size": size, "n_bins": nb, "chromosomes": ranges}
        if name == "1Mb":
            lv["bin_start"] = starts.tolist()
            lv["bin_end"] = ends.tolist()
            lv["bin_chr"] = chroms
            for k, m in mats.items():
                lv[k] = [[None if not np.isfinite(x) else int(round(x)) for x in row] for row in m]
            den = mats["major"] + mats["minor"]
            with np.errstate(invalid="ignore", divide="ignore"):
                baf = np.where(den > 0, mats["minor"] / den, np.nan)
            lv["baf"] = [[None if not np.isfinite(x) else round(float(x), 4) for x in row] for row in baf]
        else:
            lv["files"] = {}
            for k, m in mats.items():
                fn = "cn_binned_%s.%s.u16" % (name, k)
                _u16(m).astype("<u2").tofile(os.path.join(pdir, fn))
                lv["files"][k] = fn
            lv["missing_value"] = MISSING
            lv["bin_index_note"] = "bin b of chromosome c: chromosomes[c].first_bin + (start-1)//bin_size; start = k*bin_size+1"
        out["levels"][name] = lv
    with open(os.path.join(pdir, "cn_binned.json"), "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    print(patient, len(ids), "cells", flush=True)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", default=os.path.join(F.DATA_ROOT, "data"))
    ap.add_argument("--datafiles", default=None, help="datafiles.json (default <data>/../datafiles.json)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--patients", nargs="*", default=["BWH69", "BWH70", "BWH75", "MGH285", "MGH301", "MGH302", "MGH303", "MGH304", "MGH305"])
    a = ap.parse_args()
    datafiles = F.jload(a.datafiles or os.path.join(os.path.dirname(a.data.rstrip("/")), "datafiles.json"))
    for p in a.patients:
        export_patient(a.data, a.out, p, datafiles)


if __name__ == "__main__":
    main()
