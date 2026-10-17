#!/usr/bin/env python3
"""Size of the copy-number-0 segment behind each HOMDEL event, per patient, written onto the event.

A homozygous deletion of a tumour suppressor is focal (CDKN2A here: 0.2-2 Mb). A CN-0 call over a whole chromosome or
arm (e.g. all of chr22 in many BWH75 cells) cannot be a viable homozygous deletion: it is a copy-number calling problem
in those cells, so the event is flagged and gOS does not count it as a biallelic hit.

For every HOMDEL row of data/<P>/filtered.events.json: for each carrier cell, the length of the per-cell graph interval
(data/<cell>/complex.json, y == 0) overlapping the gene. Adds
  cn0_segment_mb_median, cn0_segment_mb_max, cn0_cells_measured, broad_cn0 (median >= BROAD_MB).
Backs up filtered.events.json as .bak_<yymmdd>_homdelqc. Also writes _cohort/broad_cn0.json: per patient and gene,
the carrier cells whose CN-0 call is broad, for gOS to drop from the per-cell allelic HOMDEL calls.

usage: gos_sc_homdel_qc.py [--data <gOS data dir>] [--patients P ...] [--dry-run]
"""
import argparse
import json
import os
import re
import shutil
import statistics
import time

DATA = "/gpfs/commons/projects/imielinski_web/sclarke/gOS_GBM_sc/data"
BROAD_MB = 20.0


def cn0_length(intervals, chrom, start, end):
    """Length in Mb of the CN-0 interval overlapping [start, end] on chrom, or None."""
    for i in intervals:
        if str(i.get("chromosome")).replace("chr", "") == chrom and i.get("y") == 0 and i["startPoint"] <= end and i["endPoint"] >= start:
            return (i["endPoint"] - i["startPoint"]) / 1e6
    return None


def annotate(event, load_intervals, broad_cells=None):
    chrom, rng = str(event["Genome_Location"]).split(":")
    start, end = (int(x) for x in rng.split("-"))
    lens = []
    for c in [c for c in str(event.get("cell_ids") or "").split(",") if c]:
        iv = load_intervals(c)
        if iv is None:
            continue
        x = cn0_length(iv, chrom.replace("chr", ""), start, end)
        if x is not None:
            lens.append(x)
            if broad_cells is not None and x >= BROAD_MB:
                broad_cells.append(c)
    med = round(statistics.median(lens), 3) if lens else None
    event.update(cn0_segment_mb_median=med, cn0_segment_mb_max=round(max(lens), 3) if lens else None, cn0_cells_measured=len(lens),
                 broad_cn0=bool(med is not None and med >= BROAD_MB))
    return event


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=DATA)
    ap.add_argument("--patients", nargs="*")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)
    cache = {}

    def load(cell):
        if cell not in cache:
            f = os.path.join(a.data, cell, "complex.json")
            cache[cell] = json.load(open(f)).get("intervals", []) if os.path.isfile(f) else None
        return cache[cell]

    cohort = {"format": "gos-sc-broad-cn0/1", "created": time.strftime("%Y-%m-%dT%H:%M:%S"), "broad_mb": BROAD_MB,
              "method": "per patient and gene: carrier cells of a HOMDEL event whose CN-0 interval over the gene is >= broad_mb (chromosome / arm scale)",
              "patients": {}}
    pts = a.patients or sorted(p for p in os.listdir(a.data) if re.fullmatch(r"[A-Z]+[0-9]+", p))
    for p in pts:
        f = os.path.join(a.data, p, "filtered.events.json")
        if not os.path.isfile(f):
            continue
        ev = json.load(open(f))
        n = nb = 0
        for e in ev:
            if str(e.get("vartype")).upper() == "HOMDEL" and e.get("Genome_Location") and e.get("cell_ids"):
                bc = []
                annotate(e, load, bc)
                if bc:
                    cohort["patients"].setdefault(p, {})[e.get("gene")] = {"mb": e["cn0_segment_mb_median"], "cells": bc}
                n += 1
                nb += e["broad_cn0"]
                if e["broad_cn0"] and (e.get("cell_fraction") or 0) >= 0.1:
                    print(f"  {p} {e.get('gene')}: CN-0 over {e['cn0_segment_mb_median']} Mb in {e.get('n_cells')} cells -> broad_cn0")
        print(p, "HOMDEL events", n, "broad", nb)
        if not a.dry_run and n:
            bak = f"{f}.bak_{time.strftime('%y%m%d')}_homdelqc"
            if not os.path.exists(bak):
                shutil.copy2(f, bak)
            with open(f + ".tmp", "w") as fh:
                json.dump(ev, fh)
            os.replace(f + ".tmp", f)
        cache.clear()
    out = os.path.join(a.data, "_cohort", "broad_cn0.json")
    if not a.dry_run:
        if os.path.exists(out):
            shutil.copy2(out, f"{out}.bak_{time.strftime('%y%m%d')}")
        with open(out, "w") as fh:
            json.dump(cohort, fh, separators=(",", ":"))
        print("wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
