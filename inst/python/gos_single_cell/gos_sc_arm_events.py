#!/usr/bin/env python3
"""Arm-level and chromosome-level copy-number gains and losses per patient (data/<P>/arm_events.json).

Input: data/<P>/cn_binned.json (gos-sc-cn-binned/1, the 1 Mb level: integer total copy number per cell per bin) and
a cytoband table (chromosome, startPoint, endPoint, title) that gives each chromosome's p / q arm.
The only cell label used is tumour vs normal (cn_binned `normal` flag or clone "Normal"); normal cells are left out.
No tree, no earlier event calls.

Rule (the same one the gOS convergence view uses): per cell, the median copy number over the 1 Mb bins of each arm
(bin centre inside the arm); baseline = round(median over the cell's autosomal arm medians);
gain if arm median - baseline >= +0.75, loss if <= -0.75. A chromosome-level event is a cell whose two arms both
carry the same call. Arms with fewer than MIN_BINS bins (acrocentric p arms) and chromosomes X / Y / M are not assessed.

Output events[]: {id, chrom, arm ("p" | "q" | "whole"), level ("arm" | "chromosome"), state ("gain" | "loss"), locus,
n_cells, n_tumour_cells, cell_fraction, median_cn_delta, cells[]}.

usage: gos_sc_arm_events.py --data <gOS data dir> --out <dir> [--patients P ...] [--cytobands hg38.tsv]
Writes <out>/<P>/arm_events.json (existing file backed up as .bak_<yymmdd>). --out may be a staging dir.
"""
import argparse
import csv
import json
import os
import re
import shutil
import statistics
import time

DATA = "/gpfs/commons/projects/imielinski_web/sclarke/gOS_GBM_sc/data"
CYTOBANDS = "/gpfs/commons/projects/imielinski_web/sclarke/gOS_dev/cytobands/hg38.tsv"
GAIN, LOSS = 0.75, -0.75
MIN_BINS = 3
AUTOSOMES = [str(i) for i in range(1, 23)]


def bare(c):
    return re.sub(r"^chr", "", str(c))


def read_arms(path):
    """{chrom: {"p": (start, end), "q": (start, end)}} from a cytoband table (p arm ends where the last p band ends)."""
    bands = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            c = bare(r["chromosome"])
            if c not in AUTOSOMES:
                continue
            arm = str(r["title"])[:1]
            if arm not in ("p", "q"):
                continue
            s, e = int(r["startPoint"]), int(r["endPoint"])
            lo, hi = bands.setdefault(c, {}).get(arm, (s, e))
            bands[c][arm] = (min(lo, s), max(hi, e))
    return bands


def median(v):
    v = [x for x in v if x is not None]
    return statistics.median(v) if v else None


def arm_events(cn, arms, normal_ids=()):
    """Events of one patient from a cn_binned dict (1 Mb level) and read_arms() output. Returns (events, n_tumour, notes)."""
    lvl = cn["levels"]["1Mb"]
    cells = cn["cells"]
    bin_chr = [bare(c) for c in lvl["bin_chr"]]
    bin_mid = [(s + e) / 2 for s, e in zip(lvl["bin_start"], lvl["bin_end"])]
    arm_bins = {}
    for c, parts in arms.items():
        for arm, (lo, hi) in parts.items():
            idx = [i for i in range(len(bin_chr)) if bin_chr[i] == c and lo <= bin_mid[i] <= hi]
            if len(idx) >= MIN_BINS:
                arm_bins[(c, arm)] = idx
    flags = cn.get("normal") or [False] * len(cells)
    clones = cn.get("clone_id") or [None] * len(cells)
    tumour = [k for k in range(len(cells)) if not flags[k] and clones[k] != "Normal" and cells[k] not in set(normal_ids)]
    calls = {}  # (chrom, arm, state) -> {"cells": [], "deltas": []}
    chrom_calls = {}
    for k in tumour:
        row = lvl["cn"][k]
        meds = {key: median([row[i] for i in idx]) for key, idx in arm_bins.items()}
        base_v = median([m for m in meds.values() if m is not None])
        if base_v is None:
            continue
        base = round(base_v)
        per_chrom = {}
        for (c, arm), m in meds.items():
            if m is None:
                continue
            d = m - base
            st = "gain" if d >= GAIN else "loss" if d <= LOSS else None
            per_chrom.setdefault(c, {})[arm] = (st, d)
            if st:
                e = calls.setdefault((c, arm, st), {"cells": [], "deltas": []})
                e["cells"].append(cells[k])
                e["deltas"].append(d)
        for c, parts in per_chrom.items():
            if len(parts) == 2 and parts["p"][0] and parts["p"][0] == parts["q"][0]:
                e = chrom_calls.setdefault((c, parts["p"][0]), {"cells": [], "deltas": []})
                e["cells"].append(cells[k])
                e["deltas"].append((parts["p"][1] + parts["q"][1]) / 2)
    n = len(tumour)
    out = []

    def add(c, arm, level, st, e):
        if level == "arm":
            lo, hi = arms[c][arm]
        else:
            lo, hi = min(v[0] for v in arms[c].values()), max(v[1] for v in arms[c].values())
        out.append({"id": f"{c}{arm if level == 'arm' else ''}_{st}" if level == "arm" else f"{c}_{st}", "chrom": c, "arm": arm, "level": level, "state": st,
                    "locus": f"chr{c}:{lo}-{hi}", "n_cells": len(e["cells"]), "n_tumour_cells": n, "cell_fraction": round(len(e["cells"]) / n, 4) if n else None,
                    "median_cn_delta": round(statistics.median(e["deltas"]), 2), "cells": e["cells"]})

    for (c, arm, st), e in sorted(calls.items(), key=lambda x: (int(x[0][0]), x[0][1], x[0][2])):
        add(c, arm, "arm", st, e)
    for (c, st), e in sorted(chrom_calls.items(), key=lambda x: (int(x[0][0]), x[0][1])):
        add(c, "whole", "chromosome", st, e)
    notes = "chromosomes X, Y, M and arms with < %d bins (acrocentric p) not assessed; normal cells excluded" % MIN_BINS
    return out, n, notes


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=DATA)
    ap.add_argument("--out", required=True, help="output root; writes <out>/<P>/arm_events.json")
    ap.add_argument("--patients", nargs="*")
    ap.add_argument("--cytobands", default=CYTOBANDS)
    a = ap.parse_args(argv)
    arms = read_arms(a.cytobands)
    pts = a.patients or sorted(p for p in os.listdir(a.data) if re.fullmatch(r"[A-Z]+[0-9]+", p))
    for p in pts:
        f = os.path.join(a.data, p, "cn_binned.json")
        if not os.path.isfile(f):
            print(p, "no cn_binned.json: skipped")
            continue
        cn = json.load(open(f))
        if cn.get("genome") not in (None, "hg38"):
            print(p, "genome", cn.get("genome"), "is not hg38: skipped")
            continue
        ev, n, notes = arm_events(cn, arms)
        doc = {"format": "gos-sc-arm-events/1", "patient": p, "genome": "hg38", "created": time.strftime("%Y-%m-%dT%H:%M:%S"), "source": "cn_binned.json 1 Mb level",
               "method": "arm median CN vs cell baseline (round of median arm median): gain >= +%.2f, loss <= %.2f; chromosome = both arms same call" % (GAIN, LOSS),
               "n_tumour_cells": n, "notes": notes, "events": ev}
        d = os.path.join(a.out, p)
        os.makedirs(d, exist_ok=True)
        dest = os.path.join(d, "arm_events.json")
        if os.path.exists(dest):
            bak = f"{dest}.bak_{time.strftime('%y%m%d')}"
            if not os.path.exists(bak):
                shutil.copy2(dest, bak)
        with open(dest + ".tmp", "w") as fh:
            json.dump(doc, fh, separators=(",", ":"))
        os.replace(dest + ".tmp", dest)
        print(p, "tumour cells", n, "events", len(ev), "->", dest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
