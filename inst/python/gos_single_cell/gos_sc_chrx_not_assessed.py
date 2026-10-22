#!/usr/bin/env python3
"""Mark chrX copy-number calls of dryclean-artefact cells as "not assessed" (per-cell filtered.events.json), into a staging folder.

usage: gos_sc_chrx_not_assessed.py --data <gOS data> --out <staging data dir> --flags cells_chrx_flags.tsv [--patients A,B]
Writes <out>/<cell>/filtered.events.json only for flagged cells that have chrX HOMDEL rows (all other files are left to the caller) and
<out>/_chrx_not_assessed_summary.json (rows changed per patient). Mutation rows are never touched. Reads --data only.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chrx_flags as cf  # noqa: E402


def run(data, out, flags, patients=None):
    summ = {}
    for patient, cells in sorted(flags.items()):
        if patients and patient not in patients:
            continue
        s = {"flagged_cells": len(cells), "cells_with_event_file": 0, "cells_changed": 0, "rows_changed": 0, "genes": []}
        for cell in sorted(cells):
            f = os.path.join(data, cell, "filtered.events.json")
            if not os.path.exists(f):
                continue
            s["cells_with_event_file"] += 1
            with open(f) as fh:
                ev = json.load(fh)
            new, n = cf.not_assess_events(ev)
            if not n:
                continue
            s["cells_changed"] += 1
            s["rows_changed"] += n
            s["genes"] += sorted({e["gene"] for e, o in zip(new, ev) if e is not o and e.get("gene")})
            os.makedirs(os.path.join(out, cell), exist_ok=True)
            with open(os.path.join(out, cell, "filtered.events.json"), "w") as fh:
                json.dump(new, fh)
        s["genes"] = sorted(set(s["genes"]))
        summ[patient] = s
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "_chrx_not_assessed_summary.json"), "w") as fh:
        json.dump(summ, fh, indent=1)
    return summ


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--flags", required=True)
    ap.add_argument("--patients", default="")
    a = ap.parse_args(argv)
    if os.path.abspath(a.out) == os.path.abspath(a.data):
        sys.exit("--out must differ from --data (staging only)")
    summ = run(a.data, a.out, cf.load_flagged(a.flags), set(x for x in a.patients.split(",") if x))
    for p, s in summ.items():
        print(p, s["flagged_cells"], "flagged;", s["cells_changed"], "cells /", s["rows_changed"], "rows changed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
