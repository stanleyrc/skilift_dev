#!/usr/bin/env python3
"""Mark complex SV events that the gGnome classifier reported twice under different names.

The same junction set can be called both as a templated insertion chain (TIC) and as a reciprocal deletion-duplication
(rDelDup / rDup): identical carrier cells, the same calling cell, footprints that differ by a base. gOS shows one entry per
twin group with both names. Rows of data/<P>/filtered.events.json with type "Complex SV" get
  twin_group  "<P>:<n>"             on every member of a group (same cell_ids, footprint ends within MAX_SHIFT bp)
  twin_names  ["TIC (16,19)", "rDelDup (16,19)"]
Rows that are not twins get nothing. Backs up filtered.events.json as .bak_<yymmdd>_twins.

usage: gos_sc_complex_twins.py [--data <gOS data dir>] [--patients P ...] [--dry-run]
"""
import argparse
import json
import os
import re
import shutil
import time

DATA = "/gpfs/commons/projects/imielinski_web/sclarke/gOS_GBM_sc/data"
MAX_SHIFT = 1000


def footprint(loc):
    out = []
    for part in str(loc or "").split(","):
        m = re.fullmatch(r"\s*(?:chr)?([0-9XYM]+):(\d+)-(\d+)\s*", part)
        if m:
            out.append((m.group(1), int(m.group(2)), int(m.group(3))))
    return sorted(out)


def same_footprint(a, b, max_shift=MAX_SHIFT):
    fa, fb = footprint(a), footprint(b)
    return bool(fa) and len(fa) == len(fb) and all(x[0] == y[0] and abs(x[1] - y[1]) <= max_shift and abs(x[2] - y[2]) <= max_shift for x, y in zip(fa, fb))


def cells_of(e):
    return frozenset(c for c in str(e.get("cell_ids") or "").split(",") if c)


def twin_groups(events, patient):
    """[[row index, ...], ...]: complex rows with identical carrier cells and the same footprint, under different names."""
    rows = [i for i, e in enumerate(events) if str(e.get("type")) == "Complex SV" and cells_of(e)]
    used, groups = set(), []
    for i in rows:
        if i in used:
            continue
        g = [i]
        for j in rows:
            if j != i and j not in used and cells_of(events[j]) == cells_of(events[i]) and same_footprint(events[i].get("Genome_Location"), events[j].get("Genome_Location")) \
                    and str(events[j].get("Variant")) != str(events[i].get("Variant")):
                g.append(j)
        if len(g) > 1:
            used.update(g)
            groups.append(g)
    return groups


def annotate(events, patient):
    for e in events:
        e.pop("twin_group", None)
        e.pop("twin_names", None)
    groups = twin_groups(events, patient)
    for n, g in enumerate(groups, 1):
        names = [str(events[i].get("Variant") or events[i].get("vartype")) for i in g]
        for i in g:
            events[i]["twin_group"] = f"{patient}:{n}"
            events[i]["twin_names"] = names
    return groups


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=DATA)
    ap.add_argument("--patients", nargs="*")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)
    pts = a.patients or sorted(p for p in os.listdir(a.data) if re.fullmatch(r"[A-Z]+[0-9]+", p))
    for p in pts:
        f = os.path.join(a.data, p, "filtered.events.json")
        if not os.path.isfile(f):
            continue
        ev = json.load(open(f))
        groups = annotate(ev, p)
        for g in groups:
            print(f"  {p}: {' = '.join(str(ev[i].get('Variant')) for i in g)} ({ev[g[0]].get('n_cells')} cells)")
        print(p, "twin groups", len(groups))
        if a.dry_run:
            continue
        bak = f"{f}.bak_{time.strftime('%y%m%d')}_twins"
        if not os.path.exists(bak):
            shutil.copy2(f, bak)
        with open(f + ".tmp", "w") as fh:
            json.dump(ev, fh)
        os.replace(f + ".tmp", f)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
