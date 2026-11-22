#!/usr/bin/env python3
"""Mask chrY as "not available" in the gOS single-cell files of female patients.

A female patient has no chrY: its chrY copy number is not a measurement, and low-pass graph fits can put arbitrary
copy numbers there (seen in the MGH303 recurrence cells, chrY CN 1-1000 with no chrY reads). For every report whose
patient is female in the cohort clinical metadata (_cohort/clinical_assays.json -> clinical_meta.patients.<P>.sex), this:
  - drops chrY intervals (and the connections that touch them) from each cell's complex.json / allelic.json, so the
    gOS heat map and tracks show chrY as missing instead of a copy number;
  - sets chrY bins of the report's cn_binned_*.u16 matrices (all levels) to the missing value;
  - records the mask: complex.json "masked_chromosomes" and <report>/sex_mask.json.
Male or unknown sex: nothing changes. Each overwritten file is first copied to <file><backup suffix> (once).
Idempotent: a file that already carries the mask is skipped.

usage: gos_sc_sex_mask.py --data-dir DATA --datafiles datafiles.json --reports R1,R2 [--backup-suffix .bak_261121] [--dry-run]
"""
import argparse
import gzip
import json
import os
import re
import shutil
import sys

REASON = "not available (female patient: no chrY)"
Y_NAMES = {"Y", "chrY"}


def patient_of(report):
    m = re.match(r"^([A-Za-z]+\d+)", report)
    return m.group(1) if m else report


def sex_of(clinical, patient):
    pts = (clinical or {}).get("clinical_meta", {}).get("patients", {})
    s = str((pts.get(patient) or {}).get("sex") or "").strip().lower()
    return "female" if s in ("f", "female") else "male" if s in ("m", "male") else "unknown"


def mask_graph(g):
    """Drop chrY intervals and the connections touching them. Returns (graph, n_intervals_dropped, n_connections_dropped)
    or (graph, None, None) when the graph already carries the mask."""
    if "Y" in (g.get("masked_chromosomes") or {}):
        return g, None, None
    ivs = g.get("intervals") or []
    drop = {abs(int(i["iid"])) for i in ivs if str(i.get("chromosome")) in Y_NAMES and i.get("iid") is not None}
    keep = [i for i in ivs if str(i.get("chromosome")) not in Y_NAMES]
    def touches(c):
        return any(c.get(k) is not None and abs(int(c[k])) in drop for k in ("source", "sink"))
    cons = g.get("connections") or []
    keep_c = [c for c in cons if not touches(c)]
    g = dict(g)
    g["intervals"] = keep
    if "connections" in g:
        g["connections"] = keep_c
    g["masked_chromosomes"] = {"Y": REASON}
    return g, len(ivs) - len(keep), len(cons) - len(keep_c)


def backup(path, suffix):
    b = path + suffix
    if not os.path.exists(b):
        shutil.copy2(path, b)


def rewrite_json(path, obj, suffix, dry):
    if dry:
        return
    for p in (path, path + ".gz"):
        if os.path.exists(p):
            backup(p, suffix)
    tmp = path + ".tmp_sexmask"
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=2)
    os.replace(tmp, path)
    if os.path.exists(path + ".gz"):
        with open(path, "rb") as fi, gzip.open(path + ".gz", "wb") as fo:
            shutil.copyfileobj(fi, fo)


def mask_cn_binned(rdir, suffix, dry):
    """Set chrY bins of every level's cn_binned_*.u16 matrices (cell-major uint16, cn_binned.json) to the missing value.
    Returns the number of chrY bins masked per cell, summed over levels."""
    import array
    meta_f = os.path.join(rdir, "cn_binned.json")
    if not os.path.exists(meta_f):
        return 0
    meta = json.load(open(meta_f))
    if "Y" in (meta.get("masked_chromosomes") or {}):
        return 0
    total = 0
    for name, lv in (meta.get("levels") or {}).items():
        n_bins = int(lv.get("n_bins") or 0)
        missing = int(lv.get("missing_value", meta.get("missing_value", 65535)))
        yr = [(int(c["first_bin"]), int(c["n_bins"])) for c in lv.get("chromosomes") or [] if str(c.get("chr")) in Y_NAMES]
        if not yr or not n_bins:
            continue
        first, n = yr[0]
        for key, fn in (lv.get("files") or {}).items():
            f = os.path.join(rdir, fn)
            if not os.path.exists(f):
                continue
            a = array.array("H")
            with open(f, "rb") as fh:
                a.frombytes(fh.read())
            if sys.byteorder != "little":
                a.byteswap()
            for r in range(len(a) // n_bins):
                b = r * n_bins + first
                a[b:b + n] = array.array("H", [missing] * n)
            if not dry:
                backup(f, suffix)
                if sys.byteorder != "little":
                    a.byteswap()
                with open(f + ".tmp_sexmask", "wb") as fh:
                    fh.write(a.tobytes())
                os.replace(f + ".tmp_sexmask", f)
        total += n
    meta["masked_chromosomes"] = {"Y": REASON}
    rewrite_json(meta_f, meta, suffix, dry)
    return total


def run(data_dir, datafiles, reports, clinical, suffix=".bak_sexmask", dry=False, log=print):
    cells_of = {}
    for e in datafiles:
        if e.get("entry_type") == "cell":
            cells_of.setdefault(e.get("patient_id"), []).append(e["pair"])
    out = {}
    for rep in reports:
        p = patient_of(rep)
        sex = sex_of(clinical, p)
        if sex != "female":
            log("%s: patient %s sex %s -> unchanged" % (rep, p, sex))
            out[rep] = {"sex": sex, "cells": 0}
            continue
        n_cells = n_iv = n_con = n_skip = 0
        for cell in cells_of.get(rep, []):
            for fname in ("complex.json", "allelic.json"):
                f = os.path.join(data_dir, cell, fname)
                if not os.path.exists(f):
                    continue
                g, a, b = mask_graph(json.load(open(f)))
                if a is None:
                    n_skip += 1
                    continue
                rewrite_json(f, g, suffix, dry)
                n_iv += a
                n_con += b
                if fname == "complex.json":
                    n_cells += 1
        yb = mask_cn_binned(os.path.join(data_dir, rep), suffix, dry)
        note = {"masked_chromosomes": {"Y": REASON}, "sex": sex, "sex_source": "_cohort/clinical_assays.json",
                "cells_masked": n_cells, "intervals_dropped": n_iv, "connections_dropped": n_con, "cn_binned_y_bins": yb}
        if not dry and os.path.isdir(os.path.join(data_dir, rep)):
            with open(os.path.join(data_dir, rep, "sex_mask.json"), "w") as fh:
                json.dump(note, fh, indent=2)
        log("%s: female; %d cells masked (%d chrY intervals, %d connections dropped; %d files already masked); cn_binned chrY bins %d"
            % (rep, n_cells, n_iv, n_con, n_skip, yb))
        out[rep] = dict(note, already=n_skip)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--datafiles", required=True)
    ap.add_argument("--clinical", help="default: <data-dir>/_cohort/clinical_assays.json")
    ap.add_argument("--reports", required=True, help="comma-separated report ids; 'female' = every registered report of a female patient")
    ap.add_argument("--backup-suffix", default=".bak_sexmask")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)
    clinical = json.load(open(a.clinical or os.path.join(a.data_dir, "_cohort", "clinical_assays.json")))
    datafiles = json.load(open(a.datafiles))
    reports = [r for r in a.reports.split(",") if r]
    if reports == ["female"]:
        reports = [e["pair"] for e in datafiles if e.get("entry_type") == "patient" and sex_of(clinical, patient_of(e["pair"])) == "female"]
    run(a.data_dir, datafiles, reports, clinical, a.backup_suffix, a.dry_run)


if __name__ == "__main__":
    main()
