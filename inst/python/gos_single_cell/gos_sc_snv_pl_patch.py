#!/usr/bin/env python3
"""Add the PL confidence byte to an existing columnar SNV matrix (snv_matrix.cells.bin) WITHOUT re-running the upload.

Input: <data>/<P>/snv_matrix.{sites.json,cells.bin} and <pl>/<P>/pl_bytes.u8 (srctools inst/gos_single_cell/gos_sc_pl_dump.R:
row-major sites x cells, uint8: 1..250 = 1 + round(249 P(ALT)), 255 = flat/missing PL, 0 = not in the pool matrices).
Output (--out, never in place unless --out == --data, then a .bak_<tag> copy is kept): patched snv_matrix.cells.bin(.gz) and
snv_matrix.sites.json with matrix.pl_byte = true. Entries of the CSR matrix keep cell/ref/alt/gt; only the 8th byte changes.
Cell/site pairs without a CSR entry have no reads and no hard call and stay "no evidence" in gOS.

  python3 gos_sc_snv_pl_patch.py --patients MGH303[,report,...] --data DATA --pl PLDIR --out STAGE [--name snv_matrix] [--tag 261101]
"""
import argparse
import gzip
import hashlib
import json
import os
import shutil
import sys

import numpy as np


def patch(buf, n_sites, n_cells, pl, head):
    """buf: cells.bin bytes; pl: uint8 array (n_sites * n_cells). -> new bytes, n entries patched."""
    offsets = np.frombuffer(buf, dtype="<u4", count=n_sites + 1)
    n_entries = int(offsets[-1])
    arr = np.frombuffer(buf, dtype=np.dtype([("c", "<u2"), ("r", "<u2"), ("a", "<u2"), ("g", "i1"), ("p", "u1")]), count=n_entries, offset=head).copy()
    site_of = np.repeat(np.arange(n_sites), np.diff(offsets.astype(np.int64)))
    vals = pl[site_of * n_cells + arr["c"].astype(np.int64)]
    arr["p"] = vals
    return buf[:head] + arr.tobytes(), int((vals > 0).sum()), n_entries


def run(patient, data, pl_dir, out, name="snv_matrix", tag="261101", min_frac=0.995):
    """patient = report folder name. Returns a status string (ok / skipped: reason)."""
    src = os.path.join(data, patient)
    dst = os.path.join(out, patient)
    if os.path.islink(os.path.join(src, f"{name}.cells.bin")):
        return f"skipped: {patient} cells.bin is a symlink to {os.path.realpath(os.path.join(src, name + '.cells.bin'))} (patch that report)"
    tab = json.load(open(os.path.join(src, f"{name}.sites.json")))
    n_sites, n_cells = tab["n_sites"], len(tab["cells"])
    pl = np.fromfile(os.path.join(pl_dir, patient, "pl_bytes.u8"), dtype=np.uint8)
    if pl.size != n_sites * n_cells:
        sys.exit(f"{patient}: pl_bytes.u8 has {pl.size} values, expected {n_sites} x {n_cells}")
    buf = open(os.path.join(src, tab["matrix"]["file"]), "rb").read()
    new, n_pl, n_entries = patch(buf, n_sites, n_cells, pl, tab["matrix"]["header_bytes"])
    if n_entries and n_pl / n_entries < min_frac:
        return f"skipped: {patient} only {n_pl}/{n_entries} matrix entries have PL in the pool matrices (< {min_frac})"
    os.makedirs(dst, exist_ok=True)
    with open(os.path.join(dst, "source_cells_bin.md5"), "w") as fh:   # lets the deploy refuse to overwrite a live file that changed since
        fh.write(hashlib.md5(buf).hexdigest() + "\n")
    if os.path.abspath(dst) == os.path.abspath(src):
        for f in (f"{name}.cells.bin", f"{name}.cells.bin.gz", f"{name}.sites.json"):
            if os.path.exists(os.path.join(src, f)):
                shutil.copy2(os.path.join(src, f), os.path.join(src, f"{f}.bak_{tag}"))
    for ext, op in (("cells.bin", open), ("cells.bin.gz", gzip.open)):
        p = os.path.join(dst, f"{name}.{ext}")
        if os.path.islink(p):
            os.unlink(p)
        with op(p + ".tmp", "wb") as fh:
            fh.write(new)
        os.replace(p + ".tmp", p)
    tab["matrix"]["pl_byte"] = True
    tab["matrix"]["encoding"] = "csr-u16x3-i8-u8 (8-byte entries: cell, ref, alt, gt, pl byte; Uint32 site offsets, 8-byte aligned)"
    p = os.path.join(dst, f"{name}.sites.json")
    if os.path.islink(p):
        os.unlink(p)
    with open(p + ".tmp", "w") as fh:
        json.dump(tab, fh, separators=(",", ":"))
    os.replace(p + ".tmp", p)
    gz = os.path.join(dst, f"{name}.sites.json.gz")
    if os.path.exists(gz) or os.path.islink(gz):
        os.unlink(gz)
    with gzip.open(gz, "wb") as fh:
        fh.write(open(p, "rb").read())
    flat = int((np.frombuffer(new, dtype=np.uint8)[tab["matrix"]["header_bytes"] + 7::8] == 255).sum())
    return f"ok: {patient}: {n_entries} entries, {n_pl} with a PL byte ({flat} flat) -> {dst}"


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", required=True)
    ap.add_argument("--data", required=True)
    ap.add_argument("--pl", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--name", default="snv_matrix")
    ap.add_argument("--tag", default="261101")
    a = ap.parse_args()
    for P in a.patients.split(","):
        print(run(P, a.data, a.pl, a.out, a.name, a.tag))
