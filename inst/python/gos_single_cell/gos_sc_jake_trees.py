#!/usr/bin/env python3
"""gos_sc_jake_trees: a collaborator's TBE-mapped phylogenies + per-cell SNV table as extra alternatives of the gOS "Compare Trees" tab.

Input (read only), per patient P in SRC:
  P.mapped.TBE.nwk, P.mapped.tumour.TBE.nwk   trees whose internal node labels are TBE support (0..1)
  P.per_cell_long_table.tsv.gz                one row per (mutation, alt-carrying cell): class, node, clade_size, gene, func, ...
  P.counts.alt.tsv.gz, P.counts.depth.tsv.gz  full per-cell alt / total read matrices (rows "chr1:822354_C>T")
A tip / cell name is matched to the gOS cell id by --alias (TSV alias, gos_id; case-insensitive, "-" == "_", "MR_5" == "MR5";
names like 439538_MGH305_MR_2_pl3_7h_<bc>_<bc>_merged are reduced to MGH305_MR_2_pl3_7h first). Unmatched tips are dropped.

Output under OUT/P (real files; the live tree.nwk and snv_matrix.json are never touched):
  trees/<name>.nwk            tips renamed to gOS ids (name = jake_tbe_261026, jake_tbe_tumour_261026)
  trees/index.json            the live index (copied from --live) + the new entries
  trees/support.json          the live support.json + entries {splits: [{key, size, fbp: null, tbe}]} (split keys as gos_sc_tree_support.py)
  snv_matrix.jake_261026.json sparse matrix of the long table's mutations, ref/alt reads of every covered cell from the counts matrices
  snv_sets.json               the live one (or just "full") + the new set; n_variants left null so the default set stays the full matrix
  jake_trees_report.json      matched / unmatched cells, sites, checks

python3 gos_sc_jake_trees.py --patient BWH70 --src DIR --alias BWH70.alias.tsv --live DATA_DIR --out OUT_DATA_DIR [--no-matrix]
"""
import argparse
import csv
import gzip
import json
import os
import re
import shutil
import sys

TREES = [
    ("jake_tbe_261026", "mapped.TBE", "Jake 26 Oct TBE (all cells)"),
    ("jake_tbe_tumour_261026", "mapped.tumour.TBE", "Jake 26 Oct TBE (tumour cells only)"),
]
SET_NAME = "jake_261026"
SET_FILE = "snv_matrix.jake_261026.json"
CATEGORY = {"truncal": "truncal", "subclonal": "subclonal", "singleton_tumour": "private", "germline": "germline"}


def key(x):
    """Cell-name key: lower case, '-' == '_', 'mr_5' == 'mr5'."""
    return re.sub(r"mr_([0-9])", r"mr\1", str(x).lower().replace("-", "_"))


def reduce_name(n):
    """439538_MGH305_MR_2_pl3_7h_CAATAGCATCGAGAT_CAATAGCATCGAGAT_merged -> MGH305_MR_2_pl3_7h"""
    m = re.match(r"^\d+_(.+?)_[ACGT]{10,}_[ACGT]{10,}_merged$", n)
    return m.group(1) if m else n


def read_alias(path):
    out = {}
    with open(path) as f:
        r = csv.reader(f, delimiter="\t")
        next(r)
        for a, g in r:
            out.setdefault(key(a), g)
    return out


# ---------------------------------------------------------------- newick
class Node:
    __slots__ = ("name", "label", "length", "children")

    def __init__(self):
        self.name, self.label, self.length, self.children = None, None, None, []


def parse_newick(text):
    s = text.strip().rstrip(";")
    pos = 0

    def token():
        nonlocal pos
        m = re.match(r"[^(),:;]*", s[pos:])
        pos += m.end()
        return m.group(0)

    def node():
        nonlocal pos
        n = Node()
        if s[pos] == "(":
            pos += 1
            while True:
                n.children.append(node())
                if s[pos] == ",":
                    pos += 1
                    continue
                break
            if s[pos] != ")":
                raise ValueError("unbalanced newick")
            pos += 1
            n.label = token().strip() or None
        else:
            n.name = token().strip().strip("'")
        if pos < len(s) and s[pos] == ":":
            pos += 1
            n.length = token().strip()
        return n

    root = node()
    if pos != len(s):
        raise ValueError("trailing characters in newick")
    return root


def leaves(n):
    if not n.children:
        return [n]
    out = []
    for c in n.children:
        out.extend(leaves(c))
    return out


def write_newick(n):
    def rec(x):
        if x.children:
            t = "(" + ",".join(rec(c) for c in x.children) + ")" + (x.label or "")
        else:
            t = x.name
        return t + (":" + x.length if x.length is not None else "")

    return rec(n) + ";"


def prune(n, keep):
    """Keep only the leaves whose name is in keep; unary nodes are merged into their child (branch lengths added; the child's label wins).
    Returns the new node or None."""
    if not n.children:
        return n if n.name in keep else None
    kids = [k for k in (prune(c, keep) for c in n.children) if k is not None]
    if not kids:
        return None
    if len(kids) == 1:
        k = kids[0]
        if n.length is not None:
            k.length = repr(float(n.length) + float(k.length)) if k.length is not None else n.length
        return k
    n.children = kids
    return n


def rename_tree(root, alias):
    """Tips renamed to gOS ids; unmatched and duplicate tips dropped. Returns (root | None, n_tips, matched names, unmatched names)."""
    ls = leaves(root)
    seen, keep, unmatched = set(), set(), []
    for lf in ls:
        g = alias.get(key(reduce_name(lf.name)))
        if g is None or g in seen:
            unmatched.append(lf.name)
            lf.name = "__drop__" + lf.name
            continue
        seen.add(g)
        lf.name = g
        keep.add(g)
    new = prune(root, keep) if keep else None
    return new, len(ls), sorted(keep), unmatched


# ---------------------------------------------------------------- support (same keys as gos_sc_tree_support.py / helpers/singleCell/treeSupport.js)
def fnv1a32(text):
    h = 0x811C9DC5
    for b in text.encode("utf-8"):
        h ^= b
        h = (h * 0x01000193) & 0xFFFFFFFF
    return format(h, "08x")


def tip_key(s):
    return str(s).lower().replace("-", "_")


def split_key(names):
    nm = sorted(names)
    return "%d:%s" % (len(nm), fnv1a32("|".join(nm)))


def support_entry(root):
    """{n_tree_leaves, n_leaves_shared, excluded_leaves, ref_leaf, splits} with the internal labels as tbe; fbp is not available (null)."""
    shared = sorted(tip_key(l.name) for l in leaves(root))
    n = len(shared)
    best = {}

    def rec(x, top):
        if not x.children:
            return [tip_key(x.name)]
        acc = []
        for c in x.children:
            acc.extend(rec(c, False))
        if not top and x.label not in (None, ""):
            try:
                v = float(x.label)
            except ValueError:
                v = None
            if v is not None:
                inside = set(acc)
                side = [t for t in shared if t in inside]
                if 2 <= len(side) <= n - 2:
                    if shared[0] in inside:
                        side = [t for t in shared if t not in inside]
                    k = split_key(side)
                    if k not in best or v > best[k]["tbe"]:
                        best[k] = {"key": k, "size": len(side), "fbp": None, "tbe": v}
        return acc

    rec(root, True)
    return {"n_tree_leaves": n, "n_leaves_shared": n, "excluded_leaves": [], "ref_leaf": shared[0], "splits": list(best.values())}


# ---------------------------------------------------------------- SNV matrix
def variant_id(chrom, pos, ref, alt):
    return "%s_%s_%s_%s" % (chrom, pos, ref, alt)


def chrom_order(c):
    m = re.match(r"chr(\d+)$", c)
    if m:
        return (0, int(m.group(1)))
    return (1, {"chrX": 23, "chrY": 24, "chrM": 25}.get(c, 99))


def read_long_sites(path, alias):
    """Distinct mutations of the long table: id -> record; also the cells it names."""
    sites, cells = {}, set()
    with gzip.open(path, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            vid = variant_id(r["chrom"], r["pos"], r["ref"], r["alt"])
            cells.add(r["cell"])
            if vid not in sites:
                func = r.get("exonic_func", ".")
                sites[vid] = {
                    "id": vid, "category": CATEGORY.get(r["class"], r["class"]), "cellphy_input": None,
                    "clade_cells": int(r["clade_size"]) if r["clade_size"].isdigit() else None,
                    "gene": r.get("gene") if r.get("gene") not in (None, ".") else None,
                    "consequence": (func if func not in (".", "") else r.get("func")) or None,
                    "protein": r["aachange"] if r.get("aachange") not in (None, ".", "") else "",
                    "alt_cells": 0, "_k": (chrom_order(r["chrom"]), int(r["pos"])),
                }
            sites[vid]["alt_cells"] += 1
    return sites, cells


def read_counts(path_alt, path_depth, wanted, alias):
    """cells: {gos_id: [[idx, ref, alt, None], ...]} for the sites in wanted (id -> index) ; returns (cells, header cells, unmatched cells, missing sites)."""
    out, dropped = {}, []
    with gzip.open(path_alt, "rt") as fa, gzip.open(path_depth, "rt") as fd:
        ha, hd = fa.readline().rstrip("\n").split("\t"), fd.readline().rstrip("\n").split("\t")
        if ha != hd:
            raise ValueError("alt and depth matrices have different columns")
        col = []
        for c in ha[1:]:
            g = alias.get(key(reduce_name(c)))
            col.append(g)
            if g is None:
                dropped.append(c)
            elif g not in out:
                out[g] = []
        seen = set()
        for la, ld in zip(fa, fd):
            a, d = la.rstrip("\n").split("\t"), ld.rstrip("\n").split("\t")
            if a[0] != d[0]:
                raise ValueError("alt and depth rows are out of step at " + a[0])
            m = re.match(r"^(chr[^:]+):(\d+)_([ACGTN]+)>([ACGTN]+)$", a[0])
            if not m:
                continue
            vid = variant_id(*m.groups())
            idx = wanted.get(vid)
            if idx is None:
                continue
            seen.add(vid)
            for g, x, y in zip(col, a[1:], d[1:]):
                if g is None:
                    continue
                alt, dp = int(float(x)), int(float(y))
                if dp > 0:
                    out[g].append([idx, max(dp - alt, 0), alt, None])
    return out, dropped, [v for v in wanted if v not in seen]


def build_matrix(src, P, alias, gos_cells):
    sites, long_cells = read_long_sites(os.path.join(src, P + ".per_cell_long_table.tsv.gz"), alias)
    ordered = sorted(sites.values(), key=lambda s: s["_k"])
    wanted = {s["id"]: i for i, s in enumerate(ordered)}
    cells, dropped, missing = read_counts(os.path.join(src, P + ".counts.alt.tsv.gz"), os.path.join(src, P + ".counts.depth.tsv.gz"), wanted, alias)
    variants = [{k: v for k, v in s.items() if k != "_k"} for s in ordered]
    keep_cells = {g: v for g, v in cells.items() if g in gos_cells}
    for v in variants:
        v["cellphy_input"] = None
    return {"schemaVersion": 1, "format": "compact", "variants": variants, "cells": keep_cells}, {
        "n_sites": len(variants), "n_sites_without_counts": len(missing), "counts_cells_unmatched": sorted(dropped),
        "long_table_cells": len(long_cells), "long_table_cells_unmatched": sorted(c for c in long_cells if alias.get(key(reduce_name(c))) is None),
    }


# ---------------------------------------------------------------- driver
def load_json(path, default=None):
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return default


def write_json(path, obj):
    if os.path.islink(path):
        raise RuntimeError("refusing to write through a symlink: " + path)
    with open(path, "w") as f:
        json.dump(obj, f, separators=(",", ":"))


def upsert(lst, entry):
    return [e for e in lst if e.get("name") != entry["name"]] + [entry]


def real_dir(path):
    p = path
    while p and p != os.sep:
        if os.path.islink(p):
            raise RuntimeError("staged path is (inside) a symlink: " + p)
        p = os.path.dirname(p)
    os.makedirs(path, exist_ok=True)


def run(P, src, alias_path, live, out, matrix=True):
    alias = read_alias(alias_path)
    gos_tips = leaves(parse_newick(open(os.path.join(live, P, "tree.nwk")).read()))
    gos_cells = {l.name for l in gos_tips}
    tdir = os.path.join(out, P, "trees")
    real_dir(tdir)
    live_t = os.path.join(live, P, "trees")
    for fn in ("index.json", "support.json"):  # live files are copied (never linked) and then extended
        if os.path.exists(os.path.join(live_t, fn)) and not os.path.exists(os.path.join(tdir, fn)):
            shutil.copyfile(os.path.join(live_t, fn), os.path.join(tdir, fn))
    for fn in os.listdir(live_t) if os.path.isdir(live_t) else []:
        if fn.endswith(".nwk") and not os.path.exists(os.path.join(tdir, fn)):
            shutil.copyfile(os.path.join(live_t, fn), os.path.join(tdir, fn))
    index = load_json(os.path.join(tdir, "index.json"), None) or [{"name": "gos", "label": "Pipeline (current)", "source_path": os.path.join(live, P, "tree.nwk"), "n_tips": len(gos_cells), "n_matched": len(gos_cells)}]
    support = load_json(os.path.join(tdir, "support.json"), None) or {"format": "gos-tree-support-v1", "available": True, "trees": {}}
    support.setdefault("trees", {})
    report = {"patient": P, "trees": {}}
    for name, suffix, label in TREES:
        path = os.path.join(src, "%s.%s.nwk" % (P, suffix))
        root, n_tips, matched, unmatched = rename_tree(parse_newick(open(path).read()), alias)
        if root is None or len(matched) < 4:
            report["trees"][name] = {"error": "fewer than 4 tips match gOS cells"}
            continue
        with open(os.path.join(tdir, name + ".nwk"), "w") as f:
            f.write(write_newick(root) + "\n")
        index = upsert(index, {"name": name, "label": label, "source_path": path, "n_tips": n_tips, "n_matched": len(matched),
                               "note": "TBE support only (no Felsenstein bootstrap)"})
        entry = support_entry(root)
        support["trees"][name] = entry
        report["trees"][name] = {"tips": n_tips, "matched": len(matched), "unmatched": unmatched, "gos_cells_absent": len(gos_cells - set(matched)),
                                 "internal_with_tbe": len(entry["splits"])}
    write_json(os.path.join(tdir, "index.json"), index)
    write_json(os.path.join(tdir, "support.json"), support)
    if matrix:
        pdir = os.path.join(out, P)
        mat, info = build_matrix(src, P, alias, gos_cells)
        write_json(os.path.join(pdir, SET_FILE), mat)
        sets = load_json(os.path.join(live, P, "snv_sets.json"), None) or [{"name": "full", "label": "Full matrix (current)", "file": "snv_matrix.json", "n_variants": None}]
        sets = upsert(sets, {"name": SET_NAME, "label": "Jake 26 Oct mapped sites (%d sites, slow to load)" % info["n_sites"], "file": SET_FILE, "n_variants": None})
        write_json(os.path.join(pdir, "snv_sets.json"), sets)
        info["cells_in_matrix"] = len(mat["cells"])
        info["gos_cells_without_data"] = len(gos_cells - set(mat["cells"]))
        report["matrix"] = info
    write_json(os.path.join(out, P, "jake_trees_report.json"), report)
    return report


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--patient", required=True)
    ap.add_argument("--src", required=True, help="dir with <P>.mapped.TBE.nwk ... (read only)")
    ap.add_argument("--alias", required=True)
    ap.add_argument("--live", required=True, help="live gOS data dir (read only): <P>/tree.nwk, trees/")
    ap.add_argument("--out", required=True, help="staged data dir (real files, no symlinks)")
    ap.add_argument("--no-matrix", action="store_true")
    a = ap.parse_args(argv)
    rep = run(a.patient, os.path.join(a.src, a.patient), a.alias, a.live, a.out, not a.no_matrix)
    json.dump(rep, sys.stdout, indent=1)
    print()


if __name__ == "__main__":
    main()
