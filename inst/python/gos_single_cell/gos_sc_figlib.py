"""Shared helpers of the paper-figure data exporters (gos_sc_cn_binned, gos_sc_walk_segment_cn, gos_sc_phylo_signal,
gos_sc_tree_site_fit): per-cell copy-number segments, gOS chromosome coordinates, Newick tree layout, patient cell lists.

The functions are ports of gOS_dev src/helpers/singleCell/{newick,matrix,cellFiles}.js so that exports reproduce what the
browser computes (see tests/fixtures/*golden*.json, produced by running the JS itself on a synthetic fixture).
"""
import json
import math
import os
import re

import numpy as np

DATA_ROOT = "/gpfs/commons/projects/imielinski_web/sclarke/gOS_GBM_sc"
# hg38 coordinate set of gOS public/settings.json (chromoBins: startPlace = cumulative length + 1)
CHROMS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
CHROM_LEN = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422,
             135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167,
             46709983, 50818468, 156040895, 57227415]
LEN = dict(zip(CHROMS, CHROM_LEN))
OFFSET = {}
_b = 0
for _c in CHROMS:
    OFFSET[_c] = _b
    _b += LEN[_c]
GENOME_LENGTH = _b


def jload(path):
    with open(path) as fh:
        return json.load(fh)


def norm_chr(c):
    c = str(c).strip()
    c = re.sub(r"^chr", "", c, flags=re.I)
    return c if c in LEN else None


def is_normal_clone(clone):
    return bool(re.search("normal", "" if clone is None else str(clone), re.I))


def patient_cells(datafiles, patient):
    """datafiles.json cell records of a patient: [{cell_id, clone_id, ...}] (gOS cellsForPatient)."""
    seen = set()
    out = []
    for r in datafiles:
        if r.get("entry_type") != "cell" or str(r.get("patient_id", "")).strip() != patient:
            continue
        cid = str(r.get("pair") or r.get("caseReportId") or "")
        if not cid or cid in seen:
            continue
        seen.add(cid)
        clone = r.get("clone_id")
        out.append(dict(r, cell_id=cid, clone_id=None if clone in (None, "") else str(clone)))
    out.sort(key=lambda c: natural_key(c["cell_id"]))
    return out


def natural_key(s):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)]


class SegRow:
    """Segments of one cell on gOS global coordinates; value lookup as gOS binAt (last start <= g, g <= end)."""

    def __init__(self, chroms, starts, ends, values):
        keep = [i for i, c in enumerate(chroms) if norm_chr(c)]
        g0 = np.array([OFFSET[norm_chr(chroms[i])] + 1 + float(starts[i]) for i in keep])
        g1 = np.array([OFFSET[norm_chr(chroms[i])] + 1 + float(ends[i]) for i in keep])
        v = np.array([float(values[i]) for i in keep])
        o = np.argsort(g0, kind="stable")
        self.g0, self.g1, self.v = g0[o], g1[o], v[o]

    def at(self, g):
        """values at global positions g (array); NaN where no segment covers."""
        g = np.asarray(g, dtype=np.float64)
        if not len(self.g0):
            return np.full(g.shape, np.nan)
        k = np.searchsorted(self.g0, g, side="right") - 1
        ok = k >= 0
        kk = np.where(ok, k, 0)
        ok &= g <= self.g1[kk]
        return np.where(ok, self.v[kk], np.nan)

    def mean_over(self, a, b):
        """length-weighted mean value over the global interval [a, b] (NaN when nothing overlaps)."""
        if not len(self.g0):
            return float("nan")
        ov = np.minimum(self.g1, b) - np.maximum(self.g0, a) + 1
        m = (ov > 0) & np.isfinite(self.v)
        return float(np.sum(ov[m] * self.v[m]) / np.sum(ov[m])) if m.any() else float("nan")


def total_cn_row(genome):
    """complex.json -> SegRow of total copy number (gOS cnRowFromGenome)."""
    iv = [d for d in (genome or {}).get("intervals", []) if (d.get("type") in (None, "interval")) and _finite(d.get("y")) and norm_chr(d.get("chromosome"))]
    return SegRow([d["chromosome"] for d in iv], [d["startPoint"] for d in iv], [d["endPoint"] for d in iv], [d["y"] for d in iv])


def allelic_rows(allelic):
    """allelic.json -> (SegRow major, SegRow minor) (gOS allelicRowFromAllelic: each locus has two y values, larger = major)."""
    loci = {}
    for d in (allelic or {}).get("intervals", []):
        if d.get("type") not in (None, "interval") or not norm_chr(d.get("chromosome")):
            continue
        key = (d["chromosome"], d["startPoint"], d["endPoint"])
        loc = loci.setdefault(key, {"values": [], "major": None, "minor": None})
        if d.get("majorCn") is not None or d.get("minorCn") is not None:
            loc["major"], loc["minor"] = _num(d.get("majorCn")), _num(d.get("minorCn"))
        else:
            loc["values"].append(_num(d.get("y")))
    keys, maj, mnr = [], [], []
    for key, loc in loci.items():
        if loc["major"] is None and len(loc["values"]) == 2 and all(math.isfinite(x) for x in loc["values"]):
            loc["major"], loc["minor"] = max(loc["values"]), min(loc["values"])
        keys.append(key)
        maj.append(loc["major"] if loc["major"] is not None and math.isfinite(loc["major"]) else float("nan"))
        mnr.append(loc["minor"] if loc["minor"] is not None and math.isfinite(loc["minor"]) else float("nan"))
    cols = list(zip(*keys)) if keys else ([], [], [])
    mk = lambda vals: SegRow(list(cols[0]), list(cols[1]), list(cols[2]), vals)
    return mk(maj), mk(mnr)


def _finite(x):
    try:
        return math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


# --------------------------------------------------------------------------------------------------------------
# Newick (gOS newick.js parseNewick / layoutTree)
# --------------------------------------------------------------------------------------------------------------
def parse_newick(text):
    s = text.strip()
    if not s:
        raise ValueError("Empty Newick string")
    i = 0
    n = len(s)

    def skip():
        nonlocal i
        while True:
            while i < n and s[i] in " \n\r\t":
                i += 1
            if i < n and s[i] == "[":
                close = s.find("]", i)
                if close < 0:
                    raise ValueError("Unterminated Newick comment")
                i = close + 1
            else:
                return

    def read_label():
        nonlocal i
        skip()
        if i < n and s[i] in "'\"":
            q = s[i]
            i += 1
            out = ""
            while i < n:
                if s[i] == q:
                    if i + 1 < n and s[i + 1] == q:
                        out += q
                        i += 2
                        continue
                    i += 1
                    return out
                out += s[i]
                i += 1
            raise ValueError("Unterminated quoted Newick label")
        st = i
        while i < n and s[i] not in "(),:;[":
            i += 1
        raw = s[st:i].strip()
        return raw if raw else None

    def read_length():
        nonlocal i
        skip()
        if i >= n or s[i] != ":":
            return None
        i += 1
        skip()
        st = i
        while i < n and re.match(r"[0-9eE+\-.]", s[i]):
            i += 1
        return float(s[st:i])

    def read_node():
        nonlocal i
        skip()
        node = {"name": None, "length": None, "children": []}
        if i < n and s[i] == "(":
            i += 1
            while True:
                node["children"].append(read_node())
                skip()
                if i < n and s[i] == ",":
                    i += 1
                    continue
                if i < n and s[i] == ")":
                    i += 1
                    break
                raise ValueError("Unexpected token in Newick at %d" % i)
        node["name"] = read_label()
        node["length"] = read_length()
        return node

    import sys
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 20000))
    root = read_node()
    skip()
    if i < n and s[i] == ";":
        i += 1
    skip()
    if i < n:
        raise ValueError("Unexpected trailing Newick content")
    return root


def layout_tree(root):
    """Leaves in DFS order with first/last leaf ranges and x = cumulative branch length (gOS layoutTree).
    Returns {nodes: [{x, firstLeaf, lastLeaf, isLeaf, children, name}], leaves: [names], x per leaf}."""
    has_lengths = False
    stack = [root]
    while stack:
        nd = stack.pop()
        if nd is not root and nd["length"] is not None and nd["length"] > 0:
            has_lengths = True
        stack.extend(nd["children"])
    heights = {}
    order = []
    stack = [root]
    while stack:
        nd = stack.pop()
        order.append(nd)
        stack.extend(nd["children"])
    for nd in reversed(order):
        kids = nd["children"]
        heights[id(nd)] = 0 if not kids else 1 + max(heights[id(c)] for c in kids)
    root_h = heights[id(root)]
    nodes = []
    leaves = []

    def rec(nd, px, is_root):
        nx = (px + (0.0 if is_root else max(0.0, nd["length"] or 0.0))) if has_lengths else root_h - heights[id(nd)]
        me = {"x": nx, "children": [], "name": nd["name"], "isLeaf": not nd["children"]}
        my_id = len(nodes)
        nodes.append(me)
        if not nd["children"]:
            me["firstLeaf"] = me["lastLeaf"] = len(leaves)
            leaves.append(nd["name"])
        else:
            for c in nd["children"]:
                me["children"].append(rec(c, nx, False))
            me["firstLeaf"] = nodes[me["children"][0]]["firstLeaf"]
            me["lastLeaf"] = nodes[me["children"][-1]]["lastLeaf"]
        return my_id

    rec(root, 0.0, True)
    return {"nodes": nodes, "leaves": leaves}


def prune_tree_leaves(layout, keep):
    """leaf subset helper: returns the leaf names in layout order that are in keep."""
    k = set(keep)
    return [l for l in layout["leaves"] if l in k]
