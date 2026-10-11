#!/usr/bin/env python3
"""Cohort evolution exports for the gOS single-cell cohort views (#19 aligned trees, #21 associations).

Per patient (data/<P>/tree.nwk pruned to tumour cells; per-cell CN from data/<cell>/complex.json;
driver SNVs from filtered.events.json; ecDNA from walks.json; RNA state / phase from datafiles.json):

1. Events. Arm and whole-chromosome gains / losses (length-weighted median CN of the arm vs the cell's
   genome-wide baseline, +-0.75), focal amplification (CN >= max(baseline + 3, 2 x baseline)) and
   homozygous deletion (CN 0) at a fixed driver-gene list, driver SNVs (tier <= 2), and ecDNA (curated
   circular walks, >= 4 copies in the cell). Carriers = tumour cells with the call.
2. Placement. Each event goes on the clade that best matches its carriers (F1 of carriers vs clade
   leaves; the same score as gOS cladeFit). Placed when F1 >= 0.6 and >= max(3, 2% of tumour cells)
   carry it. Truncal when that clade is the tumour root or >= 85% of tumour cells carry it.
   Relative evolutionary time of a clade = 1 - log(size) / log(N tumour cells): 0 at the trunk,
   1 at single cells, monotone along every root-to-tip path, comparable across patients.
3. Ordering. For two placed events in one patient: same clade = same branch; one clade an ancestor of
   the other = that event came first; otherwise parallel lineages. Across patients: counts per event
   pair and a two-sided sign test on A-first vs B-first.
4. Associations. Per patient: lineage x RNA state (clone labels where >= 2 clones have >= 5 RNA cells,
   else the major tree clades), permutation chi-square and Cramer's V; per-state phylogenetic signal
   (share of state-X cells whose nearest tree relative is also X, vs 2000 permutations; also permuted
   only within lineages); event-anchored state associations (carriers vs non-carriers of each subclonal
   event) pooled across patients with Mantel-Haenszel; meta-program signal from
   analysis/story/out/heritability.json. Concordance: Stouffer Z, k of n significant, Cochran Q of the
   per-patient z.

Writes data/_cohort/aligned_trees.json and data/_cohort/associations.json.

usage: gos_sc_cohort_evolution.py [--data DIR] [--out DIR] [--patients P ...]
"""
import argparse
import datetime
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict
from itertools import combinations

import numpy as np
from scipy import stats

ROOT = "/gpfs/commons/projects/imielinski_web/sclarke/gOS_GBM_sc"
CYTOBANDS = "/gpfs/commons/projects/imielinski_web/sclarke/gOS_dev/cytobands/hg38.tsv"
HERITABILITY = "/gpfs/commons/groups/imielinski_lab/projects/single_cell_GBM/analysis/story/out/heritability.json"
MODULES = "/gpfs/commons/groups/imielinski_lab/projects/single_cell_GBM/analysis/story/out/modules.json"
PATIENTS = ["BWH69", "BWH70", "BWH75", "MGH285", "MGH301", "MGH302", "MGH303", "MGH304", "MGH305"]

# hg38 gene spans (GENCODE) for the focal calls; a config list, not GBM-specific code
FOCAL_GENES = {
    "EGFR": ("7", 55019017, 55211628, "amp"),
    "PDGFRA": ("4", 54229280, 54298245, "amp"),
    "MET": ("7", 116672196, 116798377, "amp"),
    "CDK4": ("12", 57747727, 57752310, "amp"),
    "CDK6": ("7", 92604921, 92836594, "amp"),
    "MDM2": ("12", 68808177, 68850686, "amp"),
    "MDM4": ("1", 204516379, 204558120, "amp"),
    "MYC": ("8", 127735434, 127742951, "amp"),
    "MYCN": ("2", 15940550, 15947004, "amp"),
    "CDKN2A": ("9", 21967752, 21995301, "del"),
    "PTEN": ("10", 87863438, 87971930, "del"),
    "RB1": ("13", 48303744, 48481890, "del"),
    "NF1": ("17", 31094927, 31382116, "del"),
    "TP53": ("17", 7661779, 7687538, "del"),
}
ARM_DELTA = 0.75
MIN_F1 = 0.6
TRUNCAL_FRAC = 0.85
ROOT_FRAC = 0.7
NPERM = 2000
MINUS = "\u2212"
STATES = ["MES", "AC", "NPC", "OPC"]
rng = np.random.default_rng(19)


def log(*a):
    print(*a, file=sys.stderr, flush=True)


# ---------------------------------------------------------------- tree

def parse_newick(text):
    """-> nodes [{children, name, length}], root id. Internal labels (bootstrap) are ignored."""
    s = text.strip().rstrip(";").strip()
    nodes = []

    def new():
        nodes.append({"children": [], "name": None, "length": 0.0, "parent": -1})
        return len(nodes) - 1

    root = new()
    stack, cur, i = [], root, 0
    label = re.compile(r"\s*('[^']*'|[^(),:;\s]+)?\s*(?::\s*([0-9.eE+-]+))?")
    started = False
    while i < len(s):
        ch = s[i]
        if ch == "(":
            if not started:
                started = True
                stack.append(root)
                child = new()
                nodes[child]["parent"] = root
                nodes[root]["children"].append(child)
                cur = child
            else:
                stack.append(cur)
                child = new()
                nodes[child]["parent"] = cur
                nodes[cur]["children"].append(child)
                cur = child
            i += 1
        elif ch == ",":
            parent = stack[-1]
            child = new()
            nodes[child]["parent"] = parent
            nodes[parent]["children"].append(child)
            cur = child
            i += 1
        elif ch == ")":
            cur = stack.pop()
            i += 1
            m = label.match(s, i)
            if m.group(2):
                nodes[cur]["length"] = float(m.group(2))
            i = max(m.end(), i)
        else:
            m = label.match(s, i)
            if m.group(1):
                nodes[cur]["name"] = m.group(1).strip("'")
            if m.group(2):
                nodes[cur]["length"] = float(m.group(2))
            i = max(m.end(), i + 1)
    return nodes, root


def prune(nodes, root, keep):
    """Pruned tree with only `keep` leaves, unary nodes collapsed. -> (children dict, root, leaf name of id)."""
    def rec(v):
        n = nodes[v]
        if not n["children"]:
            return v if n["name"] in keep else None
        kids = [k for k in (rec(c) for c in n["children"]) if k is not None]
        if not kids:
            return None
        if len(kids) == 1:
            return kids[0]
        n["kept"] = kids
        return v

    sys.setrecursionlimit(100000)
    r = rec(root)
    return r


class Tree:
    """Pruned tumour tree in leaf order: nodes with parent, firstLeaf, lastLeaf, size, depth."""

    def __init__(self, nodes, root):
        self.parent, self.first, self.last, self.depth, self.children = [], [], [], [], []
        self.leaves = []
        order = []

        def rec(v, p, d):
            idx = len(order)
            order.append(v)
            self.parent.append(p)
            self.depth.append(d)
            self.children.append([])
            self.first.append(len(self.leaves))
            self.last.append(-1)
            if p >= 0:
                self.children[p].append(idx)
            kids = nodes[v].get("kept") if nodes[v]["children"] else None
            if not kids:
                self.leaves.append(nodes[v]["name"])
            else:
                for k in kids:
                    rec(k, idx, d + 1)
            self.last[idx] = len(self.leaves) - 1
            return idx

        rec(root, -1, 0)
        self.n = len(self.leaves)
        self.size = [b - a + 1 for a, b in zip(self.first, self.last)]
        self.leaf_node = {}
        for i, (a, b) in enumerate(zip(self.first, self.last)):
            if a == b and not self.children[i]:
                self.leaf_node[a] = i

    def time(self, v):
        if self.n <= 1:
            return 0.0
        return 1.0 - math.log(self.size[v]) / math.log(self.n)

    def is_ancestor(self, a, b):
        """a is a strict ancestor of b."""
        return a != b and self.first[a] <= self.first[b] and self.last[b] <= self.last[a]

    def best_clade(self, rows):
        """Best F1 clade for carrier leaf rows. -> (node, f1, inClade)."""
        if not rows:
            return None, 0.0, 0
        P = np.zeros(self.n + 1, dtype=np.int64)
        flag = np.zeros(self.n, dtype=np.int64)
        flag[list(rows)] = 1
        P[1:] = np.cumsum(flag)
        first = np.array(self.first)
        last = np.array(self.last)
        inc = P[last + 1] - P[first]
        size = last - first + 1
        f1 = 2 * inc / (size + len(rows))
        v = int(np.argmax(f1))
        return v, float(f1[v]), int(inc[v])

    def leaf_distance_matrix(self):
        """Edges between leaves, with patristic-free tie-break by leaf order (nearest relative = sister first)."""
        n = self.n
        D = np.zeros((n, n), dtype=np.float32)
        leaf_depth = np.zeros(n)
        for r, v in self.leaf_node.items():
            leaf_depth[r] = self.depth[v]
        for v, kids in enumerate(self.children):
            if len(kids) < 2:
                continue
            for i in range(len(kids)):
                for j in range(i + 1, len(kids)):
                    a = np.arange(self.first[kids[i]], self.last[kids[i]] + 1)
                    b = np.arange(self.first[kids[j]], self.last[kids[j]] + 1)
                    block = leaf_depth[a][:, None] + leaf_depth[b][None, :] - 2 * self.depth[v]
                    D[np.ix_(a, b)] = block
                    D[np.ix_(b, a)] = block.T
        return D


# ---------------------------------------------------------------- copy number

def load_arms():
    by = defaultdict(lambda: {"acen": [], "end": 0})
    for line in open(CYTOBANDS):
        f = line.rstrip("\n").split("\t")
        if f[0] == "chromosome":
            continue
        c = f[0].replace("chr", "")
        if not re.fullmatch(r"[0-9]{1,2}|X", c):
            continue
        by[c]["end"] = max(by[c]["end"], int(f[2]))
        if f[4] == "acen":
            by[c]["acen"] += [int(f[1]), int(f[2])]
    arms = []
    for c, e in by.items():
        c0, c1 = min(e["acen"]), max(e["acen"])
        if c not in ("13", "14", "15", "21", "22"):
            arms.append((f"{c}p", c, 0, c0))
        arms.append((f"{c}q", c, c1, e["end"]))
    return arms


def wmedian(vals, wts):
    if not len(vals):
        return float("nan")
    o = np.argsort(vals)
    v, w = np.asarray(vals)[o], np.asarray(wts)[o]
    c = np.cumsum(w)
    return float(v[np.searchsorted(c, c[-1] / 2)])


def cell_cn_calls(path, arms):
    """-> dict event -> True for one cell, or None when the profile is missing."""
    try:
        g = json.load(open(path))
    except Exception:
        return None
    iv = [(str(x["chromosome"]).replace("chr", ""), int(x["startPoint"]), int(x["endPoint"]), x.get("y")) for x in g.get("intervals", [])]
    iv = [x for x in iv if x[3] is not None]
    if not iv:
        return None
    auto = [x for x in iv if re.fullmatch(r"[0-9]{1,2}", x[0])]
    base = wmedian([x[3] for x in auto], [x[2] - x[1] + 1 for x in auto])
    if not np.isfinite(base) or base <= 0:
        return None
    bych = defaultdict(list)
    for x in iv:
        bych[x[0]].append(x)
    calls = {}
    armstate = {}
    for name, c, s, e in arms:
        vals, wts = [], []
        for _, a, b, y in bych.get(c, []):
            lo, hi = max(a, s), min(b, e)
            if hi > lo:
                vals.append(y)
                wts.append(hi - lo)
        if sum(wts) < 0.3 * (e - s):
            continue
        m = wmedian(vals, wts)
        st = 1 if m - base >= ARM_DELTA else -1 if m - base <= -ARM_DELTA else 0
        armstate[name] = st
        if st:
            calls[f"{'+' if st > 0 else '-'}{name}"] = True
    # whole chromosome: both arms (or the q arm of an acrocentric) share the call
    chroms = defaultdict(list)
    for name, c, _, _ in arms:
        chroms[c].append(name)
    for c, names in chroms.items():
        sts = [armstate.get(n) for n in names]
        if any(s is None for s in sts):
            continue
        for sign in (1, -1):
            if all(s == sign for s in sts):
                calls[f"{'+' if sign > 0 else '-'}{c}"] = True
    for gene, (c, s, e, kind) in FOCAL_GENES.items():
        ys = [y for _, a, b, y in bych.get(c, []) if b >= s and a <= e]
        if not ys:
            continue
        if kind == "amp" and max(ys) >= max(base + 3, 2 * base):
            calls[f"{gene} amp"] = True
        if kind == "del" and min(ys) == 0 and base >= 1:
            calls[f"{gene} del"] = True
    return calls


# ---------------------------------------------------------------- stats

def sign_test(a, b):
    n = a + b
    if n == 0:
        return None
    return float(min(1.0, stats.binomtest(max(a, b), n, 0.5).pvalue)) if hasattr(stats, "binomtest") else float(stats.binom_test(max(a, b), n, 0.5))


def bh(ps):
    ps = np.asarray(ps, dtype=float)
    out = np.full(len(ps), np.nan)
    ok = np.isfinite(ps)
    if not ok.any():
        return out
    p = ps[ok]
    o = np.argsort(p)
    q = p[o] * len(p) / (np.arange(len(p)) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    r = np.empty(len(p))
    r[o] = np.minimum(q, 1)
    out[ok] = r
    return out


def chi2_stat(tab):
    tab = np.asarray(tab, dtype=float)
    tab = tab[tab.sum(1) > 0][:, tab.sum(0) > 0]
    if tab.shape[0] < 2 or tab.shape[1] < 2:
        return float("nan"), tab
    e = tab.sum(1, keepdims=True) * tab.sum(0, keepdims=True) / tab.sum()
    return float(((tab - e) ** 2 / e).sum()), tab


def perm_chi2(groups, values, nperm=NPERM):
    """Permutation chi-square of a groups x values table; -> (chi2, cramer V, p)."""
    g_lv, gi = np.unique(groups, return_inverse=True)
    v_lv, vi = np.unique(values, return_inverse=True)

    def table(v):
        t = np.zeros((len(g_lv), len(v_lv)))
        np.add.at(t, (gi, v), 1)
        return t

    obs, tab = chi2_stat(table(vi))
    if not np.isfinite(obs):
        return None
    k = min(tab.shape) - 1
    V = math.sqrt(obs / (tab.sum() * k)) if k > 0 else float("nan")
    ge = 0
    for _ in range(nperm):
        s, _ = chi2_stat(table(rng.permutation(vi)))
        ge += s >= obs - 1e-9
    return {"chi2": obs, "cramer_v": V, "p": (1 + ge) / (1 + nperm), "dof": int((tab.shape[0] - 1) * (tab.shape[1] - 1))}


def nn_index(D):
    D = D.copy()
    np.fill_diagonal(D, np.inf)
    # nearest relative by edges; ties -> the closest in leaf order (sister side)
    n = D.shape[0]
    pos = np.arange(n)
    tie = np.abs(pos[:, None] - pos[None, :]) / (10.0 * n)
    return np.argmin(D + tie, axis=1)


def state_signal(is_x, nn, groups=None, nperm=NPERM):
    """Share of X cells whose nearest relative is X, vs permutations (globally or within groups)."""
    is_x = np.asarray(is_x, dtype=bool)
    nx = int(is_x.sum())
    if nx < 3 or nx > len(is_x) - 3:
        return None
    obs = float(np.mean(is_x[nn][is_x]))
    null = np.empty(nperm)
    for k in range(nperm):
        if groups is None:
            w = rng.permutation(is_x)
        else:
            w = is_x.copy()
            for g in np.unique(groups):
                idx = np.where(groups == g)[0]
                w[idx] = rng.permutation(is_x[idx])
        null[k] = np.mean(w[nn][w])
    sd = float(null.std()) or 1e-9
    return {
        "n_state": nx,
        "observed": obs,
        "expected": float(null.mean()),
        "enrichment": obs / max(float(null.mean()), 1e-9),
        "z": (obs - float(null.mean())) / sd,
        "p": float((1 + (null >= obs - 1e-12).sum()) / (1 + nperm)),
    }


def stouffer(zs, ws=None):
    zs = np.asarray(zs, dtype=float)
    ws = np.ones(len(zs)) if ws is None else np.asarray(ws, dtype=float)
    Z = float((ws * zs).sum() / math.sqrt((ws ** 2).sum()))
    return Z, float(stats.norm.sf(Z))


def cochran_q(zs):
    zs = np.asarray(zs, dtype=float)
    if len(zs) < 2:
        return None
    Q = float(((zs - zs.mean()) ** 2).sum())
    return {"Q": Q, "df": len(zs) - 1, "p": float(stats.chi2.sf(Q, len(zs) - 1))}


def mantel_haenszel(tables):
    """2x2 tables [[a,b],[c,d]] (rows carrier / non-carrier, cols state / other).
    -> pooled MH odds ratio, Robins-Breslow-Greenland 95% CI, MH chi-square p, Cochran Q heterogeneity of the log ORs."""
    R = S = PR = PSQR = QS = 0.0
    sa = se = sv = 0.0
    lors, wts = [], []
    for (a, b), (c, d) in tables:
        n = a + b + c + d
        if n == 0:
            continue
        r_, s_ = a * d / n, b * c / n
        p_, q_ = (a + d) / n, (b + c) / n
        R += r_
        S += s_
        PR += p_ * r_
        PSQR += p_ * s_ + q_ * r_
        QS += q_ * s_
        m1, n1, m0 = a + c, a + b, b + d
        sa += a
        se += n1 * m1 / n
        if n > 1:
            sv += n1 * (n - n1) * m1 * m0 / (n * n * (n - 1))
        aa, bb, cc, dd = a + 0.5, b + 0.5, c + 0.5, d + 0.5
        lors.append(math.log(aa * dd / (bb * cc)))
        wts.append(1 / (1 / aa + 1 / bb + 1 / cc + 1 / dd))
    if R == 0 or S == 0:
        return None
    OR = R / S
    var = PR / (2 * R * R) + PSQR / (2 * R * S) + QS / (2 * S * S)
    se_l = math.sqrt(var) if var > 0 else float("nan")
    chi = (abs(sa - se) - 0.5) ** 2 / sv if sv > 0 else float("nan")
    lors, wts = np.array(lors), np.array(wts)
    pooled = (wts * lors).sum() / wts.sum()
    Qh = float((wts * (lors - pooled) ** 2).sum())
    return {
        "or": OR,
        "lo": math.exp(math.log(OR) - 1.96 * se_l) if np.isfinite(se_l) else None,
        "hi": math.exp(math.log(OR) + 1.96 * se_l) if np.isfinite(se_l) else None,
        "p": float(stats.chi2.sf(chi, 1)) if np.isfinite(chi) else None,
        "het_q": Qh,
        "het_p": float(stats.chi2.sf(Qh, len(lors) - 1)) if len(lors) > 1 else None,
    }


# ---------------------------------------------------------------- per patient

def event_class(label):
    if label.startswith("ecDNA"):
        return "ecdna"
    if label.endswith(" amp"):
        return "amp"
    if label.endswith(" del"):
        return "homdel"
    if label.endswith(" mut"):
        return "snv"
    return "gain" if label.startswith("+") else "loss"


def pretty(label):
    m = re.fullmatch(r"([+-])(\d+|X)([pq]?)", label)
    if m:
        sign, c, arm = m.groups()
        return f"chr{c}{arm} {'gain' if sign == '+' else 'loss'}"
    if label.endswith(" amp"):
        return label[:-4] + " amplification"
    if label.endswith(" del"):
        return label[:-4] + " deletion"
    if label.endswith(" mut"):
        return label[:-4] + " mutation"
    return label


def short(label):
    m = re.fullmatch(r"([+-])(\d+|X)([pq]?)", label)
    if m:
        sign, c, arm = m.groups()
        return f"{c}{arm}{'+' if sign == '+' else '−'}"
    return label


def patient_export(P, data, recs, arms):
    cells = [r for r in recs if r.get("entry_type") == "cell" and r.get("patient_id") == P]
    tumour = {r["pair"]: r for r in cells if str(r.get("clone_id") or "").lower() != "normal"}
    nodes, root = parse_newick(open(f"{data}/{P}/tree.nwk").read())
    r = prune(nodes, root, set(tumour))
    if r is None:
        return None
    T = Tree(nodes, r)
    row = {c: i for i, c in enumerate(T.leaves)}
    N = T.n
    log(f"== {P}: {N} tumour leaves")
    # --- per-cell CN calls
    carriers = defaultdict(set)
    n_cn = 0
    for c in T.leaves:
        calls = cell_cn_calls(f"{data}/{c}/complex.json", arms)
        if calls is None:
            continue
        n_cn += 1
        for k in calls:
            carriers[k].add(row[c])
    # --- driver SNVs (tier <= 2)
    try:
        ev = json.load(open(f"{data}/{P}/filtered.events.json"))
    except Exception:
        ev = []
    for e in ev:
        if e.get("type") in ("Missense", "Trunc", "Splice") and (e.get("Tier") or 9) <= 2:
            ids = [x for x in str(e.get("cell_ids") or "").split(",") if x in row]
            if ids:
                carriers[f"{e['gene']} mut"] |= {row[x] for x in ids}
    # --- ecDNA (curated circular walks with a focal driver gene)
    ecdna_cn = defaultdict(dict)
    try:
        w = json.load(open(f"{data}/{P}/walks.json"))
        for wk in w.get("walks", []):
            if not wk.get("circular") or wk.get("curated") is False:
                continue
            drivers = [g for g in re.split(r"[,; ]+", str(wk.get("driver_genes") or "")) if g in FOCAL_GENES]
            if not drivers:
                drivers = [g for g in wk.get("genes", []) if g in FOCAL_GENES]
            if not drivers:
                continue
            key = f"ecDNA {drivers[0]}"
            for c, cn in (wk.get("cells") or {}).items():
                if c in row and cn is not None:
                    ecdna_cn[key][c] = max(ecdna_cn[key].get(c, 0), cn)
    except Exception:
        pass
    for key, m in ecdna_cn.items():
        carriers[key] |= {row[c] for c, cn in m.items() if cn >= 4}
    # --- arm events redundant with their whole-chromosome event are dropped
    for k in list(carriers):
        m = re.fullmatch(r"([+-])(\d+|X)([pq])", k)
        if m:
            whole = f"{m.group(1)}{m.group(2)}"
            if whole in carriers:
                a, b = carriers[k], carriers[whole]
                if len(a & b) / max(1, len(a | b)) >= 0.8:
                    del carriers[k]
    # --- place
    min_c = max(3, int(math.ceil(0.02 * N)))
    events = []
    for k, rows in carriers.items():
        if len(rows) < min_c:
            continue
        v, f1, inc = T.best_clade(sorted(rows))
        frac = len(rows) / N
        # truncal: nearly all tumour cells, or the root is the best clade and most cells carry it;
        # a root-best event carried by few cells is scattered (not tree-coherent)
        truncal = frac >= TRUNCAL_FRAC or (v == 0 and frac >= ROOT_FRAC)
        if truncal:
            v = 0
            f1 = 2 * len(rows) / (N + len(rows))
        placed = truncal or (v != 0 and f1 >= MIN_F1)
        pv = T.parent[v] if v >= 0 else -1
        events.append({
            "id": k,
            "label": pretty(k),
            "short": short(k),
            "class": event_class(k),
            "n": len(rows),
            "frac": round(frac, 4),
            "node": int(v),
            "clade": T.size[v],
            "f1": round(f1, 3),
            "placed": bool(placed),
            "truncal": bool(truncal),
            "t0": 0.0 if truncal else round(T.time(pv) if pv >= 0 else 0.0, 4),
            "t1": 0.0 if truncal else round(T.time(v), 4),
            "carriers": sorted(int(x) for x in rows),
        })
    events.sort(key=lambda e: (not e["placed"], e["t1"], -e["n"]))
    # --- pairwise ordering
    placed = [e for e in events if e["placed"]]
    order = {}
    for a, b in combinations(placed, 2):
        if a["node"] == b["node"]:
            rel = "same"
        elif T.is_ancestor(a["node"], b["node"]):
            rel = "a_first"
        elif T.is_ancestor(b["node"], a["node"]):
            rel = "b_first"
        else:
            rel = "parallel"
        order[(a["id"], b["id"])] = rel
    tree_out = {
        "nodes": [[int(T.parent[i]), int(T.first[i]), int(T.last[i])] for i in range(len(T.parent))],
        "leaves": T.leaves,
        "time": [round(T.time(i), 4) for i in range(len(T.parent))],
    }
    return T, row, tumour, events, order, tree_out, n_cn


def lineages(T, row, tumour, rna_rows):
    """Lineage label per leaf row: clone labels if >= 2 clones with >= 5 RNA cells, else major tree clades."""
    clone = {row[c]: str(r.get("clone_id") or "") for c, r in tumour.items() if c in row}
    cnt = Counter(clone[i] for i in rna_rows if clone.get(i) and clone[i].lower() not in ("unassigned", "nan", "none", ""))
    good = [k for k, n in cnt.items() if n >= 5]
    if len(good) >= 2:
        return {i: (clone[i] if clone.get(i) in good else "other") for i in rna_rows}, "clone"
    # major tree clades: repeatedly replace the largest clade by its children of >= 10% of cells
    # (descending through caterpillar steps), until 3-4 clades or none can be split
    min_size = max(5, int(0.1 * T.n))
    groups, final = [0], set()
    for _ in range(500):
        open_ = [v for v in groups if v not in final]
        if not open_ or len(groups) >= 4:
            break
        v = max(open_, key=lambda v: T.size[v])
        kids = [k for k in T.children[v] if T.size[k] >= min_size]
        groups.remove(v)
        if not kids:
            groups.append(v)
            final.add(v)
        else:
            groups.extend(kids)
    if len(groups) < 2:
        return {i: "clade A" for i in rna_rows}, "tree"
    lab = {}
    ordered = sorted(groups, key=lambda v: T.first[v])
    for i in rna_rows:
        lab[i] = next((f"clade {chr(65 + k)}" for k, v in enumerate(ordered) if T.first[v] <= i <= T.last[v]), "other")
    return lab, "tree"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=f"{ROOT}/data")
    ap.add_argument("--datafiles", default=f"{ROOT}/datafiles.json")
    ap.add_argument("--out", default=f"{ROOT}/data/_cohort")
    ap.add_argument("--patients", nargs="*", default=PATIENTS)
    args = ap.parse_args()
    recs = json.load(open(args.datafiles))
    arms = load_arms()
    created = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    per = {}
    for P in args.patients:
        if not os.path.exists(f"{args.data}/{P}/tree.nwk"):
            continue
        res = patient_export(P, args.data, recs, arms)
        if res:
            per[P] = res

    # ============================================================ aligned trees
    old = {P: not os.path.exists(f"{args.data}/{P}/snv_matrix.json") for P in per}
    pres = defaultdict(dict)  # event -> patient -> event dict
    for P, (T, row, tumour, events, order, tree_out, n_cn) in per.items():
        for e in events:
            pres[e["id"]][P] = e
    recurrent = [k for k, m in pres.items() if sum(e["placed"] for e in m.values()) >= 2]

    def pair_stats(a, b, patients):
        c = Counter()
        detail = []
        for P in patients:
            m = per[P][4]
            ea, eb = pres[a].get(P), pres[b].get(P)
            if not ea or not eb:
                continue
            if not (ea["placed"] and eb["placed"]):
                c["unplaced"] += 1
                detail.append({"patient": P, "rel": "unplaced"})
                continue
            rel = m.get((a, b))
            if rel is None:
                rel = {"a_first": "b_first", "b_first": "a_first"}.get(m.get((b, a)), m.get((b, a)))
            c[rel] += 1
            detail.append({"patient": P, "rel": rel})
        return c, detail

    pairs = []
    rec_sorted = sorted(recurrent)
    for a, b in combinations(rec_sorted, 2):
        if a.split()[0] == b.split()[0] and a != b and event_class(a) == event_class(b):
            continue
        for scope, pts in (("all", list(per)), ("current", [P for P in per if not old[P]])):
            c, detail = pair_stats(a, b, pts)
            both = c["a_first"] + c["b_first"] + c["same"] + c["parallel"]
            if both < 2:
                continue
            # orient so that A is the event that comes first more often
            A, B, af, bf = (a, b, c["a_first"], c["b_first"]) if c["a_first"] >= c["b_first"] else (b, a, c["b_first"], c["a_first"])
            if A != a:
                detail = [{**d, "rel": {"a_first": "b_first", "b_first": "a_first"}.get(d["rel"], d["rel"])} for d in detail]
            pairs.append({
                "scope": scope, "a": A, "b": B, "a_first": af, "b_first": bf, "same": c["same"], "parallel": c["parallel"],
                "unplaced": c["unplaced"], "both_placed": both, "p": sign_test(af, bf), "patients": detail,
            })

    def label(k):
        return pretty(k)

    # per-event summary: truncal / subclonal counts and median time
    ev_rows = []
    for k in recurrent:
        m = pres[k]
        placed = [e for e in m.values() if e["placed"]]
        ev_rows.append({
            "id": k, "label": label(k), "short": short(k), "class": event_class(k),
            "n_patients": len(m), "n_placed": len(placed),
            "n_truncal": sum(e["truncal"] for e in placed),
            "n_subclonal": sum(not e["truncal"] for e in placed),
            "median_time": float(np.median([e["t1"] for e in placed])) if placed else None,
        })
    ev_rows.sort(key=lambda r: (r["median_time"] if r["median_time"] is not None else 9, -r["n_placed"]))

    # headline findings (plain language with numbers)
    findings = []
    all_pairs = [p for p in pairs if p["scope"] == "all"]

    def ptxt(p):
        return "" if p is None else (f"sign test p = {p:.2g}" if p >= 0.001 else "sign test p < 0.001")

    ordered = [p for p in all_pairs if p["a_first"] + p["b_first"] >= 2]
    ordered.sort(key=lambda p: (-(p["a_first"] - p["b_first"]), p["p"] if p["p"] is not None else 1, -p["both_placed"]))
    for p in ordered[:4]:
        n_ord = p["a_first"] + p["b_first"]
        txt = f"{label(p['a'])} precedes {label(p['b'])} in {p['a_first']}/{n_ord} patients where the two are ordered"
        extra = []
        if p["same"]:
            extra.append(f"{p['same']} more on the same branch")
        if p["parallel"]:
            extra.append(f"{p['parallel']} on parallel lineages")
        if extra:
            txt += f" ({'; '.join(extra)})"
        findings.append({"kind": "order", "text": txt + ".", "stat": ptxt(p["p"]), "a": p["a"], "b": p["b"],
                         "patients": [d["patient"] for d in p["patients"] if d["rel"] != "unplaced"]})
    # the canonical GBM question, always answered when the events exist
    for a, b in (("+7", "EGFR amp"), ("-10", "EGFR amp"), ("+7", "-10"), ("CDKN2A del", "EGFR amp")):
        hit = next((p for p in all_pairs if {p["a"], p["b"]} == {a, b}), None)
        if not hit or any(f.get("a") in (a, b) and f.get("b") in (a, b) for f in findings):
            continue
        A, B = hit["a"], hit["b"]
        txt = (f"{label(A)} and {label(B)}: {label(A)} first in {hit['a_first']}, {label(B)} first in {hit['b_first']}, "
               f"same branch in {hit['same']}, parallel in {hit['parallel']} of {hit['both_placed']} patients with both placed on the tree.")
        findings.append({"kind": "order", "text": txt, "stat": ptxt(hit["p"]), "a": A, "b": B,
                         "patients": [d["patient"] for d in hit["patients"] if d["rel"] != "unplaced"]})
    trunc = sorted([r for r in ev_rows if r["n_truncal"] >= 2], key=lambda r: -r["n_truncal"])[:5]
    if trunc:
        findings.insert(0, {
            "kind": "truncal",
            "text": "Most often truncal (on the trunk of the tree): " + ", ".join(f"{r['label']} ({r['n_truncal']}/{r['n_placed']})" for r in trunc) + ".",
            "events": [r["id"] for r in trunc],
            "patients": sorted({P for r in trunc for P, e in pres[r["id"]].items() if e["truncal"]}),
        })
        sub = sorted([r for r in ev_rows if r["n_subclonal"] >= 2], key=lambda r: -r["n_subclonal"])[:4]
        if sub:
            findings.insert(1, {
                "kind": "subclonal",
                "text": "Recurrently subclonal (arose after the trunk): " + ", ".join(f"{r['label']} ({r['n_subclonal']}/{r['n_placed']})" for r in sub) + ".",
                "events": [r["id"] for r in sub],
                "patients": sorted({P for r in sub for P, e in pres[r["id"]].items() if e["placed"] and not e["truncal"]}),
            })

    aligned = {
        "format": "gos-sc-aligned-trees/1",
        "created": created,
        "method": {
            "time": "relative evolutionary time of a clade = 1 - log(clade size) / log(tumour cells): 0 = trunk, 1 = single cells",
            "placement": f"best-F1 clade of the carrier cells; placed when F1 >= {MIN_F1} and >= max(3, 2% of tumour cells) carry it; truncal when the clade is the tumour root or >= {int(TRUNCAL_FRAC * 100)}% carry it",
            "arm_calls": f"length-weighted median CN of the arm vs the cell's genome-wide median, +-{ARM_DELTA}; chromosome = both arms",
            "focal": "amplification CN >= max(baseline + 3, 2 x baseline); deletion CN 0; genes: " + ", ".join(FOCAL_GENES),
            "ecdna": "curated circular walks (walks.json) with a focal driver gene, >= 4 copies in the cell",
            "snv": "driver SNVs of tier <= 2 in filtered.events.json (carrier cells as listed)",
            "test": "two-sided sign test of A-first vs B-first over patients where both are placed and nested",
        },
        "patients": {},
        "events": ev_rows,
        "pairs": pairs,
        "findings": findings,
    }
    for P, (T, row, tumour, events, order, tree_out, n_cn) in per.items():
        aligned["patients"][P] = {
            "n_tumour": T.n, "n_cn": n_cn, "old_tree": old[P],
            "has_clones": len({str(r.get("clone_id")) for r in tumour.values()} - {"Unassigned", "None", "nan"}) >= 2,
            "tree": tree_out, "events": events,
        }

    # ============================================================ associations
    her = {}
    try:
        her = json.load(open(HERITABILITY))
    except Exception:
        pass
    try:
        mods = {m["id"]: m for m in json.load(open(MODULES))["modules"]}
    except Exception:
        mods = {}
    assoc = {"format": "gos-sc-associations/1", "created": created, "patients": {}, "concordance": [], "event_state": [], "findings": []}
    signal_rows = defaultdict(dict)
    es_tables = defaultdict(dict)
    for P, (T, row, tumour, events, order, tree_out, n_cn) in per.items():
        st = {}
        ph = {}
        for c, r in tumour.items():
            if c in row and r.get("state") in STATES:
                st[row[c]] = r["state"]
                if r.get("Phase"):
                    ph[row[c]] = r["Phase"] in ("S", "G2M")
        rna_rows = sorted(st)
        if len(rna_rows) < 20:
            assoc["patients"][P] = {"n_rna": len(rna_rows), "skipped": "fewer than 20 tumour cells with an RNA state on the tree"}
            continue
        lin, lin_src = lineages(T, row, tumour, rna_rows)
        groups = np.array([lin[i] for i in rna_rows])
        states = np.array([st[i] for i in rna_rows])
        ct = perm_chi2(groups, states)
        # composition per lineage
        comp = []
        for g in sorted(set(groups), key=lambda g: (g == "other", g)):
            m = groups == g
            comp.append({"lineage": g, "n": int(m.sum()), "counts": {s: int((states[m] == s).sum()) for s in STATES}})
        # per-state enrichment in each lineage (Fisher)
        enrich = []
        for g in set(groups):
            if g == "other":
                continue
            for s in STATES:
                a = int(((groups == g) & (states == s)).sum())
                b = int(((groups == g) & (states != s)).sum())
                c_ = int(((groups != g) & (states == s)).sum())
                d = int(((groups != g) & (states != s)).sum())
                if a + c_ < 3:
                    continue
                orr, p = stats.fisher_exact([[a, b], [c_, d]])
                enrich.append({"lineage": g, "state": s, "a": a, "b": b, "c": c_, "d": d,
                               "or": float((a + 0.5) * (d + 0.5) / ((b + 0.5) * (c_ + 0.5))), "p": float(p)})
        q = bh([e["p"] for e in enrich])
        for e, qq in zip(enrich, q):
            e["q"] = float(qq)
        enrich.sort(key=lambda e: e["p"])
        # phylogenetic signal per state (+ cycling)
        D = T.leaf_distance_matrix()
        idx = np.array(rna_rows)
        nn = nn_index(D[np.ix_(idx, idx)])
        sig = {}
        for s in STATES:
            g = state_signal(states == s, nn)
            if g:
                w = state_signal(states == s, nn, groups=groups) if len(set(groups)) > 1 else None
                g["within_lineage"] = w
                sig[s] = g
                signal_rows[s][P] = g
        cyc_rows = [i for i in rna_rows if i in ph]
        if len(cyc_rows) >= 20:
            sub = np.array([rna_rows.index(i) for i in cyc_rows])
            nn2 = nn_index(D[np.ix_(idx[sub], idx[sub])])
            g = state_signal(np.array([ph[i] for i in cyc_rows]), nn2)
            if g:
                sig["cycling"] = g
                signal_rows["cycling"][P] = g
        # event-anchored tables (subclonal placed events)
        for e in events:
            if not e["placed"] or e["truncal"]:
                continue
            carr = set(e["carriers"])
            isc = np.array([i in carr for i in rna_rows])
            if isc.sum() < 5 or (~isc).sum() < 5:
                continue
            for s in STATES:
                x = states == s
                es_tables[(e["id"], s)][P] = [[int((isc & x).sum()), int((isc & ~x).sum())], [int((~isc & x).sum()), int((~isc & ~x).sum())]]
        assoc["patients"][P] = {
            "n_rna": len(rna_rows), "lineage_source": lin_src, "old_tree": old[P],
            "clone_state": ct, "composition": comp, "enrichment": enrich[:12], "signal": sig,
            # leaf rows (into aligned_trees tree.leaves) with state and lineage, for the drill-down strips
            "cells": [[int(i), st[i], lin[i]] for i in rna_rows],
        }
        log(f"  {P}: {len(rna_rows)} RNA cells, lineages={lin_src} {Counter(groups)}; clone x state p={ct and round(ct['p'], 4)} V={ct and round(ct['cramer_v'], 2)}")

    # meta-program signal from the story heritability (7 patients)
    for P, h in her.items():
        for t in h.get("tests", []):
            ph = t.get("phenotype", "")
            if not re.fullmatch(r"MP\d+", ph) or not mods.get(ph, {}).get("meta"):
                continue
            sd = None
            z = None
            if t.get("p") is not None:
                z = float(stats.norm.isf(min(max(t["p"], 1e-6), 1 - 1e-6)))
            signal_rows[f"{ph}"][P] = {
                "observed": t["observed"], "expected": t["null_mean"], "z": z, "p": t["p"], "n_state": t["n"],
                "within_lineage": t.get("within_clone"), "source": "heritability.json", "label": t.get("label"),
            }
            _ = sd

    for ph, m in signal_rows.items():
        pts = sorted(m)
        zs = [m[P]["z"] for P in pts if m[P].get("z") is not None and np.isfinite(m[P]["z"])]
        if len(zs) < 2:
            continue
        ns = [m[P].get("n_state") or 1 for P in pts if m[P].get("z") is not None and np.isfinite(m[P]["z"])]
        Z, p = stouffer(np.clip(zs, -8, 8), np.sqrt(ns))
        k_sig = sum(1 for P in pts if (m[P].get("p") or 1) < 0.05)
        k_pos = sum(1 for z in zs if z > 0)
        within = [m[P]["within_lineage"]["p"] for P in pts if m[P].get("within_lineage")]
        lab = ph if not ph.startswith("MP") else (m[pts[0]].get("label") or ph)
        assoc["concordance"].append({
            "phenotype": ph, "label": lab if ph != "cycling" else "Cycling (S/G2M)",
            "kind": "state" if ph in STATES else "cycling" if ph == "cycling" else "program",
            "n_patients": len(zs), "k_significant": k_sig, "k_positive": k_pos,
            "stouffer_z": Z, "p": p, "sign_p": sign_test(k_pos, len(zs) - k_pos),
            "heterogeneity": cochran_q(np.clip(zs, -8, 8)),
            "k_within_lineage": sum(1 for x in within if x < 0.05), "n_within_lineage": len(within),
            "per_patient": {P: {k: m[P].get(k) for k in ("observed", "expected", "z", "p", "n_state")} | {"within_p": (m[P].get("within_lineage") or {}).get("p")} for P in pts},
        })
    qs = bh([c["p"] for c in assoc["concordance"]])
    for c, q in zip(assoc["concordance"], qs):
        c["q"] = float(q)
    assoc["concordance"].sort(key=lambda c: (c["kind"] != "state", c["p"]))

    for (eid, s), m in es_tables.items():
        if len(m) < 2:
            continue
        mh = mantel_haenszel(list(m.values()))
        if not mh:
            continue
        per_p = {}
        for P, ((a, b), (c, d)) in m.items():
            per_p[P] = {"table": [[a, b], [c, d]], "or": (a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5)), "p": float(stats.fisher_exact([[a, b], [c, d]])[1])}
        k_up = sum(1 for v in per_p.values() if v["or"] > 1)
        assoc["event_state"].append({"event": eid, "label": pretty(eid), "state": s, "n_patients": len(m), "k_up": k_up, **mh, "per_patient": per_p})
    qs = bh([r["p"] if r["p"] is not None else np.nan for r in assoc["event_state"]])
    for r, q in zip(assoc["event_state"], qs):
        r["q"] = float(q) if np.isfinite(q) else None
    assoc["event_state"].sort(key=lambda r: (r["p"] if r["p"] is not None else 1))

    # association findings
    F = assoc["findings"]
    pts_ct = {P: v["clone_state"] for P, v in assoc["patients"].items() if v.get("clone_state")}
    k = sum(1 for v in pts_ct.values() if v["p"] < 0.05)
    if pts_ct:
        med_v = float(np.median([v["cramer_v"] for v in pts_ct.values()]))
        F.append({"kind": "clone_state", "text": f"RNA state differs between lineages in {k}/{len(pts_ct)} patients (permutation chi-square p < 0.05; median Cramer's V {med_v:.2f}).",
                  "patients": sorted(P for P, v in pts_ct.items() if v["p"] < 0.05)})
    for c in assoc["concordance"]:
        if c["kind"] not in ("state", "cycling"):
            continue
        name = c["label"] + ("-like state" if c["kind"] == "state" else "")
        if c["k_significant"] >= 1 or c["p"] < 0.05:
            F.append({"kind": "heritability", "phenotype": c["phenotype"],
                      "text": f"{name} clusters on the tree in {c['k_significant']}/{c['n_patients']} patients (nearest-relative test p < 0.05); "
                              f"combined across patients Stouffer p = {c['p']:.2g}, heterogeneity p = {c['heterogeneity']['p']:.2g}"
                              + (f"; beyond lineage labels in {c['k_within_lineage']}/{c['n_within_lineage']}." if c["n_within_lineage"] else "."),
                      "patients": sorted(P for P, v in c["per_patient"].items() if (v.get("p") or 1) < 0.05)})
    none_states = [c["label"] for c in assoc["concordance"] if c["kind"] == "state" and c["k_significant"] == 0 and c["p"] >= 0.05]
    if none_states:
        F.append({"kind": "heritability_none", "text": "No consistent tree signal for: " + ", ".join(none_states) + ".", "patients": []})
    for r in assoc["event_state"][:3]:
        if r["p"] is None or r["p"] >= 0.05:
            continue
        direction = "more" if r["or"] > 1 else "less"
        F.append({"kind": "event_state", "event": r["event"], "state": r["state"],
                  "text": f"Cells carrying {r['label']} are {direction} often {r['state']}-like (Mantel-Haenszel OR {r['or']:.2f}, "
                          f"95% CI {r['lo'] or float('nan'):.2f}-{r['hi'] or float('nan'):.2f}, p = {r['p']:.2g}, q = {r['q']:.2g}; same direction in {r['k_up'] if r['or'] > 1 else r['n_patients'] - r['k_up']}/{r['n_patients']} patients).",
                  "patients": sorted(r["per_patient"])})
    if not any(f["kind"] == "event_state" for f in F) and assoc["event_state"]:
        r = assoc["event_state"][0]
        F.append({"kind": "event_state_none", "text": f"No genomic event has a significant cross-patient link to a cell state (best: {r['label']} with {r['state']}, MH OR {r['or']:.2f}, p = {r['p']:.2g}, {r['n_patients']} patients).",
                  "event": r["event"], "state": r["state"], "patients": sorted(r["per_patient"])})
    assoc["caveats"] = [
        f"Small n: {len(per)} patients; at most {max((c['n_patients'] for c in assoc['concordance']), default=0)} contribute to any concordance test, so a null result is weak evidence of absence.",
        "BWH69, BWH75, MGH285 and MGH301 are on older trees (no SNV-anchored placement); BWH75, MGH301 and MGH305 have no clone table, so their lineages are the major clades of the tree.",
        "Meta-program rows come from analysis/story heritability.json (7 patients, nearest-relative Spearman with 2000 permutations).",
        "Single-cell CN calls are noisy and RNA states come from the Seurat annotation of the matched cells.",
    ]
    os.makedirs(args.out, exist_ok=True)
    for name, obj in (("aligned_trees.json", aligned), ("associations.json", assoc)):
        path = f"{args.out}/{name}"
        tmp = path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(obj, fh, separators=(",", ":"), default=lambda x: float(x) if isinstance(x, np.floating) else int(x))
        os.replace(tmp, path)
        log(f"wrote {path} ({os.path.getsize(path) / 1e3:.0f} KB)")
    for f in findings + F:
        log("FINDING", f["text"], f.get("stat", ""))


if __name__ == "__main__":
    main()
