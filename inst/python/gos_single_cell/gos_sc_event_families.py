#!/usr/bin/env python3
"""Junction identity, per-cell read support and junction-keyed event families (stage 1) for the gOS single-cell export.

Inputs (per patient P)
  analysis/<P>/junction_identity_261018/junction_master.tsv, junction_support.tsv   (skilift sc_export_junction_tables on the
      re-collected junction set analysis/<P>/db/junctions_10kb/)
  data/<cell>/complex.json     per-cell balanced JaBbA graph (ALT connections)
  data/<P>/filtered.events.json  rows with type "Complex SV"

Outputs
  data/<P>/junction_identity.json  {junctions, cells, unmatched, support, support_measure, match_rate, ...}
  data/<P>/event_families.json     {families: [...]}
  data/<P>/filtered.events.json    complex rows get family_id (and junction_ids/n_unmatched_alts); backup .bak_<yymmdd>_families

Breakpoint convention (gGnome $json: source = left snode id, sink = -right snode id): for BOTH ends a positive signed id is the right end
of the interval with strand '+', a negative id the left end with strand '-'. (mut_junctions.py flips the sink end; that is not the
gGnome convention.) The master junctions come from a gGnome Junction object whose strand letters are the opposite of the graph's
(graph DEL-like = +-, master = -+), so master strands are flipped before comparing; ids keep the master's own letters.
"""
import argparse
import csv
import json
import os
import re
import shutil
import sys
import time
from collections import Counter, defaultdict

DATA = "/gpfs/commons/projects/imielinski_web/sclarke/gOS_GBM_sc/data"
ANALYSIS = "/gpfs/commons/groups/imielinski_lab/projects/single_cell_GBM/analysis"
TOL = 1000
FOOT_SLACK = 1000
BIG = 10
JACCARD = 0.5
CARRIER_FRAC = 0.5
MIN_BRIDGE_CELLS = 3  # bridging guard; 0 disables
BIG_RULE = "containment"  # events with > BIG ids: overlap coefficient |A&B|/min(|A|,|B|) >= 0.5 (default) or "jaccard"
CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]
SUPPORT_MEASURE = ("distinct supporting read names per cell (split reads, plus discordant pairs when the junction-calling run counted them), "
                   "taken as the maximum over the region tests merged into the junction (merged.counts.dt max_reads); a cell absent from "
                   "the table has no supporting read observed. The table does not separate split from discordant reads.")


def nk(s):
    return re.sub(r"[-_]", "", s).lower()


def norm_chrom(c):
    return re.sub(r"^chr", "", str(c))


# ------------------------------------------------------------------ per-cell graph
def cell_alts(graph):
    """complex.json dict -> [(cid, (chrom,pos,strand), (chrom,pos,strand), title, weight)] for ALT connections."""
    iv = {abs(int(d["iid"])): d for d in graph.get("intervals", [])}
    out = []
    for e in graph.get("connections", []):
        t = str(e.get("type", "")).upper()
        if t not in ("ALT", "") or e.get("source") is None or e.get("sink") is None:
            continue
        s, k = iv.get(abs(int(e["source"]))), iv.get(abs(int(e["sink"])))
        if not s or not k:
            continue
        b1 = (norm_chrom(s["chromosome"]), int(s["endPoint"] if e["source"] > 0 else s["startPoint"]), "+" if e["source"] > 0 else "-")
        b2 = (norm_chrom(k["chromosome"]), int(k["endPoint"] if e["sink"] > 0 else k["startPoint"]), "+" if e["sink"] > 0 else "-")
        w = e.get("weight")
        out.append((e.get("cid"), b1, b2, e.get("title"), float(w) if w is not None else 1.0))
    return out


# ------------------------------------------------------------------ master + matching
def read_master(path):
    m = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            m[r["id"]] = dict(chrom1=r["chrom1"], pos1=int(r["pos1"]), strand1=r["strand1"],
                              chrom2=r["chrom2"], pos2=int(r["pos2"]), strand2=r["strand2"])
    return m


def flip(s):
    return "+" if s == "-" else "-"


def build_index(master):
    """(chrom, graph strand) -> sorted [(pos, id, end)]. The master (gGnome Junction object) carries the opposite strand
    letter to the graph json (verified: a graph DEL-like is +- here and -+ in the master), so strands are flipped on entry."""
    idx = defaultdict(list)
    for jid, j in master.items():
        idx[(norm_chrom(j["chrom1"]), flip(j["strand1"]))].append((j["pos1"], jid, 1))
        idx[(norm_chrom(j["chrom2"]), flip(j["strand2"]))].append((j["pos2"], jid, 2))
    for v in idx.values():
        v.sort()
    return idx


def match_junction(b1, b2, master, idx, tol=TOL):
    """Best master id whose two breakpoints are within tol of b1 and b2 with the same strands (either end order), or None."""
    cands = {}
    for b in (b1, b2):
        for pos, jid, _ in idx.get((b[0], b[2]), ()):
            if abs(pos - b[1]) <= tol:
                cands[jid] = 1
    best, bd = None, None
    for jid in cands:
        j = master[jid]
        m1 = (norm_chrom(j["chrom1"]), j["pos1"], flip(j["strand1"]))
        m2 = (norm_chrom(j["chrom2"]), j["pos2"], flip(j["strand2"]))
        for p, q in ((m1, m2), (m2, m1)):
            if p[0] == b1[0] and q[0] == b2[0] and p[2] == b1[2] and q[2] == b2[2]:
                d1, d2 = abs(p[1] - b1[1]), abs(q[1] - b2[1])
                if d1 <= tol and d2 <= tol and (bd is None or d1 + d2 < bd):
                    best, bd = jid, d1 + d2
    return best


def _gkey(b):
    return (CHROMS.index(b[0]) if b[0] in CHROMS else 99, b[1])


def graph_junction_id(b1, b2):
    """Id (master strand letters = flipped graph strands) of a graph ALT connection, ends in genome order."""
    o = (b1, b2) if _gkey(b1) <= _gkey(b2) else (b2, b1)
    return f"{o[0][0]}:{o[0][1]}{flip(o[0][2])}|{o[1][0]}:{o[1][1]}{flip(o[1][2])}", o


def match_cells(cell_graphs, master, tol=TOL, add_graph_ids=True):
    """cell_graphs {cell: [(cid,b1,b2,title,w)]} ->
    (cells {cell:{cid:id}}, unmatched {cell:[cid]}, titles {id:Counter}, graph_master {id: junction dict}).
    ALT connections are matched to the re-collected master first. Those that match nothing (the re-collected set is filtered
    by read support and span) get a graph-derived id, clustered across the patient's cells at tol (greedy, first-seen
    coordinates), so that events can still be linked; they are listed in `unmatched` and flagged in junctions[id].source."""
    idx = build_index(master)
    cells, unmatched, titles = {}, {}, defaultdict(Counter)
    gmaster, gidx = {}, defaultdict(list)
    for c, alts in cell_graphs.items():
        cm, um = {}, []
        for cid, b1, b2, title, w in alts:
            if w <= 0:
                continue
            jid = match_junction(b1, b2, master, idx, tol)
            if jid is None:
                um.append(cid)
                if add_graph_ids:
                    jid = match_junction(b1, b2, gmaster, gidx, tol)
                    if jid is None:
                        jid, o = graph_junction_id(b1, b2)
                        if jid in master or jid in gmaster:
                            jid = jid + "~g"
                        gmaster[jid] = dict(chrom1="chr" + o[0][0], pos1=o[0][1], strand1=flip(o[0][2]),
                                            chrom2="chr" + o[1][0], pos2=o[1][1], strand2=flip(o[1][2]))
                        gidx[(o[0][0], o[0][2])].append((o[0][1], jid, 1))
                        gidx[(o[1][0], o[1][2])].append((o[1][1], jid, 2))
            if jid is not None:
                cm[str(cid)] = jid
                titles[jid][title or ""] += 1
        cells[c] = cm
        if um:
            unmatched[c] = um
    return cells, unmatched, titles, gmaster


# ------------------------------------------------------------------ events
def parse_loc(loc):
    out = []
    for part in str(loc or "").split(","):
        m = re.match(r"^\s*(?:chr)?([0-9XY]+):(\d+)-(\d+)\s*$", part)
        if m:
            out.append((m.group(1), int(m.group(2)), int(m.group(3))))
    return out


def inside(b, foot, slack=FOOT_SLACK):
    return any(b[0] == c and s - slack <= b[1] <= e + slack for c, s, e in foot)


def cell_event_junctions(alts, jmap, foot):
    """ids (matched) and number of unmatched ALT connections of one cell with both breakpoints inside the footprint."""
    ids, nun = set(), 0
    for cid, b1, b2, title, w in alts:
        if w <= 0 or not (inside(b1, foot) and inside(b2, foot)):
            continue
        j = jmap.get(str(cid))
        if j is None:
            nun += 1
        else:
            ids.add(j)
    return ids, nun


def event_junctions(ev, cell_alts_by_cell, cells_map):
    """Junction ids of one complex event: the calling cell's ALT ids inside the footprint plus ids found inside the footprint in >= 50% of carrier graphs."""
    foot = parse_loc(ev.get("Genome_Location"))
    carriers = [c for c in str(ev.get("cell_ids") or "").split(",") if c]
    ids, nun = set(), 0
    s = ev.get("sample")
    if s in cell_alts_by_cell:
        ids, nun = cell_event_junctions(cell_alts_by_cell[s], cells_map.get(s, {}), foot)
    cnt = Counter()
    ncar = 0
    for c in carriers:
        if c not in cell_alts_by_cell:
            continue
        ncar += 1
        cid_ids, _ = cell_event_junctions(cell_alts_by_cell[c], cells_map.get(c, {}), foot)
        cnt.update(cid_ids)
    shared = {j for j, n in cnt.items() if ncar and n / ncar >= CARRIER_FRAC}
    return sorted(ids | shared), nun, sorted(shared)


def link(a, b):
    """True when two junction-id sets belong to one family: share >=1 id; when either has > BIG ids, Jaccard >= JACCARD."""
    if not a or not b:
        return False
    inter = len(a & b)
    if not inter:
        return False
    if len(a) > BIG or len(b) > BIG:
        if BIG_RULE == "containment":
            return inter / min(len(a), len(b)) >= JACCARD
        return inter / len(a | b) >= JACCARD
    return True


def families(jsets, twin_groups=None):
    """Union-find over event indices. jsets list of sets; twin_groups list of index lists forced together. -> list of index lists."""
    n = len(jsets)
    par = list(range(n))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            par[ry] = rx

    for i in range(n):
        for j in range(i + 1, n):
            if link(jsets[i], jsets[j]):
                union(i, j)
    for g in twin_groups or ():
        for k in g[1:]:
            union(g[0], k)
    comp = defaultdict(list)
    for i in range(n):
        comp[find(i)].append(i)
    return list(comp.values())


def guarded_families(jsets, carriers, cell_has, twin_groups=None, min_cells=None):
    """Bridging-guarded families. jsets: list of id sets; carriers: list of carrier-cell sets; cell_has: {cell: set(ids in its graph)}.
    (a) two events are linked only if link() holds AND some shared id is present in the graphs of >= min_cells cells that carry BOTH events;
    (b) events with < min_cells carrier cells never create links; each is attached to the single family it shares most ids with
        (>= 1 shared id) after the strong events are grouped, without merging families;
    twin_groups are always together. Returns (components [list of index lists], links [{a, b, shared_ids, n_cells, kind}]) with a, b = list indices."""
    mc = MIN_BRIDGE_CELLS if min_cells is None else min_cells
    n = len(jsets)
    par = list(range(n))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            par[ry] = rx

    strong = [len(carriers[i]) >= mc for i in range(n)]
    links = []

    def evidence(i, k):
        both = carriers[i] & carriers[k]
        sh = sorted(jsets[i] & jsets[k])
        best = 0
        for j in sh:
            best = max(best, sum(1 for c in both if j in cell_has.get(c, ())))
        return sh, best

    for i in range(n):
        if not strong[i]:
            continue
        for k in range(i + 1, n):
            if strong[k] and link(jsets[i], jsets[k]):
                sh, nc = evidence(i, k)
                if nc >= mc:
                    union(i, k)
                    links.append(dict(a=i, b=k, shared_ids=sh, n_cells=nc, kind="shared_junction"))
    for g in twin_groups or ():
        for k in g[1:]:
            union(g[0], k)
            links.append(dict(a=g[0], b=k, shared_ids=sorted(jsets[g[0]] & jsets[k]), n_cells=len(carriers[g[0]] & carriers[k]), kind="twin"))
    # attach weak events to the family (of strong events) they share most ids with
    for i in range(n):
        if strong[i] or any(l["a"] == i or l["b"] == i for l in links):
            continue
        best, bk = 0, None
        for k in range(n):
            if k != i and strong[k] and jsets[i] & jsets[k]:
                c = len(jsets[i] & jsets[k])
                if c > best:
                    best, bk = c, k
        if bk is not None:
            union(i, bk)
            sh, nc = evidence(i, bk)
            links.append(dict(a=i, b=bk, shared_ids=sh, n_cells=nc, kind="weak_attach"))
    comp = defaultdict(list)
    for i in range(n):
        comp[find(i)].append(i)
    return list(comp.values()), links


def merge_footprint(foots):
    by = defaultdict(list)
    for foot in foots:
        for c, s, e in foot:
            by[c].append((s, e))
    out = []
    for c in sorted(by, key=lambda x: (CHROMS.index(x) if x in CHROMS else 99)):
        iv = sorted(by[c])
        cur = list(iv[0])
        for s, e in iv[1:]:
            if s <= cur[1]:
                cur[1] = max(cur[1], e)
            else:
                out.append((c, cur[0], cur[1]))
                cur = [s, e]
        out.append((c, cur[0], cur[1]))
    return out


def rule_name():
    return BIG_RULE + (f"+bridge_guard{MIN_BRIDGE_CELLS}" if MIN_BRIDGE_CELLS else "")


def build_families(patient, events, ev_idx, ev_juncs, ev_unmatched, cell_has=None):
    """events full list; ev_idx indices of complex rows; ev_juncs {idx: [ids]}. Returns list of family dicts (sorted by n_cells desc)."""
    jsets = [set(ev_juncs[i]) for i in ev_idx]
    tw = defaultdict(list)
    for k, i in enumerate(ev_idx):
        g = events[i].get("twin_group")
        if g:
            tw[g].append(k)
    carriers = [{c for c in str(events[i].get("cell_ids") or "").split(",") if c} for i in ev_idx]
    if MIN_BRIDGE_CELLS and cell_has is not None:
        comps, links = guarded_families(jsets, carriers, cell_has, list(tw.values()))
    else:
        comps, links = families(jsets, list(tw.values())), []
    fams = []
    for comp in comps:
        mem = [ev_idx[k] for k in comp]
        cells = set()
        for i in mem:
            cells |= {c for c in str(events[i].get("cell_ids") or "").split(",") if c}
        jids = sorted(set().union(*[set(ev_juncs[i]) for i in mem]))
        cls = Counter()
        for i in mem:
            cls[events[i].get("Variant") or events[i].get("vartype") or "?"] += int(events[i].get("n_cells") or 0)
        foot = merge_footprint([parse_loc(events[i].get("Genome_Location")) for i in mem])
        fams.append(dict(junction_ids=jids, footprint=",".join(f"{c}:{s}-{e}" for c, s, e in foot), cells=sorted(cells), n_cells=len(cells),
                         members=[dict(event_index=i, vartype=events[i].get("vartype"), Variant=events[i].get("Variant"),
                                       n_cells=int(events[i].get("n_cells") or 0), cell_ids=[c for c in str(events[i].get("cell_ids") or "").split(",") if c],
                                       n_junctions=len(ev_juncs[i]), n_unmatched_alts=ev_unmatched.get(i, 0)) for i in mem],
                         classifications=dict(cls.most_common()),
                         resolved=bool(jids),
                         links=[dict(a=ev_idx[l["a"]], b=ev_idx[l["b"]], shared_ids=l["shared_ids"], n_cells=l["n_cells"], kind=l["kind"])
                                for l in links if l["a"] in comp]))
    fams.sort(key=lambda f: (-f["n_cells"], f["members"][0]["event_index"]))
    for n, f in enumerate(fams, 1):
        f["family_id"] = f"{patient}:F{n}"
    return [dict(family_id=f.pop("family_id"), **f) for f in fams]


# ------------------------------------------------------------------ patient driver
def read_support(path, dirs_by_norm):
    sup = defaultdict(dict)
    miss = set()
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            c = dirs_by_norm.get(nk(r["cell"]))
            if c is None:
                miss.add(r["cell"])
                continue
            sup[r["id"]][c] = int(float(r["reads"]))
    return sup, miss


def run_patient(p, data=DATA, analysis=ANALYSIS, dry=False):
    pdir = f"{data}/{p}"
    cdirs = sorted(d for d in os.listdir(data) if d.startswith(p + "_") and os.path.exists(f"{data}/{d}/complex.json"))
    dirs_by_norm = {nk(d): d for d in cdirs}
    idd = f"{analysis}/{p}/junction_identity_261018"
    master = read_master(f"{idd}/junction_master.tsv")
    graphs = {c: cell_alts(json.load(open(f"{data}/{c}/complex.json"))) for c in cdirs}
    cells_map, unmatched, titles, gmaster = match_cells(graphs, master)
    n_alt = sum(1 for a in graphs.values() for x in a if x[4] > 0)
    n_unmatched = sum(len(v) for v in unmatched.values())
    n_match = n_alt - n_unmatched
    used = {j for v in cells_map.values() for j in v.values()}
    ncg = Counter(j for v in cells_map.values() for j in set(v.values()))
    jun = {}
    allj = dict(master)
    allj.update(gmaster)
    for jid, j in allj.items():
        if jid not in master and jid not in used:
            continue
        d = dict(j)
        d["source"] = "recollected" if jid in master else "graph"
        d["cls"] = titles[jid].most_common(1)[0][0] if jid in titles else None
        d["n_cells_graph"] = ncg[jid]
        jun[jid] = d
    sup, miss = read_support(f"{idd}/junction_support.tsv", dirs_by_norm)
    ident = dict(patient=p, generated=time.strftime("%Y-%m-%d %H:%M"), tolerance_bp=TOL,
                 master_source=f"analysis/{p}/db/junctions_10kb (re-collected, min_span 1e4)",
                 junctions=jun, cells={c: v for c, v in cells_map.items()}, unmatched=unmatched,  # ALT cids with no re-collected match; they carry a graph-derived id in `cells` (source "graph")
                 match=dict(n_alt_connections=n_alt, n_matched_recollected=n_match, rate=(n_match / n_alt if n_alt else None),
                            n_graph_only_ids=len(gmaster), n_master=len(master), n_master_seen_in_graphs=len(used & set(master))),
                 support={k: v for k, v in sup.items()}, support_measure=SUPPORT_MEASURE,
                 support_cells_not_in_gos=len(miss))
    ef = f"{pdir}/filtered.events.json"
    events = json.load(open(ef))
    ev_idx = [i for i, e in enumerate(events) if e.get("type") == "Complex SV"]
    ev_juncs, ev_un, ev_shared = {}, {}, {}
    for i in ev_idx:
        ids, nun, sh = event_junctions(events[i], graphs, cells_map)
        ev_juncs[i], ev_un[i], ev_shared[i] = ids, nun, sh
    cell_has = {c: set(v.values()) for c, v in cells_map.items()}
    fams = build_families(p, events, ev_idx, ev_juncs, ev_un, cell_has)
    fam_of = {m["event_index"]: f["family_id"] for f in fams for m in f["members"]}
    for i in ev_idx:
        events[i]["family_id"] = fam_of[i]
        events[i]["junction_ids"] = ev_juncs[i]
        events[i]["n_unmatched_alts"] = ev_un[i]
    summ = dict(patient=p, n_events=len(ev_idx), n_families=len(fams), match_rate=ident["match"]["rate"],
                n_unresolved_events=sum(1 for i in ev_idx if not ev_juncs[i]))
    if not dry:
        stamp = time.strftime("%y%m%d")
        for fn, obj in ((f"{pdir}/junction_identity.json", ident), (f"{pdir}/event_families.json", dict(patient=p, generated=ident["generated"],
                         family_rule=rule_name(),
                         rules=dict(footprint_slack_bp=FOOT_SLACK, link=f"share >= 1 junction id carried by >= {MIN_BRIDGE_CELLS} cells that carry both events (0 = guard off); events with < {MIN_BRIDGE_CELLS} carrier cells never bridge, they attach to the family sharing most ids; when either event has > 10 ids, {BIG_RULE} >= 0.5 ({'overlap coefficient |A&B|/min' if BIG_RULE == 'containment' else 'Jaccard'}); twin_group forced together",
                                    carrier_fraction=CARRIER_FRAC), families=fams))):
            if os.path.exists(fn) and not os.path.exists(f"{fn}.bak_{stamp}"):
                shutil.copy2(fn, f"{fn}.bak_{stamp}")
            json.dump(obj, open(fn, "w"), separators=(",", ":"))
        if not os.path.exists(f"{ef}.bak_{stamp}_families"):
            shutil.copy2(ef, f"{ef}.bak_{stamp}_families")
        json.dump(events, open(ef, "w"))
    return summ, fams


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", nargs="+", default=["BWH69", "BWH70", "BWH75", "MGH285", "MGH301", "MGH302", "MGH303", "MGH304", "MGH305"])
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--min-bridge-cells", type=int, default=3, help="bridging guard (0 = off)")
    ap.add_argument("--big-rule", choices=["jaccard", "containment"], default="containment")
    a = ap.parse_args(argv)
    global BIG_RULE, MIN_BRIDGE_CELLS
    BIG_RULE = a.big_rule
    MIN_BRIDGE_CELLS = a.min_bridge_cells
    for p in a.patients:
        s, _ = run_patient(p, dry=a.dry_run)
        print(json.dumps(s), flush=True)


if __name__ == "__main__":
    main()
