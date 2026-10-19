#!/usr/bin/env python3
"""Multi-modal patient portrait + cohort "canonical GBM biology vs novel" export for the gOS single-cell cohort.

Reads only exports that already exist under data/ (nothing is recomputed upstream, nothing is imputed):
  _cohort/aligned_trees.json   per-patient events (arm CN, focal amp, homdel, ecDNA, SNV) with carrier cells
                               (tree leaf order), truncal / subclonal placement and relative time
  _cohort/associations.json    lineage x RNA state, tree signal, event x state (arm events)
  _cohort/tsg_loss.json        per-cell tumour-suppressor loss codes (allelic CN)
  _cohort/hotspot_rescue.json  read-level evidence at TERT / IDH1 / H3-3A / PTEN / NF1 / TP53 ... hotspots
  _cohort/subclonal_drivers.json  QC'd driver SNVs with carrier fraction, tree category, read-level flags,
                               clade signature shifts
  _cohort/signatures_joint.json   subclonal / private vs truncal signature contrasts
  _cohort/ecdna_plasticity.json   per-cell ecDNA copy distributions
  <P>/rna/cells.json, <P>/filtered.events.json, <P>/clades.json, <P>/walks.json, <P>/metadata.json

Writes (new files; an existing file is first copied to <file>.bak_<yymmdd>):
  _cohort/portrait.json   canonical-vs-novel matrix, novel list, caveats, references, all patient portraits
  <P>/portrait.json       the patient's portrait card

Checks against canonical GBM biology (Neftel 2019 Cell; Brennan 2013 Cell; Snuderl 2011 Cancer Cell; Turner 2017
Nature). Verdicts: consistent / inconsistent / not_testable (data missing or underpowered) / not_applicable
(the feature is absent in this patient, nothing to test). The state <-> genotype links are carriers vs
non-carriers inside each patient (Mantel-Haenszel over informative patients, also stratified by lineage);
a patient is informative with >= MIN_GRP carriers and >= MIN_GRP non-carriers that have an RNA state. The four
Neftel links are the only pre-specified tests; nothing is selected by p value.

usage: gos_sc_portrait.py [--data DIR] [--patients P ...] [--dry]
"""
import argparse
import datetime
import importlib.util
import json
import math
import os
import shutil
import sys
from collections import Counter, defaultdict

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = "/gpfs/commons/projects/imielinski_web/sclarke/gOS_GBM_sc"
PATIENTS = ["BWH69", "BWH70", "BWH75", "MGH285", "MGH301", "MGH302", "MGH303", "MGH304", "MGH305"]
MIN_GRP = 5          # carriers and non-carriers with an RNA state for a patient to be informative
TRUNCAL_FRAC = 0.85  # same threshold as aligned_trees.json
MIN_CARRIERS = 5     # an event is "present" with at least this many carrier cells
MIN_EVENT_CELLS = 10  # and >= 10% of tumour cells, for loss / homdel to count as an alteration
LOSS_CODES = (2, 4)  # tsg_loss: LOH by loss, homozygous deletion
STATES = ["MES", "AC", "OPC", "NPC"]
CANON_GENES = {"TP53", "PTEN", "NF1", "PIK3CA", "PIK3R1", "EGFR", "TERT", "RB1", "ATRX", "IDH1", "IDH2", "H3-3A", "H3F3A", "BRAF", "PDGFRA", "CDK4", "MDM2", "SMARCA4", "PTPN11", "KRAS"}
RTK = ["EGFR", "PDGFRA", "MET"]
TODAY = datetime.date.today().strftime("%y%m%d")

# Neftel 2019: (gene event, state expected enriched, free-text)
NEFTEL = [
    ("egfr_ac", "EGFR amp", "AC", "EGFR amplification / mutation is associated with the AC-like state"),
    ("pdgfra_opc", "PDGFRA amp", "OPC", "PDGFRA amplification is associated with the OPC-like state"),
    ("cdk4_npc", "CDK4 amp", "NPC", "CDK4 amplification is associated with the NPC-like state"),
    ("nf1_mes", "NF1 loss", "MES", "NF1 alteration is associated with the MES-like state"),
]

FEATURES = [
    ("idh_wt", "IDH1/2, H3-3A wild type", "IDH1 R132, H3-3A K27M/G34 absent from read evidence (IDH-wild-type GBM)"),
    ("tert", "TERT promoter activation", "C228T / C250T mutation in about 80% of IDH-wild-type GBM (Brennan 2013), or a TERT promoter rearrangement (Peifer 2015)"),
    ("chr7", "chr7 gain, truncal", "chr7 gain on the trunk (Brennan 2013)"),
    ("chr10", "chr10 loss, truncal", "chr10 loss on the trunk (Brennan 2013)"),
    ("cdkn2a", "CDKN2A homdel, truncal", "CDKN2A/B homozygous deletion in about 50-60% of GBM"),
    ("egfr_ecdna", "EGFR on ecDNA, copy heterogeneity", "EGFR amplified on ecDNA with cell-to-cell copy variation (Turner 2017)"),
    ("rtk_mosaic", "Mosaic RTK amplification", "different RTKs (EGFR, PDGFRA, MET) amplified in different cells (Snuderl 2011)"),
    ("egfr_ac", "EGFR amp - AC-like", NEFTEL[0][3] + " (Neftel 2019)"),
    ("pdgfra_opc", "PDGFRA amp - OPC-like", NEFTEL[1][3] + " (Neftel 2019)"),
    ("cdk4_npc", "CDK4 amp - NPC-like", NEFTEL[2][3] + " (Neftel 2019)"),
    ("nf1_mes", "NF1 loss - MES-like", NEFTEL[3][3] + " (Neftel 2019)"),
    ("pten", "PTEN altered", "PTEN loss / mutation in about 40% of GBM"),
    ("nf1", "NF1 altered", "NF1 loss / mutation in about 15% of GBM"),
    ("tp53", "TP53 altered", "TP53 mutation / loss in about 30% of GBM"),
]

REFERENCES = [
    "Hadi et al. 2020 Cell: Distinct classes of complex structural variation uncovered across thousands of cancer genome graphs (JaBbA: the genome graphs behind every copy-number and junction call here).",
    "Choo et al. 2023 Nature Genetics: Most large structural variants in cancer genomes can be detected without long reads (short-read SV detection and the junction classes used here).",
    "Neftel et al. 2019 Cell: An integrative model of cellular states, plasticity and genetics for glioblastoma (MES/AC/OPC/NPC and their genetic associations).",
    "Brennan et al. 2013 Cell: The somatic genomic landscape of glioblastoma (TCGA).",
    "Snuderl et al. 2011 Cancer Cell: Mosaic amplification of multiple receptor tyrosine kinase genes in glioblastoma.",
    "Turner et al. 2017 Nature: Extrachromosomal oncogene amplification drives tumour evolution and genetic heterogeneity.",
]


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def load_cohort_evolution():
    spec = importlib.util.spec_from_file_location("gos_sc_cohort_evolution", os.path.join(HERE, "gos_sc_cohort_evolution.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def jload(path):
    with open(path) as f:
        return json.load(f)


def fnum(x, nd=3):
    return None if x is None or (isinstance(x, float) and not math.isfinite(x)) else round(float(x), nd)


def pfmt(p):
    if p is None:
        return "p = n/a"
    return "p < 0.001" if p < 0.001 else f"p = {p:.2g}" if p < 0.01 else f"p = {p:.2f}"


def pct(x):
    return f"{100 * x:.0f}%"


# ---------------------------------------------------------------- inputs

class Inputs:
    def __init__(self, data, patients):
        self.data = data
        self.patients = patients
        c = f"{data}/_cohort"
        self.aligned = jload(f"{c}/aligned_trees.json")
        self.assoc = jload(f"{c}/associations.json")
        self.tsg = jload(f"{c}/tsg_loss.json")
        self.hot = jload(f"{c}/hotspot_rescue.json")
        self.sub = jload(f"{c}/subclonal_drivers.json")
        self.sigj = jload(f"{c}/signatures_joint.json")
        self.ecpl = jload(f"{c}/ecdna_plasticity.json")
        self.evassoc = jload(f"{c}/event_associations.json")
        self.tmb = jload(f"{c}/tmb.json")
        self.rna = {}
        self.events = {}
        self.clades = {}
        self.walks = {}
        for p in patients:
            for name, store, rel in (("rna", self.rna, "rna/cells.json"), ("events", self.events, "filtered.events.json"), ("clades", self.clades, "clades.json"), ("walks", self.walks, "walks.json")):
                try:
                    store[p] = jload(f"{data}/{p}/{rel}")
                except Exception:
                    store[p] = None

    def leaves(self, p):
        return self.aligned["patients"][p]["tree"]["leaves"]

    def event(self, p, eid):
        for e in self.aligned["patients"][p]["events"]:
            if e["id"] == eid:
                return e
        return None

    def carriers(self, p, eid):
        e = self.event(p, eid)
        if not e:
            return set()
        lv = self.leaves(p)
        return {lv[i] for i in e["carriers"]}

    def rna_state(self, p):
        """cell id -> state for tumour cells with an RNA state and a DNA match."""
        out = {}
        d = self.rna.get(p)
        if not d:
            return out
        for c in d["cells"]:
            if c.get("state") in STATES and c.get("cell_id"):
                out[c["cell_id"]] = c["state"]
        return out

    def lineage(self, p):
        """cell id -> lineage label (associations.json, leaf rows into aligned_trees leaves)."""
        a = self.assoc["patients"].get(p) or {}
        lv = self.leaves(p)
        return {lv[i]: lin for i, st, lin in a.get("cells", []) if i < len(lv)}

    def tsg_codes(self, p, gene):
        t = self.tsg["patients"].get(p)
        if not t or gene not in t["genes"]:
            return None
        normal = t.get("normal") or [0] * len(t["cells"])
        return {c: code for c, code, nrm in zip(t["cells"], t["genes"][gene], normal) if not nrm}


# ---------------------------------------------------------------- state <-> genotype tests

def table(carrier_ids, states, state, strata=None):
    """2x2 [[a, b], [c, d]] = carriers (state, other), non-carriers (state, other), over cells with an RNA state.
    With strata (cell -> label) returns one table per stratum."""
    def one(cells):
        a = b = c = d = 0
        for cell in cells:
            s = states[cell]
            if cell in carrier_ids:
                a += s == state
                b += s != state
            else:
                c += s == state
                d += s != state
        return [[a, b], [c, d]]
    if strata is None:
        return one(states.keys())
    groups = defaultdict(list)
    for cell in states:
        groups[strata.get(cell, "other")].append(cell)
    return [(k, one(v)) for k, v in sorted(groups.items())]


def pooled_mh(mh, tabs):
    """MH estimate; with a zero cell that makes it inestimable, add 0.5 to every cell (Haldane-Anscombe) and flag it."""
    r = mh(tabs)
    if r and r.get("lo") is not None and r["or"] > 0 and math.isfinite(r["or"]):
        return r, False
    tabs2 = [((a + .5, b + .5), (c + .5, d + .5)) for (a, b), (c, d) in tabs]
    r = mh(tabs2)
    return r, True


def group_sizes(t):
    return t[0][0] + t[0][1], t[1][0] + t[1][1]


def link_test(inp, mh, event_carriers, state):
    """event_carriers: patient -> set of carrier cell ids. -> dict with per-patient tables and pooled estimates."""
    per, tabs, strat_tabs = [], [], []
    for p in inp.patients:
        car = event_carriers.get(p)
        states = inp.rna_state(p)
        if car is None or not states:
            continue
        t = table(car, states, state)
        nc, nn = group_sizes(t)
        row = {"patient": p, "table": t, "n_carriers": nc, "n_noncarriers": nn, "informative": nc >= MIN_GRP and nn >= MIN_GRP,
               "state_counts_carriers": dict(Counter(st for cell, st in states.items() if cell in car)),
               "state_counts_noncarriers": dict(Counter(st for cell, st in states.items() if cell not in car))}
        if nc >= 1 and nn >= 1:
            (a, b), (c, d) = t
            orr, pv = stats.fisher_exact(t)
            row["or"], row["p"] = fnum(orr), fnum(pv, 4)
            row["frac_carriers"] = fnum(a / nc) if nc else None
            row["frac_noncarriers"] = fnum(c / nn) if nn else None
        per.append(row)
        if row["informative"]:
            tabs.append(tuple(map(tuple, t)))
            lin = inp.lineage(p)
            for _, st in table(car, states, state, strata=lin):
                strat_tabs.append(tuple(map(tuple, st)))
    out = {"state": state, "per_patient": per, "n_informative": len(tabs)}
    if tabs:
        r, corrected = pooled_mh(mh, tabs)
        out["pooled"] = ({k: fnum(v) for k, v in r.items()} | {"continuity_corrected": corrected}) if r else None
        # within lineage: strata with both carriers and non-carriers only
        mixed = [t for t in strat_tabs if (t[0][0] + t[0][1]) > 0 and (t[1][0] + t[1][1]) > 0]
        if mixed:
            rl, corr_l = pooled_mh(mh, mixed)
            out["within_lineage"] = ({k: fnum(v) for k, v in rl.items()} | {"n_strata": len(mixed), "continuity_corrected": corr_l}) if rl else None
    return out


def link_verdict(res, state, label):
    """(verdict, evidence text, caveat or None). A link whose pooled MH test is inconclusive but whose within-lineage estimate
    excludes 1 in the Neftel direction is 'consistent' with a caveat (the UI marks it with an asterisk: not conclusive)."""
    inf = [r for r in res["per_patient"] if r["informative"]]
    if not res["per_patient"]:
        return "not_applicable", "no carriers of the event in this patient.", None
    if not inf:
        small = ", ".join(f"{r['patient']} {r['n_carriers']} carriers / {r['n_noncarriers']} non-carriers" for r in res["per_patient"][:6])
        return "not_testable", f"underpowered: need >= {MIN_GRP} carriers and >= {MIN_GRP} non-carriers with an RNA state ({small}).", None
    pooled = res.get("pooled")
    if not pooled or pooled.get("lo") is None:
        return "not_testable", "pooled estimate not estimable.", None
    txt = f"{state}-like share {pct(inf[0]['frac_carriers'])} in carriers vs {pct(inf[0]['frac_noncarriers'])} in non-carriers" if len(inf) == 1 else f"{len(inf)} informative patients"
    txt += f"; MH OR {pooled['or']:.2f} (95% CI {pooled['lo']:.2f}-{pooled['hi']:.2f}, {pfmt(pooled['p'])}{', 0.5 continuity correction' if pooled.get('continuity_corrected') else ''})"
    wl = res.get("within_lineage")
    wl_ok = wl and wl.get("lo") is not None
    if wl_ok:
        txt += f"; within lineage OR {wl['or']:.2f} ({wl['lo']:.2f}-{wl['hi']:.2f}, {wl['n_strata']} strata)"
    if pooled["lo"] > 1:
        return "consistent", txt + ("" if not wl_ok or wl["lo"] > 1 else "; within lineages the CI includes 1, so lineage confounding is not excluded"), None
    if pooled["hi"] < 1:
        return "inconsistent", txt + " (opposite direction)" + ("" if not wl_ok or wl["hi"] < 1 else "; within lineages the CI includes 1, so lineage confounding is not excluded"), None
    note = f" - no association detected (CI excludes OR above {pooled['hi']:.1f}), underpowered for smaller effects"
    if wl_ok and wl["lo"] > 1:
        # the Neftel direction holds within lineages (CI excludes 1) although the pooled test is not significant: call it
        # consistent, flagged as not conclusive
        return "consistent", txt + " - within lineages the effect points the Neftel way", \
            f"not conclusive: the overall (pooled) test is not significant ({pfmt(pooled['p'])}); the association holds within lineages only ({pfmt(wl['p'])}, {wl['n_strata']} strata)"
    elif wl_ok and wl["hi"] < 1:
        note = " - inconclusive overall, within lineages the effect is opposite"
    return "not_testable", txt + note, None


# ---------------------------------------------------------------- per-patient pieces

def tumour_ids(inp, p):
    return set(inp.leaves(p))


def event_summary(e):
    return {
        "id": e["id"], "label": e["label"], "cls": e["class"], "n": e["n"], "frac": fnum(e["frac"]),
        "truncal": bool(e["truncal"]), "placed": bool(e["placed"]), "t0": fnum(e.get("t0")), "t1": fnum(e.get("t1")),
    }


def pooled_tert(inp, p):
    """TERT promoter hotspot row written by the pooled driver scan import (gos_sc_pooled_drivers), or None."""
    return next((e for e in (inp.events.get(p) or []) if e.get("gene") == "TERT" and e.get("source") == "pooled driver scan" and str(e.get("type")) == "Promoter"), None)


def snv_list(inp, p):
    """Driver SNVs of one patient from filtered.events.json, subclonal_drivers.json and hotspot_rescue.json, de-duplicated."""
    n_t = inp.aligned["patients"][p]["n_tumour"]
    seen = {}
    for e in inp.events.get(p) or []:
        if e["type"] in ("SCNA", "Fusion", "Complex SV") or not e.get("gene"):
            continue
        try:
            n = int(e.get("n_cells") or 0)
        except ValueError:
            n = 0
        key = (e["gene"], (e.get("Variant") or "").split(" ")[0])
        seen[key] = {"gene": e["gene"], "variant": e.get("Variant"), "class": e["type"], "tier": e.get("Tier"), "n_cells": n,
                     "frac": fnum(n / n_t) if n_t else None, "source": "filtered.events", "qc_flags": []}
    for h in inp.sub["hits"]:
        if h["patient"] != p:
            continue
        key = (h["gene"], h.get("protein") or h["site"])
        flags = (h.get("read_qc") or {}).get("flags") or []
        rec = {"gene": h["gene"], "variant": h.get("protein") or h["site"], "class": h["class"], "tier": h["tier"], "driver_class": h.get("driver_class"),
               "n_cells": h["n_carriers"], "frac": fnum(h.get("carrier_fraction")), "category": (h.get("tree") or {}).get("category"),
               "clade_f1": fnum(h.get("clade_f1")), "source": "subclonal_drivers", "qc_flags": flags, "site": h["site"], "oncokb": h.get("oncokb_effect")}
        # prefer the dropout-aware record; drop an equal filtered.events duplicate by gene+protein
        for k in [k for k in seen if k[0] == h["gene"] and k[1] and k[1] == (h.get("protein") or "")]:
            seen.pop(k)
        seen[key] = rec
    hot = inp.hot["patients"].get(p) or {}
    for s in hot.get("sites", []):
        if s.get("called") and s.get("level") in ("supported", "weak") and not s.get("discovered"):
            key = (s.get("gene") or s["name"].split()[0], s["name"])
            if not any(k[0] == key[0] and (k[1].split(" ")[0] in s["name"] or s["name"] in (k[1] or "")) for k in seen):
                seen[key] = {"gene": key[0], "variant": s["name"], "class": "hotspot rescued", "tier": s.get("tier"), "n_cells": s.get("carriers"),
                             "frac": fnum(s["carriers"] / s["tumour_cells"]) if s.get("tumour_cells") else None, "source": "hotspot_rescue",
                             "qc_flags": [], "evidence": s.get("level")}
    keep = [v for v in seen.values() if v.get("tier") in (1, 2, "1", "2") or v["gene"] in CANON_GENES or v.get("source") != "filtered.events"]
    keep.sort(key=lambda v: -(v.get("frac") or 0))
    return keep


def gene_alteration(inp, p, gene, snvs):
    """PTEN / NF1 / TP53: SNV and allelic-loss evidence. -> (altered?, evidence text, numbers)."""
    sn = [s for s in snvs if s["gene"] == gene and s.get("class") != "Promoter"]
    sn = [s for s in sn if (s.get("tier") in (1, 2, "1", "2") or s.get("driver_class") in ("likely", "possible") or s.get("source") == "hotspot_rescue") and not s["qc_flags"]]
    codes = inp.tsg_codes(p, gene)
    nloss = None
    if codes is not None and codes:
        nloss = sum(1 for v in codes.values() if v in LOSS_CODES) / len(codes)
        n_hom = sum(1 for v in codes.values() if v == 4) / len(codes)
    else:
        n_hom = None
    parts, altered = [], False
    for s in sn:
        altered = True
        parts.append(f"{s['variant']} in {pct(s['frac'])} of cells" if s.get("frac") is not None else f"{s['variant']}")
    # an event counts only when seen in >= 10% of tumour cells AND >= MIN_EVENT_CELLS cells (a homdel call in 1-5 cells is noise)
    ncodes = len(codes) if codes else 0
    n_loss_cells = round(nloss * ncodes) if nloss is not None else 0
    n_hom_cells = round(n_hom * ncodes) if n_hom is not None else 0
    if nloss is not None and nloss >= 0.10 and n_loss_cells >= MIN_EVENT_CELLS:
        altered = True
        kind = "homozygous deletion" if n_hom and n_hom >= 0.5 * nloss else "loss / LOH (single copy; no second hit unless a mutation is listed)"
        parts.append(f"{kind} in {pct(nloss)} of tumour cells ({n_loss_cells}/{ncodes} cells; {'truncal' if nloss >= TRUNCAL_FRAC else 'subclonal'})")
        if n_hom_cells >= MIN_EVENT_CELLS and kind.startswith("loss"):
            parts.append(f"of which homozygous deletion in {n_hom_cells} cells ({pct(n_hom)})")
    num = {"gene": gene, "loss_frac": fnum(nloss), "homdel_frac": fnum(n_hom), "snvs": [(s["variant"], s["frac"]) for s in sn]}
    if codes is None and not sn:
        return None, f"{gene}: allelic CN not available (no tsg_loss entry) and no SNV recorded", num
    bi = None
    if sn and codes:
        car = set()
        for e in inp.events.get(p) or []:
            if e.get("gene") == gene and e["type"] not in ("SCNA", "Fusion", "Complex SV") and e.get("Tier") in (1, 2, "1", "2"):
                car |= {x for x in str(e.get("cell_ids") or "").split(",") if x}
        for h in inp.sub["hits"]:
            if h["patient"] == p and h["gene"] == gene and not (h.get("read_qc") or {}).get("flags"):
                car |= set(h.get("carrier_ids") or [])
        both = {c for c in car if codes.get(c) in LOSS_CODES}
        if car:
            bi = len(both) / len(codes)
            if bi >= 0.10:
                parts.append(f"mutation plus loss / LOH in the same cell for {pct(bi)} of tumour cells ({len(both)} cells)")
    biallelic = bi is not None and bi >= 0.10
    num["biallelic_frac"] = fnum(bi)
    num.update({"loss_cells": n_loss_cells, "homdel_cells": n_hom_cells, "tumour_cells": ncodes})
    return altered, "; ".join(parts) if parts else f"no {gene} SNV; loss / homozygous deletion in < 10% of tumour cells ({n_loss_cells} loss, {n_hom_cells} homdel of {ncodes} cells)", num | {"biallelic": biallelic}


def build_checks(inp, mh, p, snvs):
    checks = {}
    ap = inp.aligned["patients"][p]
    n_t = ap["n_tumour"]
    hot = inp.hot["patients"].get(p) or {}
    sites = {s["name"]: s for s in hot.get("sites", [])}
    link = lambda **k: k

    # IDH-wt
    names = ["IDH1 R132H", "IDH1 R132C", "H3-3A K27M", "H3-3A G34R", "H3-3A G34V"]
    if not hot.get("assessed"):
        checks["idh_wt"] = ("not_testable", f"site evidence was not run for {p} ({hot.get('reason') or 'no hotspot evidence'}), so IDH1 / H3-3A hotspots were not examined.", {}, link(report=p, tab=15))
    else:
        called = [n for n in names if sites.get(n, {}).get("called") or sites.get(n, {}).get("level") in ("supported", "weak")]
        cov = [f"{n} {sites[n]['carriers']}/{sites[n]['covered']} ({sites[n]['level']})" for n in names if n in sites and sites[n].get("covered")]
        other = [s for s in snvs if s["gene"] in ("IDH1", "IDH2", "H3-3A", "H3F3A")]
        note = ""
        if other:
            note = " Non-hotspot " + ", ".join(f"{s['gene']} {s['variant']} ({pct(s['frac'])} of cells, tier {s['tier']}{', ' + '; '.join(s['qc_flags']) if s['qc_flags'] else ''})" for s in other) + " is not an R132 / K27M / G34 hotspot."
        if called:
            checks["idh_wt"] = ("inconsistent", "hotspot supported: " + ", ".join(called) + "." + note, {"called": called}, link(report=p, tab=15))
        else:
            checks["idh_wt"] = ("consistent", "no IDH1 R132 or H3-3A K27M / G34 hotspot called (alt reads at C>T-context sites are classified as artefact: " + "; ".join(cov[:4]) + ")." + note, {}, link(report=p, tab=15))

    # TERT promoter: a TERT promoter rearrangement (gos_sc_tert_sv: whole TERT CDS behind another locus) counts as TERT activation
    tert_sv = [e for e in (inp.events.get(p) or []) if e.get("tert_promoter_sv") in (True, "true", "TRUE")]
    if tert_sv:
        e = tert_sv[0]
        partners = sorted({str(x.get("fusion_genes") or x.get("gene")).split("::")[0] for x in tert_sv})
        cells = e.get("cells") or f"{e.get('n_cells')} cells"
        checks["tert"] = ("consistent", f"TERT promoter rearrangement: the whole TERT coding sequence is joined to {' / '.join(partners)} (5p15.31) in {cells} tumour cells, no normal cells; TERT activation by structural rearrangement rather than a C228T / C250T mutation.",
                          {"route": "rearrangement", "partners": partners, "cells": cells}, link(report=p, tab=7, heatmap="drivers", drvf="genes:TERT"))
    elif not hot.get("assessed") and pooled_tert(inp, p):
        # no per-site evidence for this patient, but the pooled driver scan (gos_sc_pooled_drivers) called the promoter hotspot
        e = pooled_tert(inp, p)
        ev = (e.get("hotspot_evidence") or {}).get("tumor") or {}
        name = "TERT " + str(e.get("Variant") or "").split(" ")[0]
        checks["tert"] = ("consistent", f"{name} in {e.get('n_cells')} of {e.get('covered_tumour_cells')} covered tumour cells (pooled driver scan: {ev.get('alt')} of {ev.get('dp')} reads alt); "
                          f"{e.get('n_normal_cells') or 0} normal cells with alt reads. Site evidence not run for {p}.",
                          {"route": "promoter mutation", "source": "pooled driver scan"}, link(report=p, tab=7, heatmap="drivers", drvf="genes:TERT"))
    elif not hot.get("assessed"):
        checks["tert"] = ("not_testable", f"site evidence not run for {p}.", {}, link(report=p, tab=15))
    else:
        ts = [sites[n] for n in ("TERT C228T", "TERT C250T") if n in sites]
        pos = [s for s in ts if s.get("called") and s.get("level") in ("supported", "weak")]
        if pos:
            s = pos[0]
            lvl = "supported" if s["level"] == "supported" else "weak (few covered cells)"
            checks["tert"] = ("consistent", f"{s['name']} in {s['carriers']} of {s['covered']} covered tumour cells ({lvl}); no alt reads in matched normal cells." if (s.get("normal") or {}).get("cells_alt", 0) == 0 else f"{s['name']} in {s['carriers']} of {s['covered']} covered cells ({lvl}).", {"site": s["name"], "carriers": s["carriers"], "covered": s["covered"]}, link(report=p, tab=15, evsite=s["site"]))
        else:
            ok = [s for s in ts if s.get("level") in ("absent", "supported", "weak") and (s.get("covered") or 0) >= 10]
            txt = "; ".join(f"{s['name']} {s.get('carriers')}/{s.get('covered')} cells" if s.get("covered") else f"{s['name']} not genotyped" for s in ts)
            if len(ok) == 2:
                checks["tert"] = ("inconsistent", f"not detected ({txt}); the minority case (about 20% of IDH-wild-type GBM lack a TERT promoter mutation), or another TERT activation route.", {"covered": [s["covered"] for s in ok]}, link(report=p, tab=15))
            else:
                checks["tert"] = ("not_testable", f"one promoter site not evaluable or too few covered cells ({txt}).", {}, link(report=p, tab=15))

    # truncal arm / homdel events
    for key, eid, what in (("chr7", "+7", "chr7 gain"), ("chr10", "-10", "chr10 loss"), ("cdkn2a", "CDKN2A del", "CDKN2A homozygous deletion")):
        e = inp.event(p, eid)
        lk = link(scview="aligned")
        if e is None or e["n"] < MIN_CARRIERS:
            if ap["n_cn"] < 20:
                checks[key] = ("not_testable", f"per-cell CN available for only {ap['n_cn']} cells.", {}, lk)
            elif key == "cdkn2a":
                checks[key] = ("not_applicable", f"no CDKN2A homozygous deletion ({e['n'] if e else 0} of {ap['n_cn']} cells); present in only about 50-60% of GBM, so absence is not itself inconsistent.", {"frac": 0}, lk)
            else:
                checks[key] = ("inconsistent", f"{what} not detected (carriers {e['n'] if e else 0} of {ap['n_cn']} cells with CN).", {"frac": 0}, lk)
        elif e["truncal"]:
            checks[key] = ("consistent", f"{what} truncal: {e['n']} of {ap['n_cn']} cells ({pct(e['frac'])}), placed on the {'root' if e['t0'] == 0 else 'clade at t = %.2f' % e['t0']}.", {"frac": e["frac"], "t0": e["t0"]}, lk)
        else:
            where = f"placed on a clade at relative time {e['t0']:.2f}" if e["placed"] else "not on a single clade"
            checks[key] = ("inconsistent", f"{what} is subclonal here: {pct(e['frac'])} of cells ({where}), not truncal.", {"frac": e["frac"], "t0": e.get("t0")}, lk)

    # EGFR ecDNA copy heterogeneity
    ec = inp.event(p, "ecDNA EGFR")
    egfr_amp = inp.event(p, "EGFR amp")
    row = next((r for r in inp.ecpl["rows"] if r["patient"] == p and r["gene"] == "EGFR"), None)
    if ec is not None:
        wl = [x for x in (inp.walks.get(p) or {}).get("walks", []) if x.get("circular") and x.get("curated") is not False and "EGFR" in str(x.get("driver_genes"))]
        percell = defaultdict(float)
        for x in wl:
            for c, cn in (x.get("cells") or {}).items():
                if cn is not None:
                    percell[c] = max(percell[c], cn)
        cn = np.array([v for v in percell.values() if v >= 4], dtype=float)
        lab = wl[0]["label"] if wl else None
        if len(cn) >= MIN_CARRIERS:
            cv = float(cn.std() / cn.mean())
            q10, q50, q90 = np.percentile(cn, [10, 50, 90])
            het = cv >= 0.3
            txt = f"EGFR on a circular ecDNA walk in {len(cn)} cells ({pct(ec['frac'])}); copies per carrier cell median {q50:.0f}, 10-90% range {q10:.0f}-{q90:.0f} ({q90 / max(q10, 1):.1f}-fold), CV {cv:.2f}"
            checks["egfr_ecdna"] = ("consistent" if het else "inconsistent", txt + ("; copy number varies cell to cell." if het else "; narrow copy distribution."), {"cv": fnum(cv), "median": fnum(q50, 1), "p10": fnum(q10, 1), "p90": fnum(q90, 1), "n": len(cn)}, link(report=p, tab=14, walk=lab) if lab else link(scview="ecplast"))
        else:
            checks["egfr_ecdna"] = ("not_testable", "EGFR ecDNA walk has too few carrier cells for a copy distribution.", {}, link(scview="ecplast"))
    elif egfr_amp is not None and egfr_amp["n"] >= MIN_CARRIERS:
        checks["egfr_ecdna"] = ("inconsistent", f"EGFR amplified in {egfr_amp['n']} cells ({pct(egfr_amp['frac'])}) with no curated EGFR ecDNA walk: chromosomal amplification or unresolved structure.", {}, link(scview="ecplast"))
    else:
        checks["egfr_ecdna"] = ("not_applicable", "no EGFR amplification in this patient.", {}, link(scview="ecplast"))

    # RTK mosaic
    rtk = [(g, inp.carriers(p, f"{g} amp")) for g in RTK]
    rtk = [(g, c) for g, c in rtk if len(c) >= MIN_CARRIERS]
    if len(rtk) < 2:
        checks["rtk_mosaic"] = ("not_applicable", ("only " + rtk[0][0] + " amplified") if rtk else "no RTK amplification", {}, link(scview="amplicons"))
    else:
        pairs, parts, excl = [], [], 0
        for i in range(len(rtk)):
            for j in range(i + 1, len(rtk)):
                (ga, ca), (gb, cb) = rtk[i], rtk[j]
                both = len(ca & cb)
                na, nb, n = len(ca), len(cb), n_t
                tab = [[both, na - both], [nb - both, n - na - nb + both]]
                _, pv = stats.fisher_exact(tab, alternative="less")
                exp_ = na * nb / n
                pairs.append({"a": ga, "b": gb, "both": both, "n_a": na, "n_b": nb, "expected_both": fnum(exp_, 1), "p_exclusive": fnum(pv, 4)})
                parts.append(f"{ga} {na} / {gb} {nb} cells, overlap {both} (expected {exp_:.1f}; exclusivity {pfmt(pv)})")
                excl += both < exp_ and pv < 0.05
        checks["rtk_mosaic"] = ("consistent" if excl else "not_testable", "; ".join(parts) + ("." if excl else ": overlap not lower than chance, no evidence of mutually exclusive subclones."), {"pairs": pairs}, link(scview="amplicons"))

    return checks


def neftel_links(inp, mh):
    """Cohort-level tests of the four Neftel links. Carriers: aligned_trees amp calls; NF1 loss = tsg_loss code 2/4."""
    ev = {"EGFR amp": {}, "PDGFRA amp": {}, "CDK4 amp": {}, "NF1 loss": {}}
    for p in inp.patients:
        for g in ("EGFR amp", "PDGFRA amp", "CDK4 amp"):
            c = inp.carriers(p, g)
            if c:
                ev[g][p] = c
        codes = inp.tsg_codes(p, "NF1")
        if codes:
            c = {k for k, v in codes.items() if v in LOSS_CODES}
            if c:
                ev["NF1 loss"][p] = c
    res = {}
    for key, event, state, text in NEFTEL:
        # a patient without an event is still a source of non-carriers only if it has none: those contribute nothing to MH, skip
        r = link_test(inp, mh, ev[event], state)
        v, t, cav = link_verdict(r, state, event)
        res[key] = {"event": event, "state": state, "claim": text, "verdict": v, "evidence": t, "caveat": cav, "test": r}
    return res, ev


def state_block(inp, p):
    a = inp.assoc["patients"].get(p) or {}
    states = inp.rna_state(p)
    cnt = Counter(states.values())
    n = sum(cnt.values())
    d = inp.rna.get(p)
    cyc = None
    if d:
        ph = [c["Phase"] for c in d["cells"] if c.get("state") in STATES and c.get("Phase")]
        if ph:
            cyc = fnum(sum(x != "G1" for x in ph) / len(ph))
    comp = [{"lineage": r["lineage"], "n": r["n"], "counts": r["counts"]} for r in a.get("composition", [])]
    sig = {}
    for s, v in (a.get("signal") or {}).items():
        sig[s] = {"n_state": v["n_state"], "enrichment": fnum(v["enrichment"]), "p": fnum(v["p"], 4)}
        if "within_lineage" in v:
            sig[s]["within_lineage_p"] = fnum(v["within_lineage"]["p"], 4)
    cs = a.get("clone_state") or {}
    return {
        "n_rna": n, "counts": dict(cnt), "fractions": {s: fnum(cnt.get(s, 0) / n) for s in STATES} if n else {}, "cycling_frac": cyc,
        "lineage_source": a.get("lineage_source"), "composition": comp,
        "lineage_state": {"cramers_v": fnum(cs.get("cramer_v")), "p": fnum(cs.get("p"), 4), "dof": cs.get("dof")} if cs else None,
        "tree_signal": sig,
    }


def lineage_text(sb):
    ls, sig = sb.get("lineage_state"), sb.get("tree_signal") or {}
    if not ls or ls["p"] is None:
        return "not available (no lineage x state test: too few RNA cells or no lineages)"
    follows = ls["p"] < 0.05
    hits = [f"{s} (observed/expected {v['enrichment']:.2f}, {pfmt(v['p'])}" + (f"; within lineage {pfmt(v['within_lineage_p'])})" if v.get("within_lineage_p") is not None else ")") for s, v in sig.items() if s in STATES and v["p"] is not None and v["p"] < 0.05]
    t = f"State {'differs between lineages' if follows else 'does not differ detectably between lineages'} (permutation chi-square {pfmt(ls['p'])}, Cramer's V {ls['cramers_v']:.2f}, {sb['n_rna']} RNA cells)."
    if hits:
        t += " Nearest relatives share state more than chance for " + ", ".join(hits) + "."
    return t


def evolution_block(inp, p):
    ap = inp.aligned["patients"][p]
    cl = (inp.clades.get(p) or {}).get("clades", [])
    evs = [event_summary(e) for e in ap["events"] if e["placed"] and (e["class"] in ("amp", "homdel", "ecdna") or e["id"] in ("+7", "-10"))]
    evs.sort(key=lambda e: (e["t0"] if e["t0"] is not None else 9, -e["frac"]))
    shifts = [s for s in inp.sub.get("signature_shifts", []) if s["patient"] == p and s.get("q") is not None and s["q"] < 0.05]
    joint = [c for c in inp.sigj["contrasts"] if c["patient"] == p and c["contrast"] == "subclonal vs truncal" and c["q"] is not None and c["q"] < 0.05]
    tm = next((x for x in inp.tmb["patients"] if x["patient"] == p), None)
    tmb = None
    if tm and tm.get("patient_summary"):
        ps = tm["patient_summary"]
        tmb = {"n_truncal_sites": (ps.get("truncal") or {}).get("n_sites"), "n_subclonal_sites": (ps.get("subclonal") or {}).get("n_sites"), "noisy": bool(tm.get("noisy"))}
    return {
        "n_tumour_cells": ap["n_tumour"], "n_clades": len(cl), "max_clade_depth": max([c.get("depth") or 0 for c in cl] or [0]),
        "old_tree": bool(ap.get("old_tree")), "has_clones": bool(ap.get("has_clones")),
        "timing": evs,
        "signature_contrasts": [{"signature": c["signature"], "frac_subclonal": c["frac_a"], "frac_truncal": c["frac_b"], "diff": c["diff"], "q": c["q"]} for c in joint],
        "clade_signature_shifts": [{"clade": s["clade"], "signature": s["signature"], "frac": s["frac"], "frac_ref": s["frac_ref"], "diff": s["diff"], "q": s["q"], "vs": s["vs"]} for s in shifts],
        "tmb": tmb,
    }


def qc_pass_subclonal(inp, p):
    out = []
    for h in inp.sub["hits"]:
        if h["patient"] != p or (h.get("read_qc") or {}).get("flags"):
            continue
        cat = (h.get("tree") or {}).get("category")
        fr = h.get("carrier_fraction")
        if fr is None or fr >= TRUNCAL_FRAC:
            continue
        out.append({"gene": h["gene"], "variant": h.get("protein") or h["site"], "class": h["class"], "frac": fnum(fr), "category": cat, "n_carriers": h["n_carriers"],
                    "tier": h["tier"], "driver_class": h.get("driver_class"), "site": h["site"], "oncokb": h.get("oncokb_effect"), "clade_f1": fnum(h.get("clade_f1"))})
    return out


def flagged_hits(inp, p):
    return [{"gene": h["gene"], "variant": h.get("protein") or h["site"], "frac": fnum(h.get("carrier_fraction")), "flags": h["read_qc"]["flags"]}
            for h in inp.sub["hits"] if h["patient"] == p and (h.get("read_qc") or {}).get("flags")]


# ---------------------------------------------------------------- driver


def build(args):
    ce = load_cohort_evolution()
    mh = ce.mantel_haenszel
    inp = Inputs(args.data, args.patients)
    neftel, ev = neftel_links(inp, mh)

    portraits = {}
    for p in args.patients:
        log("==", p)
        snvs = snv_list(inp, p)
        checks = build_checks(inp, mh, p, snvs)
        sb = state_block(inp, p)
        # canonical gene alterations
        for key, gene in (("pten", "PTEN"), ("nf1", "NF1"), ("tp53", "TP53")):
            alt, txt, num = gene_alteration(inp, p, gene, snvs)
            if alt is None:
                checks[key] = ("not_testable", txt, num, {"scview": "drivers"})
            elif alt:
                checks[key] = ("consistent", txt, num, {"scview": "drivers"})
            else:
                checks[key] = ("not_applicable", txt + " (no alteration here)", num, {"scview": "drivers"})
        # Neftel links: this patient's slice
        for key, event, state, text in NEFTEL:
            r = neftel[key]["test"]
            mine = next((x for x in r["per_patient"] if x["patient"] == p), None)
            if mine is None:
                checks[key] = ("not_applicable", f"no {event} carriers in this patient.", {}, {"scview": "eventassoc", "evassoc": event.replace(" loss", " homdel")})
                continue
            if mine["informative"]:
                a = mine["table"]
                t = f"{mine['n_carriers']} carriers vs {mine['n_noncarriers']} non-carriers with RNA: {state}-like {pct(mine['frac_carriers'])} vs {pct(mine['frac_noncarriers'])}, Fisher OR {mine['or']:.2f} ({pfmt(mine['p'])})."
                v = "consistent" if (mine["p"] < 0.05 and mine["or"] > 1) else "inconsistent" if (mine["p"] < 0.05 and mine["or"] < 1) else "not_testable"
                if v == "not_testable":
                    t += " Not significant in this patient alone (underpowered)."
                checks[key] = (v, t, {"table": a, "or": mine["or"], "p": mine["p"]}, {"scview": "eventassoc", "evassoc": event.replace(" loss", " homdel")})
            else:
                checks[key] = ("not_testable", f"underpowered: {mine['n_carriers']} carriers vs {mine['n_noncarriers']} non-carriers with an RNA state (need >= {MIN_GRP} each).", {"n_carriers": mine["n_carriers"], "n_noncarriers": mine["n_noncarriers"]}, {"scview": "eventassoc", "evassoc": event.replace(" loss", " homdel")})
        ordered = {k: {"verdict": v[0], "evidence": v[1], "numbers": v[2], "link": v[3]} for k, v in ((f[0], checks[f[0]]) for f in FEATURES)}

        geno_events = [event_summary(e) for e in inp.aligned["patients"][p]["events"] if e["class"] in ("amp", "homdel", "ecdna", "snv") or e["id"] in ("+7", "-10", "+7q", "-13", "-9p")]
        geno_events = [e for e in geno_events if e["n"] >= MIN_CARRIERS]
        geno_events.sort(key=lambda e: -e["frac"])
        walks = []
        for w in (inp.walks.get(p) or {}).get("walks", []):
            if w.get("circular") and w.get("curated") is not False and w.get("driver_genes"):
                walks.append({"label": w["label"], "genes": w["driver_genes"], "ncells": w["ncells"], "median_cn": w["median_cn"], "max_cn": w["max_cn"], "span": w["span"]})
        sub_pass = qc_pass_subclonal(inp, p)
        portraits[p] = {
            "format": "gos-sc-portrait/1", "patient": p, "created": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
            "n_tumour_cells": inp.aligned["patients"][p]["n_tumour"],
            "genotype": {"events": geno_events, "snvs": snvs, "ecdna_walks": walks, "subclonal_qc_pass": sub_pass, "flagged_hits": flagged_hits(inp, p)},
            "evolution": evolution_block(inp, p),
            "phenotype": dict(sb, lineage_text=lineage_text(sb)),
            "checks": ordered,
            "novel": [],
        }

    # ---------------- cohort
    matrix = {"features": [{"key": k, "label": l, "canonical": c} for k, l, c in FEATURES], "patients": args.patients,
              "cells": {p: {k: {"verdict": portraits[p]["checks"][k]["verdict"], "evidence": portraits[p]["checks"][k]["evidence"], "link": portraits[p]["checks"][k]["link"]} for k, _, _ in FEATURES} for p in args.patients}}
    tally = {}
    for k, _, _ in FEATURES:
        c = Counter(portraits[p]["checks"][k]["verdict"] for p in args.patients)
        tally[k] = {v: c.get(v, 0) for v in ("consistent", "inconsistent", "not_testable", "not_applicable")}
    cohort_links = {k: {"event": v["event"], "state": v["state"], "claim": v["claim"], "verdict": v["verdict"], "evidence": v["evidence"], "caveat": v.get("caveat"),
                        "per_patient": v["test"]["per_patient"], "pooled": v["test"].get("pooled"), "within_lineage": v["test"].get("within_lineage"), "n_informative": v["test"]["n_informative"]} for k, v in neftel.items()}

    inp.tmb_caveats = None
    try:
        t = jload(f"{args.data}/_cohort/tmb_associations.json")
        inp.tmb_caveats = t["caveats"][0]
    except Exception:
        pass
    novel = cohort_novel(inp, portraits, args.patients, neftel)
    for p in args.patients:
        portraits[p]["novel"] = [n for n in novel if p in n["patients"]]
    caveats = [
        "PRKDC truncating / splice hits are flagged as artefacts by read-level QC (alt on spliced exon-exon reads, alt at read-block ends, same allele in other patients including normals); PRKDC recurrence q = 0.012 became q = 1 after removing flagged sites. Shown only as a caveat, not a finding.",
        "SOX2 amplification dosage in RNA is weak and confounded with cell cycle (S-phase cells carry inflated copy number); treat it as unconfirmed.",
        "TMB-RNA associations are explained by clade membership (the within-clone adjustment removes them) and did not survive robustness checks.",
        "Nine patients: a null state-genotype test is weak evidence of absence. Carriers of a focal amplification mostly sit on one lineage; the within-lineage estimate removes lineage but has fewer informative strata.",
        "BWH69, BWH75, MGH285 and MGH301 are on older trees; BWH75 and MGH301 have no SNV calls and no hotspot site evidence, so SNV-based checks there are not available.",
        "Focal amplification and homozygous deletion calls come from aligned_trees.json (amp = CN >= max(baseline + 3, 2 x baseline)); the broader 'EGFR amp' set in event_associations.json (filtered-event amplification) can differ, e.g. BWH75.",
        "BWH69 SNV calls are noisy (many single-cell-private sites); its subclonal driver hits are less certain.",
    ]
    cohort = {
        "format": "gos-sc-portrait-cohort/1", "created": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "method": f"Verdicts: consistent / inconsistent / not_testable / not_applicable. State links: carriers vs non-carriers inside each patient (informative if >= {MIN_GRP} of each with an RNA state), Mantel-Haenszel across patients, also stratified by lineage; consistent when the 95% CI excludes 1 in the Neftel direction, inconsistent when it excludes 1 the other way, otherwise not_testable (underpowered). Only the four Neftel links are tested; none is selected by p value. Truncal = carrier fraction >= {TRUNCAL_FRAC:.2f} or root clade (aligned_trees.json).",
        "patients": args.patients, "matrix": matrix, "tally": tally, "neftel_links": cohort_links, "novel": novel, "caveats": caveats, "caveat_items": caveat_items(inp, args.patients), "references": REFERENCES,
        "portraits": portraits,
    }
    return cohort, portraits


def cohort_novel(inp, portraits, patients, neftel):
    out = []
    shifts = inp.sub.get("signature_shifts", [])
    # APOBEC gain
    by = defaultdict(list)
    for s in shifts:
        if s["signature"] in ("SBS13", "SBS2") and s["diff"] > 0 and s["q"] < 0.05:
            by[s["patient"]].append(s)
    for p, ss in by.items():
        txt = f"{p}: " + "; ".join(f"clade {s['clade'].split('_')[-1]} {s['signature']} {pct(s['frac'])} of its SNVs vs {pct(s['frac_ref'])} in {s['vs'].replace('MGH305_', '')} (q = {s['q']:.2g})" for s in ss) + ". Canonical GBM carries little APOBEC (SBS2 / SBS13); a gain confined to late clades is unexpected."
        out.append({"id": f"apobec_{p}", "kind": "signature", "status": "candidate", "title": f"{p}: subclonal APOBEC (SBS13) gain", "text": txt, "patients": [p],
                    "caveat": "Fitted on clade-private SNVs (>= 50 sites per clade) with SBS2 / SBS13 added to the patient signature set by design; q is bootstrap-based over sites; small fractions (7-19%).",
                    "links": [{"label": "Subclonal drivers and signature shifts", "scview": "subclonal"}, {"label": f"{p} signatures", "report": p, "tab": 9}]})
    # SBS1 falls
    s1 = [s for s in shifts if s["signature"] == "SBS1"]
    s1d = [s for s in s1 if s["diff"] < 0 and s["q"] < 0.05]
    pts = sorted({s["patient"] for s in s1d})
    jt = [c for c in inp.sigj["contrasts"] if c["signature"] == "SBS1" and c["contrast"] == "subclonal vs truncal" and c["q"] is not None and c["q"] < 0.05 and c["diff"] < 0]
    if s1d:
        out.append({"id": "sbs1_drop", "kind": "signature", "status": "robust", "title": "Clock-like SBS1 drops in subclones",
                    "text": f"SBS1 is lower than in the truncal set in {len(s1d)} of {len(s1)} fitted clades ({len(pts)} patients: {', '.join(pts)}; q < 0.05), and the pooled subclonal-vs-truncal SBS1 fraction is lower in {len({c['patient'] for c in jt})} of {len({c['patient'] for c in inp.sigj['contrasts']})} patients (" + "; ".join(f"{c['patient']} {pct(c['frac_a'])} vs {pct(c['frac_b'])}" for c in jt) + "). Truncal SNVs are enriched for clock-like (age-related) mutations; later SNVs carry other processes.",
                    "patients": sorted({c["patient"] for c in jt}), "caveat": "Compares signature composition, not burden; subclonal sites are called from fewer cells and BWH69 calls are noisy.",
                    "links": [{"label": "Joint signatures (subclonal vs truncal)", "scview": "mutations"}, {"label": "Clade signature shifts", "scview": "subclonal"}]})
    # QC-pass subclonal drivers
    items = []
    for p in patients:
        for h in portraits[p]["genotype"]["subclonal_qc_pass"]:
            if (h["driver_class"] in ("likely", "possible") or h["tier"] in (1, 2)) and h["category"] in ("private", "subclonal"):
                items.append((p, h))
    if items:
        out.append({"id": "subclonal_drivers", "kind": "subclonal_driver", "status": "candidate", "title": "Subclonal driver mutations that pass read-level QC",
                    "text": "; ".join(f"{p} {h['gene']} {h['variant']} ({pct(h['frac'])} of cells)" for p, h in items) + ". The recurrent canonical GBM drivers (TP53, PTEN, PIK3CA, EGFR) are clonal in these patients; these chromatin, DNA-repair and signalling hits are confined to subclones.",
                    "patients": sorted({p for p, _ in items}), "caveat": "Carrier fractions 2-10% from single-cell DNA; every hit has per-site read evidence in the Evidence tab. BWH69 calls are noisier.",
                    "links": [{"label": "Subclonal drivers", "scview": "subclonal"}] + [{"label": f"{p} {h['gene']}", "report": p, "tab": 15, "evsite": h["site"]} for p, h in items[:8]]})
    # PDGFRA amp vs OPC (opposite to Neftel)
    pd = neftel.get("pdgfra_opc")
    if pd and pd["verdict"] == "inconsistent":
        row = next(r for r in pd["test"]["per_patient"] if r["informative"])
        out.append({"id": "pdgfra_not_opc", "kind": "genotype_state", "status": "candidate", "title": f"{row['patient']}: PDGFRA-amplified cells are never OPC-like",
                    "text": f"{row['n_carriers']} PDGFRA-amplified cells with RNA have state counts {row['state_counts_carriers']}, none OPC-like, against {row['state_counts_noncarriers'].get('OPC', 0)} of {row['n_noncarriers']} non-carriers ({pd['evidence']}). Neftel 2019 links PDGFRA amplification to OPC-like cells; here it is absent from them.",
                    "patients": [row["patient"]], "caveat": "One informative patient; OPC-like cells cluster on the tree in this patient (nearest-relative test p < 0.05), so lineage and genotype are partly confounded (within-lineage CI includes 1).",
                    "links": [{"label": "PDGFRA amp: events x state", "scview": "eventassoc", "evassoc": "PDGFRA amp"}, {"label": f"{row['patient']} RNA", "report": row["patient"], "tab": 8}]})
    # EGFR amp vs AC within lineage
    ea = neftel.get("egfr_ac")
    if ea and ea["test"].get("within_lineage") and ea["test"]["within_lineage"].get("lo") is not None and ea["test"]["within_lineage"]["lo"] > 1 and ea["verdict"] != "consistent":
        wl = ea["test"]["within_lineage"]
        out.append({"id": "egfr_ac_lineage", "kind": "genotype_state", "status": "exploratory", "title": "EGFR amplification - AC-like, visible within lineages only",
                    "text": f"Across patients the EGFR amp - AC-like link is not significant ({ea['evidence'].split(';')[1].strip()}), but stratified by lineage the odds ratio is {wl['or']:.2f} (95% CI {wl['lo']:.2f}-{wl['hi']:.2f}, {wl['n_strata']} strata).",
                    "patients": [r["patient"] for r in ea["test"]["per_patient"] if r["informative"]], "caveat": "Secondary (lineage-stratified) analysis; the primary patient-stratified estimate includes 1. Most EGFR-amplified cells sit in near-clonal EGFR patients with few non-carriers.",
                    "links": [{"label": "EGFR amp: events x state", "scview": "eventassoc", "evassoc": "EGFR amp"}]})
    # RTK mosaic exclusivity
    for p in patients:
        c = portraits[p]["checks"]["rtk_mosaic"]
        if c["verdict"] == "consistent":
            out.append({"id": f"rtk_mosaic_{p}", "kind": "genotype", "status": "robust", "title": f"{p}: different RTKs amplified in different cells",
                        "text": c["evidence"] + " Mutually exclusive RTK amplification in separate subclones matches the mosaic RTK amplification described by Snuderl 2011.",
                        "patients": [p], "caveat": "Amplification calls are single-cell CN thresholds (CN >= max(baseline + 3, 2 x baseline)).",
                        "links": [{"label": "Amplicons", "scview": "amplicons"}, {"label": f"{p} CN", "report": p, "tab": 7, "heatmap": "cn"}]})
    # chr10 loss not truncal
    sub10 = [(p, portraits[p]["checks"]["chr10"]) for p in patients if portraits[p]["checks"]["chr10"]["verdict"] == "inconsistent" and portraits[p]["checks"]["chr10"]["numbers"].get("frac")]
    if sub10:
        out.append({"id": "chr10_subclonal", "kind": "genotype", "status": "candidate", "title": "chr10 loss is subclonal in several patients",
                    "text": "chr10 loss is present in only part of the tumour cells in " + "; ".join(f"{p} ({pct(c['numbers']['frac'])})" for p, c in sub10) + ", where chr7 gain is " + ", ".join(f"{p} {pct(portraits[p]['checks']['chr7']['numbers'].get('frac') or 0)}" for p, _ in sub10 if portraits[p]["checks"]["chr7"]["numbers"].get("frac") is not None) + ". The canonical sequence is early chr7 gain plus chr10 loss; here chr10 loss is a later or branch-specific event.",
                    "patients": [p for p, _ in sub10], "caveat": "Arm-level CN calls (+-0.75 vs cell median) are noisy in low-coverage single cells; a cell-level miss lowers the fraction.",
                    "links": [{"label": "Aligned trees (event placement)", "scview": "aligned"}]})
    # PTEN second hit subclonal
    pt = []
    for p in patients:
        n = portraits[p]["checks"]["pten"]["numbers"]
        sn = [x for x in n.get("snvs", []) if x[1] is not None and x[1] >= TRUNCAL_FRAC]
        if sn and n.get("biallelic_frac") is not None and n["biallelic_frac"] < TRUNCAL_FRAC:
            pt.append((p, sn[0], n["biallelic_frac"]))
    if pt:
        out.append({"id": "pten_second_hit", "kind": "genotype", "status": "candidate", "title": "PTEN: clonal mutation, subclonal second hit",
                    "text": "; ".join(f"{p} PTEN {v[0]} in {pct(v[1])} of cells, with loss / LOH of the other allele in the same cell for only {pct(b)}" for p, v, b in pt) + ". Loss of the second PTEN allele follows the clonal point mutation in a subset of cells.",
                    "patients": [p for p, _, _ in pt], "caveat": "LOH from allele-specific single-cell CN; cells lacking an LOH call may be undercalled.",
                    "links": [{"label": "Drivers (two-hit)", "scview": "drivers"}] + [{"label": f"{p} drivers", "report": p, "tab": 7, "heatmap": "drivers"} for p, _, _ in pt]})
    # arm event x state, exploratory
    for r in inp.assoc["event_state"]:
        if r["q"] is not None and r["q"] < 0.05 and r.get("het_p") is not None and r["het_p"] > 0.05 and r["n_patients"] >= 2 and (r["k_up"] == r["n_patients"] or r["k_up"] == 0):
            up = r["or"] > 1
            out.append({"id": f"evst_{r['event']}_{r['state']}", "kind": "event_state", "status": "exploratory", "title": f"{r['label']} - {r['state']}-like",
                        "text": f"Cells with {r['label']} are {'more' if up else 'less'} often {r['state']}-like (MH OR {r['or']:.2f}, 95% CI {r['lo']:.2f}-{r['hi']:.2f}, q = {r['q']:.2g}, same direction in {r['n_patients']}/{r['n_patients']} patients, heterogeneity {pfmt(r['het_p'])}). Not one of the Neftel genotype-state links.",
                        "patients": list(r["per_patient"]), "caveat": "Exploratory: BH over 32 tests, 2-3 patients, carriers cluster on lineages, arm-level CN calls are noisy.",
                        "links": [{"label": "Lineage x state and events", "scview": "associations"}, {"label": "Aligned trees", "scview": "aligned"}]})
    # state follows lineage
    fnd = inp.assoc["findings"]
    cs = next((f for f in fnd if f["kind"] == "clone_state"), None)
    if cs:
        out.append({"id": "state_lineage", "kind": "state_lineage", "status": "robust", "title": "RNA state follows lineage in most patients; MES does not",
                    "text": " ".join([cs["text"]] + [f["text"] for f in fnd if f["kind"] == "heritability_none"] + [f["text"] for f in fnd if f["kind"] == "heritability"][:3]),
                    "patients": cs["patients"], "caveat": "Nine patients; older-tree patients use major clades rather than SNV-anchored lineages.",
                    "links": [{"label": "Lineage x state", "scview": "associations"}]})
    return out


def caveat_items(inp, patients):
    prk = {}
    for st in inp.sub.get("recurrence", []):
        for r in st["rows"]:
            if r["gene"] == "PRKDC" and st["set"] == "trunc":
                prk = r
    sox = {r["outcome"]: r for r in inp.evassoc["results"] if r["event"] == "SOX2 amp"}
    out = []
    if prk:
        out.append({"id": "prkdc", "title": "PRKDC recurrence is an artefact", "text": f"PRKDC truncating / splice hits look recurrent ({prk['hits']} hits in {prk['patients']} patients, q = {prk['q_patients']:.3f}) but every site is flagged by read-level QC (alt on spliced exon-exon reads, alt at read-block ends, same allele in other patients and normals); after removing flagged sites {prk['patients_read_qc_pass']} patients remain (q = {prk['q_patients_read_qc_pass']:.0f}).", "links": [{"label": "Subclonal drivers", "scview": "subclonal"}]})
    c, ss = sox.get("cycling (not G1)"), sox.get("S score")
    if c and ss:
        out.append({"id": "sox2", "title": "SOX2 amplification dosage tracks the cell cycle", "text": f"SOX2-amplified cells are more often cycling (log OR {c['effect']:.2f}, p = {c['p']:.1e}, {c['k']} patients) and have higher S-phase scores (effect {ss['effect']:.2f}, p = {ss['p']:.1e}); S-phase cells carry inflated copy number, so SOX2 RNA dosage is confounded.", "links": [{"label": "SOX2 amp events x state (dosage)", "scview": "eventassoc", "evassoc": "SOX2 amp", "evrank": "dosage"}]})
    tc = inp.__dict__.get("tmb_caveats")
    out.append({"id": "tmb", "title": "TMB-RNA associations are explained by clade", "text": "From the TMB association export: " + (tc or "Total and subclonal TMB mostly encode clade membership; the within-clone adjustment separates TMB from lineage."), "links": [{"label": "TMB associations", "scview": "mutations"}]})
    return out


def write_json(path, obj, dry=False):
    if dry:
        log("dry:", path, os.path.getsize(path) if os.path.exists(path) else "new")
        return
    if os.path.exists(path):
        shutil.copy2(path, f"{path}.bak_{TODAY}")
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f, separators=(",", ":"), allow_nan=False)
    os.replace(tmp, path)
    log("wrote", path, os.path.getsize(path))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=f"{ROOT}/data")
    ap.add_argument("--patients", nargs="*", default=PATIENTS)
    ap.add_argument("--dry", action="store_true")
    args = ap.parse_args()
    cohort, portraits = build(args)
    write_json(f"{args.data}/_cohort/portrait.json", cohort, args.dry)
    for p, d in portraits.items():
        write_json(f"{args.data}/{p}/portrait.json", d, args.dry)


if __name__ == "__main__":
    main()
