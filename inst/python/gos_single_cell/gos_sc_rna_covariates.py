#!/usr/bin/env python3
"""Extra per-cell covariates for the cross-patient RNA differential (data/<P>/rna/covariates.json).

The Seurat object only carries state / phase / region / QC columns. This adds, for every RNA cell that has a
matched DNA cell, variables the cohort RNA view can group or stratify by:

  Clade           smallest called clade (data/<P>/clades.json, >= 5 cells, not the root) holding the cell
  DNA clone       clone_id of datafiles.json
  Ploidy          numeric + "Ploidy bin"
  Plate           sequencing plate of the cell id
  Driver: <ev>    carrier / non-carrier of each recurrent event (events of _cohort/event_associations.json
                  with >= 3 carriers in the patient; carriers = filtered.events.json cell_ids)
  CN: <gene>      total copy number at the gene midpoint (data/<cell>/complex.json), numeric, plus
                  "CN level: <gene>" bins (0 / 1 / 2 / 3-4 / 5-9 / 10+); genes: GBM drivers and amp / homdel
                  genes of the recurrent events, only when >= 10% of cells are off the modal copy number
  ecDNA copies: <walk>  copies of each curated circular amplicon (0 when the DNA cell lacks it), numeric + bins

File format: {"format": "gos-sc-rna-covariates/1", "fields": [{"name", "group", "numeric"}], "cells": {rna_id: {field: value}}}

usage: gos_sc_rna_covariates.py [--data DIR] [--patients P ...]
"""
import argparse
import datetime
import gzip
import importlib.util
import json
import math
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("evassoc", os.path.join(HERE, "gos_sc_cohort_event_assoc.py"))
ev = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ev)

ROOT = ev.ROOT
GTF = "/gpfs/commons/groups/imielinski_lab/projects/single_cell_GBM/db/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz"
DATAFILES = os.path.join(os.path.dirname(ROOT.rstrip("/")), "datafiles.json")
DATAFILES = "/gpfs/commons/projects/imielinski_web/sclarke/gOS_GBM_sc/datafiles.json"


def cn_bin(v):
    if v is None:
        return None
    v = round(v)
    return "0" if v <= 0 else "1" if v == 1 else "2" if v == 2 else "3-4" if v <= 4 else "5-9" if v <= 9 else "10+"


def copy_bin(v):
    return "0" if v <= 0 else "1-9" if v < 10 else "10-49" if v < 50 else "50+"


def plate_of(cell_id):
    parts = str(cell_id).split("_")
    return "_".join(parts[1:-1]) if len(parts) > 2 else ""


def gene_midpoints(wanted):
    out = {}
    with gzip.open(GTF, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 9 or c[2] != "gene":
                continue
            attrs = c[8]
            i = attrs.find('gene_name "')
            if i < 0:
                continue
            name = attrs[i + 11: attrs.find('"', i + 11)]
            if name in wanted and name not in out:
                out[name] = (c[0].replace("chr", ""), (int(c[3]) + int(c[4])) // 2)
    return out


def cn_at(intervals, chrom, pos):
    for i in intervals.get(chrom, ()):
        if i[0] <= pos <= i[1]:
            return i[2]
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=ROOT)
    ap.add_argument("--patients", nargs="*", default=ev.PATIENTS)
    a = ap.parse_args()
    assoc = json.load(open(os.path.join(a.data, "_cohort", "event_associations.json")))
    recurrent = {e["event"] for e in assoc["events"]}
    cn_genes = set(ev.GBM_GENES) | {e["gene"] for e in assoc["events"] if e["cls"] in ("amp", "homdel") and e["gene"]}
    coords = gene_midpoints(cn_genes)
    print(f"{len(cn_genes)} CN genes, {len(coords)} with coordinates", flush=True)
    meta = {d["pair"]: d for d in json.load(open(DATAFILES)) if d.get("entry_type") == "cell"}
    for p in a.patients:
        d = os.path.join(a.data, p)
        try:
            rna = json.load(open(os.path.join(d, "rna", "cells.json")))["cells"]
            clades = json.load(open(os.path.join(d, "clades.json")))["clades"]
            events = json.load(open(os.path.join(d, "filtered.events.json")))
        except FileNotFoundError as err:
            print(f"[{p}] skipped: {err}")
            continue
        try:
            walks = json.load(open(os.path.join(d, "walks.json")))["walks"]
        except FileNotFoundError:
            walks = []
        matched = [c for c in rna if c.get("cell_id")]
        dna = [c["cell_id"] for c in matched]
        fields, cols = {}, defaultdict(dict)

        def put(cell, name, value, group, numeric=False):
            if value is None or (isinstance(value, float) and not math.isfinite(value)):
                return
            fields.setdefault(name, {"name": name, "group": group, "numeric": numeric})
            cols[cell["rna_id"]][name] = value

        # clades: smallest clade (>= 5 cells, not a root) containing the cell
        usable = sorted([c for c in clades if c.get("parent") not in (None, "") and c["n_cells"] >= 5], key=lambda c: c["n_cells"])
        of_cell = {}
        for c in usable:
            for x in c["cells"]:
                of_cell.setdefault(x, c["id"])
        for c in matched:
            put(c, "Clade", of_cell.get(c["cell_id"]), "Lineage")
            m = meta.get(c["cell_id"])
            if m:
                put(c, "DNA clone", m.get("clone_id"), "Lineage")
                pl = m.get("ploidy")
                if pl is not None:
                    put(c, "Ploidy", float(pl), "Genome", True)
                    put(c, "Ploidy bin", "<= 2.5" if pl <= 2.5 else "2.5 - 3.5" if pl <= 3.5 else "> 3.5", "Genome")
            put(c, "Plate", plate_of(c["cell_id"]), "Sample")
        # driver carrier status
        pos = set(dna)
        carriers = defaultdict(set)
        for e in events:
            cls = ev.event_class(e)
            if not cls or not ev.tier_ok(e, cls) or str(e.get("upstream_only")).lower() == "true":
                continue
            eid, label, gene = ev.event_id(e, cls)
            if eid in recurrent:
                carriers[eid] |= {x for x in str(e.get("cell_ids", "")).split(",") if x}
        for eid, cs in carriers.items():
            if len(cs & pos) < 3 or len(pos - cs) < 3:
                continue
            for c in matched:
                put(c, f"Driver: {eid}", "carrier" if c["cell_id"] in cs else "non-carrier", "Driver event")
        # gene copy number
        ints = {}
        for cid in dna:
            f = os.path.join(a.data, cid, "complex.json")
            if os.path.exists(f):
                by = defaultdict(list)
                for i in json.load(open(f)).get("intervals", []):
                    by[str(i["chromosome"]).replace("chr", "")].append((i["startPoint"], i["endPoint"], i["y"]))
                ints[cid] = by
        for g, (chrom, mid) in coords.items():
            vals = {cid: cn_at(iv, chrom, mid) for cid, iv in ints.items()}
            v = [x for x in vals.values() if x is not None]
            if len(v) < 10:
                continue
            modal = Counter(round(x) for x in v).most_common(1)[0][0]
            if sum(1 for x in v if round(x) != modal) < 0.1 * len(v):
                continue
            for c in matched:
                x = vals.get(c["cell_id"])
                put(c, f"CN: {g}", float(x) if x is not None else None, "Copy number", True)
                put(c, f"CN level: {g}", cn_bin(x), "Copy number")
        # ecDNA copies
        for w in walks:
            if not w.get("circular") or not w.get("curated") or w.get("ncells", 0) < 3:
                continue
            name = w.get("label") or w.get("name")
            wc = w.get("cells", {})
            for c in matched:
                if c["cell_id"] in ints or c["cell_id"] in meta:
                    n = float(wc.get(c["cell_id"], 0))
                    put(c, f"ecDNA copies: {name}", n, "ecDNA", True)
                    put(c, f"ecDNA level: {name}", copy_bin(n), "ecDNA")
        out = {"format": "gos-sc-rna-covariates/1", "created": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), "patient": p,
               "fields": sorted(fields.values(), key=lambda f: (f["group"], f["name"])), "cells": cols}
        path = os.path.join(d, "rna", "covariates.json")
        if os.path.exists(path) and not os.path.exists(path + ".bak_261014"):
            os.replace(path, path + ".bak_261014")
        with open(path, "w") as f:
            json.dump(out, f, allow_nan=False)
        print(f"[{p}] {len(matched)} matched cells, {len(fields)} fields -> {path}", flush=True)


if __name__ == "__main__":
    main()
