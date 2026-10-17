#!/usr/bin/env python3
"""Where each splice alteration comes from: DNA (a splice-site SNV in filtered.events.json), RNA (aberrant splicing in
rna/splicing.json) or both, per patient and gene. Read by the gOS cohort "Altered genes across patients" panel.

DNA: filtered.events rows with type Splice (splice-site SNVs from the SNV calls, rescue or the pooled driver scan).
RNA: splicing.json findings (known variants such as EGFRvIII, patient-specific and clone-differential junctions) and
unannotated junction events used in >= MIN_RNA_CELLS cells.
Both: an RNA finding / event in the same gene with a junction end within MAX_DIST bp of the DNA splice-site SNV.
RNA-only entries are kept only for findings (the curated top junctions per patient), with their reads and cells.

usage: gos_sc_splice_support.py [--data <gOS data dir>] [--dry-run]
writes <data>/_cohort/splice_support.json (backs up an existing file as .bak_<yymmdd>)
"""
import argparse
import json
import os
import re
import shutil
import time

DATA = "/gpfs/commons/projects/imielinski_web/sclarke/gOS_GBM_sc/data"
MAX_DIST = 1000
MIN_RNA_CELLS = 2


def dna_splice(events):
    out = []
    for e in events:
        if "splice" not in str(e.get("type", "")).lower() or not e.get("gene"):
            continue
        try:
            pos = int(str(e.get("Genome_Location", "")).split(":")[1].split("-")[0])
        except (IndexError, ValueError):
            pos = None
        out.append({"gene": e["gene"], "variant": e.get("Variant") or "", "chrom": str(e.get("seqnames") or "").replace("chr", ""), "pos": pos,
                    "n_cells": e.get("n_cells"), "source": e.get("source") or "SNV calls"})
    return out


def junction_ends(j):
    return [x for x in ((j or {}).get("start"), (j or {}).get("end")) if x is not None]


def rna_splice(doc):
    """RNA findings (kept for display) and unannotated events (used only to confirm a DNA site)."""
    finds, evs = [], []
    for f in doc.get("findings") or []:
        n = (f.get("numbers") or {})
        j = f.get("junction") or {}
        finds.append({"gene": f.get("gene"), "kind": f.get("kind"), "severity": f.get("severity"), "text": f.get("text"),
                      "chrom": str(j.get("chromosome") or "").replace("chr", ""), "ends": junction_ends(j),
                      "n_cells": n.get("n_cells") if n.get("n_cells") is not None else len(f.get("cells") or []), "reads": n.get("reads"), "total": n.get("total")})
    for e in doc.get("events") or []:
        if e.get("annotated") or (e.get("cells_using") or 0) < MIN_RNA_CELLS:
            continue
        evs.append({"gene": e.get("gene"), "chrom": str(e.get("chromosome") or "").replace("chr", ""), "ends": [e.get("start"), e.get("end")],
                    "n_cells": e.get("cells_using"), "reads": e.get("reads"), "text": e.get("summary")})
    return finds, evs


def near(site, rna):
    return site["pos"] is not None and rna["gene"] == site["gene"] and (not site["chrom"] or not rna["chrom"] or rna["chrom"] == site["chrom"]) \
        and any(x is not None and abs(int(x) - site["pos"]) <= MAX_DIST for x in rna["ends"])


def patient_support(events, splicing):
    """gene -> {support: DNA | RNA | DNA+RNA, dna: [...], rna: [...]}"""
    dna = dna_splice(events)
    finds, evs = rna_splice(splicing) if splicing else ([], [])
    out = {}
    for s in dna:
        hits = [r for r in finds + evs if near(s, r)]
        g = out.setdefault(s["gene"], {"dna": [], "rna": []})
        g["dna"].append(s)
        for r in hits:
            if r not in g["rna"]:
                g["rna"].append(r)
    for f in finds:
        if not f["gene"]:
            continue
        g = out.setdefault(f["gene"], {"dna": [], "rna": []})
        if f not in g["rna"]:
            g["rna"].append(f)
    for g in out.values():
        confirmed = any(near(s, r) for s in g["dna"] for r in g["rna"])
        g["support"] = "DNA+RNA" if g["dna"] and confirmed else ("DNA" if g["dna"] else "RNA")
        for r in g["rna"]:
            r.pop("ends", None)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=DATA)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)
    doc = {"format": "gos-sc-splice-support/1", "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "params": {"max_dist_bp": MAX_DIST, "min_rna_cells": MIN_RNA_CELLS},
           "method": __doc__.split("usage:")[0].strip(), "patients": {}, "rna_available": {}}
    for p in sorted(os.listdir(a.data)):
        ev_f = os.path.join(a.data, p, "filtered.events.json")
        if not re.fullmatch(r"[A-Z]+[0-9]+", p) or not os.path.isfile(ev_f):  # patient dirs only (not cells / pools / _cohort)
            continue
        sp_f = os.path.join(a.data, p, "rna", "splicing.json")
        splicing = json.load(open(sp_f)) if os.path.isfile(sp_f) else None
        doc["rna_available"][p] = splicing is not None
        doc["patients"][p] = patient_support(json.load(open(ev_f)), splicing)
        tally = {}
        for g in doc["patients"][p].values():
            tally[g["support"]] = tally.get(g["support"], 0) + 1
        print(p, tally, "rna" if splicing else "no splicing.json")
    out = os.path.join(a.data, "_cohort", "splice_support.json")
    if a.dry_run:
        return 0
    if os.path.exists(out):
        shutil.copy2(out, f"{out}.bak_{time.strftime('%y%m%d')}")
    with open(out + ".tmp", "w") as fh:
        json.dump(doc, fh, separators=(",", ":"))
    os.replace(out + ".tmp", out)
    print("wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
