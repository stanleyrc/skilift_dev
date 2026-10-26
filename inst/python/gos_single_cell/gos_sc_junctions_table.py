#!/usr/bin/env python3
"""Junction site table per patient (data/<P>/junctions_table.json) from the called_v3 junction calls.

Input: <v3-dir>/called_v3_<P>.tsv.gz (junction_calling_allcells_261023: merged junctions with per-row tumour / normal counts and
the called_v3 flags) and gene_spans.json (gos-gene-spans/1, hg38). Rows kept: tierlike == TRUE (the tiered backbone: merge-called, >= 3 tumour cells,
normal rule L1; 2k-56k rows per patient). The full called_v3 column (backbone + hotspot lift + single-cell path) has 55k-265k rows
per patient, too many for a browser table (70 MB for MGH285); --scope called_v3 exports it anyway.
No phylogeny and no earlier graph is used.

Output junctions[]: {id, chrA, posA, strA, chrB, posB, strB, class, span, tier (path that called it), n_cells_tumour,
n_reads_tumour, n_cells_normal, n_reads_normal, genes (genes spanning either breakend), splice_mimic}.

usage: gos_sc_junctions_table.py --out <dir> [--v3-dir <dir>] [--gene-spans <json>] [--patients P ...] [--scope backbone|called_v3]
"""
import argparse
import csv
import gzip
import json
import os
import re
import shutil
import time

V3 = "/gpfs/commons/groups/imielinski_lab/projects/single_cell_GBM/analysis/junction_calling_allcells_261023/v3"
GENES = "/gpfs/commons/projects/imielinski_web/sclarke/gOS_GBM_sc/data/_cohort/gene_spans.json"
BIN = 100_000
TRUE = {"TRUE", "True", "true", "T", "1"}


def bare(c):
    return re.sub(r"^chr", "", str(c))


def gene_index(genes):
    idx = {}
    for g in genes:
        for b in range(int(g["start"]) // BIN, int(g["end"]) // BIN + 1):
            idx.setdefault((bare(g["chrom"]), b), []).append(g)
    return idx


def genes_at(idx, chrom, pos):
    hit = [g["gene"] for g in idx.get((bare(chrom), pos // BIN), []) if g["start"] <= pos <= g["end"]]
    return sorted(set(hit))


def num(v, cast=int):
    try:
        return cast(float(v))
    except (TypeError, ValueError, OverflowError):
        return None


def rows(fh, idx, scope="backbone"):
    key = "tierlike" if scope == "backbone" else "called_v3"
    out = []
    for r in csv.DictReader(fh, delimiter="\t"):
        if r.get(key) not in TRUE:
            continue
        pa, pb = num(r["posA"]), num(r["posB"])
        ga, gb = genes_at(idx, r["chrA"], pa), genes_at(idx, r["chrB"], pb)
        out.append({"id": f"{r['chrA']}:{pa}{r['strA']}|{r['chrB']}:{pb}{r['strB']}", "chrA": r["chrA"], "posA": pa, "strA": r["strA"], "chrB": r["chrB"], "posB": pb,
                    "strB": r["strB"], "class": r["class"], "span": num(r["span"]), "tier": r.get("path") or None,
                    "n_cells_tumour": num(r["n_cells_tumour"]), "n_reads_tumour": num(r["n_reads_tumour"]),
                    "n_cells_normal": num(r["n_cells_normal"]), "n_reads_normal": num(r["n_reads_normal"]),
                    "genes": ",".join(sorted(set(ga + gb))), "splice_mimic": r.get("splice_mimic") in TRUE})
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--v3-dir", default=V3)
    ap.add_argument("--gene-spans", default=GENES)
    ap.add_argument("--patients", nargs="*")
    ap.add_argument("--scope", default="backbone", choices=["backbone", "called_v3"])
    a = ap.parse_args(argv)
    idx = gene_index(json.load(open(a.gene_spans))["genes"])
    pts = a.patients or sorted(re.findall(r"called_v3_([A-Z]+[0-9]+)\.tsv\.gz", " ".join(os.listdir(a.v3_dir))))
    for p in pts:
        f = os.path.join(a.v3_dir, f"called_v3_{p}.tsv.gz")
        if not os.path.isfile(f):
            print(p, "no", f)
            continue
        with gzip.open(f, "rt") as fh:
            js = rows(fh, idx, a.scope)
        doc = {"format": "gos-sc-junctions-table/1", "patient": p, "genome": "hg38", "created": time.strftime("%Y-%m-%dT%H:%M:%S"), "source": os.path.basename(f),
               "scope": "called_v3 tiered backbone (tierlike)" if a.scope == "backbone" else "called_v3 (all)", "n_rows": len(js), "junctions": js}
        d = os.path.join(a.out, p)
        os.makedirs(d, exist_ok=True)
        dest = os.path.join(d, "junctions_table.json")
        if os.path.exists(dest):
            bak = f"{dest}.bak_{time.strftime('%y%m%d')}"
            if not os.path.exists(bak):
                shutil.copy2(dest, bak)
        with open(dest + ".tmp", "w") as fh:
            json.dump(doc, fh, separators=(",", ":"))
        os.replace(dest + ".tmp", dest)
        print(p, len(js), "junctions ->", dest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
