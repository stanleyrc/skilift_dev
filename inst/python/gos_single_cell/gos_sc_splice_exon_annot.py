#!/usr/bin/env python3
"""Adds gene exon context to splicing.json files (patient + _cohort) so single-exon / intergenic calls can be filtered.
Per event / finding / cluster (and the cohort's findings / clusters) with a gene symbol:
  n_exons          exons of the gene's Ensembl-canonical transcript (largest transcript when none is tagged)
  intron_overlap   the junction (start-end) overlaps an annotated intron of that gene (any transcript)
Usage: gos_sc_splice_exon_annot.py <data dir with <P>/rna/splicing.json and _cohort/rna/splicing.json> [GTF]
Edits the files in place (back them up first, <file>.bak_<yymmdd>). Idempotent."""
import gzip, json, os, re, sys, glob
from collections import defaultdict

DATA = sys.argv[1]
GTF = sys.argv[2] if len(sys.argv) > 2 else "/gpfs/commons/groups/landau_lab/jzinno/ref/refdata-gex-GRCh38-2024-A/genes/genes.gtf"


def load_gtf(path):
    tx_exons = defaultdict(list)  # transcript -> [(start, end)]
    tx_gene, canonical = {}, set()
    gene_chr = {}
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.split("\t")
            if c[2] not in ("exon", "transcript"):
                continue
            a = c[8]
            tid = re.search(r'transcript_id "([^"]+)"', a).group(1)
            gm = re.search(r'gene_name "([^"]+)"', a) or re.search(r'gene_id "([^"]+)"', a)
            if c[2] == "transcript":
                if 'tag "Ensembl_canonical"' in a:
                    canonical.add(tid)
                continue
            tx_exons[tid].append((int(c[3]), int(c[4])))
            tx_gene[tid] = gm.group(1)
            gene_chr[gm.group(1)] = c[0].replace("chr", "")
    n_exons, introns = {}, defaultdict(list)
    best = {}
    for tid, ex in tx_exons.items():
        g = tx_gene[tid]
        ex.sort()
        for (a, b), (c2, d) in zip(ex[:-1], ex[1:]):
            if b + 1 <= c2 - 1:
                introns[g].append((b + 1, c2 - 1))
        if tid in canonical:
            n_exons[g] = len(ex)
        elif g not in n_exons:
            best[g] = max(best.get(g, 0), len(ex))
    for g, n in best.items():
        n_exons.setdefault(g, n)
    return n_exons, introns


def annotate(row, n_exons, introns, start, end):
    g = row.get("gene")
    if not g or g not in n_exons:
        row["n_exons"] = None
        row["intron_overlap"] = False
        return
    row["n_exons"] = n_exons[g]
    row["intron_overlap"] = bool(start is not None and end is not None and any(s <= end and e >= start for s, e in introns.get(g, ())))


def patch(path, n_exons, introns):
    d = json.load(open(path))
    n = 0
    for e in d.get("events") or []:
        annotate(e, n_exons, introns, e.get("start"), e.get("end")); n += 1
    for f in d.get("findings") or []:
        j = f.get("junction") or {}
        annotate(f, n_exons, introns, j.get("start"), j.get("end")); n += 1
    for c in d.get("clusters") or []:
        # a cluster spans its junctions: overlap with any junction of the cluster
        js = c.get("junctions") or [{"start": c.get("start"), "end": c.get("end")}]
        g = c.get("gene")
        c["n_exons"] = n_exons.get(g)
        c["intron_overlap"] = bool(g in n_exons and any(s <= j["end"] and e >= j["start"] for j in js if j.get("start") is not None for s, e in introns.get(g, ())))
        n += 1
    json.dump(d, open(path, "w"), separators=(",", ":"), allow_nan=False)
    return n


if __name__ == "__main__":
    n_exons, introns = load_gtf(GTF)
    print("genes with exon counts:", len(n_exons))
    for p in sorted(glob.glob(os.path.join(DATA, "*", "rna", "splicing.json"))):
        print(p, patch(p, n_exons, introns), "rows")
