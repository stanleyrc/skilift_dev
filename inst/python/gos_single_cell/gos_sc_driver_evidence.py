#!/usr/bin/env python3
"""RNA evidence for the DNA fusions of single-cell patients, for the gOS driver deep-dive.

For every DNA fusion in data/<P>/filtered.events.json (fusion_genes "A::B"), per RNA cell:
  - Arriba reads joining the two genes (split reads + discordant mates), from fusions.tsv
    and, separately, from fusions.discarded.tsv (filtered-out calls still show real chimeric reads);
  - read counts on every exon of both partner genes (one samtools view -L pass over the per-cell STAR BAM),
    so 5' vs 3' exon usage of a partner can be compared between carriers and non-carriers
    (e.g. a breakpoint inside MLLT3 between exons 4 and 5).
Exons follow each gene's Ensembl canonical transcript (else its longest) in the 10x GRCh38-2024-A GTF.

Writes data/<P>/drivers/evidence.json:
  { format, patient, created, transcript_source,
    fusions: [{ label, genes: [{ gene, chrom, strand, transcript, exons: [{ n, start, end }] }],
                cells: { <rna cell_id>: { arriba, discarded, exons: { <gene>: [counts...] } } } }] }

usage: gos_sc_driver_evidence.py PATIENT [PATIENT ...] [--cores N]
"""
import argparse
import datetime
import gzip
import json
import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

DATA = "/gpfs/commons/projects/imielinski_web/sclarke/gOS_GBM_sc/data"
FUS = "/gpfs/commons/groups/imielinski_lab/projects/single_cell_GBM/Flow/rna_fusion_2026_10_08"
GTF = "/gpfs/commons/groups/landau_lab/jzinno/ref/refdata-gex-GRCh38-2024-A/genes/genes.gtf"
SAMTOOLS = "/nfs/sw/easybuild/software/samtools/1.21/bin/samtools"


def attr(s, key):
    m = re.search(key + r' "([^"]+)"', s)
    return m.group(1) if m else None


def load_exons(genes):
    """gene -> {chrom, strand, transcript, exons [(start, end)] in transcript order}."""
    tx = {}  # transcript -> {gene, chrom, strand, exons, canonical}
    opener = gzip.open if GTF.endswith(".gz") else open
    with opener(GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            g = attr(f[8], "gene_name")
            if g not in genes:
                continue
            t = attr(f[8], "transcript_id")
            r = tx.setdefault(t, {"gene": g, "chrom": f[0], "strand": f[6], "exons": [], "canonical": False})
            r["exons"].append((int(f[3]), int(f[4])))
            if 'tag "Ensembl_canonical"' in f[8]:
                r["canonical"] = True
    out = {}
    for t, r in tx.items():
        cur = out.get(r["gene"])
        length = sum(e - s + 1 for s, e in r["exons"])
        key = (r["canonical"], length)
        if cur is None or key > cur["_key"]:
            ex = sorted(r["exons"], reverse=(r["strand"] == "-"))
            out[r["gene"]] = {"chrom": r["chrom"], "strand": r["strand"], "transcript": t, "exons": ex, "_key": key}
    for v in out.values():
        v.pop("_key")
    return out


def arriba_counts(path, pairs):
    """{(a, b) unordered: reads} for the gene pairs from one Arriba table."""
    counts = {}
    if not os.path.exists(path):
        return counts
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 12:
                continue
            g1 = set(f[0].split(","))
            g2 = set(f[1].split(","))
            for a, b in pairs:
                if (a in g1 and b in g2) or (b in g1 and a in g2):
                    try:
                        n = int(f[9]) + int(f[10]) + int(f[11])
                    except ValueError:
                        n = 0
                    counts[(a, b)] = counts.get((a, b), 0) + n
    return counts


CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def exon_counts(bam, bed, regions):
    """Reads overlapping each region (aligned blocks only, a read counted once per region), from one
    `samtools view -L` pass over the BAM (much faster than bedcov's per-region queries for hundreds of exons).
    regions: [(chrom, start1, end1)] in output order."""
    by_chrom = {}
    for k, (c, s, e) in enumerate(regions):
        by_chrom.setdefault(c, []).append((s, e, k))
    for v in by_chrom.values():
        v.sort()
    counts = [0] * len(regions)
    proc = subprocess.Popen([SAMTOOLS, "view", "-F", "0x904", "-L", bed, bam], stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        f = line.split("\t", 7)
        ivs = by_chrom.get(f[2])
        if not ivs:
            continue
        pos = int(f[3])
        hit = set()
        for n, op in CIGAR.findall(f[5]):
            n = int(n)
            if op in "M=X":
                bs, be = pos, pos + n - 1
                for s, e, k in ivs:
                    if s > be:
                        break
                    if e >= bs:
                        hit.add(k)
            if op in "MDN=X":
                pos += n
        for k in hit:
            counts[k] += 1
    if proc.wait() != 0:
        return None
    return counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("patients", nargs="+")
    ap.add_argument("--cores", type=int, default=8)
    args = ap.parse_args()
    for P in args.patients:
        ev_path = os.path.join(DATA, P, "filtered.events.json")
        cells_path = os.path.join(DATA, P, "rna", "cells.json")
        if not (os.path.exists(ev_path) and os.path.exists(cells_path)):
            print(P, ": no filtered.events.json or rna/cells.json, skipped", file=sys.stderr)
            continue
        events = json.load(open(ev_path))
        labels = sorted({e["fusion_genes"] for e in events if e.get("fusion_genes") and "::" in e["fusion_genes"]
                   and not e.get("upstream_only") and "upstream of " not in e["fusion_genes"]})
        pairs = [tuple(l.split("::")[:2]) for l in labels]
        genes = sorted({g for p in pairs for g in p})
        ex = load_exons(set(genes))
        rna_cells = json.load(open(cells_path))["cells"]
        # per-cell folders may spell the gOS ids differently (BWH69_..._10f vs _10F, MGH285_MR_3_pl1 vs MR3_pl1,
        # MGH305_MR4_pl_5 vs MR_4_pl5): match on lowercase letters and digits only
        norm = lambda x: re.sub(r"[^a-z0-9]", "", str(x).lower())
        pdir = os.path.join(FUS, "per_cell", P)
        dir_of = {norm(d): os.path.join(pdir, d) for d in (os.listdir(pdir) if os.path.isdir(pdir) else [])}
        cell_dirs = {}
        for c in rna_cells:
            for k in ("rna_id", "cell_id"):
                d = dir_of.get(norm(c.get(k) or ""))
                if c.get(k) and d:
                    cell_dirs[c.get("cell_id") or c.get("rna_id")] = d
                    break
        # one BED with every exon of every partner gene (bedcov keeps line order)
        bed_rows = [(g, i) for g in genes if g in ex for i in range(len(ex[g]["exons"]))]
        tmp = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False)
        for g, i in bed_rows:
            s, e = ex[g]["exons"][i]
            tmp.write(f"{ex[g]['chrom']}\t{s - 1}\t{e}\n")
        tmp.close()
        regions = [(ex[g]["chrom"], ex[g]["exons"][i][0], ex[g]["exons"][i][1]) for g, i in bed_rows]

        def one(item):
            cid, d = item
            bam = os.path.join(d, "Aligned.sorted.bam")
            cov = exon_counts(bam, tmp.name, regions) if os.path.exists(bam) and bed_rows else None
            kept = arriba_counts(os.path.join(d, "fusions.tsv"), pairs)
            disc = arriba_counts(os.path.join(d, "fusions.discarded.tsv"), pairs)
            return cid, cov, kept, disc

        with ThreadPoolExecutor(args.cores) as pool:
            results = list(pool.map(one, sorted(cell_dirs.items())))
        os.unlink(tmp.name)

        fusions = []
        for label, (a, b) in zip(labels, pairs):
            gl = [g for g in (a, b) if g in ex]
            entry = {
                "label": label,
                "genes": [{"gene": g, "chrom": ex[g]["chrom"], "strand": ex[g]["strand"], "transcript": ex[g]["transcript"],
                           "exons": [{"n": i + 1, "start": s, "end": e} for i, (s, e) in enumerate(ex[g]["exons"])]} for g in gl],
                "cells": {},
            }
            for cid, cov, kept, disc in results:
                rec = {"arriba": kept.get((a, b), 0), "discarded": disc.get((a, b), 0)}
                if cov is not None:
                    rec["exons"] = {g: [cov[k] for k, (gg, _) in enumerate(bed_rows) if gg == g] for g in gl}
                entry["cells"][cid] = rec
            fusions.append(entry)
        out = {"format": "gos-sc-driver-evidence/1", "patient": P,
               "created": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
               "transcript_source": "Ensembl canonical (else longest) transcript, 10x GRCh38-2024-A genes.gtf",
               "fusions": fusions}
        od = os.path.join(DATA, P, "drivers")
        os.makedirs(od, exist_ok=True)
        tmp_out = os.path.join(od, ".evidence.json.tmp")
        with open(tmp_out, "w") as fh:
            json.dump(out, fh, separators=(",", ":"))
        os.replace(tmp_out, os.path.join(od, "evidence.json"))
        os.chmod(os.path.join(od, "evidence.json"), 0o664)
        n_rna = sum(1 for f in fusions for c in f["cells"].values() if c["arriba"] + c["discarded"] > 0)
        print(f"{P}: {len(fusions)} fusions, {len(results)} RNA cells, {n_rna} cell-fusion pairs with chimeric reads", flush=True)


if __name__ == "__main__":
    main()
