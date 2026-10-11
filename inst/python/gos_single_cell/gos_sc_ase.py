#!/usr/bin/env python3
"""Allele-specific expression (ASE) against allelic copy number in the same single cells (G&T-seq: DNA + RNA).

For every patient with per-cell RNA BAMs and per-cell phased germline hets:
  1. Phased germline het SNVs (Eagle, analysis/<P>/eagle_phasing/hets_with_phase_*) inside exons of
     protein-coding / lncRNA genes (10x GRCh38-2024-A GTF; SNPs in more than one gene dropped).
  2. RNA: REF / ALT reads at those SNPs in each cell's STAR BAM (samtools mpileup, MAPQ >= 20,
     BQ >= 20), turned into haplotype A / B reads with the Eagle phase (A = first haplotype of
     PHASED_GT) and summed per gene and cell.
  3. DNA: the cell's allelic copy number at the gene (major / minor intervals of its balanced allelic
     JaBbA graph, data/<cell>/allelic.json), oriented onto haplotypes A / B with the cell's own phased
     DNA het read counts within 100 kb of the gene (widened to 1 Mb if sparse; inside the CN segment) (the haplotype with more reads carries the larger copy
     number); a cell with too few DNA reads in the segment borrows its clone's pooled orientation.
     Expected RNA fraction of haplotype A in that cell = CN_A / (CN_A + CN_B).
  4. Per gene, over cells with RNA reads: pooled RNA fraction of haplotype A vs the DNA expectation and
     vs 1/2. The test treats cells as units (a weighted mean of per-cell deviations, overdispersed
     binomial floor on its variance) because transcriptional bursting makes single cells monoallelic
     at random; BH across genes. Per clone pooled fractions, and allele switches between clones.
  5. Findings (ranked, plain-language sentences): monoallelic expression with both DNA copies present,
     RNA imbalance beyond DNA, expression from the minority DNA allele, imbalance explained by dosage
     (LOH / allelic gain), dosage buffered in RNA, and the expressed allele switching between clones.

Writes
  data/<P>/rna/ase.json            per-patient genes, clones, per-cell values for the top genes, findings
  data/<P>/rna/ase_reads/<id>.bam  RNA reads at the finding SNPs, per cell (IGV evidence)
  data/_cohort/rna/ase.json        per-patient summary, cohort findings, genes recurrent across patients

usage: gos_sc_ase.py [PATIENT ...] [--cores N] [--pairs-map all_pairs_ase_map.tsv] [--no-slices]
"""
import argparse
import bisect
import datetime
import json
import math
import os
import pickle
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from scipy import stats

DATA = "/gpfs/commons/projects/imielinski_web/sclarke/gOS_GBM_sc/data"
PROJ = "/gpfs/commons/groups/imielinski_lab/projects/single_cell_GBM"
FUS = f"{PROJ}/Flow/rna_fusion_2026_10_08/per_cell"
WORK = f"{PROJ}/Flow/ase_2026_10_10"
GTF = "/gpfs/commons/groups/landau_lab/jzinno/ref/refdata-gex-GRCh38-2024-A/genes/genes.gtf"
FASTA = "/gpfs/commons/groups/landau_lab/jzinno/ref/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
SAMTOOLS = "/nfs/sw/easybuild/software/samtools/1.21/bin/samtools"
CGC = "/gpfs/commons/groups/imielinski_lab/DB/COSMIC/v97_GRCh38/cancer_gene_census.csv"
PATIENTS = ["BWH69", "BWH70", "BWH75", "MGH285", "MGH301", "MGH302", "MGH303", "MGH304", "MGH305"]
GENE_TYPES = {"protein_coding", "lncRNA"}

PARAMS = {
    "min_mapq": 20,
    "min_baseq": 20,
    "rho": 0.2,  # per-cell overdispersion (bursting) for the variance floor and cell weights
    "min_cells": 4,  # cells with RNA reads at the gene
    "min_reads": 20,  # RNA reads at the gene, all cells
    "min_dna_reads_orient": 20,  # DNA het reads near the gene to orient a cell on its own
    "orient_flanks": [100_000, 300_000, 1_000_000],  # bp around the gene (within its CN segment) for that orientation
    "q_max": 0.05,
    "delta_min": 0.2,  # |RNA - DNA| for an RNA-beyond-DNA call
    "mono": 0.1,  # RNA fraction of the weaker allele below this = monoallelic
    "dna_imbalanced": 0.2,  # |DNA expectation - 1/2| for dosage-driven calls
    "switch_hi": 0.75,  # clone switch: one clone >= this, another <= 1 - this
    "switch_min_cells": 4,
    "switch_agree": 0.75,  # share of a clone's cells on its side
    "switch_min_reads": 10,
    "top_genes_cells": 300,  # genes with per-cell values in the JSON
    "slice_findings": 40,  # findings whose SNPs get read slices
}

# Known human imprinted genes (geneimprint / Otago "imprinted" with human evidence): monoallelic by design.
IMPRINTED = set(
    """H19 IGF2 MEG3 MEG8 DLK1 RTL1 PEG3 ZIM2 PEG10 SGCE PEG13 MEST MESTIT1 SNRPN SNURF NDN MAGEL2 MKRN3 NPAP1
    ZDBF2 PLAGL1 HYMAI KCNQ1OT1 CDKN1C PHLDA2 SLC22A18 KCNQ1 GRB10 NNAT BLCAP ZNF597 NAA60 L3MBTL1 DIRAS3
    NAP1L5 INPP5F GNAS GNAS-AS1 KCNK9 PPP1R9A DLX5 TP73 ZNF331 MIMT1 UBE3A IPW SNHG14 PWAR6 NLRP2 IGF2-AS
    ZFAT LIN28B MIR296 MIR298 FAM50B PXDC1 RB1 LRRTM1 GPR1 ZNF264 AIM1 ANO1 CPA4 TFPI2 DCN""".split()
)
FLAG_GENE = re.compile(r"^(HLA-|KIR|IG[HKL][VDJC]|TR[ABDG][VDJC]|MT-|RPL|RPS)")


def log(*a):
    print(datetime.datetime.now().strftime("%H:%M:%S"), *a, flush=True)


def norm(x):
    return re.sub(r"[^a-z0-9]", "", str(x).lower())


def attr(s, key):
    m = re.search(key + r' "([^"]+)"', s)
    return m.group(1) if m else None


# ---------------------------------------------------------------- reference


def load_genes():
    """Exons (per chrom, sorted, merged per gene) and gene spans of protein-coding / lncRNA genes. Cached."""
    cache = os.path.join(WORK, "cache", "gtf_exons.pkl")
    if os.path.exists(cache):
        with open(cache, "rb") as fh:
            return pickle.load(fh)
    exons = defaultdict(list)
    genes = {}
    with open(GTF) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("exon", "gene"):
                continue
            if attr(f[8], "gene_type") not in GENE_TYPES:
                continue
            g = attr(f[8], "gene_name")
            if f[2] == "gene":
                genes[g] = {"chrom": f[0], "start": int(f[3]), "end": int(f[4]), "strand": f[6], "type": attr(f[8], "gene_type")}
            else:
                exons[f[0]].append((int(f[3]), int(f[4]), g))
    merged = {}
    for chrom, ex in exons.items():
        ex.sort()
        out = []
        for s, e, g in ex:
            if out and out[-1][2] == g and s <= out[-1][1] + 1:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e, g])
        merged[chrom] = out
    res = {"exons": merged, "genes": genes}
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    with open(cache, "wb") as fh:
        pickle.dump(res, fh)
    return res


class ExonIndex:
    """Genes whose exons cover a position (exons may overlap between genes)."""

    def __init__(self, exons):
        self.by = {}
        self.union = {}
        for chrom, ex in exons.items():
            starts = [x[0] for x in ex]
            maxlen = max((x[1] - x[0] for x in ex), default=0)
            self.by[chrom] = (starts, ex, maxlen)
            u = []
            for s0, e0, _ in ex:  # ex is sorted by start
                if u and s0 <= u[-1][1] + 1:
                    u[-1][1] = max(u[-1][1], e0)
                else:
                    u.append([s0, e0])
            self.union[chrom] = (np.array([x[0] for x in u]), np.array([x[1] for x in u]))

    def exonic(self, chrom, pos):
        """Vectorised: positions (np array) inside any exon of the chromosome."""
        if chrom not in self.union:
            return np.zeros(len(pos), bool)
        us, ue = self.union[chrom]
        i = np.searchsorted(us, pos, "right") - 1
        ok = i >= 0
        out = np.zeros(len(pos), bool)
        out[ok] = pos[ok] <= ue[i[ok]]
        return out

    def genes_at(self, chrom, pos):
        if chrom not in self.by:
            return set()
        starts, ex, maxlen = self.by[chrom]
        i = bisect.bisect_right(starts, pos) - 1
        out = set()
        while i >= 0 and ex[i][0] >= pos - maxlen:
            if ex[i][0] <= pos <= ex[i][1]:
                out.add(ex[i][2])
            i -= 1
        return out


def cancer_genes():
    try:
        return set(pd.read_csv(CGC, usecols=["Gene Symbol"])["Gene Symbol"].astype(str))
    except Exception:  # noqa: BLE001
        return set()


# ---------------------------------------------------------------- inputs per patient


def read_hets(path, usecols=None):
    cols = usecols or ["seqnames", "start", "Reference_Allele", "Tumor_Seq_Allele1", "ref.count.t", "alt.count.t", "PHASED_GT"]
    df = pd.read_csv(path, sep="\t", usecols=cols, dtype={"seqnames": str, "PHASED_GT": str}, low_memory=False)
    return df


def patient_cells(P, pairs_map):
    """DNA-linked RNA cells of a patient with their BAM, phased-het file, allelic graph and clone."""
    rna_cells_path = os.path.join(DATA, P, "rna", "cells.json")
    if not os.path.exists(rna_cells_path):
        return []
    rc = json.load(open(rna_cells_path))
    rc = rc["cells"] if isinstance(rc, dict) else rc
    pdir = os.path.join(FUS, P)
    dir_of = {norm(d): os.path.join(pdir, d) for d in (os.listdir(pdir) if os.path.isdir(pdir) else [])}
    # phased hets: all_pairs map first (MGH285 files are numbered), else the folder by name
    het_of = {}
    if pairs_map is not None:
        sub = pairs_map[(pairs_map["patient"] == P) & pairs_map["het_pileups_wgs_phase"].notna()]
        for _, r in sub.iterrows():
            for k in ("pair", "cell_id"):
                if isinstance(r[k], str) and r[k]:
                    het_of.setdefault(norm(r[k]), r["het_pileups_wgs_phase"])
    eg = os.path.join(PROJ, "analysis", P, "eagle_phasing")
    hp = sorted(d for d in os.listdir(eg) if d.startswith("hets_with_phase")) if os.path.isdir(eg) else []
    if hp:
        hdir = os.path.join(eg, hp[-1])
        for f in os.listdir(hdir):
            het_of.setdefault(norm(f.replace("_hets_with_phase.txt", "")), os.path.join(hdir, f))
    out = []
    for c in rc:
        cid = c.get("cell_id")
        if not cid:
            continue
        bamdir = dir_of.get(norm(c.get("rna_id"))) or dir_of.get(norm(cid))
        bam = os.path.join(bamdir, "Aligned.sorted.bam") if bamdir else None
        het = het_of.get(norm(cid))
        allelic = os.path.join(DATA, cid, "allelic.json")
        meta = os.path.join(DATA, cid, "metadata.json")
        clone = None
        if os.path.exists(meta):
            try:
                m = json.load(open(meta))
                m = m[0] if isinstance(m, list) else m
                clone = m.get("clone_id")
            except Exception:  # noqa: BLE001
                pass
        out.append(
            {
                "cell_id": cid,
                "rna_id": c.get("rna_id"),
                "clone": clone,
                "bam": bam if bam and os.path.exists(bam) and os.path.exists(bam + ".bai") else None,
                "hets": het if het and os.path.exists(het) else None,
                "allelic": allelic if os.path.exists(allelic) else None,
            }
        )
    return out


def snp_table(cells, eidx):
    """Phased exonic het SNVs of the patient (from one cell's phased-het file; phase is patient-wide)."""
    for c in cells:
        if c["hets"]:
            df = read_hets(c["hets"], ["seqnames", "start", "Reference_Allele", "Tumor_Seq_Allele1", "PHASED_GT"])
            break
    else:
        return None
    df = df[df["PHASED_GT"].isin(["0|1", "1|0"])]
    df = df[(df["Reference_Allele"].str.len() == 1) & (df["Tumor_Seq_Allele1"].str.len() == 1)]
    df = df.assign(chrom="chr" + df["seqnames"].str.replace("^chr", "", regex=True))
    keep = np.zeros(len(df), bool)
    chroms = df["chrom"].values
    posv = df["start"].values.astype(np.int64)
    for ch in pd.unique(chroms):
        m = chroms == ch
        keep[m] = eidx.exonic(ch, posv[m])
    df = df[keep]
    genes = []
    for chrom, pos in zip(df["chrom"].values, df["start"].values):
        g = eidx.genes_at(chrom, int(pos))
        genes.append(next(iter(g)) if len(g) == 1 else None)
    df = df.assign(gene=genes)
    df = df[df["gene"].notna()].copy()
    df["alt_on_a"] = df["PHASED_GT"] == "1|0"
    df = df.rename(columns={"start": "pos", "Reference_Allele": "ref", "Tumor_Seq_Allele1": "alt"})
    return df[["chrom", "pos", "ref", "alt", "alt_on_a", "gene"]].drop_duplicates(["chrom", "pos"]).reset_index(drop=True)


# ---------------------------------------------------------------- RNA pileup

_INDEL = re.compile(r"[+-](\d+)")


def count_bases(s):
    """Upper-cased base counts of an mpileup read-base column ('.'/',' = ref; indels and markers skipped)."""
    out = defaultdict(int)
    i = 0
    n = len(s)
    while i < n:
        ch = s[i]
        if ch == "^":
            i += 2
            continue
        if ch == "$":
            i += 1
            continue
        if ch in "+-":
            m = _INDEL.match(s, i)
            if m:
                i = m.end() + int(m.group(1))
                continue
        if ch in ".,":
            out["REF"] += 1
        elif ch.upper() in "ACGT":
            out[ch.upper()] += 1
        i += 1
    return out


def rna_pileup(bam, sites_file, cfg):
    """{(chrom, pos): (ref_reads, alt_reads)} at the patient's sites, for one cell."""
    cmd = [
        SAMTOOLS, "mpileup", "-B", "-d", "0", "-q", str(cfg["min_mapq"]), "-Q", str(cfg["min_baseq"]),
        "-l", sites_file, "-f", FASTA, "--ff", "UNMAP,SECONDARY,QCFAIL,DUP", bam,
    ]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        log("mpileup failed", bam, p.stderr[-300:])
        return None
    out = {}
    for line in p.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 5 or f[3] == "0":
            continue
        out[(f[0], int(f[1]))] = count_bases(f[4])
    return out


# ---------------------------------------------------------------- DNA allelic CN


def allelic_segments(path):
    """Per chrom: sorted [(start, end, y)] of the 'major' (red) and 'minor' (blue) intervals."""
    d = json.load(open(path))
    red = defaultdict(list)
    blue = defaultdict(list)
    for iv in d.get("intervals", []):
        col = (iv.get("metadata") or {}).get("color", "")
        chrom = "chr" + str(iv["chromosome"]).replace("chr", "")
        rec = (int(iv["startPoint"]), int(iv["endPoint"]), float(iv.get("y") or 0))
        (red if col.upper().startswith("#FF0000") else blue)[chrom].append(rec)
    for dd in (red, blue):
        for k in dd:
            dd[k].sort()
    return red, blue


def seg_at(segs, chrom, pos):
    lst = segs.get(chrom)
    if not lst:
        return None
    i = bisect.bisect_right([s[0] for s in lst], pos) - 1
    if i >= 0 and lst[i][0] <= pos <= lst[i][1]:
        return lst[i]
    return None


def dna_hap_counts(path):
    """Per chrom: (sorted positions, cumulative hapA reads, cumulative hapB reads) from the cell's phased hets."""
    df = read_hets(path, ["seqnames", "start", "ref.count.t", "alt.count.t", "PHASED_GT"])
    df = df[df["PHASED_GT"].isin(["0|1", "1|0"])]
    alt_a = (df["PHASED_GT"] == "1|0").values
    ref = df["ref.count.t"].fillna(0).values.astype(np.int64)
    alt = df["alt.count.t"].fillna(0).values.astype(np.int64)
    a = np.where(alt_a, alt, ref)
    b = np.where(alt_a, ref, alt)
    chroms = ("chr" + df["seqnames"].astype(str).str.replace("^chr", "", regex=True)).values
    pos = df["start"].values.astype(np.int64)
    out = {}
    for ch in pd.unique(chroms):
        m = chroms == ch
        o = np.argsort(pos[m], kind="stable")
        out[ch] = (pos[m][o], np.concatenate([[0], np.cumsum(a[m][o])]), np.concatenate([[0], np.cumsum(b[m][o])]))
    return out


def window_counts(cum, chrom, start, end):
    if chrom not in cum:
        return 0, 0
    p, ca, cb = cum[chrom]
    i = np.searchsorted(p, start, "left")
    j = np.searchsorted(p, end, "right")
    return int(ca[j] - ca[i]), int(cb[j] - cb[i])


# ---------------------------------------------------------------- statistics


def gene_test(h_a, n, e, rho):
    """Cells as units: weighted mean deviation of RNA hapA fraction from e (per cell), t-test with a binomial floor."""
    h_a = np.asarray(h_a, float)
    n = np.asarray(n, float)
    e = np.asarray(e, float)
    k = len(n)
    if k < 2:
        return math.nan, math.nan, math.nan
    w = n / (1 + (n - 1) * rho)
    d = h_a / n - e
    W = w.sum()
    D = float((w * d).sum() / W)
    emp = float((w**2 * (d - D) ** 2).sum() / W**2 * k / (k - 1))
    ee = np.clip(e, 0.05, 0.95)
    floor = float((w**2 * ee * (1 - ee) * (1 + (n - 1) * rho) / n).sum() / W**2)
    se = math.sqrt(max(emp, floor))
    t = D / se if se > 0 else math.nan
    p = float(2 * stats.t.sf(abs(t), df=k - 1)) if math.isfinite(t) else math.nan
    return D, t, p


def bh(p):
    p = np.asarray(p, float)
    q = np.full(len(p), np.nan)
    ok = np.where(np.isfinite(p))[0]
    if not len(ok):
        return q
    o = ok[np.argsort(p[ok])]
    m = len(ok)
    r = p[o] * m / np.arange(1, m + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    q[o] = np.minimum(r, 1)
    return q


def r3(x):
    return None if x is None or not math.isfinite(x) else round(float(x), 3)


def pct(x):
    return f"{round(100 * x)}%"


# ---------------------------------------------------------------- per patient


def run_patient(P, ref, cgc, pairs_map, cores, slices):
    t0 = datetime.datetime.now()
    cells = patient_cells(P, pairs_map)
    usable = [c for c in cells if c["bam"] and c["hets"] and c["allelic"]]
    log(P, f"{len(cells)} DNA-linked RNA cells; {len(usable)} with RNA BAM + phased hets + allelic graph",
        f"(no BAM {sum(not c['bam'] for c in cells)}, no hets {sum(not c['hets'] for c in cells)}, no graph {sum(not c['allelic'] for c in cells)})")
    if len(usable) < PARAMS["min_cells"]:
        return {"patient": P, "status": "skipped", "reason": f"{len(usable)} usable cells (RNA BAM, phased hets and allelic graph)", "n_cells_linked": len(cells)}
    eidx = ExonIndex(ref["exons"])
    snps = snp_table(usable, eidx)
    log(P, f"{len(snps)} phased exonic het SNVs in {snps['gene'].nunique()} genes")
    tmp = tempfile.NamedTemporaryFile("w", suffix=".sites", delete=False, dir=os.path.join(WORK, "cache"))
    for chrom, pos in sorted(zip(snps["chrom"], snps["pos"])):
        tmp.write(f"{chrom}\t{pos}\n")
    tmp.close()
    snp_key = {(c, int(p)): i for i, (c, p) in enumerate(zip(snps["chrom"], snps["pos"]))}
    s_ref = snps["ref"].values
    s_alt = snps["alt"].values
    s_alt_a = snps["alt_on_a"].values
    s_gene = snps["gene"].values
    genes = ref["genes"]

    def one(c):
        pile = rna_pileup(c["bam"], tmp.name, PARAMS)
        if pile is None:
            return c, None, None, None
        rows = []  # (snp index, hapA reads, hapB reads)
        for key, cnt in pile.items():
            i = snp_key.get(key)
            if i is None:
                continue
            r = cnt.get("REF", 0) + cnt.get(s_ref[i], 0)
            a = cnt.get(s_alt[i], 0)
            if r + a == 0:
                continue
            rows.append((i, a, r) if s_alt_a[i] else (i, r, a))
        # DNA allelic CN at the genes this cell has RNA reads for (the per-SNP arrays are dropped here)
        cum = dna_hap_counts(c["hets"])
        red, blue = allelic_segments(c["allelic"])
        dna_g = {}
        for g in {s_gene[i] for i, _, _ in rows}:
            gi = genes.get(g)
            if not gi:
                continue
            mid = (gi["start"] + gi["end"]) // 2
            sr = seg_at(red, gi["chrom"], mid)
            sb = seg_at(blue, gi["chrom"], mid)
            if not sr and not sb:
                continue
            ys = sorted([sr[2] if sr else 0.0, sb[2] if sb else 0.0], reverse=True)
            lo = max(x[0] for x in (sr, sb) if x)
            hi = min(x[1] for x in (sr, sb) if x)
            # orient with DNA hets close to the gene: Eagle switch errors scramble phase over Mb-scale segments
            da = db = 0
            for flank in PARAMS["orient_flanks"]:
                da, db = window_counts(cum, gi["chrom"], max(lo, gi["start"] - flank), min(hi, gi["end"] + flank))
                if da + db >= PARAMS["min_dna_reads_orient"]:
                    break
            dna_g[g] = [ys[0], ys[1], da, db]
        return c, rows, dna_g, None

    with ThreadPoolExecutor(max_workers=cores) as ex:
        res = list(ex.map(one, usable))
    os.unlink(tmp.name)
    res = [r for r in res if r[1] is not None]
    log(P, f"pileups done for {len(res)} cells ({(datetime.datetime.now() - t0).seconds}s)")

    # per gene x cell RNA counts; per SNP pooled counts
    gene_cells = defaultdict(dict)  # gene -> cell index -> [a, b]
    snp_tot = defaultdict(lambda: [0, 0, 0])  # snp -> [a, b, cells]
    snp_cell = defaultdict(dict)  # snp -> cell index -> (a, b)
    cell_meta = []
    for ci, (c, rows, _, _) in enumerate(res):
        reads = 0
        for i, a, b in rows:
            g = s_gene[i]
            gc = gene_cells[g].setdefault(ci, [0, 0])
            gc[0] += a
            gc[1] += b
            st = snp_tot[i]
            st[0] += a
            st[1] += b
            st[2] += 1
            snp_cell[i][ci] = (a, b)
            reads += a + b
        cell_meta.append({"cell_id": c["cell_id"], "rna_id": c["rna_id"], "clone": c["clone"], "snps": len(rows), "reads": reads})

    # DNA: allelic CN at each gene for each cell with RNA reads there (computed per cell above)
    dna = {}  # (gene, ci) -> [cnMaj, cnMin, dnaA, dnaB]
    for ci, (_, _, dna_g, _) in enumerate(res):
        for g, v in dna_g.items():
            dna[(g, ci)] = v

    # clone-pooled orientation for cells with few DNA reads in the segment
    clone_orient = defaultdict(lambda: [0, 0])
    for (g, ci), v in dna.items():
        cl = cell_meta[ci]["clone"]
        clone_orient[(g, cl)][0] += v[2]
        clone_orient[(g, cl)][1] += v[3]

    def expected(g, ci):
        v = dna.get((g, ci))
        if v is None:
            return None
        cmaj, cmin, da, db = v
        if cmaj + cmin <= 0:
            return None
        if abs(cmaj - cmin) < 1e-6:
            return 0.5, cmaj, cmin, "balanced"
        if da + db >= PARAMS["min_dna_reads_orient"]:
            a_major = da >= db
            how = "cell"
        else:
            pa, pb = clone_orient[(g, cell_meta[ci]["clone"])]
            if pa + pb < PARAMS["min_dna_reads_orient"]:
                return None
            a_major = pa >= pb
            how = "clone"
        cn_a, cn_b = (cmaj, cmin) if a_major else (cmin, cmaj)
        return cn_a / (cn_a + cn_b), cn_a, cn_b, how

    # per gene statistics
    rows = []
    for g, cmap in gene_cells.items():
        cis = sorted(cmap)
        n = np.array([cmap[ci][0] + cmap[ci][1] for ci in cis], float)
        if len(cis) < PARAMS["min_cells"] or n.sum() < PARAMS["min_reads"]:
            continue
        ha = np.array([cmap[ci][0] for ci in cis], float)
        ex_ = [expected(g, ci) for ci in cis]
        has = np.array([x is not None for x in ex_])
        e = np.array([x[0] if x is not None else np.nan for x in ex_])
        D_half, t_half, p_half = gene_test(ha, n, np.full(len(n), 0.5), PARAMS["rho"])
        if has.sum() >= PARAMS["min_cells"]:
            D_dna, t_dna, p_dna = gene_test(ha[has], n[has], e[has], PARAMS["rho"])
            exp_pooled = float((n[has] * e[has]).sum() / n[has].sum())
        else:
            D_dna = t_dna = p_dna = math.nan
            exp_pooled = math.nan
        obs = float(ha.sum() / n.sum())
        # per cell: monoallelic in either direction (bursting check)
        frac = ha / n
        n2 = n >= 2
        cells_mono_a = int(((frac >= 0.9) & n2).sum())
        cells_mono_b = int(((frac <= 0.1) & n2).sum())
        cells_both = int(((frac > 0.1) & (frac < 0.9) & n2).sum())
        loh = int(sum(1 for x in ex_ if x is not None and min(x[1], x[2]) <= 0.05))
        # clones
        clones = defaultdict(lambda: {"cells": 0, "a": 0.0, "b": 0.0, "en": 0.0, "n_e": 0.0})
        for j, ci in enumerate(cis):
            cl = cell_meta[ci]["clone"] or "Unassigned"
            r = clones[cl]
            r["cells"] += 1
            r["a"] += ha[j]
            r["b"] += n[j] - ha[j]
            if has[j]:
                r["en"] += n[j] * e[j]
                r["n_e"] += n[j]
        clone_rows = []
        for cl, r in sorted(clones.items()):
            tot = r["a"] + r["b"]
            clone_rows.append({"clone": cl, "cells": r["cells"], "a": int(r["a"]), "b": int(r["b"]), "rna": r3(r["a"] / tot) if tot else None, "dna": r3(r["en"] / r["n_e"]) if r["n_e"] else None})
        gi = genes.get(g, {})
        rows.append({
            "gene": g, "chrom": gi.get("chrom"), "start": gi.get("start"), "end": gi.get("end"),
            "cells": len(cis), "reads": int(n.sum()), "a": int(ha.sum()), "b": int(n.sum() - ha.sum()),
            "rna": obs, "dna": exp_pooled, "delta": D_dna, "t": t_dna, "p": p_dna,
            "imb": D_half, "p_imb": p_half, "cells_dna": int(has.sum()), "loh_cells": loh,
            "mono_a": cells_mono_a, "mono_b": cells_mono_b, "both": cells_both,
            "clones": clone_rows, "_cis": cis, "_ha": ha, "_n": n, "_ex": ex_,
        })
    if not rows:
        return {"patient": P, "status": "skipped", "reason": "no informative genes", "n_cells": len(res)}
    q = bh([r["p"] for r in rows])
    q_imb = bh([r["p_imb"] for r in rows])
    snps_of = defaultdict(list)
    for i, st in snp_tot.items():
        snps_of[s_gene[i]].append(i)
    for r, qq, qi in zip(rows, q, q_imb):
        r["q"] = float(qq)
        r["q_imb"] = float(qi)
        r["snps"] = len(snps_of[r["gene"]])
    log(P, f"{len(rows)} informative genes")

    # ---- clone switches (pre-pass: BH across genes before the findings)
    for r in rows:
        # clone switch: expressed allele differs between clones, judged on cells (bursting makes pooled reads
        # follow one or two high-expression cells): >= switch_min_cells cells with >= 2 reads per clone, most
        # cells of each clone on their side, Mann-Whitney on per-cell shares (BH across genes below)
        sw = None
        per_clone = defaultdict(list)
        for j, ci in enumerate(r["_cis"]):
            if r["_n"][j] >= 2:
                per_clone[cell_meta[ci]["clone"] or "Unassigned"].append(r["_ha"][j] / r["_n"][j])
        cl_rows = {c["clone"]: c for c in r["clones"]}
        cand = []
        for cl, fr in per_clone.items():
            if cl in ("Unassigned",) or len(fr) < PARAMS["switch_min_cells"]:
                continue
            c = cl_rows.get(cl)
            if not c or c["rna"] is None or c["a"] + c["b"] < PARAMS["switch_min_reads"]:
                continue
            fr = np.array(fr)
            if c["rna"] >= PARAMS["switch_hi"] and (fr > 0.5).mean() >= PARAMS["switch_agree"]:
                cand.append(("hi", cl, fr, c))
            elif c["rna"] <= 1 - PARAMS["switch_hi"] and (fr < 0.5).mean() >= PARAMS["switch_agree"]:
                cand.append(("lo", cl, fr, c))
        his = [x for x in cand if x[0] == "hi"]
        los = [x for x in cand if x[0] == "lo"]
        if his and los:
            hi = max(his, key=lambda x: len(x[2]))
            lo = max(los, key=lambda x: len(x[2]))
            p_sw = float(stats.mannwhitneyu(hi[2], lo[2], alternative="two-sided").pvalue)
            dna_sw = hi[3]["dna"] is not None and lo[3]["dna"] is not None and hi[3]["dna"] > 0.6 and lo[3]["dna"] < 0.4
            sw = {"hi": hi[1], "lo": lo[1], "rna_hi": hi[3]["rna"], "rna_lo": lo[3]["rna"], "cells_hi": len(hi[2]), "cells_lo": len(lo[2]),
                  "dna_hi": hi[3]["dna"], "dna_lo": lo[3]["dna"], "p": p_sw, "dna_explained": bool(dna_sw)}
        r["switch"] = sw
    sw_rows = [r for r in rows if r["switch"]]
    for r, qq in zip(sw_rows, bh([r["switch"]["p"] for r in sw_rows])):
        r["switch"]["q"] = float(qq)
        if qq >= PARAMS["q_max"]:
            r["switch"] = None

    # ---- classify + findings
    findings = []
    for r in rows:
        g = r["gene"]
        flags = []
        if g in IMPRINTED:
            flags.append("imprinted")
        if FLAG_GENE.match(g):
            flags.append("mapping-prone")
        if r["snps"] == 1:
            flags.append("single SNP")
        cancer = g in cgc
        rna, dnae, delta, qv = r["rna"], r["dna"], r["delta"], r["q"]
        weak = min(rna, 1 - rna)
        cls = None
        if math.isfinite(dnae) and math.isfinite(qv):
            dna_bal = abs(dnae - 0.5) < PARAMS["dna_imbalanced"]
            if qv < PARAMS["q_max"] and weak <= PARAMS["mono"] and dna_bal:
                cls = "monoallelic"
            elif qv < PARAMS["q_max"] and abs(dnae - 0.5) >= PARAMS["dna_imbalanced"] and (rna - 0.5) * (dnae - 0.5) < 0 and abs(rna - 0.5) >= 0.15:
                cls = "minority-allele"
            elif qv < PARAMS["q_max"] and abs(delta) >= PARAMS["delta_min"] and abs(rna - 0.5) > abs(dnae - 0.5):
                cls = "beyond-dna"
            elif qv < PARAMS["q_max"] and abs(dnae - 0.5) >= PARAMS["dna_imbalanced"] and abs(rna - 0.5) < 0.1:
                cls = "buffered"
            elif abs(dnae - 0.5) >= PARAMS["dna_imbalanced"] and abs(delta) < 0.15 and abs(rna - 0.5) >= PARAMS["dna_imbalanced"] and r["q_imb"] < PARAMS["q_max"]:
                cls = "dosage"
        elif r["q_imb"] < PARAMS["q_max"] and weak <= PARAMS["mono"]:
            cls = "monoallelic-nodna"
        # expressed from an allele the DNA says is absent (LOH in most cells): local phase / orientation conflict
        if cls == "minority-allele" and r["loh_cells"] >= max(1, r["cells_dna"]) / 2:
            cls = "phase-conflict"
            flags.append("phase conflict: RNA from the allele lost in DNA")
        r["class"] = cls
        r["flags"] = flags
        r["cancer"] = cancer
        sw = r["switch"]

        def add(kind, text, sev, score, extra=None):
            if cancer:
                score += 2.0
            if "imprinted" in flags:
                score -= 1.5
            if "mapping-prone" in flags:
                score -= 3.0
            if "single SNP" in flags:
                score -= 1.0
            findings.append({"id": f"{P}:{g}:{kind}", "patient": P, "gene": g, "kind": kind, "severity": sev, "score": round(score, 2), "text": text, "flags": flags + (["cancer gene"] if cancer else []), "chrom": r["chrom"], "start": r["start"], "end": r["end"]})

        sig = -math.log10(max(qv, 1e-12)) if math.isfinite(qv) else -math.log10(max(r["q_imb"], 1e-12))
        nn = f"{r['reads']} reads, {r['cells']} cells"
        if cls == "monoallelic":
            main = "A" if rna > 0.5 else "B"
            sev = "high" if (cancer or r["cells"] >= 10) and "mapping-prone" not in flags else "moderate"
            add("monoallelic", f"{g} is expressed from one allele only (haplotype {main}: {pct(max(rna, 1 - rna))} of {nn}) although the DNA carries both copies (expected {pct(dnae)} A) - monoallelic expression not explained by copy number{' (known imprinted gene)' if 'imprinted' in flags else ''}.", sev, 4 + sig * 0.5 + abs(delta) * 6)
        elif cls == "minority-allele":
            add("minority-allele", f"{g} is expressed mainly from the allele with fewer DNA copies: RNA {pct(rna)} haplotype A vs {pct(dnae)} expected from allelic copy number ({nn}).", "high" if cancer else "moderate", 3.5 + sig * 0.5 + abs(delta) * 6)
        elif cls == "beyond-dna":
            add("beyond-dna", f"{g} allelic imbalance in RNA exceeds the DNA: {pct(rna)} haplotype A in RNA vs {pct(dnae)} expected from allelic copy number ({nn}) - cis-regulatory imbalance.", "moderate", 2 + sig * 0.4 + abs(delta) * 5)
        elif cls == "buffered":
            add("buffered", f"{g} is expressed equally from both alleles ({pct(rna)} A, {nn}) despite allelic copy imbalance (DNA {pct(dnae)} A) - dosage buffered.", "moderate" if cancer else "low", 1.5 + sig * 0.3 + abs(delta) * 3)
        elif cls == "dosage":
            what = f"LOH in {r['loh_cells']}/{r['cells_dna']} cells" if r["loh_cells"] >= r["cells_dna"] / 2 else "allelic copy imbalance"
            add("dosage", f"{g} allelic imbalance follows copy number: RNA {pct(rna)} vs DNA {pct(dnae)} haplotype A ({what}; {nn}).", "moderate" if cancer else "low", 0.8 + sig * 0.15 + (1.5 if cancer else 0))
        elif cls == "monoallelic-nodna":
            add("monoallelic", f"{g} is expressed from one allele only ({pct(max(rna, 1 - rna))} haplotype {'A' if rna > 0.5 else 'B'}, {nn}); allelic copy number not available for these cells.", "moderate", 2.5 + sig * 0.3)
        if sw:
            txt = (f"{g}: the expressed allele switches between clones - {sw['hi']} expresses haplotype A ({pct(sw['rna_hi'])}), {sw['lo']} haplotype B ({pct(1 - sw['rna_lo'])}; {sw['cells_hi']} vs {sw['cells_lo']} cells)"
                   + (" - mirrored by the DNA (allelic copy number switches too)." if sw["dna_explained"] else " - not explained by allelic copy number."))
            add("clone-switch", txt, "high" if not sw["dna_explained"] else "moderate", (5 if not sw["dna_explained"] else 3) + min(6, -math.log10(max(sw["p"], 1e-12)) * 0.4))
    findings.sort(key=lambda f: -f["score"])
    log(P, f"{len(findings)} findings:", ", ".join(f"{f['gene']}({f['kind']})" for f in findings[:12]))

    # ---- per-cell values + SNPs for the top genes
    finding_genes = [f["gene"] for f in findings]
    by_gene = {r["gene"]: r for r in rows}
    ranked = sorted(rows, key=lambda r: (-(r["q"] < 0.05 if math.isfinite(r["q"]) else 0), -(abs(r["delta"]) if math.isfinite(r["delta"]) else abs(r["imb"]))))
    top = list(dict.fromkeys(finding_genes + [g for g in sorted(by_gene, key=lambda g: -by_gene[g]["reads"]) if g in cgc] + [r["gene"] for r in ranked]))[: max(PARAMS["top_genes_cells"], len(set(finding_genes)))]
    top_set = set(top)
    genes_out = []
    for r in rows:
        o = {k: (r3(v) if isinstance(v, float) else v) for k, v in r.items() if not k.startswith("_") and k not in ("p", "p_imb")}
        o["p"] = r["p"] if math.isfinite(r["p"]) else None
        o["q"] = r["q"] if math.isfinite(r["q"]) else None
        o["q_imb"] = r["q_imb"] if math.isfinite(r["q_imb"]) else None
        if r["gene"] in top_set:
            cells_rows = []
            for j, ci in enumerate(r["_cis"]):
                x = r["_ex"][j]
                cells_rows.append([ci, int(r["_ha"][j]), int(r["_n"][j] - r["_ha"][j]), r3(x[0]) if x else None, r3(x[1]) if x else None, r3(x[2]) if x else None])
            o["cell_values"] = cells_rows
            o["snp_list"] = [
                [snps["pos"].iat[i], s_ref[i], s_alt[i], bool(s_alt_a[i]), snp_tot[i][0], snp_tot[i][1], snp_tot[i][2]]
                for i in sorted(snps_of[r["gene"]], key=lambda i: snps["pos"].iat[i])
            ]
        genes_out.append(o)
    genes_out.sort(key=lambda o: (o["q"] if o["q"] is not None else 2, -(o["reads"])))

    # ---- read slices at finding SNPs (per cell)
    ase_reads = {}
    if slices:
        sl_genes = list(dict.fromkeys(f["gene"] for f in findings[: PARAMS["slice_findings"]]))
        regions = []
        for g in sl_genes:
            for i in snps_of[g]:
                regions.append((snps["chrom"].iat[i], int(snps["pos"].iat[i])))
        if regions:
            odir = os.path.join(DATA, P, "rna", "ase_reads")
            os.makedirs(odir, exist_ok=True)
            bed = os.path.join(WORK, "cache", f"{P}_slices.bed")
            with open(bed, "w") as fh:
                for ch, pos in sorted(set(regions)):
                    fh.write(f"{ch}\t{max(0, pos - 76)}\t{pos + 75}\n")
            need = set()
            for g in sl_genes:
                for i in snps_of[g]:
                    need.update(snp_cell[i].keys())

            def cut(ci):
                c = res[ci][0]
                out = os.path.join(odir, f"{c['rna_id']}.bam")
                p1 = subprocess.run([SAMTOOLS, "view", "-b", "-L", bed, "-o", out, c["bam"]], capture_output=True)
                if p1.returncode == 0 and subprocess.run([SAMTOOLS, "index", out], capture_output=True).returncode == 0:
                    return c["rna_id"], f"rna/ase_reads/{c['rna_id']}.bam"
                return c["rna_id"], None

            with ThreadPoolExecutor(max_workers=cores) as ex:
                for rid, path in ex.map(cut, sorted(need)):
                    if path:
                        ase_reads[rid] = path
            log(P, f"{len(ase_reads)} read slices")

    for f in findings:
        r = by_gene[f["gene"]]
        f["cells"] = [res[ci][0]["rna_id"] for ci in r["_cis"]]
        top_snps = sorted(snps_of[f["gene"]], key=lambda i: -(snp_tot[i][0] + snp_tot[i][1]))[:3]
        f["loci"] = [{"chromosome": snps["chrom"].iat[i].replace("chr", ""), "position": int(snps["pos"].iat[i])} for i in top_snps]

    counts = defaultdict(int)
    for r in rows:
        if r["class"]:
            counts[r["class"]] += 1
    summary = {
        "n_cells_linked": len(cells), "n_cells": len(res), "n_snps": int(len(snps)), "n_snps_covered": len(snp_tot),
        "n_genes": len(rows), "n_genes_dna": sum(1 for r in rows if math.isfinite(r["dna"])),
        "n_imbalanced_q05": sum(1 for r in rows if math.isfinite(r["q_imb"]) and r["q_imb"] < 0.05),
        "n_beyond_dna_q05": sum(1 for r in rows if math.isfinite(r["q"]) and r["q"] < 0.05),
        "by_class": dict(counts), "n_switch": sum(1 for r in rows if r["switch"]),
        "median_reads_per_cell": float(np.median([m["reads"] for m in cell_meta])) if cell_meta else 0,
        "median_snps_per_cell": float(np.median([m["snps"] for m in cell_meta])) if cell_meta else 0,
    }
    out = {
        "format": "gos-sc-ase-v1", "patient": P, "created": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "method": "Phased germline het SNVs (Eagle) in exons; RNA haplotype reads per gene and cell (samtools mpileup on the cell's STAR BAM); DNA expectation = haplotype A share of the cell's allelic copy number (balanced allelic JaBbA graph) oriented by the cell's phased DNA het reads; cells-as-units weighted test, BH.",
        "params": PARAMS, "summary": summary,
        "cells": cell_meta, "genes": genes_out, "findings": findings, "ase_reads": ase_reads,
    }
    path = os.path.join(DATA, P, "rna", "ase.json")
    with open(path + ".tmp", "w") as fh:
        json.dump(out, fh, separators=(",", ":"), default=lambda x: x.item() if hasattr(x, "item") else str(x))
    os.replace(path + ".tmp", path)
    log(P, f"wrote {path} ({os.path.getsize(path) // 1024} kB) in {(datetime.datetime.now() - t0).seconds}s")
    return {"patient": P, "status": "ok", "summary": summary, "findings": findings[:15],
            "genes": {r["gene"]: {"rna": r3(r["rna"]), "dna": r3(r["dna"]), "q": r["q"] if math.isfinite(r["q"]) else None, "class": r["class"], "cells": r["cells"]} for r in rows if r["class"] or (r["switch"] is not None)}}


def write_cohort(results):
    recurrent = defaultdict(list)
    for r in results:
        for g, v in (r.get("genes") or {}).items():
            if v["class"]:
                recurrent[g].append({"patient": r["patient"], **v})
    rec = [{"gene": g, "n_patients": len(v), "patients": v} for g, v in recurrent.items() if len(v) >= 2]
    rec.sort(key=lambda x: (-x["n_patients"], x["gene"]))
    findings = sorted([f for r in results for f in (r.get("findings") or [])], key=lambda f: -f["score"])
    out = {
        "format": "gos-sc-ase-cohort-v1", "created": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "patients": [{"patient": r["patient"], "status": r["status"], "reason": r.get("reason"), "summary": r.get("summary")} for r in results],
        "findings": findings[:60], "recurrent": rec[:100],
    }
    os.makedirs(os.path.join(DATA, "_cohort", "rna"), exist_ok=True)
    path = os.path.join(DATA, "_cohort", "rna", "ase.json")
    with open(path + ".tmp", "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    os.replace(path + ".tmp", path)
    log("cohort", f"wrote {path}: {len(findings)} findings, {len(rec)} recurrent genes")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("patients", nargs="*", default=PATIENTS)
    ap.add_argument("--cores", type=int, default=4)
    ap.add_argument("--pairs-map", default=os.path.join(WORK, "all_pairs_ase_map.tsv"))
    ap.add_argument("--no-slices", action="store_true")
    ap.add_argument("--no-cohort", action="store_true")
    a = ap.parse_args()
    os.makedirs(os.path.join(WORK, "cache"), exist_ok=True)
    pairs_map = pd.read_csv(a.pairs_map, sep="\t", dtype=str) if os.path.exists(a.pairs_map) else None
    if pairs_map is None:
        log("warning: no all_pairs map", a.pairs_map, "- phased-het files matched by name only")
    ref = load_genes()
    cgc = cancer_genes()
    log(f"{len(ref['genes'])} genes, {len(cgc)} CGC genes")
    results = []
    for P in a.patients:
        try:
            results.append(run_patient(P, ref, cgc, pairs_map, a.cores, not a.no_slices))
        except Exception as e:  # noqa: BLE001
            import traceback

            traceback.print_exc()
            results.append({"patient": P, "status": "error", "reason": str(e)})
    # keep the other patients' cohort rows when rerunning a subset
    if not a.no_cohort:
        prev = os.path.join(DATA, "_cohort", "rna", "ase.json")
        done = {r["patient"] for r in results}
        if os.path.exists(prev) and set(a.patients) != set(PATIENTS):
            old = json.load(open(prev))
            for p in old.get("patients", []):
                if p["patient"] not in done:
                    pj = os.path.join(DATA, p["patient"], "rna", "ase.json")
                    if os.path.exists(pj):
                        d = json.load(open(pj))
                        results.append({"patient": p["patient"], "status": "ok", "summary": d["summary"], "findings": d["findings"][:15],
                                        "genes": {g["gene"]: {"rna": g["rna"], "dna": g["dna"], "q": g["q"], "class": g["class"], "cells": g["cells"]} for g in d["genes"] if g.get("class") or g.get("switch")}})
                    else:
                        results.append(p)
        results.sort(key=lambda r: r["patient"])
        write_cohort(results)


if __name__ == "__main__":
    sys.exit(main())
