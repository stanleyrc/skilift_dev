#!/usr/bin/env python3
"""Per-cell read pileups at the candidate driver-gene SNVs of gos_sc_allsnv (read evidence for sites that are not tree sites).

Usage (mskilab_ne1 python, has pysam):
  MODE=own    PATIENT=MGH305 python3 gos_sc_allsnv_pileup.py   # the patient's candidate sites in the patient's own cells
  MODE=cohort PATIENT=MGH305 python3 gos_sc_allsnv_pileup.py   # the shortlisted sites (WORK/shortlist.tsv) in this patient's cells: cross-patient background
  CORES=8 WORK=... DATA_DIR=... PRECOMPUTE=.../gos_datasets/gOS_GBM_sc/precompute

Reads every read over each site from the cell's DNA BAM (RG-free: one BAM per cell, inputs.tsv `dna_bam`; no read filter), exactly the pileup
of precompute/steps/site_evidence.py (pile_site's base_at / flags), plus a `spliced` flag (N in the CIGAR, a >= 15 bp clip or an SA tag: RNA / cDNA carry-over or a
processed pseudogene gives alt reads on spliced / split alignments). Output WORK/pileup_<MODE>_<PATIENT>.json:
  own:    sites { id: { cells: {alt,dp} caller-filtered per cell (region_calls.json cell order), presets (site_evidence.py summarise: all / caller / strict),
                        splice: {alt, alt_spl, ref, ref_spl} tumour reads, n_reads } }
  cohort: sites { id: { alt_tumor, dp_tumor, cells_alt, cells_cov, alt_normal, dp_normal, n_mut } } (caller filter, region_calls.json-style counts)
"""
import csv
import json
import multiprocessing as mp
import os
import sys
import time

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "precompute", "steps"))
import site_evidence as se  # noqa: E402

DATA = os.environ.get("DATA_DIR", "/gpfs/commons/projects/imielinski_web/sclarke/gOS_GBM_sc/data")
WORK = os.environ.get("WORK", "/gpfs/commons/groups/imielinski_lab/projects/single_cell_GBM/gos_datasets/gOS_GBM_sc/work/allsnv")
PRE = os.environ.get("PRECOMPUTE", "/gpfs/commons/groups/imielinski_lab/projects/single_cell_GBM/gos_datasets/gOS_GBM_sc/precompute")


CLIP_MIN = 15


def pile(bam, s):
    """Read dicts over one site (se.pile_site's fields + spl)."""
    chrom, pos = s["chrom"], int(s["pos"])
    pos0 = pos - 1
    alt, ref = s["alt"].upper(), s["ref"].upper()
    out = []
    try:
        it = bam.fetch(chrom, pos0, pos)
        for a in it:
            if a.is_unmapped or a.cigartuples is None:
                continue
            qp = se.base_at(a, pos0)
            if qp == -1:
                continue
            # spliced (N) or split / clipped alignment: cDNA carry-over (or a processed pseudogene) spans exon junctions, which BWA aligns as a split read
            spl = any(op == 3 or (op in (4, 5) and ln >= CLIP_MIN) for op, ln in a.cigartuples) or a.has_tag("SA")
            if qp is None:
                al, bq, end = "D", 0, 0
            else:
                b = a.query_sequence[qp].upper() if a.query_sequence else "N"
                al = "A" if b == alt else "R" if b == ref else "O"
                q = a.query_qualities
                bq = q[qp] if q is not None else 0
                end = min(qp - a.query_alignment_start, a.query_alignment_end - 1 - qp)
            out.append(dict(cell=0, al=al, bq=int(bq), mq=int(a.mapping_quality), rev=a.is_reverse, dup=a.is_duplicate, sec=a.is_secondary,
                            supp=a.is_supplementary, qcf=a.is_qcfail, end=int(end), spl=spl))
    except (ValueError, OSError):
        return None      # contig not in this BAM / unreadable
    return out


def work_cell(args):
    i, path, sites = args
    if not path or not os.path.exists(path):
        return i, None
    try:
        bam = pysam.AlignmentFile(path)
    except Exception:
        return i, None
    res = []
    for s in sites:
        r = pile(bam, s)
        if r is None:
            r = []
        for x in r:
            x["cell"] = i
        res.append(r)
    bam.close()
    return i, res


def load_inputs(patient):
    rc = json.load(open(os.path.join(DATA, patient, "precompute", "region_calls.json")))
    bam = {}
    with open(os.path.join(PRE, patient, "inputs.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            bam[r["cell_id"]] = r["dna_bam"]
    cells = rc["cells"]
    normal = {i for i, c in enumerate(rc.get("clone") or []) if str(c).lower() == "normal"}
    return cells, normal, [bam.get(c) for c in cells]


def run(patient, sites, mode):
    cells, normal, paths = load_inputs(patient)
    n = len(cells)
    t0 = time.time()
    reads = [[] for _ in sites]
    covered_cells = 0
    with mp.Pool(int(os.environ.get("CORES", "8"))) as pool:
        for k, (i, res) in enumerate(pool.imap_unordered(work_cell, [(i, paths[i], sites) for i in range(n)])):
            if res is None:
                continue
            covered_cells += 1
            for j, r in enumerate(res):
                reads[j].extend(r)
            if k % 20 == 0:
                print(patient, mode, "cell", k, "of", n, round(time.time() - t0), "s", flush=True)
    out = {}
    for j, s in enumerate(sites):
        rs = reads[j]
        caller = se.PRESETS["caller"]
        if mode == "own":
            alt = [0] * n
            dp = [0] * n
            for r in rs:
                if r["cell"] is None or not se.passes(r, caller):
                    continue
                if r["al"] == "A":
                    alt[r["cell"]] += 1
                    dp[r["cell"]] += 1
                elif r["al"] == "R":
                    dp[r["cell"]] += 1
            tum = [r for r in rs if r["cell"] not in normal and se.passes(r, caller)]
            sp = dict(alt=sum(r["al"] == "A" for r in tum), alt_spl=sum(r["al"] == "A" and r["spl"] for r in tum),
                      ref=sum(r["al"] == "R" for r in tum), ref_spl=sum(r["al"] == "R" and r["spl"] for r in tum))
            out[s["id"]] = dict(cells=dict(alt=alt, dp=dp), splice=sp, n_reads=len(rs),
                                presets={k: se.summarise(rs, normal, p) for k, p in se.PRESETS.items()})
        else:
            pr = se.summarise(rs, normal, caller)
            t, nm = pr["tumor"], pr["normal"]
            cell_alt = {}
            for r in rs:
                if r["al"] == "A" and se.passes(r, caller) and r["cell"] not in normal:
                    cell_alt[r["cell"]] = cell_alt.get(r["cell"], 0) + 1
            out[s["id"]] = dict(alt_tumor=t["alt_f"] + t["alt_r"], dp_tumor=t["alt_f"] + t["alt_r"] + t["ref_f"] + t["ref_r"], cells_alt=t["cells_alt"],
                                cells_cov=t["cells_cov"], alt_normal=nm["alt_f"] + nm["alt_r"], dp_normal=nm["alt_f"] + nm["alt_r"] + nm["ref_f"] + nm["ref_r"],
                                n_mut=len(cell_alt))
    return dict(format="gos-sc-allsnv-pileup/1", mode=mode, patient=patient, created=time.strftime("%Y-%m-%dT%H:%M:%S"), cells_with_bam=covered_cells,
                n_cells=n, sites=out)


def read_sites(path):
    out = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out.append(dict(id=r["id"], chrom=r["chrom"], pos=int(r["pos"]), ref=r["ref"], alt=r["alt"]))
    return out


def main():
    patient, mode = os.environ["PATIENT"], os.environ.get("MODE", "own")
    sites = read_sites(os.environ.get("SITES_FILE") or os.path.join(WORK, "pileup_sites_%s.tsv" % patient if mode == "own" else "shortlist.tsv"))
    if mode == "cohort":
        own = {s["id"] for s in read_sites(os.path.join(WORK, "pileup_sites_%s.tsv" % patient))} if os.path.exists(os.path.join(WORK, "pileup_sites_%s.tsv" % patient)) else set()
    print(patient, mode, len(sites), "sites", flush=True)
    doc = run(patient, sites, mode)
    path = os.path.join(WORK, "pileup_%s_%s.json" % (os.environ.get("OUT_MODE", mode), patient))
    with open(path + ".tmp", "w") as fh:
        json.dump(doc, fh, separators=(",", ":"))
    os.replace(path + ".tmp", path)
    print("wrote", path, flush=True)


if __name__ == "__main__":
    main()
