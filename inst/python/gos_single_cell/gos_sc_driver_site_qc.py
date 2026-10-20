#!/usr/bin/env python3
"""Site-level QC for subclonal driver hits -> driver_site_qc.json {"<patient>|<site>": {...}} (needs pysam + minimap2).

  paralog   ref flank (+-FLANK bp) aligned back to the genome with minimap2 (up to 20 secondary loci); the base at the homologous
            position of each other locus with identity >= 0.85 / coverage >= 0.8; `paralog` when it equals the hit's alt allele.
  segdup    site inside genomicSuperDups (hg38 BED).
  umap_k100 site inside Umap k100 uniquely-mappable regions (False = low mappability).
  mapq      mean MAPQ of alt vs ref reads at the site in the patient's carrier cells (CRAMs read-only; `mapq` subcommand).

Usage:
  gos_sc_driver_site_qc.py sites --hits subclonal_drivers.json --ref hg38.fa --segdup bed --umap bed.gz --minimap2 minimap2 --out out.json
  gos_sc_driver_site_qc.py mapq --hits subclonal_drivers.json --crams cell_cram_table.tsv --ref hg38.fa --into out.json [--max-cells 60]
Patient data stay local; nothing is sent anywhere.
"""
import argparse, csv, gzip, json, os, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import driver_filters as F

FLANK = 50


def load_hits(path):
    d = json.load(open(path))
    return d["hits"] if isinstance(d, dict) else d


def read_bed(path, usecols=3):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for ln in fh:
            if ln.startswith(("track", "#", "browser")):
                continue
            t = ln.split("\t")
            if len(t) >= 3:
                yield t[0], int(t[1]), int(t[2])


def parse_sam(text, flank_idx):
    """-> {qname: [records]}; record = dict(chrom,pos0,cigar,flag,nm,qlen,seq)"""
    out = {}
    for ln in text.splitlines():
        if not ln or ln[0] == "@":
            continue
        t = ln.split("\t")
        if int(t[1]) & 4:
            continue
        nm = next((int(x[5:]) for x in t[11:] if x.startswith("NM:i:")), 0)
        cig = F.parse_cigar(t[5])
        rev = bool(int(t[1]) & 16)
        qlen = sum(n for op, n in cig if op in "MIS=X")
        qidx = qlen - 1 - flank_idx if rev else flank_idx
        out.setdefault(t[0], []).append(dict(chrom=t[2], pos0=int(t[3]) - 1, cigar=cig, nm=nm, qlen=qlen, qidx=qidx, rev=rev, mapq=int(t[4])))
    return out


def comp(b):
    return {"A": "T", "C": "G", "G": "C", "T": "A"}.get(b.upper(), "N")


def sites_cmd(a):
    import pysam
    hits = load_hits(a.hits)
    ref = pysam.FastaFile(a.ref)
    seg = F.IntervalIndex(read_bed(a.segdup)) if a.segdup else None
    umap = F.IntervalIndex(read_bed(a.umap)) if a.umap else None
    keys, fa = {}, []
    for h in hits:
        k = f"{h['patient']}|{h['site']}"
        if k in keys:
            continue
        c, p = h["chrom"], int(h["pos"])
        if c not in ref.references or p - FLANK < 1:
            continue
        seq = ref.fetch(c, p - 1 - FLANK, p + FLANK).upper()
        keys[k] = (h, seq)
        fa.append(f">{len(keys) - 1}\n{seq}\n")
    names = list(keys)
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as tf:
        tf.write("".join(fa))
    sam = subprocess.run([a.minimap2, "-ax", "sr", "-N", "20", "-p", "0.2", "-t", str(a.threads), "--eqx", "--secondary=yes", a.ref, tf.name],
                         capture_output=True, text=True, check=True).stdout
    os.unlink(tf.name)
    recs = parse_sam(sam, FLANK)
    res = {}
    for i, k in enumerate(names):
        h, seq = keys[k]
        c, p = h["chrom"], int(h["pos"])
        homs = []
        for r in recs.get(str(i), []):
            tpos = F.query_to_ref(r["cigar"], r["pos0"], r["qidx"])
            if tpos is None:
                continue
            if r["chrom"] == c and abs(tpos - (p - 1)) <= 5:
                continue  # the site itself
            tb = ref.fetch(r["chrom"], tpos, tpos + 1).upper()
            base = comp(tb) if r["rev"] else tb  # homologous base in the orientation of the hit's own strand
            ident, cov = F.aln_identity(r["cigar"], r["nm"], r["qlen"])
            homs.append(dict(chrom=r["chrom"], pos=tpos + 1, base=base, strand="-" if r["rev"] else "+", identity=round(ident, 3), coverage=round(cov, 3)))
        pc = F.paralog_call(h["alt"], h["ref"], homs)
        res[k] = dict(patient=h["patient"], gene=h["gene"], site=h["site"], paralog=dict(pc, n_homologs=len(homs), homologs=homs[:6]),
                      segdup=bool(seg.contains(c, p)) if seg else None, umap_k100=bool(umap.contains(c, p)) if umap else None)
    json.dump(res, open(a.out, "w"), indent=1)
    print("sites", len(res), "paralog", sum(1 for v in res.values() if v["paralog"]["paralog"]),
          "segdup", sum(1 for v in res.values() if v["segdup"]), "low_map", sum(1 for v in res.values() if v["umap_k100"] is False))


def mapq_cell(args):
    import pysam
    cram, ref, sites = args
    out = []
    f = pysam.AlignmentFile(cram, "rc", reference_filename=ref)
    for k, c, p, r_, a_ in sites:
        for pc in f.pileup(c, p - 1, p, truncate=True, stepper="nofilter", min_base_quality=0, max_depth=10000, ignore_overlaps=False):
            for pr in pc.pileups:
                rd = pr.alignment
                if pr.is_del or pr.query_position is None or rd.is_secondary or rd.is_supplementary or rd.is_duplicate:
                    continue
                b = rd.query_sequence[pr.query_position]
                if b in (r_, a_):
                    out.append((k, b == a_, rd.mapping_quality))
    return out


def mapq_cmd(a):
    import multiprocessing as mp
    hits = load_hits(a.hits)
    tab = {r["pair"]: r["cram"] for r in csv.DictReader(open(a.crams), delimiter="\t")}
    res = json.load(open(a.into))
    jobs = {}
    for h in hits:
        k = f"{h['patient']}|{h['site']}"
        for cid in (h.get("carrier_ids") or [])[: a.max_cells]:
            if cid in tab:
                jobs.setdefault(cid, []).append((k, h["chrom"], int(h["pos"]), h["ref"], h["alt"]))
    acc = {}
    with mp.Pool(a.threads) as pool:
        for part in pool.imap_unordered(mapq_cell, [(tab[c], a.ref, s) for c, s in jobs.items()], chunksize=2):
            for k, isalt, mq in part:
                acc.setdefault(k, {True: [], False: []})[isalt].append(mq)
    for k, v in res.items():
        d = acc.get(k)
        if not d:
            v["mapq"] = {"available": False}
            continue
        m = lambda x: round(sum(x) / len(x), 2) if x else None
        v["mapq"] = {"available": True, "alt_reads": len(d[True]), "ref_reads": len(d[False]), "mean_alt": m(d[True]), "mean_ref": m(d[False]),
                     "alt_mapq_lt30": sum(1 for x in d[True] if x < 30), "carrier_cells_scanned": min(a.max_cells, len(next(h for h in hits if f"{h['patient']}|{h['site']}" == k).get("carrier_ids") or []))}
    json.dump(res, open(a.into, "w"), indent=1)
    print("mapq done", sum(1 for v in res.values() if v["mapq"].get("available")), "of", len(res))


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("sites")
    s.add_argument("--hits", required=True); s.add_argument("--ref", required=True); s.add_argument("--segdup"); s.add_argument("--umap")
    s.add_argument("--minimap2", default="minimap2"); s.add_argument("--threads", type=int, default=4); s.add_argument("--out", required=True)
    m = sp.add_parser("mapq")
    m.add_argument("--hits", required=True); m.add_argument("--crams", required=True); m.add_argument("--ref", required=True)
    m.add_argument("--into", required=True); m.add_argument("--max-cells", type=int, default=60); m.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    sites_cmd(a) if a.cmd == "sites" else mapq_cmd(a)


if __name__ == "__main__":
    main()
