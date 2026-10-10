#!/usr/bin/env python3
"""Exon models of the genes in the single-cell splicing exports, so gOS can draw exons under
sashimi plots and name events (e.g. "PTBP2 exon 10 inclusion") instead of bare junctions.

Reads data/_cohort/rna/splicing.json and every data/<P>/rna/splicing.json; for each cluster gene
keeps, from the 10x GRCh38-2024-A GTF, the Ensembl canonical transcript plus the transcripts whose
introns match at least one cluster junction (up to MAX_TX per gene, most matching junctions first).
Junctions in the exports are introns, 1-based inclusive (first..last intronic base).

Writes data/_cohort/rna/splice_exons.json:
  { format, created, source, genes: { <gene>: { chrom, strand,
      transcripts: [{ id, canonical, matched, exons: [[start, end], ...] }] } } }
exons are 1-based inclusive, sorted by genomic position; exon numbers follow the transcript's direction.
"""
import datetime
import glob
import json
import os
import re

DATA = "/gpfs/commons/projects/imielinski_web/sclarke/gOS_GBM_sc/data"
GTF = "/gpfs/commons/groups/landau_lab/jzinno/ref/refdata-gex-GRCh38-2024-A/genes/genes.gtf"
MAX_TX = 6


def attr(s, key):
    m = re.search(key + r' "([^"]+)"', s)
    return m.group(1) if m else None


def main():
    files = [os.path.join(DATA, "_cohort", "rna", "splicing.json")] + sorted(glob.glob(os.path.join(DATA, "*", "rna", "splicing.json")))
    junctions = {}  # gene -> set of (start, end)
    for f in files:
        if not os.path.exists(f) or "/_cohort/" in f and f != files[0]:
            continue
        for c in json.load(open(f)).get("clusters", []):
            g = c.get("gene")
            if not g:
                continue
            s = junctions.setdefault(g, set())
            for j in c.get("junctions", []):
                s.add((int(j["start"]), int(j["end"])))
    tx = {}
    with open(GTF) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            g = attr(f[8], "gene_name")
            if g not in junctions:
                continue
            t = attr(f[8], "transcript_id")
            r = tx.setdefault(t, {"gene": g, "chrom": f[0], "strand": f[6], "exons": [], "canonical": False})
            r["exons"].append((int(f[3]), int(f[4])))
            if 'tag "Ensembl_canonical"' in f[8]:
                r["canonical"] = True
    genes = {}
    for t, r in tx.items():
        ex = sorted(r["exons"])
        introns = {(ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1)}
        matched = len(introns & junctions[r["gene"]])
        if not (r["canonical"] or matched):
            continue
        e = genes.setdefault(r["gene"], {"chrom": r["chrom"].replace("chr", ""), "strand": r["strand"], "transcripts": []})
        e["transcripts"].append({"id": t, "canonical": r["canonical"], "matched": matched, "exons": [list(x) for x in ex]})
    for e in genes.values():
        e["transcripts"].sort(key=lambda x: (not x["canonical"], -x["matched"], x["id"]))
        canon = [x for x in e["transcripts"] if x["canonical"]][:1]
        rest = [x for x in e["transcripts"] if not x["canonical"] and x["matched"]][: MAX_TX - len(canon)]
        e["transcripts"] = canon + rest
    out = {"format": "gos-sc-splice-exons/1", "created": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
           "source": "10x GRCh38-2024-A genes.gtf: Ensembl canonical + transcripts with introns matching cluster junctions",
           "genes": genes}
    od = os.path.join(DATA, "_cohort", "rna")
    tmp = os.path.join(od, ".splice_exons.json.tmp")
    with open(tmp, "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    os.replace(tmp, os.path.join(od, "splice_exons.json"))
    os.chmod(os.path.join(od, "splice_exons.json"), 0o664)
    print(f"{len(genes)} genes, {sum(len(e['transcripts']) for e in genes.values())} transcripts", flush=True)


if __name__ == "__main__":
    main()
