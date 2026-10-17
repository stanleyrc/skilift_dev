#!/usr/bin/env python3
"""HGVS cDNA (c.) -> hg38 genomic coordinates for the clinical-assay variants, to look them up in the single-cell data.

Uses the 10x GRCh38-2024-A GTF (protein-coding transcripts of the gene, Ensembl canonical first) and genome.fa. An SNV is
accepted on the first transcript whose reference base equals the report's reference base (so a wrong transcript is
rejected, never guessed); indels/splice/5'UTR variants are given as a genomic region (and flagged ref-unchecked).
Reads data/_cohort/clinical_assays.json; writes data/_cohort/clinical_sites.json:
  { format, patients: { <P>: [ { gene, raw, cdna, ok, chrom, start, end (1-based inclusive), ref, alt, kind,
                                  transcript, strand, ref_check, reason } ] } }
Run with the mskilab_ne1 python (pysam).   usage: gos_sc_clinical_sites.py
"""
import json
import os
import re
import sys

DATA = "/gpfs/commons/projects/imielinski_web/sclarke/gOS_GBM_sc/data"
GTF = "/gpfs/commons/groups/landau_lab/jzinno/ref/refdata-gex-GRCh38-2024-A/genes/genes.gtf"
FA = "/gpfs/commons/groups/landau_lab/jzinno/ref/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
COMP = str.maketrans("ACGTacgt", "TGCAtgca")
POS = r"(-?\d+)([+-]\d+)?"
SNV = re.compile(r"^c\.%s([ACGT])>([ACGT])$" % POS)
RNG = re.compile(r"^c\.%s(?:_%s)?(del|delins([ACGT]+)|dup)$" % (POS, POS))


def attr(s, k):
    m = re.search(k + r' "([^"]+)"', s)
    return m.group(1) if m else None


def load_gtf(genes, gtf=GTF):
    """-> {gene: [transcript dict(id, chrom, strand, exons[(s,e)] in transcript order, start_codon_g, canonical)]}"""
    tx = {}
    with open(gtf) as fh:
        for line in fh:
            if not any(('gene_name "%s"' % g) in line for g in genes):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] not in ("exon", "start_codon"):
                continue
            g = attr(f[8], "gene_name")
            if g not in genes or attr(f[8], "transcript_type") != "protein_coding":
                continue
            t = tx.setdefault(g, {}).setdefault(attr(f[8], "transcript_id"), {
                "id": attr(f[8], "transcript_id"), "chrom": f[0], "strand": f[6], "exons": [], "start": [],
                "canonical": 'tag "Ensembl_canonical"' in f[8]})
            (t["exons"] if f[2] == "exon" else t["start"]).append((int(f[3]), int(f[4])))
    out = {}
    for g, d in tx.items():
        ts = []
        for t in d.values():
            if not t["start"]:
                continue
            plus = t["strand"] == "+"
            t["exons"].sort(reverse=not plus)
            t["start_g"] = min(s for s, _ in t["start"]) if plus else max(e for _, e in t["start"])
            ts.append(t)
        ts.sort(key=lambda t: (not t["canonical"], -sum(e - s for s, e in t["exons"])))
        out[g] = ts
    return out


def g2t(t, g):
    """genomic -> 0-based transcript position (None if intronic)"""
    off, plus = 0, t["strand"] == "+"
    for s, e in t["exons"]:
        if s <= g <= e:
            return off + ((g - s) if plus else (e - g))
        off += e - s + 1
    return None


def t2g(t, p):
    off, plus = 0, t["strand"] == "+"
    for s, e in t["exons"]:
        ln = e - s + 1
        if p < off + ln:
            return (s + p - off) if plus else (e - (p - off))
        off += ln
    return None


def pos_g(t, num, offs):
    """c. position (+ intronic offset) -> genomic"""
    n = int(num)
    cds0 = g2t(t, t["start_g"])
    tp = cds0 + n - 1 if n > 0 else cds0 + n  # c.-1 is the base before the A of ATG
    base = t2g(t, tp)
    if base is None:
        return None
    if offs:
        k = int(offs)
        base += k if t["strand"] == "+" else -k
    return base


BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AA[i] for i, (a, b, c) in enumerate((a, b, c) for a in BASES for b in BASES for c in BASES)}
PROT = re.compile(r"^p\.([A-Z])(\d+)([A-Z*])$")
ONE = {}


def resolve_protein(gene, protein, txs, fa):
    """p.R130* style (no c. in the report): the codon's genomic bases and every single-base change giving that residue."""
    base = {"gene": gene, "cdna": None, "protein": protein}
    m = PROT.match(protein or "")
    if not m or gene not in txs:
        return {**base, "ok": False, "reason": "no c. notation and protein change not a single substitution"}
    ref_aa, n, alt_aa = m.group(1), int(m.group(2)), m.group(3)
    for k, t in enumerate(txs[gene]):
        cds0 = g2t(t, t["start_g"])
        gs = [t2g(t, cds0 + 3 * (n - 1) + i) for i in range(3)]
        if None in gs:
            continue
        cod = "".join(fa.fetch(t["chrom"], g - 1, g).upper() for g in gs)
        if t["strand"] == "-":
            cod = cod.translate(COMP)
        if CODON.get(cod) != ref_aa:
            continue
        opts = []
        for i in range(3):
            for b in "ACGT":
                if b != cod[i]:
                    new = cod[:i] + b + cod[i + 1:]
                    if CODON[new] == alt_aa:
                        gb = b if t["strand"] == "+" else b.translate(COMP)
                        gr = fa.fetch(t["chrom"], gs[i] - 1, gs[i]).upper()
                        opts.append({"pos": gs[i], "ref": gr, "alt": gb})
        return {**base, "ok": True, "chrom": t["chrom"], "start": min(gs), "end": max(gs), "ref": None, "alt": None,
                "kind": "SNV_PROTEIN", "options": opts, "transcript": t["id"], "strand": t["strand"], "tx_rank": k,
                "ref_check": "reference codon %s encodes %s" % (cod, ref_aa)}
    return {**base, "ok": False, "reason": "no protein-coding transcript of %s has %s at codon %d" % (gene, ref_aa, n)}


def predicted_protein(t, num, offs, alt_coding, fa):
    """protein change (e.g. R175H) of a c. SNV on transcript t (None outside the CDS / for intronic offsets)."""
    n = int(num)
    if n <= 0 or offs:
        return None
    cds0 = g2t(t, t["start_g"])
    i = (n - 1) // 3
    gs = [t2g(t, cds0 + 3 * i + k) for k in range(3)]
    if None in gs:
        return None
    cod = "".join(fa.fetch(t["chrom"], x - 1, x).upper() for x in gs)
    alt = alt_coding
    if t["strand"] == "-":
        cod = cod.translate(COMP)
    k = (n - 1) % 3
    new = cod[:k] + alt + cod[k + 1:]
    return "%s%d%s" % (CODON[cod], i + 1, CODON[new])


def resolve(gene, cdna, txs, fa, protein=None):
    """SNV on EVERY protein-coding transcript: the reference base must match and, when the report gives a protein change,
    the codon translation on that transcript must reproduce it (a base match alone proves nothing: 1 in 4 by chance)."""
    base = {"gene": gene, "cdna": cdna}
    if not cdna and protein:
        return resolve_protein(gene, protein, txs, fa)
    if not cdna or gene not in txs:
        return {**base, "ok": False, "reason": "no c. notation in the report" if not cdna else "gene not in GTF protein-coding set"}
    ms = SNV.match(cdna)
    mr = None if ms else RNG.match(cdna)
    if not ms and not mr:
        return {**base, "ok": False, "reason": "c. notation not parsed"}
    want = re.sub(r"^p\.", "", protein or "")
    if ms:
        cands, mism = [], []
        for k, t in enumerate(txs[gene]):
            g = pos_g(t, ms.group(1), ms.group(2))
            if g is None:
                continue
            ref = fa.fetch(t["chrom"], g - 1, g).upper()
            rep, alt = ms.group(3), ms.group(4)
            if t["strand"] == "-":
                rep, alt = rep.translate(COMP), alt.translate(COMP)
            if ref != rep:
                continue
            pp = predicted_protein(t, ms.group(1), ms.group(2), ms.group(4), fa)
            rec = {"transcript": t["id"], "canonical": t["canonical"], "pos": g, "ref": ref, "alt": alt, "predicted_protein": pp, "tx_rank": k}
            (cands if (not want or pp is None or pp == want) else mism).append(rec)
        if not cands:
            return {**base, "ok": False, "mismatches": mism[:6],
                    "reason": "no transcript reproduces the report: %s" % (
                        "reference base matches on %d transcripts but the translation gives %s, not %s" % (len(mism), sorted({m["predicted_protein"] for m in mism}), want) if mism
                        else "no transcript has the reported reference base at that c. position")}
        best = cands[0]
        positions = sorted({c["pos"] for c in cands})
        verified = bool(want) and best["predicted_protein"] == want
        return {**base, "ok": True, "chrom": txs[gene][0]["chrom"], "start": best["pos"], "end": best["pos"], "ref": best["ref"], "alt": best["alt"],
                "kind": "SNV", "transcript": best["transcript"], "strand": txs[gene][best["tx_rank"]]["strand"], "tx_rank": best["tx_rank"],
                "alt_positions": [c for c in positions if c != best["pos"]],
                "alt_candidates": [c for c in cands if c["pos"] != best["pos"]],
                "ref_check": ("reference base and translation (%s) reproduce the report on %s" % (want, best["transcript"])) if verified else
                             "reference base matches on %s; no protein change to verify (UTR / splice / non-coding)" % best["transcript"],
                "protein_verified": verified}
    for k, t in enumerate(txs[gene]):
        g1 = pos_g(t, mr.group(1), mr.group(2))
        g2 = pos_g(t, mr.group(3), mr.group(4)) if mr.group(3) else g1
        if g1 is None or g2 is None:
            continue
        kind = "DEL" if mr.group(5) == "del" else "DUP" if mr.group(5) == "dup" else "DELINS"
        return {**base, "ok": True, "chrom": t["chrom"], "start": min(g1, g2), "end": max(g1, g2), "ref": None,
                "alt": mr.group(6), "kind": kind, "transcript": t["id"], "strand": t["strand"],
                "ref_check": "not checkable (indel/range): canonical transcript used", "tx_rank": k, "protein_verified": False}
    return {**base, "ok": False, "reason": "no transcript maps the interval"}


def main():
    import pysam
    cl = json.load(open(os.path.join(DATA, "_cohort", "clinical_assays.json")))
    genes = {v["gene"] for r in cl["patients"].values() for a in r["assays"] for v in a["variants"]}
    txs = load_gtf(genes)
    fa = pysam.FastaFile(FA)
    out = {}
    for p, r in cl["patients"].items():
        rows = []
        for a in r["assays"]:
            for v in a["variants"]:
                rows.append({**resolve(v["gene"], v["cdna"], txs, fa, v["protein"]), "assay": a["assay"], "raw": v["raw"], "protein": v["protein"]})
        out[p] = rows
    dest = os.path.join(DATA, "_cohort", "clinical_sites.json")
    json.dump({"format": "gos-sc-clinical-sites/1", "reference": "GRCh38 (10x 2024-A)", "patients": out}, open(dest, "w"), indent=1)
    for p, rows in out.items():
        for x in rows:
            print(p, x["gene"], x["cdna"], x.get("chrom"), x.get("start"), x.get("end"), x.get("kind"), x.get("reason", ""))


def write_regions(out_dir, pad=5):
    """regions.tsv (ds_pileup.py format: chrom,start0,end,gene,exon,strand,tx) around every resolved clinical site."""
    d = json.load(open(os.path.join(DATA, "_cohort", "clinical_sites.json")))["patients"]
    seen, rows = set(), []
    for p, xs in d.items():
        for x in xs:
            if x.get("ok"):
                for ap in x.get("alt_positions", []):
                    rows.append((x["chrom"], max(0, ap - 1 - pad), ap + pad, x["gene"], 0, x["strand"], x["transcript"]))
                k = (x["chrom"], x["start"], x["end"])
                if k not in seen:
                    seen.add(k)
                    rows.append((x["chrom"], max(0, x["start"] - 1 - pad), x["end"] + pad, x["gene"], 0, x["strand"], x["transcript"]))
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "regions.tsv"), "w") as fh:
        fh.write("chrom\tstart0\tend\tgene\texon\tstrand\ttx\n")
        for r in sorted(rows):
            fh.write("\t".join(map(str, r)) + "\n")
    return len(rows)


if __name__ == "__main__":
    main()
    if len(sys.argv) > 1:
        print("regions:", write_regions(sys.argv[1]))
