#!/usr/bin/env python3
"""Columnar SNV site sets for gOS ("gos-sc-snv-cols/1"): a site table that loads on its own and a per-cell read matrix fetched on demand.

Why: snv_matrix.json stores every (cell, site) read pair as a JSON array (MGH285: 385 MB, 25 M arrays), which the browser must parse in one go
before the patient view can draw. The columnar layout keeps the same content in two files:
  <name>.sites.json   site table. Every column is either dictionary-encoded ({"type": "dict", "dict": [...], "codes": [...]}), numeric
                      ({"type": "num", "values": [...]}, null = missing), boolean ({"type": "bool", "values": [1, 0, null]}) or plain
                      ({"type": "json", "values": [...]}). Sites are sorted by `group_key` (the branch anchor) so each branch's sites are
                      contiguous; `matrix.groups` gives the byte range of every branch in the cell matrix (HTTP Range fetch per branch).
  <name>.cells.bin    (and the same bytes gzipped as <name>.cells.bin.gz, for whole-file fetches: the web server does not compress) sparse per-cell reads, CSR by site, little-endian: Uint32 offsets[n_sites + 1] (entry index; padded to 8 bytes), then
                      8-byte entries {u16 cell index into sites.json "cells", u16 ref reads, u16 alt reads, i8 genotype call (1 alt, 0 ref,
                      -1 none), u8 0}. Counts above 65534 are clipped; ref = alt = 65535 marks a genotype call without read counts.
Sets:
  tree       data/<P>/snv_matrix.json converted as is (category, anchors, annotation, per-cell reads and PL calls) -> snv_matrix.sites.json/.cells.bin
  cn_mapped  sites placed on the current tree by the CN-aware mapping (srctools inst/analysis/cn_tree_mapping, staged per patient as
             final/cn_mapped.tsv.gz + final/cn_mapped_cells.tsv.gz + map/cell_index.json) -> cn_mapped_snvs.sites.json/.cells.bin, every site
             flagged cn_mapped = true and category "cn_mapped" so the UI can tell them from tree-building sites.
Drivers (--drivers, cn_mapped only): CN-mapped protein-altering / splice-site SNVs in a cancer gene (--cancer-genes TSV with a `gene` column)
become filtered.events.json rows with evidence "CN-aware mapped" (Tier 3: not tiered by OncoKB); a site that is already an event (e.g. a
Tier K known-driver rescue) is not duplicated, it gets the branch as `cn_mapped_branch`.

usage: gos_sc_snv_columnar.py tree --data <gOS data> --out <dir> --patients A,B
       gos_sc_snv_columnar.py cn_mapped --data <gOS data> --work <staging work> --out <dir> --patients A,B [--drivers --cancer-genes genes.tsv]
  --out: data/<P>/ files are written under <out>/<P>/ (a staging copy; --out == --data writes in place after a .<tag> backup).
"""
import argparse
import csv
import datetime
import gzip
import json
import os
import shutil
import struct
import sys

FORMAT = "gos-sc-snv-cols/1"
EVIDENCE = "CN-aware mapped"
SOURCE = "CN-aware tree mapping (srctools cn_tree_mapping)"
PROTEIN_ALTERING = {"nonsynonymous_SNV", "stopgain", "stoploss", "startloss"}
# ANNOVAR refGene classes -> the SnpEff (Sequence Ontology) terms of snv_matrix.json, so the gOS region / effect tracks read them alike
ANNOVAR_SO = {"nonsynonymous_SNV": "missense_variant", "synonymous_SNV": "synonymous_variant", "stopgain": "stop_gained", "stoploss": "stop_lost",
              "startloss": "start_lost", "splicing": "splice_region_variant", "intronic": "intron_variant", "intergenic": "intergenic_variant",
              "UTR3": "3_prime_UTR_variant", "UTR5": "5_prime_UTR_variant", "upstream": "upstream_gene_variant", "downstream": "downstream_gene_variant",
              "upstream;downstream": "upstream_gene_variant", "ncRNA_exonic": "non_coding_transcript_exon_variant", "ncRNA_intronic": "intron_variant",
              "ncRNA_splicing": "splice_region_variant", "UTR5;UTR3": "5_prime_UTR_variant", "exonic;splicing": "splice_region_variant"}


def so_term(func, exonic_func):
    if exonic_func and exonic_func not in (".", "unknown"):
        return ANNOVAR_SO.get(exonic_func, exonic_func)
    return ANNOVAR_SO.get(func, func or None)
MAX_U16 = 65535
MAX_COUNT = 65534


# ---------------------------------------------------------------- encoding

def encode_column(values):
    """Pick a compact JSON encoding for one column (list of Python values, None = missing)."""
    present = [v for v in values if v is not None]
    if all(isinstance(v, bool) for v in present):
        return {"type": "bool", "values": [None if v is None else int(v) for v in values]}
    if all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in present):
        return {"type": "num", "values": [None if v is None else (round(v, 4) if isinstance(v, float) else v) for v in values]}
    if all(isinstance(v, str) for v in present):
        uniq = {}
        codes = []
        for v in values:
            if v is None:
                codes.append(-1)
                continue
            if v not in uniq:
                uniq[v] = len(uniq)
            codes.append(uniq[v])
        if len(uniq) <= max(1, len(values) // 2):
            return {"type": "dict", "dict": list(uniq), "codes": codes}
        return {"type": "str", "values": values}
    return {"type": "json", "values": values}


def decode_column(col):
    t = col["type"]
    if t == "dict":
        d = col["dict"]
        return [None if c < 0 else d[c] for c in col["codes"]]
    if t == "bool":
        return [None if v is None else bool(v) for v in col["values"]]
    return list(col["values"])


def encode_matrix(entries):
    """entries: per site a list of (cell index, ref, alt, gt[, pl]). -> (bytes, offsets list, header bytes).

    pl is the PL confidence byte stored in the 8th byte of the entry (see pl_byte): 0 = no PL data, 1..250 = P(ALT) quantised,
    255 = PL present but flat (no read evidence). Entries without a 5th element get 0."""
    n = len(entries)
    offsets = [0]
    for e in entries:
        offsets.append(offsets[-1] + len(e))
    head = struct.pack(f"<{n + 1}I", *offsets)
    pad = (-len(head)) % 8
    head += b"\0" * pad
    try:
        import numpy as np
        dt = np.dtype([("c", "<u2"), ("r", "<u2"), ("a", "<u2"), ("g", "i1"), ("p", "u1")])
        arr = np.zeros(offsets[-1], dtype=dt)
        flat = [(tuple(x) + (0,))[:5] for e in entries for x in e]
        if flat:
            m = np.array(flat, dtype=np.int64).reshape(-1, 5)
            arr["c"] = m[:, 0]
            both = (m[:, 1] == MAX_U16) & (m[:, 2] == MAX_U16)
            arr["r"] = np.where(both, MAX_U16, np.minimum(m[:, 1], MAX_COUNT))
            arr["a"] = np.where(both, MAX_U16, np.minimum(m[:, 2], MAX_COUNT))
            arr["g"] = np.where((m[:, 3] == 0) | (m[:, 3] == 1), m[:, 3], -1)
            arr["p"] = np.where((m[:, 4] >= 0) & (m[:, 4] <= 255), m[:, 4], 0)
        body = arr.tobytes()
    except ImportError:
        body = bytearray(8 * offsets[-1])
        k = 0
        for e in entries:
            for ci, r, al, gt, *rest in e:
                pl = int(rest[0]) if rest and 0 <= int(rest[0]) <= 255 else 0
                both = r == MAX_U16 and al == MAX_U16
                struct.pack_into("<HHHbB", body, 8 * k, ci, MAX_U16 if both else min(int(r), MAX_COUNT), MAX_U16 if both else min(int(al), MAX_COUNT),
                                 int(gt) if gt in (0, 1) else -1, pl)
                k += 1
        body = bytes(body)
    return head + body, offsets, len(head)


PL_NONE = 0       # no PL data for this cell/site (older exports)
PL_FLAT = 255     # PL present but flat (all genotypes equally likely): no read evidence, P(ALT) meaningless


def pl_byte(pl00, pl01, pl11, q=None):
    """Phred PLs of genotypes 0/0, 0/1, 1/1 -> PL confidence byte (same maths as srctools gos_sc_upload.R pl_conf_byte).

    L = 10^(-PL/10). q = per-site prior probability of an ALT genotype (fraction of cells called ALT at the site, split equally
    between het and hom-alt); None = flat prior over the three genotypes. P(ALT) = q/2 (L01 + L11) / ((1-q) L00 + q/2 (L01 + L11)),
    byte = 1 + round(249 P) (1..250). Missing PL (None/NaN) or a flat PL (all three equal: no reads) -> PL_FLAT.
    Real PLs only: nothing is derived from the hard genotype call."""
    vals = (pl00, pl01, pl11)
    try:
        if any(v is None or v != v for v in vals):
            return PL_FLAT
        if pl00 == pl01 == pl11:
            return PL_FLAT
        l = [10.0 ** (-float(v) / 10.0) for v in vals]
    except (TypeError, ValueError):
        return PL_FLAT
    if q is None:
        p = (l[1] + l[2]) / (l[0] + l[1] + l[2])
    else:
        alt = q / 2.0 * (l[1] + l[2])
        p = alt / ((1.0 - q) * l[0] + alt)
    return 1 + int(round(249 * p))


def decode_matrix(buf, n_sites, with_pl=False):
    """Inverse of encode_matrix (tests / checks): per site a list of (cell, ref, alt, gt[, pl])."""
    offsets = struct.unpack_from(f"<{n_sites + 1}I", buf, 0)
    head = 4 * (n_sites + 1)
    head += (-head) % 8
    out = []
    for s in range(n_sites):
        out.append([struct.unpack_from("<HHHbB", buf, head + 8 * k)[:5 if with_pl else 4] for k in range(offsets[s], offsets[s + 1])])
    return out


def write_set(out_dir, name, meta, variants, entries, cells, group_key="anchor"):
    """Sort sites by group, write <name>.sites.json + <name>.cells.bin. variants: list of dicts (same keys); entries: per site cell reads."""
    order = sorted(range(len(variants)), key=lambda i: (str(variants[i].get(group_key) or "~"), i))
    variants = [variants[i] for i in order]
    entries = [entries[i] for i in order]
    blob, offsets, head = encode_matrix(entries)
    groups = []
    for i, v in enumerate(variants):
        key = v.get(group_key)
        if not groups or groups[-1]["key"] != key:
            groups.append({"key": key, "node": v.get("node"), "first_site": i, "n_sites": 0})
        groups[-1]["n_sites"] += 1
    for g in groups:
        g["byte_start"] = head + 8 * offsets[g["first_site"]]
        g["byte_end"] = head + 8 * offsets[g["first_site"] + g["n_sites"]]
    keys = []
    for v in variants:
        for k in v:
            if k not in keys:
                keys.append(k)
    table = dict(meta, format=FORMAT, n_sites=len(variants), cells=list(cells),
                 columns={k: encode_column([v.get(k) for v in variants]) for k in keys},
                 matrix={"file": f"{name}.cells.bin", "file_gz": f"{name}.cells.bin.gz", "encoding": "csr-u16x3-i8-u8 (8-byte entries: cell, ref, alt, gt, pl byte; Uint32 site offsets, 8-byte aligned)", "pl_byte": any(len(x) > 4 and x[4] for e in entries for x in e),
                         "header_bytes": head, "bytes": len(blob), "n_entries": offsets[-1], "group_key": group_key, "groups": groups})
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, f"{name}.cells.bin.tmp"), "wb") as fh:
        fh.write(blob)
    with gzip.open(os.path.join(out_dir, f"{name}.cells.bin.gz.tmp"), "wb", compresslevel=6) as fh:
        fh.write(blob)
    os.replace(os.path.join(out_dir, f"{name}.cells.bin.gz.tmp"), os.path.join(out_dir, f"{name}.cells.bin.gz"))
    with open(os.path.join(out_dir, f"{name}.sites.json.tmp"), "w") as fh:
        json.dump(table, fh, separators=(",", ":"))
    os.replace(os.path.join(out_dir, f"{name}.cells.bin.tmp"), os.path.join(out_dir, f"{name}.cells.bin"))
    os.replace(os.path.join(out_dir, f"{name}.sites.json.tmp"), os.path.join(out_dir, f"{name}.sites.json"))
    return table


# ---------------------------------------------------------------- sets

def tree_set(snv):
    """snv_matrix.json (compact) -> variants, entries, cells."""
    cells = list((snv.get("cells") or {}).keys())
    ci = {c: i for i, c in enumerate(cells)}
    variants = [dict(v) for v in snv.get("variants") or []]
    entries = [[] for _ in variants]
    for c, rows in (snv.get("cells") or {}).items():
        for x in rows:
            vi = x[0]
            r, al = x[1], x[2]
            gt = x[3] if len(x) > 3 and x[3] in (0, 1) else -1
            pl = x[4] if len(x) > 4 and isinstance(x[4], int) and 0 <= x[4] <= 255 else PL_NONE   # 5th element: PL confidence byte
            if r is None or al is None:
                if gt in (0, 1) or pl:      # a genotype call / PL without read counts: ref = alt = 65535
                    entries[vi].append((ci[c], MAX_U16, MAX_U16, gt) + ((pl,) if pl else ()))
                continue
            entries[vi].append((ci[c], r, al, gt) + ((pl,) if pl else ()))
    for e in entries:
        e.sort()
    return variants, entries, cells


def _num(x, cast=float):
    try:
        return cast(x)
    except (TypeError, ValueError):
        return None


def cn_variant(r, cancer_genes=None):
    """One row of final/cn_mapped.tsv.gz -> a site-table variant."""
    ef = r.get("exonic_func") or ""
    func = r.get("func") or ""
    consequence = so_term(func, ef)
    gene = r.get("gene") or ""
    driver = bool(cancer_genes is not None and gene in cancer_genes and (ef in PROTEIN_ALTERING or func.startswith("splicing")))
    status = r.get("status")
    return {
        "id": r["mutation"], "category": "cn_mapped", "cn_mapped": True, "status": status,
        "gene": gene, "consequence": consequence, "annovar_func": func, "annovar_exonic_func": ef, "protein": r.get("protein") or "", "driver": driver,
        "node": r.get("node"), "anchor": r.get("anchor"), "clade_cells": _num(r.get("clade_cells"), int),
        "map_confidence": _num(r.get("posterior")), "posterior": _num(r.get("posterior")), "confidence": r.get("confidence"),
        "multiplicity": r.get("multiplicity"), "lost_in": r.get("lost_in") or "", "lost_anchors": r.get("lost_anchors") or "",
        "lost_cells": _num(r.get("lost_cells"), int), "modal_cn": _num(r.get("modal_cn"), int), "cn_counts": r.get("cn_counts") or "",
        "alt_cells": _num(r.get("t_alt_cells"), int), "alt_in_clade": _num(r.get("alt_in_clade"), int), "alt_outside_clade": _num(r.get("alt_out_clade"), int),
        "cov_in_clade": _num(r.get("cov_in_clade"), int), "cov_outside_clade": _num(r.get("cov_out_clade"), int),
        "d_iid": _num(r.get("d_iid")), "d_art": _num(r.get("d_art")), "germline_llr": _num(r.get("germline_llr")),
        "normal_reads": _num(r.get("rq_normal_reads"), int), "alt_mapq_mean": _num(r.get("rq_alt_mapq_mean")), "alt_bq_mean": _num(r.get("rq_alt_bq_mean")),
        "alt_rev_frac": _num(r.get("rq_alt_rev_frac")), "cosmic_occurrences": _num(r.get("cosmic_occurrences"), int),
    }


def read_tsv_gz(path):
    with gzip.open(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def cn_set(work, patient, cancer_genes=None):
    W = os.path.join(work, patient)
    rows = read_tsv_gz(os.path.join(W, "final", "cn_mapped.tsv.gz"))
    cells = json.load(open(os.path.join(W, "map", "cell_index.json")))["cells"]
    reads = {}
    with gzip.open(os.path.join(W, "final", "cn_mapped_cells.tsv.gz"), "rt") as fh:
        fh.readline()
        for ln in fh:
            mid, s = ln.rstrip("\n").split("\t")
            reads[mid] = sorted((int(i), int(r), int(al), -1) for i, r, al in (e.split(":") for e in s.split(",")))
    variants, entries = [], []
    for r in rows:
        if r["mutation"] not in reads:
            continue
        variants.append(cn_variant(r, cancer_genes))
        entries.append(reads[r["mutation"]])
    return variants, entries, cells, rows


# ---------------------------------------------------------------- drivers

def driver_rows(patient, variants, entries, cells, tumour):
    """filtered.events.json rows for the CN-mapped driver sites (variant['driver'])."""
    out = []
    for v, e in zip(variants, entries):
        if not v.get("driver"):
            continue
        chrom, pos, ref, alt = v["id"].split("_")
        chrom = chrom.replace("chr", "")
        car = [cells[i] for i, r, al, _ in e if al > 0 and cells[i] in tumour]
        cov = [cells[i] for i, r, al, _ in e if r + al > 0 and cells[i] in tumour]
        ta = sum(al for i, r, al, _ in e if cells[i] in tumour)
        tr = sum(r for i, r, al, _ in e if cells[i] in tumour)
        trunk = v.get("status") == "trunk"
        where = "the trunk" if trunk else f"branch {v.get('node')} ({v.get('clade_cells')} cells)"
        loss = f"; mutant allele lost in {v['lost_in']} ({v.get('lost_cells')} cells, LOH)" if v.get("lost_in") else ""
        stype = {"stop_gained": "Nonsense", "stop_lost": "Nonstop", "start_lost": "Start loss", "missense_variant": "Missense"}.get(v.get("consequence"), "Splice site")
        out.append({
            "gene": v["gene"], "fusion_genes": "", "vartype": "SNV", "type": stype, "Variant": v.get("protein") or f"{ref}>{alt}",
            "gene_summary": f"{v['gene']} (cancer gene) SNV in a copy-number-altered region, placed on the tree by the CN-aware mapping (not a tree-building site).",
            "Variant_g": f"{chrom}:{pos}-{pos} {ref}>{alt}", "Tier": 3, "tier_basis": "CN-aware mapped site in a cancer gene (OncoKB cancer gene list); not tiered by OncoKB",
            "variant_summary": f"Mapped to {where} with posterior {v.get('posterior')}{loss}; expected per-cell VAF from the cell's copy number ({v.get('multiplicity')}).",
            "effect": None, "is_multi_hit_per_gene": False, "Genome_Location": f"{chrom}:{pos}-{pos}", "seqnames": chrom, "start": pos, "end": pos,
            "sample": car[0] if car else None, "role": None, "id": patient, "alt": ta, "ref": tr, "n_cells": len(car), "n_normal_cells": 0,
            "cell_ids": ",".join(car), "VAF": round(ta / max(1, ta + tr), 4), "cells": f"{len(car)}/{len(cov)}",
            "cell_fraction": round(len(car) / len(cov), 3) if cov else None, "tumour_cells": len(tumour), "covered_tumour_cells": len(cov),
            "covered_ref_cell_ids": ",".join(c for c in cov if c not in set(car)), "coverage_fraction": round(len(cov) / len(tumour), 3) if tumour else None,
            "hotspot_evidence": {"preset": "genotyping VCF AD_BAM (per-cell CRAM recount)", "cell_source": "genotyping VCF AD_BAM",
                                 "verdict": {"level": "supported", "label": "read QC pass",
                                             "summary": f"0 alt reads in {v.get('normal_reads') or 0} normal-cell reads; alt reads mean MAPQ {v.get('alt_mapq_mean')}, "
                                                        f"BQ {v.get('alt_bq_mean')}; reverse-strand fraction {v.get('alt_rev_frac')}", "reasons": []}},
            "source": SOURCE, "evidence": EVIDENCE, "rescued": True, "cn_mapped": True,
            "clonality": "clonal" if trunk else "subclonal", "clonality_method": "CN-aware tree mapping",
            "cn_mapped_branch": {k: v.get(k) for k in ("node", "anchor", "clade_cells", "posterior", "confidence", "multiplicity", "lost_in", "lost_anchors", "lost_cells", "modal_cn", "cn_counts")},
        })
    return out


def merge_driver_rows(events, rows):
    """Append new rows; a site already present (same seqnames / start / alt base) gets cn_mapped_branch instead. -> (events, n_added, n_annotated)."""
    def key(e):
        g = str(e.get("Variant_g") or "")
        return (str(e.get("seqnames", "")).replace("chr", ""), str(e.get("start")), g[-1:] if ">" in g else None)
    have = {}
    for i, e in enumerate(events):
        if e.get("source") == SOURCE:
            continue
        have.setdefault(key(e), i)
    kept = [e for e in events if e.get("source") != SOURCE]
    added = annotated = 0
    for r in rows:
        k = key(r)
        if k in have:
            i = have[k]
            target = events[i]
            target["cn_mapped_branch"] = r["cn_mapped_branch"]
            annotated += 1
        else:
            kept.append(r)
            added += 1
    return kept, added, annotated


# ---------------------------------------------------------------- main

def backup_write(path, obj, tag):
    if os.path.exists(path) and tag and not os.path.exists(f"{path}.{tag}"):
        shutil.copy2(path, f"{path}.{tag}")
    with open(path + ".tmp", "w") as fh:
        json.dump(obj, fh)
    os.replace(path + ".tmp", path)


def tumour_cells(data, patient):
    p = os.path.join(data, patient, "precompute", "region_calls.json")
    if not os.path.exists(p):
        return set()
    rc = json.load(open(p))
    clone = rc.get("clone") or []
    return {c for i, c in enumerate(rc.get("cells") or []) if not (i < len(clone) and str(clone[i]).lower() == "normal")}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("set", choices=("tree", "cn_mapped"))
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--patients", required=True)
    ap.add_argument("--work")
    ap.add_argument("--drivers", action="store_true")
    ap.add_argument("--cancer-genes")
    ap.add_argument("--tag", default=None)
    a = ap.parse_args(argv)
    now = datetime.datetime.now().isoformat(timespec="seconds")
    genes = None
    if a.cancer_genes:
        with open(a.cancer_genes) as fh:
            genes = {r["gene"] for r in csv.DictReader(fh, delimiter="\t")}
    summary = {}
    for P in a.patients.split(","):
        od = os.path.join(a.out, P)
        if a.set == "tree":
            p = os.path.join(a.data, P, "snv_matrix.json")
            if not os.path.exists(p):
                summary[P] = "not available (no snv_matrix.json)"
                continue
            snv = json.load(open(p))
            variants, entries, cells = tree_set(snv)
            t = write_set(od, "snv_matrix", {"patient": P, "set": "tree", "created": now, "source": "snv_matrix.json (same content, columnar)"}, variants, entries, cells)
            summary[P] = {"sites": t["n_sites"], "entries": t["matrix"]["n_entries"], "bytes": t["matrix"]["bytes"]}
        else:
            W = os.path.join(a.work or "", P, "final", "cn_mapped.tsv.gz")
            if not os.path.exists(W):
                summary[P] = "not available (no CN-aware mapping for this patient)"
                continue
            variants, entries, cells, rows = cn_set(a.work, P, genes)
            vs = os.path.join(a.work, P, "final", "summary.json")
            val = json.load(open(vs)) if os.path.exists(vs) else {}
            meta = {"patient": P, "set": "cn_mapped", "created": now, "source": SOURCE,
                    "method": "sites in copy-number-altered regions (dropped by the diploid tree-building filter) placed on the current tree with a copy-number-aware "
                              "likelihood: expected per-cell VAF = mutant copies / total CN, one copy or all copies of the major allele, loss of the mutant allele "
                              "allowed in LOH clades, PTA dropout; 0 normal alt reads and read-level QC (MAPQ, BQ, strand, read end, split reads, other patients)",
                    "tree": val.get("tree"), "validation": {k: val.get(k) for k in ("counts", "tree_site_consistency", "holdout_normal_B", "d_iid_threshold", "rule", "model")}}
            t = write_set(od, "cn_mapped_snvs", meta, variants, entries, cells)
            summary[P] = {"sites": t["n_sites"], "entries": t["matrix"]["n_entries"], "bytes": t["matrix"]["bytes"], "drivers": sum(1 for v in variants if v["driver"])}
            if a.drivers:
                ev_in = os.path.join(a.out, P, "filtered.events.json")
                ev_src = ev_in if os.path.exists(ev_in) and not os.path.islink(ev_in) else os.path.join(a.data, P, "filtered.events.json")
                events = json.load(open(ev_src)) if os.path.exists(ev_src) else []
                rows_d = driver_rows(P, variants, entries, cells, tumour_cells(a.data, P))
                events, n_add, n_ann = merge_driver_rows(events, rows_d)
                dest = os.path.join(od, "filtered.events.json")
                if os.path.islink(dest):
                    os.remove(dest)
                backup_write(dest, events, a.tag)
                summary[P].update(driver_rows_added=n_add, driver_rows_annotated=n_ann, driver_genes=sorted({r["gene"] for r in rows_d}))
        print(P, json.dumps(summary[P]), flush=True)
    if a.set == "cn_mapped":
        cohort = {"format": "gos-sc-cn-mapped-summary/1", "created": now, "source": SOURCE, "patients": {}}
        for P in a.patients.split(","):
            vs = os.path.join(a.work, P, "final", "summary.json")
            v = json.load(open(vs)) if os.path.exists(vs) else None
            cohort["patients"][P] = {"export": summary.get(P), **({k: v.get(k) for k in ("counts", "tree_site_consistency", "holdout_normal_B", "spectrum_sbs6",
                                                                                         "d_iid_threshold", "positives", "read_flag_counts", "cn_altered_fraction_autosomal")} if v else
                                                                  {"status": "not available (no CN-aware mapping for this patient)"})}
        os.makedirs(os.path.join(a.out, "_cohort"), exist_ok=True)
        backup_write(os.path.join(a.out, "_cohort", "cn_mapped_summary.json"), cohort, a.tag)
    return summary


if __name__ == "__main__":
    main(sys.argv[1:])
