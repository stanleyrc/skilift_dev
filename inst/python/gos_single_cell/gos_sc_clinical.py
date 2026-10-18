#!/usr/bin/env python3
"""Parse the clinical assay reports (MSK-style Oncopanel and Genexus free text) of the single-cell patients.

Input: one text file with one block per patient (header line "<PATIENT> ..." then the report text). Encoding is
Mac-Roman with CR line endings in the file we received; any ASCII/UTF-8/latin text works.

Writes data/_cohort/clinical_assays.json:
  { format, created, source, patients: { <P>: { raw, assays: [ {
        assay: "oncopanel"|"genexus", status: "ok"|"failed"|"nothing detected", raw,
        tmb_per_mb, variants:[{gene, cdna, protein, exon, tier, vaf, reads, raw, ...}],
        svs:[{kind, tier, left, right, raw}], cnvs:[{location, scope, type, genes, copies, raw}],
        negatives:[...], tert_promoter:{C228T, C250T}, notes:[...] } ] } },
    patients_without_clinical_data: [..] }
Every parsed entry keeps the verbatim report string in `raw`; nothing is inferred beyond the text.

usage: gos_sc_clinical.py [--text FILE] [--data DIR]   (also: import parse_report / parse_text)
"""
import argparse
import datetime
import glob
import json
import os
import re

DATA = "/gpfs/commons/projects/imielinski_web/sclarke/gOS_GBM_sc/data"
TEXT = "/gpfs/commons/groups/imielinski_lab/projects/single_cell_GBM/data/2026_10_15_figures/GBM_PTAGBM_PTA_cohort_mutation_from_oncopanel_genexus.txt"
META = "/gpfs/commons/groups/imielinski_lab/projects/single_cell_GBM/data/metadata/clinical/Patient.Metadata.30Mar26.txt"
# primary vs recurrent (from the user, 2026-10-18; not in Patient.Metadata)
STATUS = "/gpfs/commons/groups/imielinski_lab/projects/single_cell_GBM/data/metadata/clinical/primary_recurrent_261018.tsv"
MANIFEST = "/gpfs/commons/groups/imielinski_lab/projects/single_cell_GBM/analysis/site_filter/patient_manifest.tsv"
FORMAT = "gos-sc-clinical/1"

HEADER = re.compile(r"^\s*((?:BWH|MGH|ZZZ)\d+)\b(.*)$")
ASSAY_MARK = re.compile(r"^\s*=+\s*(ONCOPANEL|GENEXUS)\s*=+\s*$", re.I)
TIER = re.compile(r"^\s*Tier\s+(\d)\s+variants?:\s*(.*)$", re.I)
VARIANT = re.compile(
    r"^\s*([A-Za-z0-9][A-Za-z0-9\-]*)\s+(c\.\S*)?\s*\(([^)]*)\)\s*(?:,\s*exon\s+(\d+))?\s*-\s*in\s*(\d+)?\s*%\s*of\s*(\d+)?\s*reads")
CNV_TYPES = r"Amplification|Gain|Loss|Two copy deletion|Deletion|Homozygous deletion"
LOC = re.compile(r"^(?:\d{1,2}|X|Y)(?:[pq][\dpq.\-]*)?(?:\s+(?:Chromosomal level|Arm level|Focal event))?$")
SV_START = re.compile(r"^\s*(Rearrangement\w*|Deletion|Duplication|Insertion|Inversion|Translocation)\s+-\s+")
SV_FULL = re.compile(r"^(\S+)\s+-\s+(.+?)\s+\((chr\w+):(\d+)\)\s*::\s*(.+?)\s+\((chr\w+):(\d+)\)\s*$")
GX_VARIANT = re.compile(r"^\s*([A-Z][A-Za-z0-9\-]+)\s+(p\.\S+)\s+([\d.]+)%\s*(?:\(([^)]+)\))?")
GX_FUSION = re.compile(r"rearrangement was identified involving\s+(\S+)\s+Exon\s+(\d+)\s+and\s+(\S+)\s+Exon\s+(\d+)", re.I)
GX_CNV = re.compile(r"^\s*([A-Z][A-Za-z0-9\-]+)\s+(gain|loss|amplification|amp|deletion)\s*$", re.I)
GX_NEG = re.compile(r"^\s*No\s+([A-Z][A-Za-z0-9\-]+)\s+(amp\w*|loss|deletion|gain)\s*$", re.I)
GX_CNV_OF = re.compile(r"^\s*((?:Homozygous|Heterozygous)\s+)?(loss|gain|deletion|amplification)\s+of\s+([A-Z][A-Za-z0-9\-]+)\s*$", re.I)
TERT = re.compile(r"TERT\s+(C228T|C250T)\s*\([^)]*\):\s*(present|absent)", re.I)
STOP = re.compile(r"^\s*(Mutations?:|Structural Variants:|Copy Number Variants:|Investigational Variants|ADDITIONAL FINDINGS|"
                  r"Investigational Mutational Signatures|Tier\s+\d|Cytoband/Size|====|Mutational Burden|Tumor Mutational)", re.I)


def decode(raw_bytes):
    """Report bytes -> text with \\n line ends. The file is Mac-Roman (0xD0 is an en dash) with CR line ends."""
    for enc in ("utf-8", "mac_roman"):
        try:
            text = raw_bytes.decode(enc)
            break
        except UnicodeDecodeError:
            continue
    return text.replace("\r\n", "\n").replace("\r", "\n")


def split_patients(text):
    """-> ordered {patient: block text}; a block starts at a line beginning with a patient id."""
    blocks, cur = {}, None
    for line in text.split("\n"):
        m = HEADER.match(line)
        if m and (cur is None or not line.startswith(" ")):
            cur = m.group(1)
            blocks.setdefault(cur, []).append(line)
        elif cur:
            blocks[cur].append(line)
    return {p: "\n".join(v).strip("\n") for p, v in blocks.items()}


def split_assays(block):
    """-> [(assay, text)] by the ==== ONCOPANEL/GENEXUS ==== markers; header-only mentions give failed/none."""
    lines = block.split("\n")
    segs, cur = [], None
    head = []
    for line in lines:
        m = ASSAY_MARK.match(line)
        if m:
            cur = [m.group(1).lower(), []]
            segs.append(cur)
        elif cur:
            cur[1].append(line)
        else:
            head.append(line)
    return "\n".join(head), [(a, "\n".join(t)) for a, t in segs]


def _num(s):
    return float(s) if s not in (None, "") else None


def parse_oncopanel(text):
    out = {"tmb_per_mb": None, "variants": [], "svs": [], "cnvs": [], "notes": [], "tert_promoter": {}, "negatives": []}
    m = re.search(r"Tumor Mutational Burden/Megabase:\s*([\d.]+)", text)
    if m:
        out["tmb_per_mb"] = float(m.group(1))
    ctx, tier, investigational = "mutations", None, False
    cur_sv = None
    cur_cnv = None

    def flush_sv():
        nonlocal cur_sv
        if cur_sv:
            joined = re.sub(r"\s+", " ", cur_sv["raw"]).strip()
            m2 = SV_FULL.match(joined)
            cur_sv["raw"] = cur_sv["raw"].strip()
            if m2:
                cur_sv.update(kind=m2.group(1), left={"feature": m2.group(2), "chrom": m2.group(3), "pos": int(m2.group(4))},
                              right={"feature": m2.group(5), "chrom": m2.group(6), "pos": int(m2.group(7))})
            else:
                cur_sv.update(kind=None, left=None, right=None, parse_error="unrecognised SV line")
            out["svs"].append(cur_sv)
            cur_sv = None

    def flush_cnv():
        nonlocal cur_cnv
        if cur_cnv:
            genes_txt = re.sub(r"\s+", " ", cur_cnv["genes_txt"]).strip()
            mc = re.search(r"\(estimated\s+(\d+)\s+copies\)", genes_txt)
            cur_cnv["copies"] = int(mc.group(1)) if mc else None
            genes_txt = re.sub(r"\(estimated\s+\d+\s+copies\)", "", genes_txt)
            cur_cnv["genes"] = [g.strip() for g in genes_txt.split(",") if g.strip()]
            del cur_cnv["genes_txt"]
            cur_cnv["raw"] = re.sub(r"[ \t]+", " ", cur_cnv["raw"]).strip()
            out["cnvs"].append(cur_cnv)
            cur_cnv = None

    for line in text.split("\n"):
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            out["notes"].append(s)
            continue
        if re.match(r"^Investigational Variants", s, re.I):
            flush_sv(); flush_cnv(); investigational, ctx, tier = True, "mutations", None
            continue
        if re.match(r"^Mutations?:", s, re.I):
            flush_sv(); flush_cnv(); ctx, tier = "mutations", None
            continue
        if re.match(r"^Structural Variants:", s, re.I):
            flush_sv(); flush_cnv(); ctx, tier = "svs", None
            continue
        if re.match(r"^Copy Number Variants:", s, re.I):
            flush_sv(); flush_cnv(); ctx, tier = "cnvs", None
            continue
        if re.match(r"^(ADDITIONAL FINDINGS|Investigational Mutational Signatures|Mutational Burden|Tumor Mutational)", s, re.I):
            flush_sv(); flush_cnv()
            if not re.match(r"^Tumor Mutational", s, re.I) and not re.match(r"^Mutational Burden", s, re.I):
                ctx = "other"
            continue
        mt = TIER.match(s)
        if mt:
            flush_sv(); flush_cnv()
            tier = int(mt.group(1))
            continue
        mtert = TERT.search(s)
        if mtert:
            out["tert_promoter"][mtert.group(1).upper()] = {"status": mtert.group(2).lower(), "raw": s}
            continue
        if ctx == "mutations" and tier:
            mv = VARIANT.match(s)
            if mv:
                gene, cdna, prot, exon, vaf, reads = mv.groups()
                out["variants"].append({
                    "gene": gene, "cdna": cdna or None, "protein": prot or None, "exon": int(exon) if exon else None,
                    "tier": tier, "investigational": tier >= 3, "vaf": (int(vaf) / 100.0 if vaf else None),
                    "reads": int(reads) if reads else None, "raw": s,
                    "note": None if vaf else "report gives no allele fraction or read count"})
            continue
        if ctx == "svs" and tier:
            if SV_START.match(s):
                flush_sv()
                cur_sv = {"tier": tier, "investigational": investigational, "raw": s}
            elif cur_sv:
                cur_sv["raw"] += " " + s
            continue
        if ctx == "cnvs":
            if s.startswith("Cytoband/Size"):
                continue
            parts = re.split(r"\s{2,}", s, maxsplit=2)
            if len(parts) == 3 and LOC.match(parts[0]) and re.match(r"^(%s)$" % CNV_TYPES, parts[1], re.I):
                flush_cnv()
                loc = parts[0]
                scope = ("chromosome" if "Chromosomal" in loc else "arm" if "Arm level" in loc else
                         "focal" if "Focal" in loc else "band")
                cur_cnv = {"location": loc, "scope": scope, "investigational": investigational, "type": parts[1], "genes_txt": parts[2], "raw": s}
            elif cur_cnv:
                cur_cnv["genes_txt"] += " " + s
                cur_cnv["raw"] += " " + s
            continue
    flush_sv()
    flush_cnv()
    return out


def parse_genexus(text):
    out = {"variants": [], "svs": [], "cnvs": [], "negatives": [], "notes": [], "tert_promoter": {}, "summary_genes": []}
    ctx = None
    for line in text.split("\n"):
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            out["notes"].append(s)
            continue
        ms = re.match(r"^Findings with Potential Clinical Significance:\s*(.*)$", s)
        if ms and ctx is None and ms.group(1):
            out["summary_genes"] = [g.strip() for g in ms.group(1).split(",") if g.strip()]
            continue
        if re.match(r"^(RESULTS:|SUMMARY:)", s, re.I):
            ctx = None if s.upper().startswith("SUMMARY") else ctx
            continue
        if re.match(r"^Sequence Variants:?", s, re.I):
            ctx = "snv"; continue
        if re.match(r"^Fusions and Translocations:?", s, re.I):
            ctx = "fusion"; continue
        if re.match(r"^Copy Number Variants:?", s, re.I):
            ctx = "cnv"; continue
        mt = TERT.search(s)
        if mt:
            out["tert_promoter"][mt.group(1).upper()] = {"status": mt.group(2).lower(), "raw": s}
            continue
        if re.match(r"^Findings with ", s) or re.match(r"^Gene\s+Variant", s, re.I):
            continue
        mo = GX_CNV_OF.match(s)
        if mo:
            kind = ((mo.group(1) or "").strip() + " " + mo.group(2).lower()).strip().capitalize()
            out["cnvs"].append({"location": None, "scope": "gene", "type": kind, "genes": [mo.group(3)], "copies": None, "raw": s})
            continue
        if ctx == "snv":
            mv = GX_VARIANT.match(s)
            if mv:
                out["variants"].append({"gene": mv.group(1), "cdna": None, "protein": mv.group(2), "exon": None, "tier": None,
                                        "investigational": False, "vaf": round(float(mv.group(3)) / 100.0, 4), "reads": None,
                                        "transcript": mv.group(4), "raw": s})
            elif s.lower() != "none":
                out["notes"].append(s)
        elif ctx == "fusion":
            mf = GX_FUSION.search(s)
            if mf:
                out["svs"].append({"kind": "Fusion", "tier": None, "raw": s,
                                   "left": {"feature": mf.group(1), "exon": int(mf.group(2))},
                                   "right": {"feature": mf.group(3), "exon": int(mf.group(4))}})
            elif s.lower() != "none":
                out["notes"].append(s)
        elif ctx == "cnv":
            mn = GX_NEG.match(s)
            mc = GX_CNV.match(s)
            if mn:
                out["negatives"].append({"gene": mn.group(1), "what": mn.group(2).lower(), "raw": s})
            elif mc:
                t = mc.group(2).lower()
                out["cnvs"].append({"location": None, "scope": "gene", "type": {"amp": "amplification"}.get(t, t).capitalize(),
                                    "genes": [mc.group(1)], "copies": None, "raw": s})
            elif s.lower() != "none":
                out["notes"].append(s)
    return out


def parse_report(patient, block):
    """One patient block -> {raw, assays:[...]}; status is only what the text says."""
    head, segs = split_assays(block)
    assays = []
    if not segs:
        low = block.lower()
        if "failed" in low:
            kind = "oncopanel" if "oncopanel" in low else "unknown"
            assays.append({"assay": kind, "status": "failed", "raw": block, "variants": [], "svs": [], "cnvs": [],
                           "negatives": [], "notes": [head.strip().split("\n")[0]], "tmb_per_mb": None, "tert_promoter": {}})
    for kind, text in segs:
        low = text.lower()
        if re.search(r"nothing was detected|no (variants|findings)", low) and not re.search(r"tier|sequence variants:\s*\n\s*\S", low):
            rec = {"assay": kind, "status": "nothing detected", "raw": text.strip("\n"), "variants": [], "svs": [], "cnvs": [],
                   "negatives": [], "notes": [], "tmb_per_mb": None, "tert_promoter": {}}
        else:
            body = parse_oncopanel(text) if kind == "oncopanel" else parse_genexus(text)
            rec = {"assay": kind, "status": "ok", "raw": text.strip("\n"), **body}
            rec.setdefault("tmb_per_mb", None)
        assays.append(rec)
    return {"raw": block, "assays": assays}


def parse_text(text):
    return {p: parse_report(p, b) for p, b in split_patients(text).items()}


def parse_meta(text):
    """Patient.Metadata TSV (header without a patient column; row names first) -> {P: {location, mgmt, age, sex}}.
    MRCA_mut (old_truncal_mutations) is obsolete (trees were rebuilt) and kept only as old_truncal_mutations, not for display. 'N/A' -> None."""
    lines = [l for l in text.split("\n") if l.strip()]
    cols = lines[0].split("\t")
    out = {}
    for l in lines[1:]:
        f = l.split("\t")
        row = dict(zip(cols, f[1:]))
        na = lambda v: None if v is None or v.strip() in ("", "N/A", "NA") else v.strip()
        out[f[0].strip()] = {"location": na(row.get("Location")), "mgmt": na(row.get("MGMT_MET")),
                             "age": int(row["Age"]) if na(row.get("Age")) else None, "sex": na(row.get("Sex")),
                             "old_truncal_mutations": int(row["MRCA_mut"]) if na(row.get("MRCA_mut")) else None}
    return out


def parse_status(text):
    """primary_recurrent TSV ('#' comments, header patient / disease_status) -> {P: "Primary" | "Recurrent"}."""
    lines = [l for l in text.split("\n") if l.strip() and not l.startswith("#")]
    cols = lines[0].split("\t")
    return {r["patient"].strip(): r["disease_status"].strip() for r in (dict(zip(cols, l.split("\t"))) for l in lines[1:]) if r.get("disease_status", "").strip()}


def mapped_truncal(patient, manifest=MANIFEST):
    """Sites of the sf9_map mutation map placed on the root/trunk branch (mapped node = the whole tree, i.e. the node
    with the maximal descendant tip count). -> (count, source file) or (None, reason)."""
    import csv
    import gzip
    if not os.path.exists(manifest):
        return None, "not available (no patient_manifest.tsv)"
    row = next((r for r in csv.DictReader(open(manifest), delimiter="\t") if r["patient"] == patient), None)
    if not row:
        return None, "not available (no mapping for current tree)"
    base = os.path.basename(row.get("tree_search") or "").replace(".raxml.bestTree", "")
    hits = glob.glob(os.path.join(row["run_dir"], "site_qc", "final", "mapping", "*%s.mutation_map.tsv.gz" % base)) if base else []
    if not hits:
        return None, "not available (no mapping for current tree)"
    n, top = 0, 0
    with gzip.open(hits[0], "rt") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            if not r["descendant_tip_count"]:
                continue
            t = int(r["descendant_tip_count"])
            if t > top:
                n, top = 0, t
            if t == top:
                n += 1
    return n, hits[0]


def tree_site_truncal(patient, data=DATA):
    f = os.path.join(data, patient, "snv_matrix.json")
    if not os.path.exists(f):
        return None
    with open(f) as fh:
        return sum(1 for v in json.load(fh)["variants"] if v.get("category") == "truncal")


def build(text_path=TEXT, data=DATA, meta_path=META, status_path=STATUS):
    with open(text_path, "rb") as fh:
        parsed = parse_text(decode(fh.read()))
    cohort = sorted(os.path.basename(os.path.dirname(f)) for f in glob.glob(os.path.join(data, "[BM]*", "tree.nwk"))
                    if re.fullmatch(r"(BWH|MGH)\d+", os.path.basename(os.path.dirname(f))))
    meta = {}
    if meta_path and os.path.exists(meta_path):
        with open(meta_path, "rb") as fh:
            meta = parse_meta(decode(fh.read()))
    status = {}
    if status_path and os.path.exists(status_path):
        with open(status_path, "rb") as fh:
            status = parse_status(decode(fh.read()))
    for pid in meta:
        meta[pid]["disease_status"] = status.get(pid)
        meta[pid]["disease_status_source"] = os.path.basename(status_path) if pid in status else None
        n, src = mapped_truncal(pid)
        meta[pid]["truncal_mutations_mapped"] = n
        meta[pid]["truncal_mutations_mapped_source"] = src if n is not None else None
        meta[pid]["truncal_mutations_mapped_note"] = None if n is not None else src
        meta[pid]["truncal_tree_sites"] = tree_site_truncal(pid, data)
    return {"format": FORMAT, "clinical_meta": {"source": os.path.basename(meta_path or ""), "patients": meta,
            "missing_fields": [] if status else ["primary/recurrent"],
            "note": "old_truncal_mutations is obsolete (trees rebuilt); not for display"}, "created": datetime.datetime.now().isoformat(timespec="seconds"),
            "source": os.path.basename(text_path), "patients": parsed,
            "patients_without_clinical_data": [p for p in cohort if p not in parsed]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--text", default=TEXT)
    ap.add_argument("--data", default=DATA)
    a = ap.parse_args()
    out = build(a.text, a.data)
    dest = os.path.join(a.data, "_cohort", "clinical_assays.json")
    if os.path.exists(dest) and not os.path.exists(dest + ".bak_261017"):
        os.replace(dest, dest + ".bak_261017")
    with open(dest, "w") as fh:
        json.dump(out, fh, indent=1)
    for p, r in out["patients"].items():
        for a_ in r["assays"]:
            print(p, a_["assay"], a_["status"], "tmb", a_.get("tmb_per_mb"), "var", len(a_["variants"]), "sv", len(a_["svs"]),
                  "cnv", len(a_["cnvs"]))
    print("no clinical data:", out["patients_without_clinical_data"])


if __name__ == "__main__":
    main()
