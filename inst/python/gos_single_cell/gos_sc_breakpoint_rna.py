#!/usr/bin/env python3
"""gos_sc_breakpoint_rna: RNA transcribed across somatic DNA breakpoints, for the gOS single-cell RNA tab.

Reads the tables of the srctools analysis splicing_amplicons_junctions (work/sv_crossing.tsv, work/sv_crossing_reads.tsv.gz,
work/sv_junctions.tsv, work/rna_only_vs_sv.tsv, tables/sv_transcribed.tsv) and writes, per patient (stdlib only):
  <out>/<P>/rna/breakpoint_transcripts.json  format breakpoint_transcripts_v1: one entry per SV with RNA reads joining both retained
      sides in more carrier than non-carrier cells (the "transcribed" SVs of sv_transcribed.tsv): carrier / non-carrier cells with such
      reads, read counts, RNA junctions (regtools linear + Arriba incl. discarded), frame (Arriba reading_frame) and retained domains,
      confidence tier, IGV loci of the two breakends.
  <out>/<P>/rna/fusions_dna_support.json  format fusions_dna_support_v1 (sidecar; fusions.json is never rewritten): the RNA-only fusions
      that match a DNA SV (rna_only_vs_sv.tsv), tied to the exact `id` of the fusions.json row (matched on the unordered pair of
      breakpoints, then the gene pair; unmatched -> fusion_id null + match_note), with class canonical_splice (both ends at annotated
      splice sites) or cryptic_end (cryptic / intergenic end next to the breakpoint).
  <out>/<P>/rna/fusion_isoforms.json  format fusion_isoforms_v1: per DNA+RNA fusion gene pair (fusion_pairs.tsv + fusion_isoforms.tsv) the
      breakpoint isoforms (joined exons, frame, cells, reads, retained domains), the matching fusions.json ids, in-frame / read-through flags.
  EGFRvIII (BWH69): the genomic deletion (intron 1 -> intron 7) is the DNA evidence; the exact exon 1 / exon 8 DNA junction is cDNA in the DNA
  libraries (splice mimic) in addition to it. work/egfrviii_del.tsv (srctools 08_egfrviii_del.R) adds the genomic-deletion entry to
  breakpoint_transcripts.json (it is not in sv_transcribed.tsv: every paired tumour cell carries it, so there is no non-carrier group), with the
  `nearby_genomic_sv` junctions found around the cDNA junction (work/nearby_genomic_sv.tsv: tier_v2 >= 2 non-mimic junctions with both breakends
  within 10 kb of the mimic's breakends; shown for cells >= 10).
Patients without rows still get a file with empty arrays so the UI can say "none found" vs "not available".

Confidence rule (levels as gOS findingConfidence.js):
  consistent  carrier cells with RNA across the breakpoint >= 5 AND carrier rate >= 5x non-carrier rate (non-carrier rate 0 passes)
              AND Fisher p < 0.05 (when computable);
  suggestive  carrier cells >= 2 AND carrier rate > non-carrier rate;
  low         otherwise.
  rate = cells with crossing reads / DNA-carrier (or non-carrier) cells with RNA of that SV (n_car_rna, n_non_rna of sv_crossing.tsv).
  Where the denominators are not available (patient without per-cell RNA runs) counts only are used (carrier >= 5 and >= 5x the
  non-carrier count) and that is said in confidence_reasons; no p value then, so such a call can at most be consistent on counts.

Usage:
  gos_sc_breakpoint_rna.py --analysis-dir <analysis/splicing_amplicons_junctions_261022> --gos-data-dir <gOS data> --out-dir <dir>
                           [--patients P1,P2]
Reads the gOS data dir only (<P>/rna/fusions.json for id matching; patients in the analysis tables without one get null fusion_id); writes only below --out-dir."""
import argparse
import csv
import gzip
import json
import math
import os
import re
import sys
import time
from collections import Counter, defaultdict

BT_FORMAT = "breakpoint_transcripts_v1"
FS_FORMAT = "fusions_dna_support_v1"
MIN_CONS, MIN_SUGG, FOLD_CONS, P_CONS = 5, 2, 5.0, 0.05
SOURCE = "srctools inst/splicing_amplicons_junctions (analysis splicing_amplicons_junctions_261022)"
RULE_TEXT = ("consistent: >= 5 carrier cells with RNA across the breakpoint, carrier rate >= 5x non-carrier rate (non-carrier rate 0 "
             "passes) and Fisher p < 0.05; suggestive: >= 2 carrier cells and carrier rate > non-carrier rate; else low. Rate = cells "
             "with crossing reads / cells of that group with RNA for the SV.")
RULE = {"text": RULE_TEXT,
        "consistent": {"min_carrier_cells": MIN_CONS, "min_fold": FOLD_CONS, "max_fisher_p": P_CONS},
        "suggestive": {"min_carrier_cells": MIN_SUGG, "rate": "carrier > non-carrier"}}
CLASS_LABEL = {"canonical_splice": "normal splicing across a DNA breakpoint",
               "cryptic_end": "cryptic / intergenic end next to the breakpoint"}


def rd(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def to_int(v, default=None):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return default


def to_float(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def fisher_two_sided(a, b, c, d):
    """Two-sided Fisher exact p for [[a,b],[c,d]] (sum of tables no more likely than the observed one)."""
    if min(a, b, c, d) < 0:
        return None
    r1, r2, c1, n = a + b, c + d, a + c, a + b + c + d
    if n == 0:
        return None
    lo, hi = max(0, c1 - r2), min(r1, c1)
    denom = math.comb(n, c1)
    p_obs = math.comb(r1, a) * math.comb(r2, c1 - a) / denom
    tot = 0.0
    for x in range(lo, hi + 1):
        p = math.comb(r1, x) * math.comb(r2, c1 - x) / denom
        if p <= p_obs * (1 + 1e-9):
            tot += p
    return min(1.0, tot)


def confidence(car, non, n_car, n_non, fisher_p=None):
    """-> (tier, reasons, car_rate, non_rate, fold, p). n_car / n_non None or 0 = denominators not available."""
    reasons = []
    have = bool(n_car) and n_non is not None and n_non > 0 and car <= n_car and non <= n_non
    if have:
        car_rate, non_rate = car / n_car, non / n_non
        fold = None if non_rate == 0 else car_rate / non_rate
        p = fisher_p if fisher_p is not None else fisher_two_sided(car, n_car - car, non, n_non - non)
        reasons.append("%d/%d carrier cells vs %d/%d non-carrier cells with crossing RNA" % (car, n_car, non, n_non))
        fold_ok = non_rate == 0 or car_rate >= FOLD_CONS * non_rate
        p_ok = p is None or p < P_CONS
        if car >= MIN_CONS and fold_ok and p_ok:
            tier = "consistent"
            reasons.append("carrier cells >= %d, rate %s, Fisher p %s" % (
                MIN_CONS, "fold %.1f" % fold if fold is not None else "non-carrier rate 0", "%.2g" % p if p is not None else "n/a"))
        elif car >= MIN_SUGG and car_rate > non_rate:
            tier = "suggestive"
            why = []
            if car < MIN_CONS:
                why.append("fewer than %d carrier cells" % MIN_CONS)
            if not fold_ok:
                why.append("fold below %g" % FOLD_CONS)
            if not p_ok:
                why.append("Fisher p %.2g >= %g" % (p, P_CONS))
            reasons.append("carrier rate above non-carrier rate but " + "; ".join(why))
        else:
            tier = "low"
            reasons.append("fewer than %d carrier cells or carrier rate not above non-carrier rate" % MIN_SUGG)
        return tier, reasons, car_rate, non_rate, fold, p
    reasons.append("denominators not available (no per-cell RNA runs for this SV's carrier / non-carrier cells): counts only, "
                   "%d carrier vs %d non-carrier cells with crossing RNA" % (car, non))
    if car >= MIN_CONS and car >= FOLD_CONS * non:
        tier = "consistent"
        reasons.append("counts: carrier cells >= %d and >= %gx non-carrier cells; no p value" % (MIN_CONS, FOLD_CONS))
    elif car >= MIN_SUGG and car > non:
        tier = "suggestive"
        reasons.append("counts: carrier cells > non-carrier cells but below the consistent thresholds")
    else:
        tier = "low"
        reasons.append("counts: fewer than %d carrier cells or not above non-carrier cells" % MIN_SUGG)
    return tier, reasons, None, None, None, None


def strip_dist(g):
    """'TNS3(19839),LINC01447(6068)' -> ['TNS3','LINC01447']"""
    out = []
    for part in (g or "").split(","):
        part = re.sub(r"\(\d+\)$", "", part.strip())
        if part and part not in (".", "NA") and part not in out:
            out.append(part)
    return out


def locus(chrom, pos, flank=100):
    c = chrom if str(chrom).startswith("chr") else "chr" + str(chrom)
    return "%s:%d-%d" % (c, max(1, pos - flank), pos + flank)


def source_name(s):
    s = (s or "").split(":")[0]
    return "regtools" if s == "linear" or not s else s


def build_transcripts(patient, trans_rows, reads_by_jid, crossing_by_jid, generated):
    out = []
    for t in trans_rows:
        jid = t["jid"]
        ch_a, ch_b = t["chrA"], t["chrB"]
        pa, pb = to_int(t["posA"]), to_int(t["posB"])
        cx = crossing_by_jid.get(jid, t)
        car, non = to_int(cx["cross_car"], 0), to_int(cx["cross_non"], 0)
        n_car, n_non = to_int(cx["n_car_rna"]), to_int(cx["n_non_rna"])
        rows = reads_by_jid.get(jid, [])
        carrier = [r for r in rows if r["status"] == "carrier"]
        reads_car = sum(to_int(r["reads"], 0) for r in carrier)
        noncar = [r for r in rows if r["status"] == "noncarrier"]
        reads_non = sum(to_int(r["reads"], 0) for r in noncar)
        reads_note = None
        if n_non is None or n_non == 0:
            reads_non, reads_note = None, "no non-carrier cells with RNA for this SV"
        # RNA junctions among carrier cells (group by exact coordinates + caller)
        grp = defaultdict(lambda: {"cells": set(), "genes": Counter(), "frames": Counter(), "domains": Counter()})
        for r in carrier:
            k = ("%s:%s-%s:%s" % (r["chr1"], r["pos1"], r["chr2"], r["pos2"]), source_name(r["source"]))
            g = grp[k]
            g["cells"].add(r["rna_id"])
            g["genes"]["%s::%s" % (",".join(strip_dist(r["gene1"])) or "intergenic", ",".join(strip_dist(r["gene2"])) or "intergenic")] += 1
            if r["reading_frame"] not in ("", ".", None):
                g["frames"][r["reading_frame"]] += 1
            if r["domains"] not in ("", ".", None):
                g["domains"][r["domains"]] += 1
        junctions = []
        for (jx, src), g in sorted(grp.items(), key=lambda kv: -len(kv[1]["cells"])):
            junctions.append({"junction": jx, "source": src, "genes": g["genes"].most_common(1)[0][0], "n_cells": len(g["cells"]),
                              "reading_frame": g["frames"].most_common(1)[0][0] if g["frames"] else None})
        top_frames, top_dom = Counter(), Counter()
        for (jx, src), g in grp.items():
            top_frames.update(g["frames"])
            top_dom.update(g["domains"])
        # frame / consequence from the most-supported Arriba junction that reports one
        frame, consequence, frame_note = None, None, None
        best = None
        for (jx, src), g in sorted(grp.items(), key=lambda kv: -len(kv[1]["cells"])):
            if g["frames"] or g["domains"]:
                best = g
                break
        if best is not None and best["frames"]:
            frame = best["frames"].most_common(1)[0][0]
        if best is not None and best["domains"]:
            consequence = "retained domains: " + best["domains"].most_common(1)[0][0].lstrip("|")
        if frame is None:
            frame_note = "reading frame not reported by Arriba for the carrier-cell junctions (non-coding, intergenic or regtools-only calls)"
        # genes by end: the end nearer to the A (chrA:posA) breakend
        ga, gb = Counter(), Counter()
        for r in carrier:
            for chrom, pos, gene in ((r["chr1"], to_int(r["pos1"]), r["gene1"]), (r["chr2"], to_int(r["pos2"]), r["gene2"])):
                if pos is None:
                    continue
                da = abs(pos - pa) if chrom == ch_a else float("inf")
                db = abs(pos - pb) if chrom == ch_b else float("inf")
                (ga if da <= db else gb).update(strip_dist(gene))
        genes_a, genes_b = [g for g, _ in ga.most_common(5)], [g for g, _ in gb.most_common(5)]
        tier, reasons, car_rate, non_rate, fold, p = confidence(car, non, n_car, n_non, to_float(cx.get("cross_fisher_p")))
        label = "%s::%s" % ("/".join(genes_a[:2]) or "intergenic", "/".join(genes_b[:2]) or "intergenic")
        if junctions and not (genes_a and genes_b):
            label = junctions[0]["genes"]
        span = to_int(t["span"]) if ch_a == ch_b else None
        out.append({
            "jid": jid, "stratum": t["stratum"], "cls": t["cls"], "chrA": ch_a, "posA": pa, "chrB": ch_b, "posB": pb, "span": span,
            "genes_a": genes_a, "genes_b": genes_b, "label": label, "n_car_rna": n_car, "n_non_rna": n_non,
            "cross_car": car, "cross_non": non, "cross_other": to_int(cx["cross_other"], 0),
            "car_rate": car_rate, "non_rate": non_rate, "fold": fold, "fisher_p": p,
            "reads_car": reads_car, "reads_non": reads_non, "reads_note": reads_note,
            "rna_junctions": junctions[:8], "n_rna_junctions": len(junctions),
            "frame": frame, "frame_note": frame_note, "consequence": consequence,
            "consequence_note": None if consequence else "no retained domains reported by Arriba (frame/consequence need an annotated Arriba call)",
            "confidence": tier, "confidence_reasons": reasons,
            "igv_locus_a": locus(ch_a, pa), "igv_locus_b": locus(ch_b, pb)})
    out.sort(key=lambda e: ({"consistent": 0, "suggestive": 1, "low": 2}[e["confidence"]], -e["cross_car"]))
    return {"format": BT_FORMAT, "patient": patient, "generated": generated, "source": SOURCE, "rule": RULE,
            "n_sv": len(out), "transcripts": out}


EGFRVIII = {"patient": "BWH69", "jid": "207007759", "cdna_jid": "207007575"}


def amplicon_overlap(analysis_dir, patient, chrom, a, b):
    path = os.path.join(analysis_dir, "work", "amplicons.tsv")
    if not os.path.isfile(path):
        return None
    for r in rd(path):
        if r.get("patient") == patient and str(r.get("chr")) == str(chrom) and to_int(r["start"]) <= b and to_int(r["end"]) >= a:
            return True
    return False


def build_egfrviii(analysis_dir, generated):
    """Entry for the BWH69 genomic EGFRvIII deletion, or None when work/egfrviii_del.tsv is absent."""
    p = os.path.join(analysis_dir, "work", "egfrviii_del.tsv")
    if not os.path.isfile(p):
        return None
    rows = [r for r in rd(p) if r["patient"] == EGFRVIII["patient"] and r["jid"] == EGFRVIII["jid"]]
    if not rows:
        return None
    r = rows[0]
    n_car, n_non = to_int(r["n_car_rna"]), to_int(r["n_non_rna"])
    car, non = to_int(r["vIII_car"], 0), to_int(r["vIII_non"], 0)
    all_v3, tum = to_int(r["vIII_cells_all_rna"], 0), to_int(r["vIII_cells_tumour_paired"], 0)
    tier, reasons, car_rate, non_rate, fold, pv = confidence(car, non, n_car, n_non)
    if not n_non:
        reasons = ["%d of %d DNA-carrier tumour cells with RNA show the EGFRvIII junction; there are no non-carrier tumour cells with RNA "
                   "(every paired tumour cell carries the deletion), so no carrier vs non-carrier comparison is possible" % (car, n_car),
                   "counts: carrier cells >= %d; no p value" % MIN_CONS if car >= MIN_CONS else "counts: fewer than %d carrier cells" % MIN_CONS]
    nb = []
    nbp = os.path.join(analysis_dir, "work", "nearby_genomic_sv.tsv")
    dna = {}
    if os.path.isfile(nbp):
        for x in rd(nbp):
            if x["patient"] == EGFRVIII["patient"] and x["cdna_jid"] == EGFRVIII["cdna_jid"]:
                if x["g_jid"] == EGFRVIII["jid"]:
                    dna = {"n_cells": to_int(x["g_cells"]), "n_reads": to_int(x["g_reads"]), "tier_v2": to_int(x["g_tier_v2"]),
                           "cdna_cells": to_int(x["cdna_cells"]), "cdna_reads": to_int(x["cdna_reads"])}
                if to_int(x["g_cells"], 0) >= 10:
                    nb.append({"jid": x["g_jid"], "cls": x["g_cls"], "posA": to_int(x["g_lo"]), "posB": to_int(x["g_hi"]),
                               "n_cells": to_int(x["g_cells"]), "n_reads": to_int(x["g_reads"]), "tier_v2": to_int(x["g_tier_v2"])})
    nb.sort(key=lambda e: -e["n_cells"])
    amp = amplicon_overlap(analysis_dir, "BWH69", "7", 55022984, 55155076)
    cons = None
    cp = os.path.join(analysis_dir, "tables", "consequences.tsv")
    if os.path.isfile(cp):
        for c in rd(cp):
            if c["gene"] == "EGFR" and c["start"] == "55019366" and c["end"] == "55155829" and "BWH69" in c["patients"]:
                cons = "EGFRvIII exon 1 joined to exon 8: %s, %s, predicted protein %s aa vs %s (canonical transcript simulation, consequences.tsv)" % (
                    c["effect"], c["frame"].replace("_", "-"), c["protein_len"], c["ref_protein_len"])
    return {
        "jid": EGFRVIII["jid"], "stratum": "amplicon" if amp else ("non-amplicon" if amp is False else None), "cls": "DEL-like",
        "chrA": "7", "posA": 55022984, "chrB": "7", "posB": 55155076, "span": 55155076 - 55022984,
        "genes_a": ["EGFR"], "genes_b": ["EGFR"], "label": "EGFRvIII (EGFR exon 1 -> exon 8 deletion)",
        "n_car_rna": n_car, "n_non_rna": n_non, "cross_car": car, "cross_non": non, "cross_other": max(0, all_v3 - tum),
        "car_rate": car_rate, "non_rate": non_rate, "fold": fold, "fisher_p": pv,
        "reads_car": to_int(r["vIII_reads_car"]), "reads_non": None if not n_non else to_int(r["vIII_reads_non"]),
        "reads_note": None if n_non else "no non-carrier cells with RNA for this SV",
        "rna_junctions": [{"junction": "7:55019366-7:55155829", "source": "regtools", "genes": "EGFR::EGFR", "n_cells": all_v3,
                           "reading_frame": "in-frame"}],
        "n_rna_junctions": 1, "frame": "in-frame", "frame_note": None, "consequence": cons,
        "consequence_note": None if cons else "consequence table row not found",
        "confidence": tier, "confidence_reasons": reasons,
        "igv_locus_a": locus("7", 55022984), "igv_locus_b": locus("7", 55155076),
        "from_sv_transcribed": False,
        "added_reason": ("Genomic EGFRvIII deletion called in BWH69 DNA (chr7:55,022,984(-) -> 55,155,076(+), intron 1 -> intron 7). It is not in "
                         "sv_transcribed.tsv because every paired tumour cell carries it (no non-carrier group); the RNA evidence is the EGFR "
                         "exon 1 -> exon 8 junction (chr7:55,019,366-55,155,829; %d cells in all, %d in tumour cells with paired DNA)." % (all_v3, tum)),
        "dna_evidence": dict(dna, junction="chr7:55022984(-) -> chr7:55155076(+)", splice_mimic=False),
        "cdna_junction": {"jid": EGFRVIII["cdna_jid"], "junction": "chr7:55019369 -> chr7:55155828 (exact exon 1 / exon 8 boundaries)",
                          "splice_mimic": True, "note": "this exact exon-boundary junction in the DNA libraries is cDNA, in addition to the genomic deletion"},
        "nearby_genomic_sv": nb[:30], "n_nearby_genomic_sv": len(nb)}


def build_fusion_isoforms(patient, pair_rows, iso_rows, fusions, generated):
    by_pair = defaultdict(list)
    for x in iso_rows:
        if x["patient"] == patient:
            by_pair[x["pair"]].append(x)
    fidx = defaultdict(list)
    for f in fusions:
        fidx[frozenset([f.get("breakpoint1"), f.get("breakpoint2")])].append(f["id"])
    pairs = []
    for pr in pair_rows:
        if pr["patient"] != patient:
            continue
        isos, ids = [], []
        for x in sorted(by_pair.get(pr["pair"], []), key=lambda y: -to_int(y["n_cells"], 0)):
            g1, g2 = ",".join(strip_dist(x["gene1"])) or "intergenic", ",".join(strip_dist(x["gene2"])) or "intergenic"
            e1 = x["exon1"] or x["site1"] or "?"
            e2 = x["exon2"] or x["site2"] or "?"
            k = frozenset(["chr%s:%s" % (x["chr1"], x["pos1"]), "chr%s:%s" % (x["chr2"], x["pos2"])])
            for fid in fidx.get(k, []):
                if fid not in ids:
                    ids.append(fid)
            fr = x["reading_frame"]
            isos.append({"exon1": x["exon1"] or None, "site1": x["site1"] or None, "exon2": x["exon2"] or None, "site2": x["site2"] or None,
                         "pos1": to_int(x["pos1"]), "pos2": to_int(x["pos2"]), "chr1": x["chr1"], "chr2": x["chr2"],
                         "label": "%s %s -> %s %s" % (g1, e1, g2, e2), "reading_frame": None if fr in ("", ".") else fr,
                         "n_cells": to_int(x["n_cells"], 0), "n_cells_dnalinked": to_int(x["n_cells_dnalinked"], 0),
                         "total_reads": to_int(x["total_reads"], 0),
                         "retained_domains": None if x["retained_domains"] in ("", ".") else x["retained_domains"].lstrip("|"),
                         "source": x["source"]})
        g = pr["pair"].split("::")
        pairs.append({"pair": pr["pair"], "gene1": g[0], "gene2": g[1] if len(g) > 1 else None, "fusion_ids": ids,
                      "n_isoforms": to_int(pr["n_isoforms"], 0), "n_isoforms_2cells": to_int(pr["n_isoforms_2cells"], 0),
                      "n_isoforms_splice_both": to_int(pr["n_isoforms_splice_both"], 0), "total_cells": to_int(pr["total_cells"], 0),
                      "any_inframe": pr["any_inframe"].upper() == "TRUE", "readthrough": pr["readthrough"].upper() == "TRUE",
                      "top_isoform": pr["top_isoform"], "domains_top": None if pr["domains_top"] in ("", ".") else pr["domains_top"].lstrip("|"),
                      "isoforms": isos})
    pairs.sort(key=lambda e: (-e["n_isoforms_2cells"], -e["total_cells"]))
    return {"format": "fusion_isoforms_v1", "patient": patient, "generated": generated, "source": SOURCE,
            "n_pairs": len(pairs), "pairs": pairs}


def match_fusion(row, fusions):
    """-> (fusion row or None, note). Match on the unordered breakpoint pair, then the gene pair, then the most cells."""
    k = frozenset(["chr%s:%s" % (row["chr1"], row["pos1"]), "chr%s:%s" % (row["chr2"], row["pos2"])])
    if not fusions:
        return None, "no fusions.json rows for this patient (file absent or empty)"
    cand = [f for f in fusions if frozenset([f.get("breakpoint1"), f.get("breakpoint2")]) == k]
    if not cand:
        return None, "no fusions.json row with breakpoints %s" % " / ".join(sorted(k))
    note = None
    if len(cand) > 1:
        genes = {row["gene1"], row["gene2"]}
        same = [f for f in cand if {f.get("gene1"), f.get("gene2")} == genes]
        note = "%d fusions.json rows share the breakpoint pair" % len(cand)
        if same:
            cand = same
        cand = sorted(cand, key=lambda f: -(f.get("n_cells") or 0))
        note += "; chose the gene-pair match with most cells" if same else "; gene pair did not match, chose most cells"
    return cand[0], note


def build_fusion_support(patient, only_rows, fusions, sv_pos, crossing_by_jid, generated):
    entries = []
    for x in only_rows:
        f, mnote = match_fusion(x, fusions)
        p1, p2 = to_int(x["pos1"]), to_int(x["pos2"])
        best = None
        for jid in [j for j in re.split(r"[,;]", x["jids"]) if j]:
            sv = sv_pos.get(jid)
            if not sv:
                continue
            ends = [(sv["chrA"], sv["posA"]), (sv["chrB"], sv["posB"])]
            for e1, e2 in (ends, ends[::-1]):
                if e1[0] != str(x["chr1"]) or e2[0] != str(x["chr2"]):
                    continue
                d1, d2 = abs(p1 - e1[1]), abs(p2 - e2[1])
                if best is None or d1 + d2 < best[0]:
                    best = (d1 + d2, jid, d1, d2)
        if best is None:
            jid, dist, d1, d2 = (re.split(r"[,;]", x["jids"])[0] or None), None, None, None
        else:
            _, jid, d1, d2 = best
            dist = max(d1, d2)
        cx = crossing_by_jid.get(jid, {})
        cc, nc = to_int(x["carrier_cells"], 0), to_int(x["noncarrier_cells"], 0)
        tier, reasons, car_rate, non_rate, fold, p = confidence(cc, nc, to_int(cx.get("n_car_rna")), to_int(cx.get("n_non_rna")))
        cls = "canonical_splice" if str(x["splice_both"]).upper() == "TRUE" else "cryptic_end"
        if x.get("cls"):
            reasons.append("analysis flag: " + x["cls"])
        entries.append({
            "fusion_id": f["id"] if f else None, "match_note": mnote,
            "gene1": x["gene1"], "gene2": x["gene2"],
            "breakpoint1": "chr%s:%s" % (x["chr1"], x["pos1"]), "breakpoint2": "chr%s:%s" % (x["chr2"], x["pos2"]),
            "sv_jid": jid, "distance_bp": dist, "distance_bp_ends": [d1, d2] if d1 is not None else None,
            "class": cls, "class_label": CLASS_LABEL[cls], "explanation": x.get("explanation") or None,
            "n_cells": to_int(x["n_cells"]), "carrier_cells": cc, "noncarrier_cells": nc, "other_cells": to_int(x["other_cells"], 0),
            "n_car_rna": to_int(cx.get("n_car_rna")), "n_non_rna": to_int(cx.get("n_non_rna")),
            "car_rate": car_rate, "non_rate": non_rate, "fold": fold, "fisher_p": p,
            "confidence": tier, "confidence_reasons": reasons})
    entries.sort(key=lambda e: ({"consistent": 0, "suggestive": 1, "low": 2}[e["confidence"]], -e["carrier_cells"]))
    return {"format": FS_FORMAT, "patient": patient, "generated": generated, "source": SOURCE, "rule": RULE,
            "n": len(entries), "entries": entries}


def run(analysis_dir, gos_data_dir, out_dir, patients=None):
    A = analysis_dir
    cross = rd(os.path.join(A, "work", "sv_crossing.tsv"))
    crossing_by_jid = {r["jid"]: r for r in cross}
    sv_pos = {r["jid"]: {"chrA": r["chrA"], "posA": to_int(r["posA"]), "chrB": r["chrB"], "posB": to_int(r["posB"])}
              for r in rd(os.path.join(A, "work", "sv_junctions.tsv"))}
    trans = rd(os.path.join(A, "tables", "sv_transcribed.tsv"))
    only = rd(os.path.join(A, "work", "rna_only_vs_sv.tsv"))
    pair_rows = rd(os.path.join(A, "tables", "fusion_pairs.tsv"))
    iso_rows = rd(os.path.join(A, "work", "fusion_isoforms.tsv"))
    egfr = build_egfrviii(A, time.strftime("%Y-%m-%d"))
    need = {t["jid"] for t in trans}
    reads_by_jid = defaultdict(list)
    for r in rd(os.path.join(A, "work", "sv_crossing_reads.tsv.gz")):
        if r["jid"] in need:
            reads_by_jid[r["jid"]].append(r)
    if not patients:
        pats = {t["patient"] for t in trans} | {x["patient"] for x in only}
        pats |= {d for d in os.listdir(gos_data_dir) if os.path.isfile(os.path.join(gos_data_dir, d, "rna", "fusions.json"))}
        patients = sorted(p for p in pats if os.path.isdir(os.path.join(gos_data_dir, p, "rna")))
    generated = time.strftime("%Y-%m-%d")
    summary = {}
    for p in patients:
        fpath = os.path.join(gos_data_dir, p, "rna", "fusions.json")
        if os.path.isfile(fpath):
            with open(fpath) as fh:
                fusions = json.load(fh).get("fusions", [])
        else:
            fusions = []
        bt = build_transcripts(p, [t for t in trans if t["patient"] == p], reads_by_jid, crossing_by_jid, generated)
        if egfr is not None and p == EGFRVIII["patient"]:
            bt["transcripts"].append(egfr)
            bt["transcripts"].sort(key=lambda e: ({"consistent": 0, "suggestive": 1, "low": 2}[e["confidence"]], -e["cross_car"]))
            bt["n_sv"] = len(bt["transcripts"])
        fi = build_fusion_isoforms(p, pair_rows, iso_rows, fusions, generated)
        fs = build_fusion_support(p, [x for x in only if x["patient"] == p], fusions, sv_pos, crossing_by_jid, generated)
        od = os.path.join(out_dir, p, "rna")
        os.makedirs(od, exist_ok=True)
        for name, obj in (("breakpoint_transcripts.json", bt), ("fusions_dna_support.json", fs), ("fusion_isoforms.json", fi)):
            with open(os.path.join(od, name), "w") as fh:
                json.dump(obj, fh, indent=1, allow_nan=False)
        summary[p] = {"n_sv": bt["n_sv"], "bt_conf": dict(Counter(e["confidence"] for e in bt["transcripts"])),
                      "n_fusions": fs["n"], "matched": sum(1 for e in fs["entries"] if e["fusion_id"]),
                      "fs_conf": dict(Counter(e["confidence"] for e in fs["entries"])),
                      "n_pairs": fi["n_pairs"], "n_pairs_ge2": sum(1 for e in fi["pairs"] if e["n_isoforms_2cells"] >= 2)}
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--analysis-dir", required=True)
    ap.add_argument("--gos-data-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--patients", default="", help="comma-separated (default: every patient with rna/fusions.json)")
    a = ap.parse_args(argv)
    s = run(a.analysis_dir, a.gos_data_dir, a.out_dir, [x for x in a.patients.split(",") if x])
    tot = Counter()
    for p, v in s.items():
        print(p, json.dumps(v, sort_keys=True))
        tot["sv"] += v["n_sv"]
        tot["fus"] += v["n_fusions"]
        tot["matched"] += v["matched"]
        tot["pairs"] += v["n_pairs"]
        tot["pairs_ge2"] += v["n_pairs_ge2"]
    print("TOTAL", dict(tot))
    return 0


if __name__ == "__main__":
    sys.exit(main())
