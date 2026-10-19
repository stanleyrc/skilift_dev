#!/usr/bin/env python3
"""Clinical-assay alterations as filtered.events-style rows, so the gOS Drivers views can show every clinical finding.

Inputs : data/_cohort/{clinical_assays,clinical_sites,clinical_recovery}.json (gos_sc_clinical*.py),
         data/<P>/filtered.events.json, per-cell <cell>/complex.json (CNV / SV carriers),
         targeted pooled pileups analysis/clinical_recovery_261017/<P>/s*/ (SNV carriers).
Output : data/<P>/clinical_events.json  { format, patient, tumour_cells, summary, events: [row, ...] }
One row per clinical alteration (SNV / indel, CNV gene, SV, TERT promoter hotspot), real counts only. Row fields follow
filtered.events.json (gene, vartype, type, Variant, Genome_Location, seqnames, start, end, n_cells, cell_ids, cells, VAF ...) plus
  source "clinical assay", clinical {assay, tier, raw, status, reason, vaf, germline_like, normal_alt, evidence ...}, germline_like,
  clinical_key (unique id), matches_event (the driverKey of the filtered event that already holds this alteration, else null).
Not-recovered / not-assessable alterations still get a row (no carriers, status says why). Negative findings (the report says
"no amplification", TERT site absent) are not alterations and get no row (counted in summary.negatives).
usage: gos_sc_clinical_events.py [--data DIR] [--patients A,B]
"""
import argparse
import datetime
import json
import os
import re
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gos_sc_clinical_recovery as rec  # noqa: E402
import gos_sc_clinical_annot as ann  # noqa: E402

FORMAT = "gos-sc-clinical-events/1"
SOURCE = "clinical assay"
COUNT_STATUS = ("recovered", "partial", "not_recovered", "not_covered", "not_assessable")


# ---------- pure helpers (unit tested) ----------
def event_key(e):
    """Same string as gOS driverKey(): label|Genome_Location."""
    fg = e.get("fusion_genes") or ""
    kind = e.get("vartype") or e.get("type") or ""
    label = fg if "::" in fg else "%s%s" % (e.get("gene") or "?", " " + kind if kind else "")
    loc = e.get("Genome_Location") or "%s:%s-%s" % (e.get("seqnames"), e.get("start"), e.get("end"))
    return "%s|%s" % (label, loc)


def snv_type(protein, cdna):
    """-> (type, vartype) in the filtered.events vocabulary: Trunc / Splice / Missense / Noncoding, SNV or INDEL."""
    prot, c = protein or "", cdna or ""
    vt = "INDEL" if re.search(r"del|ins|dup|fs", prot + " " + c) else "SNV"
    if "*" in prot or "fs" in prot or "Ter" in prot:
        return "Trunc", vt
    if prot and prot not in ("p.?", "p.="):
        return "Missense", vt
    m = re.search(r"\d([+-])(\d+)(?!\d)", c)
    if re.search(r"_\d+[+-]\d+", c) or (m and m.group(1) in "+-" and int(m.group(2)) <= 2 and not c.startswith("c.-")):
        return "Splice", vt
    return "Noncoding", vt


def location(chrom, start, end):
    ch = str(chrom).replace("chr", "")
    return "%s:%d-%d" % (ch, start, end), ch


def ids_str(ids):
    return ",".join(sorted(ids))


def base_row(p, kind, gene, vartype, typ, variant, loc, ch, start, end, carriers, covered, n_tumour, clinical, key, vaf=None):
    n = len(carriers)
    row = {"gene": gene, "fusion_genes": "", "vartype": vartype, "type": typ, "Variant": variant,
           "Genome_Location": loc, "seqnames": ch, "start": str(start) if start is not None else None, "end": str(end) if end is not None else None,
           "id": p, "alt": 0, "ref": 0, "n_cells": n, "n_normal_cells": 0, "cell_ids": ids_str(carriers),
           "cells": "%d/%d" % (n, int(covered) if covered is not None else n_tumour), "cell_fraction": round(n / int(covered), 4) if covered else None,
           "tumour_cells": n_tumour, "covered_tumour_cells": covered, "source": SOURCE, "clinical_kind": kind,
           "clinical_key": key, "matches_event": None, "germline_like": bool(clinical.get("germline_like")), "clinical": clinical}
    if vaf is not None:
        row["VAF"] = vaf
    return row


def find_snv_event(events, gene, protein, chrom, start, end):
    """Filtered event holding this SNV / indel (same protein change, or position inside the site)."""
    c = (chrom or "").replace("chr", "")
    for e in events:
        if e.get("gene") != gene or e.get("type") in ("SCNA", "Fusion", "Complex SV") or e.get("source") == SOURCE:
            continue
        if protein and e.get("Variant") == protein:
            return e
        if c and start is not None and str(e.get("Genome_Location", "")).startswith(c + ":") and e.get("vartype") in ("SNV", "INDEL", "DEL", "INS"):
            try:
                pos = int(str(e["Genome_Location"]).split(":")[1].split("-")[0])
            except ValueError:
                continue
            if start <= pos <= end:
                return e
    return None


def find_cn_event(events, gene, cls):
    """amp / gain match an AMP event of the gene, deep matches HOMDEL; loss has no filtered-event counterpart."""
    want = {"amp": "AMP", "gain": "AMP", "deep": "HOMDEL"}.get(cls)
    for e in events:
        if want and e.get("type") == "SCNA" and e.get("gene") == gene and e.get("vartype") == want:
            return e
    return None


def cn_vartype(cls):
    return "AMP" if cls in ("amp", "gain") else "HOMDEL"


def germline_text(c):
    return "also in normal cells (%s of %s reads): germline-like" % (c.get("normal_alt"), c.get("normal_dp"))


def summarize(rows):
    s = {"findings": len(rows), "germline_like": sum(1 for r in rows if r["germline_like"]), "linked_to_filtered_event": sum(1 for r in rows if r["matches_event"])}
    for k in COUNT_STATUS:
        s[k] = sum(1 for r in rows if r["clinical"]["status"] == k)
    return s


# ---------- per patient ----------
def annotate(row, site, pe, vcf, release):
    """ClinVar + gnomAD (local) and the highlight tags of one SNV / indel row; row['clinical'] gets clinvar, gnomad, lof, highlight, note."""
    clin = row['clinical']
    chrom = (site or {}).get('chrom')
    lof = ann.is_lof(row['type'], row.get('cdna'))
    clin['lof'] = lof
    if chrom and pe and pe.get('ref') and pe.get('alt') and site.get('kind') in ('SNV', 'SNV_PROTEIN'):
        pos = pe['pos']
        clin['clinvar'] = ann.clinvar_lookup(vcf, release, chrom, pos, pe['ref'], pe['alt']) if vcf else None
        clin['gnomad'] = ann.gnomad_lookup(chrom, pos, pe['ref'], pe['alt'])
        clin['alleles'] = '%s:%d %s>%s (GRCh38)' % (chrom, pos, pe['ref'], pe['alt'])
    elif chrom and site.get('kind') == 'DEL' and vcf:
        clin['clinvar'] = ann.clinvar_lookup(vcf, release, chrom, site['start'], None, None, indel_len=site['end'] - site['start'] + 1)
        clin['gnomad'] = {'found': False, 'text': 'not looked up (deletion alleles are not given in the report)', 'release': ann.GNOMAD_RELEASE}
    else:
        clin['clinvar'] = None if not vcf else ann.clinvar_none(release) if not chrom else {'found': False, 'release': release, 'text': 'not looked up (complex indel or no genomic position)'}
        clin['gnomad'] = {'found': False, 'text': 'not looked up (complex indel or no genomic position)', 'release': ann.GNOMAD_RELEASE}
    clin['highlight'] = ann.highlight_tags(clin.get('clinvar'), lof, row['germline_like'])
    note = ann.recessive_note(clin.get('clinvar'), lof)
    if note:
        clin['recessive_note'] = note


def snv_rows(p, rep, srows, events, obs, depth, n_tumour, vcf=None, release=None):
    rows = []
    for k, v in enumerate(rep["snv"]):
        site = srows[k] if k < len(srows) else {}
        clin = {kk: v.get(kk) for kk in ("assay", "tier", "clinical_vaf", "clinical_reads", "status", "reason", "outcome", "note", "germline_like",
                                       "normal_alt", "normal_dp", "sc_vaf", "sc_carriers", "sc_covered", "evidence", "source")}
        clin["raw"] = v["clinical_raw"]
        clin["kind"] = "SNV"
        typ, vt = snv_type(v.get("protein"), v.get("cdna"))
        variant = v.get("protein") or v.get("cdna") or "None"
        key = "clinical|%s|%s|%s" % (p, v["gene"], v.get("cdna") or v.get("protein") or k)
        st = v.get("site") or {}
        chrom, start, end = st.get("chrom"), st.get("resolved_pos") or st.get("start"), st.get("end") or st.get("start")
        if chrom and st.get("resolved_pos"):
            end = start + (st.get("end", start) - st.get("start", start))
        loc, ch = location(chrom, start, end) if chrom else (None, None)
        carriers, pe = [], None
        ev = find_snv_event(events, v["gene"], v.get("protein"), chrom, st.get("start"), st.get("end"))
        if site.get("ok") and obs is not None and depth:
            pe = rec.pileup_eval(site, obs, depth)
            if pe:
                carriers = pe.get("carrier_ids", [])
                if len(carriers) != v.get("sc_carriers"):
                    print("WARNING %s %s: pileup carriers %d != recovery %s" % (p, v["gene"], len(carriers), v.get("sc_carriers")), file=sys.stderr)
        if not carriers and ev is not None:
            carriers = [c for c in str(ev.get("cell_ids") or "").split(",") if c]
        covered = v.get("sc_covered")
        row = base_row(p, "snv", v["gene"], vt, typ, variant, loc, ch, start, end, carriers, covered, n_tumour, clin, key, vaf=v.get("sc_vaf"))
        if v.get("cdna"):
            row["cdna"] = v["cdna"]
        if pe is None and site.get('ok') and site.get('kind') == 'SNV':
            pe = {'pos': site['start'], 'ref': site['ref'], 'alt': site['alt']}
        annotate(row, site if site.get('ok') else {}, pe, vcf, release)
        if ev is not None:
            row["matches_event"] = event_key(ev)
            inherit_clonality(row, ev)
        rows.append(row)
    return rows


def inherit_clonality(row, ev):
    """Copy the matched filtered event's dropout-aware clonality (clonality.py) onto a clinical row; nothing when it has none."""
    if ev.get("clonality_method") == "dropout-aware":
        row.update(clonality=ev.get("clonality"), clonality_method=ev["clonality_method"], clonality_stats=ev.get("clonality_stats"),
                   clonality_inherited_from=event_key(ev))
    return row


def cnv_rows(p, rep, events, coords, tcells, n_tumour):
    """One gene-level CN row per listed gene; carriers = tumour cells whose own CN at the gene shows the clinical class."""
    rows, neg = [], 0
    for c in rep["cnv"]:
        if c.get("negative_finding"):
            neg += 1
            continue
        cls = rec.clinical_class(c["type"])
        for g in c["genes"]:
            clin = {"assay": c["assay"], "kind": "CNV", "raw": c["clinical_raw"], "cnv_type": c["type"], "scope": c["scope"], "location": c.get("location"),
                    "clinical_copies": c.get("clinical_copies"), "event_genes": len(c["genes"])}
            key = "clinical|%s|%s|%s" % (p, g, c["type"])
            carriers, vals = [], []
            if g in coords and cls and tcells:
                ch, s, e = coords[g]
                for cid, segs, base, alt in tcells:
                    cn = rec.cn_at(segs, ch, s, e)
                    if cn is None:
                        continue
                    vals.append(cn)
                    st = rec.cell_state(cn, base)
                    if st and cls in st:
                        carriers.append(cid)
                frac = len(carriers) / len(vals) if vals else None
                status = rec.cnv_status(frac)
                clin.update(status=status, reason="%d of %d tumour cells with %s at %s (median CN %s)" % (len(carriers), len(vals), cls, g, statistics.median(vals) if vals else None),
                            median_cn=statistics.median(vals) if vals else None)
                loc, chs = location(ch, s, e)
                start, end = s, e
            else:
                clin.update(status="not_assessable", reason="no tumour genome graphs" if not tcells else "gene coordinates or alteration class unresolved (not available)")
                loc = chs = start = end = None
            row = base_row(p, "cnv", g, cn_vartype(cls), "SCNA", cn_vartype(cls), loc, chs, start, end, carriers, len(vals) or None, n_tumour, clin, key)
            row["cn_kind"] = cls
            ev = find_cn_event(events, g, cls)
            if ev is not None:
                row["matches_event"] = event_key(ev)
            rows.append(row)
    return rows, neg


def sv_rows(p, rep, events, coords, tcells, n_tumour):
    rows = []
    for k, s in enumerate(rep["sv"]):
        clin = {"assay": s["assay"], "kind": "SV", "sv_kind": s["kind"], "tier": s.get("tier"), "raw": s["clinical_raw"], "status": s["status"], "reason": s["reason"],
                "evidence": [{"source": "recovery", "n_cell_junction_pairs": s.get("n_cells"), "nearest_bp_distance": s.get("nearest_bp_distance")}]}
        lf, rf = (s.get("left") or {}), (s.get("right") or {})
        g1, g2 = lf.get("feature", "").split(" ")[0], rf.get("feature", "").split(" ")[0]
        genes = [g for g in dict.fromkeys([g1, g2]) if g]
        fg = "::".join([g1, g2]) if g1 and g2 and g1 != g2 else ""
        carriers, loc, chs, start, end, ev = [], None, None, None, None, None
        if "chrom" in lf and "chrom" in rf:
            c1, p1 = lf["chrom"].replace("chr", ""), lf["pos"]
            c2, p2 = rf["chrom"].replace("chr", ""), rf["pos"]
            loc, chs = location(c1, p1, p1)
            start = end = p1
            if s["status"] in ("recovered", "partial"):
                carriers = [cid for cid, segs, base, alt in tcells if rec.junction_match(alt, (c1, p1), (c2, p2)) is not None]
        else:
            fus = [e for e in events if e.get("type") == "Fusion" and {g1, g2} <= rec.gene_tokens(e.get("fusion_genes"))]
            if fus:
                ev = fus[0]
                carriers = [c for c in str(ev.get("cell_ids") or "").split(",") if c]
            if g1 in coords:
                ch, s0, e0 = coords[g1]
                loc, chs = location(ch, s0, e0)
                start, end = s0, e0
                if not fus and g2 in coords and coords[g1][0] == coords[g2][0] and s["status"] == "partial":
                    mid = lambda g: (coords[g][1] + coords[g][2]) // 2
                    carriers = [cid for cid, segs, base, alt in tcells if rec.junction_match(alt, (ch, mid(g1)), (ch, mid(g2))) is not None]
        key = "clinical|%s|SV%d|%s" % (p, k, "-".join(genes) or s["kind"])
        row = base_row(p, "sv", fg or (genes[0] if genes else ""), "fusion", "Fusion", "%s %s" % (s["kind"], " :: ".join(genes)) if genes else s["kind"], loc, chs, start, end,
                       carriers, n_tumour if carriers or s["status"] in ("recovered", "partial", "not_recovered") else None, n_tumour, clin, key)
        row["fusion_genes"] = fg
        if ev is not None:
            row["matches_event"] = event_key(ev)
        rows.append(row)
    return rows


def tert_rows(p, rep, events, n_tumour):
    rows, neg = [], 0
    for t in rep["tert_promoter"]:
        if t["clinical"] != "present":
            neg += 1
            continue
        ev = next((e for e in events if e.get("gene") == "TERT" and e.get("vartype") == "Hotspot" and str(e.get("Variant", "")).startswith(t["site"])), None)
        clin = {"assay": t["assay"], "kind": "TERT promoter", "raw": t["clinical_raw"], "status": t["status"], "reason": t.get("reason"), "evidence": [
            {"source": t.get("source"), "carriers": t.get("carriers"), "covered": t.get("covered"), "alt_reads": t.get("alt"), "depth": t.get("dp")}]}
        carriers = [c for c in str((ev or {}).get("cell_ids") or "").split(",") if c]
        loc = (ev or {}).get("Genome_Location")
        row = base_row(p, "tert", "TERT", "Hotspot", "Promoter", t["site"], loc, (ev or {}).get("seqnames"), (ev or {}).get("start"), (ev or {}).get("end"),
                       carriers, t.get("covered"), n_tumour, clin, "clinical|%s|TERT|%s" % (p, t["site"]))
        if ev is not None:
            row["matches_event"] = event_key(ev)
            inherit_clonality(row, ev)
        rows.append(row)
    return rows, neg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=rec.DATA)
    ap.add_argument("--patients", default="")
    a = ap.parse_args()
    data = a.data
    recov = json.load(open(os.path.join(data, "_cohort", "clinical_recovery.json")))["patients"]
    sites = json.load(open(os.path.join(data, "_cohort", "clinical_sites.json")))["patients"]
    want = [x for x in a.patients.split(",") if x] or sorted(recov)
    genes = {g for p in want for c in recov[p]["cnv"] for g in c["genes"]} | {g for p in want for s in recov[p]["sv"] for g in
             [(s.get("left") or {}).get("feature", "").split(" ")[0], (s.get("right") or {}).get("feature", "").split(" ")[0]] if g}
    coords = rec.load_genes(genes)
    vcf = ann.clinvar_vcf()
    release = ann.clinvar_release(vcf) if vcf else None
    print('ClinVar', vcf, release)
    for p in want:
        rep = recov[p]
        events = json.load(open(os.path.join(data, p, "filtered.events.json")))
        cells = rec.tumour_cells(data, p)
        tum = [c for c, nrm in cells if not nrm]
        obs, depth = rec.read_pileup(p)
        tcells = []
        if rep["cnv"] or rep["sv"]:
            for c in tum:
                segs, base, alt = rec.cell_graph(os.path.join(data, c, "complex.json"))
                tcells.append((c, segs, base, alt))
        n_t = len(tum)
        rows = snv_rows(p, rep, sites.get(p, []), events, obs, depth, n_t, vcf, release)
        cr, neg_c = cnv_rows(p, rep, events, coords, tcells, n_t)
        tr, neg_t = tert_rows(p, rep, events, n_t)
        rows += cr + sv_rows(p, rep, events, coords, tcells, n_t) + tr
        out = {"format": FORMAT, "patient": p, "created": datetime.datetime.now().isoformat(timespec="seconds"), "tumour_cells": n_t,
               "clinical_status": rep["clinical_status"], "annotation": {"clinvar": {"file": vcf, "release": release}, "gnomad": ann.GNOMAD_RELEASE}, "summary": dict(summarize(rows), negatives=neg_c + neg_t), "events": rows}
        dest = os.path.join(data, p, "clinical_events.json")
        if os.path.exists(dest) and not os.path.exists(dest + ".bak_261017"):
            os.replace(dest, dest + ".bak_261017")
        json.dump(out, open(dest, "w"))
        print(p, out["summary"])


if __name__ == "__main__":
    main()
