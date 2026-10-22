#!/usr/bin/env python3
"""Confidence tier of every splicing call in the gOS single-cell files, from the cell-level re-test (differential-splicing
deep-dive, srctools inst/diff_splicing_deepdive, events_annotated.tsv: one row per cohort cluster with a verdict).

Tiers (no call is removed; every cluster / finding / event gets a tier and a reason):
  consistent  within-patient differential in the same direction in >= 3 patients (malignant vs non-malignant cells, or one
              malignant state vs the other malignant cells), a patient-specific junction explained by a somatic splice-site
              variant in cis, or a known GBM splice variant (EGFRvIII) tracked outside the cluster test.
  suggestive  tumour-to-tumour difference that survives the cell-level and library-batch tests but has no within-patient
              contrast; a within-patient contrast in one patient only; or a cluster dropped only for the BWH75 library-batch
              split whose within-patient contrast has the same direction in all of >= RESCUE_MIN_PATIENTS patients with
              q < RESCUE_MAX_Q (BIN1).
  low         dropped by the re-test (not significant once cells are the unit / explained by library batch / differs between
              the BWH75 batches / paralog-repeat gene), or not re-tested (cluster outside the re-tested cohort set).

Added fields:
  clusters[] and findings[] (patient + _cohort rna/splicing.json): confidence, confidence_reason, confidence_detail
      (verdict, contrast, effect statistics, consequence class from the deep-dive table, first_intron = alternative
      promoter / first exon).
  events[] (patient rna/splicing.json): confidence, confidence_reason (tier of the event's cluster).
  confidence_meta (each file): format gos-sc-splice-confidence/1, source table, parameters, counts per tier.
  _cohort/splice_support.json: confidence / confidence_reason on each RNA entry, rna_confidence (best tier) per gene.
  _cohort/vignettes.json: confidence / confidence_reason on the "splice" vignette items (by their finding id).

Usage:
  gos_sc_splice_confidence.py --data <gOS data dir> --deepdive <events_annotated.tsv> --out <staging dir>
  gos_sc_splice_confidence.py --data <gOS data dir> --deepdive <events_annotated.tsv> --in-place   (back up first)
With --out the data dir is only read; the annotated files are written below <out> at the same relative paths.
Idempotent: the fields are recomputed from the table each run."""
import argparse
import csv
import json
import os
import sys
import time

FORMAT = "gos-sc-splice-confidence/1"
TIERS = ("consistent", "suggestive", "low")
RANK = {t: i for i, t in enumerate(TIERS)}
RESCUE_MIN_PATIENTS = 4
RESCUE_MAX_Q = 0.1

NORMAL_LABEL = "malignant vs non-malignant cells (non-malignant: mostly oligodendrocytes, some OPC, myeloid and astrocytes)"
NOT_RETESTED = "not re-tested: outside the cohort clusters of the cell-level re-test"

# protein_consequence of the deep-dive table (03_consequence.R) -> short class
CONSEQUENCE = [
    ("known alternative protein-coding isoform", "known isoform switch"),
    ("annotated protein-coding junction, but on the canonical frame", "alternative start (annotated isoform)"),
    ("known NMD isoform", "known NMD isoform"),
    ("novel junction: frameshift / PTC", "novel, NMD-type (frameshift / PTC)"),
    ("novel junction: in-frame change", "novel in-frame (candidate new isoform)"),
    ("novel junction: altered protein that escapes", "novel, escapes NMD (altered protein)"),
    ("novel junction in UTR", "novel UTR junction"),
    ("junction of an annotated retained-intron", "retained-intron / non-coding transcript"),
    ("not resolved", "not resolved"),
]


def num(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if x == x else None  # NaN -> None


def inum(v):
    x = num(v)
    return int(x) if x is not None else None


def consequence_class(text):
    t = (text or "").strip()
    for prefix, label in CONSEQUENCE:
        if t.startswith(prefix):
            return label
    return t or None


def fmt(x, d=2, signed=False):
    if x is None:
        return "n/a"
    return f"{x:+.{d}f}" if signed else f"{x:.{d}f}"


def fmt_p(x):
    if x is None:
        return "n/a"
    if x == 0:
        return "<1e-300"
    return f"{x:.1e}" if x < 1e-3 else f"{x:.2g}"


def contrast_stats(row, pre):
    """Within-patient contrast statistics of the table (pre = 'tn' malignant vs normal, 'st' state vs other malignant)."""
    return {
        "k": inum(row.get(f"{pre}_k")),
        "dpsi": num(row.get(f"{pre}_dpsi")),
        "lo": num(row.get(f"{pre}_lo")),
        "hi": num(row.get(f"{pre}_hi")),
        "q": num(row.get(f"{pre}_q")),
        "perm_p": num(row.get(f"{pre}_perm_p")),
        "same_dir": inum(row.get(f"{pre}_same_dir")),
        "glmm_p": num(row.get(f"glmm_p_{pre}")),
        "i2": num(row.get(f"{pre}_I2")),
    }


def stats_text(s):
    return (f"ΔPSI {fmt(s['dpsi'], signed=True)} (95% CI {fmt(s['lo'])} to {fmt(s['hi'])}), same direction in "
            f"{s['same_dir'] if s['same_dir'] is not None else 'n/a'}/{s['k'] if s['k'] is not None else 'n/a'} patients, "
            f"q {fmt_p(s['q'])}")


def classify(row, rescue_min_patients=RESCUE_MIN_PATIENTS, rescue_max_q=RESCUE_MAX_Q):
    """{tier, reason, detail} of one deep-dive table row."""
    verdict = (row.get("verdict") or "").strip()
    v = verdict.lower()
    tn = contrast_stats(row, "tn")
    st = contrast_stats(row, "st")
    state = (row.get("st_best") or "").strip()
    detail = {
        "verdict": verdict,
        "contrast": None,
        "stats": None,
        "consequence": consequence_class(row.get("protein_consequence")),
        "consequence_text": (row.get("protein_consequence") or "").strip() or None,
        "annotation_class": (row.get("annotation_class") or "").strip() or None,
        "effect": (row.get("effect") or "").strip() or None,
        "frame": (row.get("frame") or "").strip().strip('"') or None,
        "nmd": {"TRUE": True, "FALSE": False}.get((row.get("nmd_50nt") or "").strip()),
        "protein_len": inum(row.get("protein_len")),
        "ref_protein_len": inum(row.get("ref_protein_len")),
        "pfam_affected": (row.get("pfam_affected") or "").strip().strip('"') or None,
        "first_intron": (row.get("terminal_intron") or "").strip() == "first",
        "cis_dna": (row.get("cis_dna") or "").strip() or None,
        "cell_q": num(row.get("patient_qb_q")),
        "batch_q": num(row.get("patient_given_batch_q")),
        "bwh75_batch_dpsi": num(row.get("bwh75_batch_dpsi")),
    }
    if v.startswith("keep"):
        if "in cis" in v:
            detail["contrast"] = "DNA-explained: somatic splice-site variant in cis"
            reason = f"patient-specific junction explained by DNA ({detail['cis_dna'] or 'splice-site variant in cis'})"
        elif "cell states" in v or "malignant cell states" in v:
            detail["contrast"] = f"{state or 'one'}-like vs other malignant cells"
            detail["stats"] = st
            reason = f"consistent within patients: {detail['contrast']}, {stats_text(st)}"
        else:
            detail["contrast"] = NORMAL_LABEL
            detail["stats"] = tn
            reason = f"consistent within patients: {NORMAL_LABEL}, {stats_text(tn)}"
        return {"tier": "consistent", "reason": reason, "detail": detail}
    if v.startswith("descriptive"):
        return {"tier": "suggestive", "detail": detail,
                "reason": "tumour-to-tumour difference: survives the cell-level and library-batch tests, no within-patient contrast"}
    if v.startswith("tentative"):
        detail["stats"] = tn if tn["k"] else st
        return {"tier": "suggestive", "detail": detail, "reason": "within-patient contrast in one patient only"}
    if v.startswith("drop"):
        if "bwh75" in v and tn["k"] is not None and tn["k"] >= rescue_min_patients and tn["same_dir"] == tn["k"] \
                and tn["q"] is not None and tn["q"] < rescue_max_q:
            detail["contrast"] = NORMAL_LABEL
            detail["stats"] = tn
            return {"tier": "suggestive", "detail": detail,
                    "reason": (f"within patients {NORMAL_LABEL}: {stats_text(tn)}; batch-sensitive (differs by "
                               f"{fmt(detail['bwh75_batch_dpsi'], signed=True)} between the BWH75 library batches)")}
        if "not significant" in v:
            why = "dropped: not significant once cells (not reads) are the unit"
        elif "explained by library batch" in v:
            why = "dropped: between-patient difference explained by library batch"
        elif "bwh75" in v:
            why = f"dropped: differs between the BWH75 library batches (ΔPSI {fmt(detail['bwh75_batch_dpsi'], signed=True)})"
        elif "paralog" in v or "repeat" in v:
            why = "dropped: paralog / repeat gene family (mapping artefact)"
        else:
            why = f"dropped: {verdict.split(':', 1)[-1].strip()}"
        return {"tier": "low", "reason": why, "detail": detail}
    return {"tier": "low", "reason": f"unrecognised verdict: {verdict or 'none'}", "detail": detail}


def read_table(path, **kw):
    """cluster_id -> classify(row)."""
    out = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            cid = (row.get("cluster_id") or "").strip()
            if cid:
                out[cid] = classify(row, **kw)
    return out


def known_variant_tier(finding):
    n = finding.get("numbers") or {}
    alt, cov = n.get("n_cells_alt"), n.get("n_cells")
    cells = f"; {alt}/{cov} covered cells carry it" if alt is not None and cov is not None else ""
    return {"tier": "consistent",
            "reason": f"known GBM splice variant ({finding.get('variant_id') or finding.get('gene')}), tracked outside the cluster test{cells}",
            "detail": {"verdict": "known variant", "contrast": "known variant", "consequence": None, "first_intron": False}}


def not_retested():
    return {"tier": "low", "reason": NOT_RETESTED, "detail": {"verdict": "not re-tested", "contrast": None, "first_intron": False}}


def tier_of_finding(f, table):
    if f.get("kind") == "known_variant":
        return known_variant_tier(f)
    return table.get(f.get("cluster_id")) or not_retested()


def put(obj, c, detail=True):
    obj["confidence"] = c["tier"]
    obj["confidence_reason"] = c["reason"]
    if detail:
        obj["confidence_detail"] = c.get("detail")
    else:
        obj.pop("confidence_detail", None)


def counts(rows):
    out = {t: 0 for t in TIERS}
    for r in rows:
        out[r.get("confidence", "low")] = out.get(r.get("confidence", "low"), 0) + 1
    return out


def meta(source, params):
    try:
        mtime = time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(os.path.getmtime(source)))
    except OSError:
        mtime = None
    return {"format": FORMAT, "source": source, "source_mtime": mtime, "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "tiers": list(TIERS), "params": params,
            "tier_help": {
                "consistent": "same direction within >= 3 patients (malignant vs non-malignant cells, or one malignant state vs the rest), DNA-explained, or a known GBM variant",
                "suggestive": "survives the cell-level and library-batch tests but no consistent within-patient contrast (or a batch-sensitive consistent one)",
                "low": "dropped by the cell-level re-test, or not re-tested"}}


def annotate_splicing(doc, table, source, params):
    """Tier the clusters, findings and events of one splicing.json (patient or cohort) in place; returns the counts."""
    for c in doc.get("clusters") or []:
        put(c, table.get(c.get("id")) or not_retested())
    for f in doc.get("findings") or []:
        put(f, tier_of_finding(f, table))
    by_finding = {f.get("id"): f for f in doc.get("findings") or []}
    for e in doc.get("events") or []:
        f = by_finding.get(e.get("finding_id"))
        if f is not None and f.get("kind") == "known_variant":
            put(e, tier_of_finding(f, table), detail=False)
        else:
            put(e, table.get(e.get("cluster_id")) or not_retested(), detail=False)
    m = meta(source, params)
    m["counts"] = {k: counts(doc.get(k) or []) for k in ("clusters", "findings", "events") if k in doc}
    doc["confidence_meta"] = m
    return m["counts"]


def annotate_support(support, patient_docs, source, params):
    """splice_support.json: tier each RNA entry by its matching finding / event (same gene and text) of the patient's splicing.json."""
    for p, genes in (support.get("patients") or {}).items():
        doc = patient_docs.get(p) or {}
        by_text = {}
        for f in doc.get("findings") or []:
            by_text[(f.get("gene"), f.get("text"))] = f
        for e in doc.get("events") or []:
            by_text.setdefault((e.get("gene"), e.get("summary")), e)
        for g, entry in genes.items():
            best = None
            for r in entry.get("rna") or []:
                hit = by_text.get((r.get("gene") or g, r.get("text")))
                if hit is None or "confidence" not in hit:
                    r["confidence"], r["confidence_reason"] = "low", ("not matched to a tiered splicing call" if doc else "patient splicing.json not available")
                else:
                    r["confidence"], r["confidence_reason"] = hit["confidence"], hit["confidence_reason"]
                if best is None or RANK[r["confidence"]] < RANK[best]:
                    best = r["confidence"]
            entry["rna_confidence"] = best
    support["confidence_meta"] = meta(source, params)
    return support


def annotate_vignettes(vig, finding_tiers):
    """vignettes.json: confidence on the 'splice' items, by the finding id of their evidence."""
    n = 0

    def walk(o):
        nonlocal n
        if isinstance(o, dict):
            if o.get("kind") == "splice" and isinstance(o.get("evidence"), list):
                fid = next((e.get("finding") for e in o["evidence"] if isinstance(e, dict) and e.get("finding")), None)
                c = finding_tiers.get(fid)
                if c:
                    o["confidence"], o["confidence_reason"] = c
                else:
                    o["confidence"], o["confidence_reason"] = "low", "splicing finding not found in rna/splicing.json"
                n += 1
            for x in o.values():
                walk(x)
        elif isinstance(o, list):
            for x in o:
                walk(x)

    walk(vig)
    return n


def load(path):
    with open(path) as f:
        return json.load(f)


def dump(doc, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(doc, f, separators=(",", ":"), allow_nan=False)
    os.replace(tmp, path)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="gOS data dir (<P>/rna/splicing.json, _cohort/...)")
    ap.add_argument("--deepdive", required=True, help="events_annotated.tsv of the differential-splicing deep-dive")
    ap.add_argument("--out", help="write the annotated files here (same relative paths); the data dir is only read")
    ap.add_argument("--in-place", action="store_true", help="overwrite the files in the data dir (back them up first)")
    ap.add_argument("--patients", default="", help="comma-separated patients (default: every <P>/rna/splicing.json)")
    ap.add_argument("--rescue-min-patients", type=int, default=RESCUE_MIN_PATIENTS)
    ap.add_argument("--rescue-max-q", type=float, default=RESCUE_MAX_Q)
    a = ap.parse_args(argv)
    if not a.out and not a.in_place:
        ap.error("give --out <staging dir> or --in-place")
    out = a.out or a.data
    params = {"rescue_min_patients": a.rescue_min_patients, "rescue_max_q": a.rescue_max_q}
    table = read_table(a.deepdive, rescue_min_patients=a.rescue_min_patients, rescue_max_q=a.rescue_max_q)
    print(f"deep-dive table: {len(table)} clusters", {t: sum(1 for c in table.values() if c['tier'] == t) for t in TIERS})
    pats = [p for p in a.patients.split(",") if p] or sorted(
        d for d in os.listdir(a.data) if not d.startswith("_") and os.path.isfile(os.path.join(a.data, d, "rna", "splicing.json")))
    written, docs, finding_tiers = [], {}, {}
    for p in pats:
        src = os.path.join(a.data, p, "rna", "splicing.json")
        doc = load(src)
        cnt = annotate_splicing(doc, table, a.deepdive, params)
        docs[p] = doc
        for f in doc.get("findings") or []:
            finding_tiers[f.get("id")] = (f["confidence"], f["confidence_reason"])
        dst = os.path.join(out, p, "rna", "splicing.json")
        dump(doc, dst)
        written.append(dst)
        print(p, cnt)
    coh = os.path.join(a.data, "_cohort", "rna", "splicing.json")
    if os.path.isfile(coh):
        doc = load(coh)
        print("_cohort", annotate_splicing(doc, table, a.deepdive, params))
        dst = os.path.join(out, "_cohort", "rna", "splicing.json")
        dump(doc, dst)
        written.append(dst)
    sup = os.path.join(a.data, "_cohort", "splice_support.json")
    if os.path.isfile(sup):
        dst = os.path.join(out, "_cohort", "splice_support.json")
        dump(annotate_support(load(sup), docs, a.deepdive, params), dst)
        written.append(dst)
    vig = os.path.join(a.data, "_cohort", "vignettes.json")
    if os.path.isfile(vig):
        v = load(vig)
        n = annotate_vignettes(v, finding_tiers)
        dst = os.path.join(out, "_cohort", "vignettes.json")
        dump(v, dst)
        written.append(dst)
        print("vignettes: splice items tiered", n)
    for w in written:
        print("wrote", w)
    return written


if __name__ == "__main__":
    main(sys.argv[1:])
