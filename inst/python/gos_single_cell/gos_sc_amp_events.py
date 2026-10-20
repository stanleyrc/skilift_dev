"""Focal amplifications as AMP events in filtered.events.json (genes the oncotable driver path missed).

The driver path (oncotable on a fixed gene list) misses amplicons in genes outside the list and calls
some listed genes in only a few cells (MGH285 FOXO1: 1 cell, SOX5 in MGH301: none). This exporter adds
AMP rows (only where the live file has no AMP row for the gene, or one with < MIN_CELLS cells) from the per-cell 100 kb copy number (cn_binned) for the genes that carry an amplicon of
<P>/ecdna_status.json (comma-separated amplicon `label`s) or are named with --genes:

  gene CN (median of the 100 kb bins overlapping the gene) >= max(5, 2 x cell ploidy, 2 x the cell's
  median CN on that chromosome) in >= MIN_CELLS (5) tumour cells.

Usage
  export: gos_sc_amp_events.py export DATA_DIR OUT_DIR [--genes G1,G2] [--min-cells 5] PATIENT ...
      reads DATA_DIR/<P>/{cn_binned.json,cn_binned_100kb.cn.u16,ecdna_status.json,filtered.events.json},
      DATA_DIR/_cohort/gene_spans.json; writes OUT_DIR/<P>/amp_events_added.json (rows + base sha256).
      Never writes filtered.events.json: promotion merges into the LIVE file (merge, below).
  merge:  gos_sc_amp_events.py merge LIVE_EVENTS.json ADDED.json OUT.json [--expect-base SHA]
      Union of carrier cells into an existing AMP row of the gene, or append a new row. Idempotent. All
      other rows are untouched, so a staging copy built on an older base cannot overwrite newer edits.
"""
import hashlib
import json
import os
import sys

import numpy as np

FORMAT = "gos-sc-amp-events/1"
MISSING = 65535
SOURCE = "gos_sc_amp_events (100 kb CN)"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def load_cn(pdir):
    meta = json.load(open(os.path.join(pdir, "cn_binned.json")))
    lv = meta["levels"]["100kb"]
    n, nb = len(meta["cells"]), lv["n_bins"]
    cn = np.fromfile(os.path.join(pdir, lv["files"]["cn"]), dtype="<u2").reshape(n, nb).astype(float)
    cn[cn == lv["missing_value"]] = np.nan
    return meta, lv, cn


def cell_stats(cn, lv):
    """Per-cell ploidy (autosomal median) and per-chromosome median."""
    chroms = {c["chr"]: c for c in lv["chromosomes"]}
    auto = [c for k, c in chroms.items() if k not in ("X", "Y", "chrX", "chrY")]
    ploidy = np.nanmedian(np.concatenate([cn[:, c["first_bin"]:c["first_bin"] + c["n_bins"]] for c in auto], axis=1), axis=1)
    chrmed = {k: np.nanmedian(cn[:, c["first_bin"]:c["first_bin"] + c["n_bins"]], axis=1) for k, c in chroms.items()}
    return chroms, ploidy, chrmed


def gene_cn(cn, lv, chroms, g):
    c = chroms.get(str(g["chrom"]).replace("chr", ""))
    if c is None:
        return None
    b = lv["bin_size"]
    lo, hi = (g["start"] - 1) // b, (g["end"] - 1) // b
    lo, hi = max(lo, 0), min(hi, c["n_bins"] - 1)
    if hi < lo:
        return None
    return np.nanmedian(cn[:, c["first_bin"] + lo:c["first_bin"] + hi + 1], axis=1)


def amplified_genes(cn, lv, normal, genes, min_cells=5):
    """genes: list of dict(gene, chrom, start, end). Returns {gene: dict(cells=[idx], segment_cn, ...)}."""
    chroms, ploidy, chrmed = cell_stats(cn, lv)
    tumour = ~np.asarray(normal, bool)
    out = {}
    for g in genes:
        v = gene_cn(cn, lv, chroms, g)
        if v is None:
            continue
        ck = str(g["chrom"]).replace("chr", "")
        thr = np.maximum.reduce([np.full(len(v), 5.0), 2 * ploidy, 2 * chrmed[ck]])
        with np.errstate(invalid="ignore"):
            hit = tumour & (v >= thr)
        idx = np.where(hit)[0]
        if len(idx) >= min_cells:
            out[g["gene"]] = dict(cells=idx.tolist(), segment_cn=float(np.nanmedian(v[idx])),
                                  ploidy=float(np.nanmedian(ploidy[idx])), g=g)
    return out


def make_row(patient, gene, info, cell_names, n_tumour, template=None):
    g = info["g"]
    ids = [cell_names[i] for i in info["cells"]]
    seg = int(round(info["segment_cn"]))
    row = dict(template) if template else dict(
        gene=gene, fusion_genes="", gene_summary=None, role=None, Tier=3,
        variant_summary="%s amplification (focal copy-number call; not reviewed by OncoKB)" % gene,
        effect="Unknown", effect_description=None)
    row.update(vartype="AMP", type="SCNA", Variant="AMP",
               Genome_Location="%s:%d-%d" % (g["chrom"], g["start"], g["end"]), seqnames=str(g["chrom"]),
               start=str(g["start"]), end=str(g["end"]), sample=ids[0], id=patient,
               estimated_altered_copies=int(round(seg - info["ploidy"])), segment_cn=seg,
               n_cells=len(ids), n_normal_cells=0, cell_ids=",".join(ids),
               cells="%d/%d" % (len(ids), n_tumour), cell_fraction=round(len(ids) / n_tumour, 3),
               amp_source=SOURCE)
    return row


def export_patient(data_dir, out_dir, patient, extra_genes=(), min_cells=5):
    pdir = os.path.join(data_dir, patient)
    meta, lv, cn = load_cn(pdir)
    spans = {g["gene"]: g for g in json.load(open(os.path.join(data_dir, "_cohort", "gene_spans.json")))["genes"]}
    names = set(extra_genes)
    ef = os.path.join(pdir, "ecdna_status.json")
    if os.path.exists(ef):
        for a in json.load(open(ef)).get("amplicons", []):
            names |= {x.strip() for x in (a.get("label") or "").split(",") if x.strip()}
    genes = [spans[n] for n in sorted(names) if n in spans]
    missing = sorted(n for n in names if n not in spans)
    normal = meta["normal"]
    n_tumour = int(sum(1 for x in normal if not x))
    amp = amplified_genes(cn, lv, normal, genes, min_cells)
    base = os.path.join(pdir, "filtered.events.json")
    ev = json.load(open(base))
    rows = []
    for gene, info in amp.items():
        live = next((r for r in ev if r.get("gene") == gene and r.get("vartype") == "AMP" and r.get("type") == "SCNA"), None)
        if live and (live.get("n_cells") or 0) >= min_cells:
            continue  # the driver path already reports this amplification: leave its row untouched
        tmpl = next((r for r in ev if r.get("gene") == gene), None)
        keep = ("gene", "fusion_genes", "gene_summary", "role", "Tier", "variant_summary", "effect", "effect_description")
        t = {k: tmpl.get(k) for k in keep} if tmpl else None
        rows.append(make_row(patient, gene, info, meta["cells"], n_tumour, t))
    res = dict(format=FORMAT, patient=patient, rule="gene CN >= max(5, 2x ploidy, 2x chromosome median) in >= %d tumour cells" % min_cells,
               base_sha256=sha256(base), n_tumour_cells=n_tumour, genes_considered=sorted(names), genes_without_span=missing,
               rows=rows)
    os.makedirs(os.path.join(out_dir, patient), exist_ok=True)
    json.dump(res, open(os.path.join(out_dir, patient, "amp_events_added.json"), "w"), indent=1)
    return res


def merge_events(events, added_rows):
    """Merge AMP rows into a list of events; returns (new_list, summary). Idempotent."""
    out = [dict(r) for r in events]
    summary = []
    for add in added_rows:
        gene = add["gene"]
        i = next((k for k, r in enumerate(out) if r.get("gene") == gene and r.get("vartype") == "AMP" and r.get("type") == "SCNA"), None)
        if i is None:
            out.append(dict(add))
            summary.append((gene, "added", add["n_cells"]))
            continue
        cur = out[i]
        ids = list(dict.fromkeys([x for x in (cur.get("cell_ids") or "").split(",") if x] + add["cell_ids"].split(",")))
        if len(ids) == cur.get("n_cells") and cur.get("amp_source") == add["amp_source"]:
            continue
        n_tum = int(add["cells"].split("/")[1])
        cur.update(cell_ids=",".join(ids), n_cells=len(ids), cells="%d/%d" % (len(ids), n_tum), cell_fraction=round(len(ids) / n_tum, 3),
                   segment_cn=max(cur.get("segment_cn") or 0, add["segment_cn"]),
                   estimated_altered_copies=max(cur.get("estimated_altered_copies") or 0, add["estimated_altered_copies"]),
                   amp_source=add["amp_source"])
        summary.append((gene, "extended", len(ids)))
    return out, summary


def merge_cli(live, added, out, expect_base=None):
    a = json.load(open(added))
    if expect_base and sha256(live) != expect_base:
        print("note: live file differs from export base (expected: rows are merged into the live file)", file=sys.stderr)
    ev = json.load(open(live))
    new, summary = merge_events(ev, a["rows"])
    json.dump(new, open(out, "w"))
    return summary


def main(argv):
    if len(argv) < 2 or argv[0] not in ("export", "merge"):
        print(__doc__)
        return 2
    if argv[0] == "merge":
        exp = argv[argv.index("--expect-base") + 1] if "--expect-base" in argv else None
        for s in merge_cli(argv[1], argv[2], argv[3], exp):
            print("merge", *s)
        return 0
    rest = argv[1:]
    genes, mc = [], 5
    for flag in ("--genes", "--min-cells"):
        if flag in rest:
            i = rest.index(flag)
            v = rest[i + 1]
            del rest[i:i + 2]
            if flag == "--genes":
                genes = v.split(",")
            else:
                mc = int(v)
    data_dir, out_dir, pats = rest[0], rest[1], rest[2:]
    for p in pats:
        r = export_patient(data_dir, out_dir, p, genes, mc)
        print(p, "AMP rows:", [(x["gene"], x["n_cells"], x["segment_cn"]) for x in r["rows"]], "no span:", r["genes_without_span"])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
