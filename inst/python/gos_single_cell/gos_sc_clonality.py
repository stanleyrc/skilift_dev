#!/usr/bin/env python3
"""Annotate read-evidence driver rows of data/<P>/filtered.events.json with dropout-aware clonality (clonality.py).

  python3 gos_sc_clonality.py --data <gOS data dir> [--patients A,B] [--tag bak_<yymmdd>] [--dry-run] [--report out.tsv]
         [--allsnv-work gos_datasets/gOS_GBM_sc/work/allsnv]

Rows with source in clonality.READ_EVIDENCE_SOURCES (hotspot rescue, pooled driver scan, called-SNV rows) get
clonality (clonal / subclonal / indeterminate), clonality_method = "dropout-aware" and clonality_stats; cell_fraction is unchanged.
Run after gos_sc_hotspot_rescue / gos_sc_pooled_drivers / gos_sc_allsnv and before gos_sc_subclonal_drivers, gos_sc_portrait and
gos_sc_clinical_recovery / gos_sc_clinical_events, which read the field. Idempotent; the first write keeps filtered.events.json.<tag>.
"""
import argparse
import json
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import clonality as cl  # noqa: E402

DEFAULT_ALLSNV = "/gpfs/commons/groups/imielinski_lab/projects/single_cell_GBM/gos_datasets/gOS_GBM_sc/work/allsnv"


def run(data, patients, tag, dry, allsnv_work, n_perm=cl.N_PERM, log=print):
    rep = []
    refs = {p: cl.patient_reference(data, p, cl.tumour_cells(data, p)) for p in patients}
    for p, r in refs.items():
        log(p, "truncal reference:", "none" if r is None else f"{r['n_sites']} sites, rate by depth {dict(zip(r['labels'], r['rate']))}")
    cohort = cl.cohort_reference(refs)
    for p in patients:
        path = os.path.join(data, p, "filtered.events.json")
        if not os.path.exists(path):
            continue
        with open(path) as fh:
            events = json.load(fh)
        if not any(cl.is_read_evidence(e) for e in events):
            continue
        ctx = cl.PatientContext(data, p, allsnv_work, reference=refs.get(p), fallback_reference=cohort)
        before = {id(e): (e.get("clonality"), e.get("clonality_method")) for e in events}
        res = cl.annotate_events(events, ctx, n_perm)
        for e, r in res:
            old = before[id(e)]
            rep.append(dict(patient=p, gene=e.get("gene"), variant=e.get("Variant"), source=e.get("source"), cells=e.get("cells"), cell_fraction=e.get("cell_fraction"),
                            old_clonality=old[0], new_clonality=r["status"], n_covered=r["n_covered"], observed=r["observed_detection"], expected=r["expected_detection"],
                            p_detection=r["p_detection"], p_clade=r["p_clade_restriction"], best_clade_f1=r["best_clade_f1"], reference=(ctx.reference or {}).get("scope"), evidence_source=r["evidence_source"], reason=r["reason"]))
            log(p, e.get("gene"), e.get("Variant"), e.get("cells"), "->", r["status"], "|", r["reason"])
        if not dry:
            bak = f"{path}.{tag}"
            if not os.path.exists(bak) and not os.path.islink(path):
                shutil.copy2(path, bak)
            elif os.path.islink(path):
                raise SystemExit(f"{path} is a symlink: refusing to write through it (stage a real file first)")
            with open(path + ".tmp", "w") as fh:
                json.dump(events, fh)
            os.replace(path + ".tmp", path)
    return rep


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True)
    ap.add_argument("--patients", default="")
    ap.add_argument("--tag", default="bak_clonality")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--report", default=None, help="TSV of every annotated row")
    ap.add_argument("--allsnv-work", default=DEFAULT_ALLSNV)
    ap.add_argument("--n-perm", type=int, default=cl.N_PERM)
    a = ap.parse_args(argv)
    pats = [p for p in a.patients.split(",") if p] or sorted(p for p in os.listdir(a.data) if p[:3] in ("BWH", "MGH") and p[3:].isdigit())
    rep = run(a.data, pats, a.tag, a.dry_run, a.allsnv_work, a.n_perm)
    if a.report:
        cols = list(rep[0]) if rep else []
        with open(a.report, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rep:
                fh.write("\t".join("" if r[c] is None else str(r[c]) for c in cols) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
