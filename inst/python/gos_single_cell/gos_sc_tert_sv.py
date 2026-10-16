#!/usr/bin/env python3
"""Mark TERT promoter rearrangements in gOS single-cell filtered.events.json.

A DNA fusion whose 3' partner is TERT and whose TERT side starts at exon 1 (p.1) or upstream keeps the whole TERT
coding sequence and puts it behind another locus's promoter / enhancers: TERT activation by rearrangement, the
structural alternative to the C228T / C250T promoter mutations (Peifer et al. 2015, Nature 526:700). The fusion
caller tiers these as an unknown in-frame fusion (tier 3); this step sets them to tier 1 like the promoter hotspots
(gos_sc_hotspot_rescue) and keeps the caller's tier in Tier_original.

Fields set on each matching event: tert_promoter_sv = true, Tier = "1", Tier_original, effect, effect_description,
variant_summary. Idempotent. Each changed file is backed up first as <file>.bak_<yymmdd>_tertsv.

usage: gos_sc_tert_sv.py --data <gOS data dir> [--patients P1,P2] [--dry-run]
"""
import argparse
import datetime
import json
import os
import re
import shutil
import sys

TERT_SIDE = re.compile(r"(?:Exon\s+(-?\d+)\s*\(p\.([^)]*)\))", re.I)
EFFECT = "Likely Gain-of-function"
DESCRIPTION = ("Rearrangement joins the whole TERT coding sequence (from exon 1) to a {partner} locus: TERT promoter "
               "rearrangement / enhancer hijacking, the structural route to TERT activation instead of a C228T / C250T "
               "promoter mutation (Peifer et al. 2015).")


def genes_of(e):
    g = e.get("fusion_genes") or e.get("gene") or ""
    return [x.strip() for x in str(g).split("::") if x.strip()]


def is_tert_promoter_sv(e):
    """Fusion with TERT as the 3' partner whose TERT side starts at exon <= 1 at p.1 or in the 5' UTR."""
    vt = f"{e.get('vartype') or ''} {e.get('type') or ''}".lower()
    if "fusion" not in vt:
        return False
    genes = genes_of(e)
    if len(genes) < 2 or genes[-1].upper() != "TERT" or any(g.upper() == "TERT" for g in genes[:-1]):
        return False
    variant = str(e.get("Variant") or "")
    tert_side = variant.split("::")[-1] if "::" in variant else ""
    m = TERT_SIDE.search(tert_side)
    if not m:
        return False
    exon, prot = int(m.group(1)), m.group(2).strip().lower()
    return exon <= 1 and (prot in ("1", "5'utr", "5utr") or exon < 1)


def mark(e):
    """Annotate one event in place; returns True when anything changed."""
    partner = "::".join(genes_of(e)[:-1])
    want = {
        "tert_promoter_sv": True,
        "Tier": "1",
        "effect": EFFECT,
        "effect_description": DESCRIPTION.format(partner=partner),
        "variant_summary": f"TERT promoter rearrangement ({partner}::TERT): TERT activation by structural rearrangement.",
    }
    changed = False
    if "Tier_original" not in e:
        e["Tier_original"] = e.get("Tier")
        changed = True
    for k, v in want.items():
        if e.get(k) != v:
            e[k] = v
            changed = True
    return changed


def process(path, dry_run=False, today=None):
    with open(path) as fh:
        events = json.load(fh)
    if not isinstance(events, list):
        return 0, []
    hits = [e for e in events if is_tert_promoter_sv(e)]
    n_changed = sum(mark(e) for e in hits)
    labels = [f"{'::'.join(genes_of(e))} ({e.get('cells') or e.get('n_cells')})" for e in hits]
    if n_changed and not dry_run:
        stamp = (today or datetime.date.today()).strftime("%y%m%d")
        bak = f"{path}.bak_{stamp}_tertsv"
        if not os.path.exists(bak):
            shutil.copy2(path, bak)
        tmp = f"{path}.tmp_tertsv"
        with open(tmp, "w") as fh:
            json.dump(events, fh)
        os.replace(tmp, path)
    return n_changed, labels


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", required=True)
    ap.add_argument("--patients", default="")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)
    pats = [p for p in a.patients.split(",") if p] or sorted(
        d for d in os.listdir(a.data) if "_" not in d and os.path.isfile(os.path.join(a.data, d, "filtered.events.json")))  # patients; cell folders are <patient>_<well>
    for p in pats:
        f = os.path.join(a.data, p, "filtered.events.json")
        if not os.path.isfile(f):
            print(f"{p}\tnot available (no filtered.events.json)")
            continue
        n, labels = process(f, dry_run=a.dry_run)
        print(f"{p}\t{len(labels)} TERT promoter SV event(s)\t{n} changed{' (dry run)' if a.dry_run else ''}\t{'; '.join(labels)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
