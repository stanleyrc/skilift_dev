"""chrX coverage-normalisation artefact: cells whose dryclean-zeroed chrX makes every chrX CN call (CN 0, HOMDEL) meaningless.

Flag table (one row per cell; columns `pair` = cell id, `patient`, `cls`): cls == "flagged" means chrX CN 0 over > 50% of chrX bp with a normal
autosomal profile (srctools analysis/chrx_artefact: chrx_artefact_count.R, cells_chrx_flags.tsv). The path is always a parameter.
Policy (cheap route, no dryclean / JaBbA rerun): CN-derived chrX calls from flagged cells are "not assessed"; mutation calls stay.
"""
import csv
import re

NOT_ASSESSED = "not assessed (chrX coverage-normalisation artefact)"
AMBIGUOUS_SEX = ("MGH304", "MGH305")   # chrY reads and X/autosome ratio disagree (README of the chrX analysis): no expected chrX CN


def load_flagged(path):
    """{patient: set(cell ids)} of cells with cls == 'flagged'."""
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if (r.get("cls") or "").strip() == "flagged":
                out.setdefault(r["patient"], set()).add(r["pair"])
    return out


def load_sex(path, ambiguous=AMBIGUOUS_SEX):
    """{patient: 'male'|'female'} from patient_sex_evidence.tsv; patients in `ambiguous` or without a clear call are left out."""
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            s = (r.get("sex") or "").lower()
            if r["patient"] not in ambiguous and (s.startswith("male") or s.startswith("female")):
                out[r["patient"]] = "male" if s.startswith("male") else "female"
    return out


def is_chrx(chrom):
    return re.sub("^chr", "", str(chrom), flags=re.I).upper() == "X"


def location_is_chrx(loc):
    return bool(re.match(r"^(chr)?X:", str(loc or ""), flags=re.I))


def is_chrx_cn_event(e):
    """A per-cell SCNA row on chrX derived from copy number (homdel / amp / loss), not a mutation or fusion."""
    return e.get("type") == "SCNA" and location_is_chrx(e.get("Genome_Location"))


def not_assess_events(events):
    """Return (new events, n changed): chrX HOMDEL rows of a flagged cell become NOT_ASSESSED. The input is not mutated.

    Only HOMDEL is rewritten: the other chrX SCNA calls in flagged cells (none exist in the current export, README) would need the same
    treatment and are counted as `other` by the caller if they appear."""
    out, n = [], 0
    for e in events:
        if is_chrx_cn_event(e) and str(e.get("vartype", e.get("Variant", ""))).upper() == "HOMDEL":
            e = dict(e, Variant=NOT_ASSESSED, vartype="NOT_ASSESSED", segment_cn=None, estimated_altered_copies=None,
                     not_assessed_reason=NOT_ASSESSED, original_call="HOMDEL")
            n += 1
        out.append(e)
    return out, n


def expected_chrx_cn(sex):
    return {"male": 1, "female": 2}.get(sex)
