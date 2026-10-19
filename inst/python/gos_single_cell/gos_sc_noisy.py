"""Data-driven 'noisy SNV' patient set shared by the exporters.

NOISY_PATIENTS (env, comma-separated) lists patients whose SNV calls are treated as noisy and reported with and without
them. Default: empty (no patient is flagged). CLEANED_SNV_PATIENTS (env, default "BWH69") lists patients whose SNVs come
from an artefact-audited cleaned site set; CLEANED_SNV_NOTE_BWH69 text below is used in caveats when they are present.
"""
import os

CLEANED_NOTES = {"BWH69": "BWH69 SNVs use a cleaned site set (5,675 sites) after an artefact audit."}


def _split(v):
    return [x.strip() for x in (v or "").split(",") if x.strip()]


def noisy_patients(env=None):
    env = os.environ if env is None else env
    return set(_split(env.get("NOISY_PATIENTS", "")))


def cleaned_patients(env=None):
    env = os.environ if env is None else env
    return set(_split(env.get("CLEANED_SNV_PATIENTS", "BWH69")))


def is_noisy(patient, env=None):
    return patient in noisy_patients(env)


def noisy_caveat(patients=None, env=None):
    """Caveat text naming noisy patients and/or the cleaned-set note, or None when neither applies."""
    nz = sorted(p for p in noisy_patients(env) if patients is None or p in patients)
    cl = sorted(p for p in cleaned_patients(env) if p in CLEANED_NOTES and (patients is None or p in patients) and p not in nz)
    parts = []
    if nz:
        parts.append(f"{', '.join(nz)} SNV calls are noisy (many single-cell-private sites); results are shown with and without them.")
    parts += [CLEANED_NOTES[p] for p in cl]
    return " ".join(parts) or None
