#!/usr/bin/env python3
"""Local ClinVar (GRCh38 VCF) and gnomAD r2.1 genomes (GRCh37) lookups for clinical SNV / indel sites. No web queries: patient
variants never leave the cluster.

ClinVar : DB/ClinVar/<download date>/clinvar.vcf.gz (NCBI vcf_GRCh38, ##fileDate recorded), matched on chrom/pos/ref/alt.
gnomAD  : DB/gnomAD/gnomad.genomes.r2.1.sites.chr*.vcf.bgz (GRCh37); the hg38 site is lifted with DB/liftover/hg38ToHg19.over.chain
          (liftOver), then matched on pos/ref/alt (strand-flipped alleles are tried when the lifted reference base disagrees).
The gnomAD .tbi files are older than the .bgz files, which tabix refuses, so fresh copies of the indexes are kept next to symlinks in
analysis/clinical_recovery_261017/gnomad_idx.
"""
import glob
import os
import re
import shutil
import subprocess
import tempfile

DB = "/gpfs/commons/groups/imielinski_lab/DB"
PROJ = "/gpfs/commons/groups/imielinski_lab/projects/single_cell_GBM"
CLINVAR_GLOB = DB + "/ClinVar/2*/clinvar.vcf.gz"
GNOMAD_DIR = DB + "/gnomAD"
GNOMAD_IDX = PROJ + "/analysis/clinical_recovery_261017/gnomad_idx"
CHAIN = DB + "/liftover/hg38ToHg19.over.chain"
LIFTOVER = "/gpfs/commons/home/sclarke/scripts/liftOver"
TABIX = "/gpfs/commons/groups/imielinski_lab/Software/speedseq-1.0/bin/tabix"
GNOMAD_RELEASE = "gnomAD r2.1 genomes (GRCh37)"
STARS = {"practice_guideline": 4, "reviewed_by_expert_panel": 3, "criteria_provided,_multiple_submitters,_no_conflicts": 2,
         "criteria_provided,_conflicting_classifications": 1, "criteria_provided,_conflicting_interpretations": 1,
         "criteria_provided,_single_submitter": 1}
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
# conditions ClinVar names as recessive (or whose name states it); used only for the optional "recessive" note
RECESSIVE_CONDITIONS = ("fanconi_anemia", "recessive", "ataxia-telangiectasia", "mutyh-associated_polyposis", "werner_syndrome", "bloom_syndrome",
                        "xeroderma_pigmentosum", "nijmegen_breakage_syndrome")


# ---------- pure helpers (unit tested) ----------
def parse_info(s):
    out = {}
    for kv in s.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
        elif kv:
            out[kv] = True
    return out


def review_stars(revstat):
    """CLNREVSTAT -> 0-4 stars."""
    return STARS.get((revstat or "").replace("|", ","), 0) if revstat else 0


def is_plp(sig):
    """Pathogenic / likely pathogenic (not conflicting, not 'risk factor')."""
    s = (sig or "").lower()
    return ("pathogenic" in s and "conflicting" not in s and "benign" not in s) or s.startswith("pathogenic")


def clinvar_from_line(line, release):
    f = line.rstrip("\n").split("\t")
    info = parse_info(f[7])
    sig = (info.get("CLNSIG") or "").replace("_", " ")
    rs = info.get("CLNREVSTAT", "")
    return {"found": True, "variation_id": f[2], "significance": sig, "review_status": rs.replace("_", " "), "stars": review_stars(rs),
            "condition": (info.get("CLNDN") or "").replace("|", "; "), "plp": is_plp(sig), "release": release,
            "link": "https://www.ncbi.nlm.nih.gov/clinvar/variation/%s/" % f[2], "source": "ClinVar GRCh38 VCF, local"}


def clinvar_none(release):
    return {"found": False, "release": release, "text": "not in ClinVar (release %s)" % release, "source": "ClinVar GRCh38 VCF, local"}


def gnomad_from_line(line, alt):
    f = line.rstrip("\n").split("\t")
    alts = f[4].split(",")
    if alt not in alts:
        return None
    k = alts.index(alt)
    info = parse_info(f[7])

    def pick(key, cast):
        v = info.get(key)
        if v is None or v is True:
            return None
        parts = v.split(",")
        try:
            return cast(parts[k] if len(parts) > k else parts[0])
        except ValueError:
            return None
    return {"found": True, "af": pick("AF", float), "ac": pick("AC", int), "an": pick("AN", int), "nhomalt": pick("nhomalt", int),
            "af_popmax": pick("AF_popmax", float), "popmax": pick("popmax", str), "rsid": f[2] if f[2] != "." else None,
            "filter": f[6], "release": GNOMAD_RELEASE, "source": "gnomAD r2.1 genomes, local"}


def gnomad_link(chrom37, pos37, ref, alt):
    return "https://gnomad.broadinstitute.org/variant/%s-%d-%s-%s?dataset=gnomad_r2_1" % (str(chrom37).replace("chr", ""), pos37, ref, alt)


def is_lof(typ, cdna):
    """Predicted loss of function from the report: nonsense / frameshift (Trunc) or canonical splice site (+-1 or 2)."""
    if typ == "Trunc":
        return True
    m = re.search(r"\d([+-])(\d+)(?!\d)", cdna or "")
    return typ == "Splice" and bool(m) and int(m.group(2)) <= 2 and not (cdna or "").startswith("c.-")


def recessive_note(clinvar, lof):
    """Only when ClinVar's own condition names a recessive disease; nothing is inferred otherwise."""
    if not (lof and clinvar and clinvar.get("found")):
        return None
    for cond in (clinvar.get("condition") or "").split("; "):
        if any(c in cond.lower() for c in RECESSIVE_CONDITIONS):
            return "one copy; recessive condition per ClinVar (%s): carriers are usually unaffected" % cond.replace("_", " ")
    return None


def highlight_tags(clinvar, lof, germline):
    tags = []
    if clinvar and clinvar.get("plp"):
        tags.append("ClinVar P/LP")
    if lof:
        tags.append("LoF, germline" if germline else "LoF")
    return tags


# ---------- I/O ----------
def clinvar_vcf():
    c = sorted(glob.glob(CLINVAR_GLOB))
    return c[-1] if c else None


def clinvar_release(path):
    import gzip
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##fileDate="):
                return line.strip().split("=", 1)[1]
            if not line.startswith("##"):
                break
    return "unknown"


def _tabix(path, region):
    r = subprocess.run([TABIX, path, region], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError("tabix %s %s: %s" % (path, region, r.stderr.strip()[:200]))
    return [l for l in r.stdout.split("\n") if l]


def clinvar_lookup(vcf, release, chrom, pos, ref, alt, indel_len=None, window=10):
    """SNV: exact chrom/pos/ref/alt. Deletion (indel_len = deleted bases): a ClinVar deletion of that length within +-window of pos."""
    c = chrom.replace("chr", "")
    if ref and alt:
        for l in _tabix(vcf, "%s:%d-%d" % (c, pos, pos)):
            f = l.split("\t")
            if f[3] == ref and alt in f[4].split(","):
                return clinvar_from_line(l, release)
    elif indel_len:
        for l in _tabix(vcf, "%s:%d-%d" % (c, pos - window, pos + window)):
            f = l.split("\t")
            if any(len(f[3]) - len(a) == indel_len and len(a) < len(f[3]) for a in f[4].split(",")):
                out = clinvar_from_line(l, release)
                out["matched_by"] = "deletion of the same length within %d bp of the reported position (alleles not given in the report)" % window
                return out
    return clinvar_none(release)


def lift(chrom, pos):
    """hg38 -> (chrom37, pos37, strand) or None when liftOver fails."""
    d = tempfile.mkdtemp()
    try:
        bed, out, un = [os.path.join(d, x) for x in ("in.bed", "out.bed", "un.bed")]
        with open(bed, "w") as fh:
            fh.write("%s\t%d\t%d\tx\t0\t+\n" % ("chr" + chrom.replace("chr", ""), pos - 1, pos))
        subprocess.run([LIFTOVER, "-bedPlus=6", bed, CHAIN, out, un], capture_output=True, text=True)
        rows = [l.split("\t") for l in open(out).read().split("\n") if l]
        if not rows:
            return None
        return rows[0][0].replace("chr", ""), int(rows[0][2]), rows[0][5] if len(rows[0]) > 5 else "+"
    finally:
        shutil.rmtree(d, ignore_errors=True)


def gnomad_path(chrom37):
    os.makedirs(GNOMAD_IDX, exist_ok=True)
    name = "gnomad.genomes.r2.1.sites.chr%s.vcf.bgz" % chrom37
    dst = os.path.join(GNOMAD_IDX, name)
    if not os.path.exists(dst + ".tbi"):
        if os.path.lexists(dst):
            os.remove(dst)
        os.symlink(os.path.join(GNOMAD_DIR, name), dst)
        shutil.copy(os.path.join(GNOMAD_DIR, name + ".tbi"), dst + ".tbi")
    return dst


def gnomad_lookup(chrom, pos, ref, alt):
    """SNV (ref/alt given): lift, match, return gnomAD fields + link. Always returns a dict with `found` and a text / reason."""
    lo = lift(chrom, pos)
    if lo is None:
        return {"found": False, "lift_failed": True, "text": "liftover from GRCh38 to GRCh37 failed (not available)", "release": GNOMAD_RELEASE}
    c37, p37, strand = lo
    if not (ref and alt):
        return {"found": False, "text": "not looked up (alleles not given in the report)", "release": GNOMAD_RELEASE}
    cands = [(ref, alt)] + ([(COMP.get(ref), COMP.get(alt))] if strand == "-" else [])
    path = gnomad_path(c37)
    for l in _tabix(path, "%s:%d-%d" % (c37, p37, p37)):
        f = l.split("\t")
        for r, a in cands:
            if f[3] == r:
                g = gnomad_from_line(l, a)
                if g:
                    g.update(link=gnomad_link(c37, p37, r, a), pos37="%s:%d" % (c37, p37))
                    return g
    return {"found": False, "text": "not in gnomAD r2.1 genomes (AF < ~1/30,000)", "release": GNOMAD_RELEASE, "pos37": "%s:%d" % (c37, p37)}
