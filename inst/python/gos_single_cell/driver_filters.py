"""Pure helpers for the subclonal-driver artefact filters (no pysam / numpy; unit-testable).

Flag codes (read QC + site QC). EXCLUDING flags remove a hit from the strict recurrence test; INFORMATIONAL flags
are only shown. Hits keep their rows either way, with `excluded_reason` listing the excluding codes.

  EXCLUDING     exon_split        alt sits on exon-exon (spliced / SA) reads: RNA or processed-copy carry-over
                split_only        alt on split (SA) reads while ref reads are not
                read_end          alt within 10 bp of the aligned-block end in >= 80% of alt reads and ref <= 40%
                low_carrier_vaf   carrier VAF <= 0.5x the patient's well-placed subclonal baseline (and binomially unlikely)
                single_read       > 50% of carriers (>= 5) have a single alt read
                paralog           alt allele equals the base at the homologous position of a close paralog
                recurrent_allele  same allele seen (>= 2 cells) in >= 3 other patients' cells (a somatic subclonal
                                  SNV is not seen in unrelated patients; the scan has no per-hit normal cells)
  INFORMATIONAL segdup            site in a genomicSuperDups region
                low_mappability     site outside Umap k100 uniquely-mappable regions
                alt_mapq_lower      mean alt-read MAPQ lower than ref-read MAPQ by >= 3 (supporting evidence only)
                same_allele_other   any other patient with the allele (legacy flag text)
"""
import bisect

EXCLUDING = ("exon_split", "split_only", "read_end", "low_carrier_vaf", "single_read", "paralog", "recurrent_allele")
INFORMATIONAL = ("segdup", "low_mappability", "alt_mapq_lower", "same_allele_other")
RECURRENT_MIN_PATIENTS = 3
RECURRENT_MIN_CELLS = 2
READ_END_REF_MAX = 0.4
MAPQ_GAP = 3.0
MIN_PARALOG_IDENTITY = 0.85
MIN_PARALOG_COVERAGE = 0.8


def read_qc_codes(counts, other_patients=None, min_alt=3):
    """counts: hit_read_qc counts dict; other_patients: {patient: {cells_with_alt, ...}}. Returns list of codes."""
    c = counts or {}
    na, nr = c.get("alt_reads", 0), c.get("ref_reads", 0)
    f = lambda k, n: (c.get(k, 0) / n) if n else None
    codes = []
    if na >= min_alt and (f("alt_exon_split", na) or 0) >= 0.5:
        codes.append("exon_split")
    elif na >= min_alt and (f("alt_sa", na) or 0) >= 0.5 and (f("ref_sa", nr) or 0) < 0.25:
        codes.append("split_only")
    if na >= min_alt and (f("alt_end10", na) or 0) >= 0.8 and (f("ref_end10", nr) or 0) <= READ_END_REF_MAX:
        codes.append("read_end")
    n_other = sum(1 for v in (other_patients or {}).values() if (v.get("cells_with_alt") or 0) >= RECURRENT_MIN_CELLS)
    if n_other >= RECURRENT_MIN_PATIENTS:
        codes.append("recurrent_allele")
    return codes


def query_to_ref(cigar, rstart, qidx):
    """cigar: [(op, len)] with op letters MIDNSHP=X; rstart 0-based target start of the alignment;
    qidx: 0-based index into the (unclipped, as stored in SAM incl. soft clips) query sequence.
    Returns the 0-based target position aligned to qidx, or None (insertion / clip / deletion)."""
    q = 0
    r = rstart
    for op, n in cigar:
        if op in "M=X":
            if q <= qidx < q + n:
                return r + (qidx - q)
            q += n
            r += n
        elif op in "IS":
            if q <= qidx < q + n:
                return None
            q += n
        elif op in "DN":
            r += n
        # H, P consume nothing
    return None


def parse_cigar(s):
    out, num = [], ""
    for ch in s:
        if ch.isdigit():
            num += ch
        else:
            out.append((ch, int(num)))
            num = ""
    return out


def aln_identity(cigar, nm, qlen):
    """(identity over aligned columns, fraction of the query aligned)."""
    aligned = sum(n for op, n in cigar if op in "M=X")
    cols = sum(n for op, n in cigar if op in "M=XID")
    if not cols or not qlen:
        return 0.0, 0.0
    return max(0.0, 1.0 - nm / cols), aligned / qlen


def paralog_call(alt, ref, homologs):
    """homologs: [{chrom,pos,base,identity,coverage}] second-locus bases for the site. A site is `paralog` when the alt allele
    equals the homologous base of a close paralog (identity and coverage thresholds) and that base differs from ref."""
    hits = [h for h in homologs
            if h["identity"] >= MIN_PARALOG_IDENTITY and h["coverage"] >= MIN_PARALOG_COVERAGE
            and (h.get("base") or "").upper() == alt.upper() and alt.upper() != ref.upper()]
    return {"paralog": bool(hits), "matching_loci": [f"{h['chrom']}:{h['pos']}" for h in hits][:5]}


class IntervalIndex:
    """0-based half-open BED intervals per chromosome, merged; point queries by bisect."""
    def __init__(self, intervals):
        by = {}
        for c, s, e in intervals:
            by.setdefault(c, []).append((int(s), int(e)))
        self.starts, self.ends = {}, {}
        for c, v in by.items():
            v.sort()
            ms = []
            for s, e in v:
                if ms and s <= ms[-1][1]:
                    ms[-1][1] = max(ms[-1][1], e)
                else:
                    ms.append([s, e])
            self.starts[c] = [m[0] for m in ms]
            self.ends[c] = [m[1] for m in ms]

    def has(self, chrom):
        return chrom in self.starts

    def contains(self, chrom, pos1):
        """pos1: 1-based position."""
        st = self.starts.get(chrom)
        if not st:
            return False
        i = bisect.bisect_right(st, pos1 - 1) - 1
        return i >= 0 and self.ends[chrom][i] > pos1 - 1


def site_codes(site_qc):
    """site_qc: one record from gos_sc_driver_site_qc.py. Returns list of codes."""
    codes = []
    if not site_qc:
        return codes
    if (site_qc.get("paralog") or {}).get("paralog"):
        codes.append("paralog")
    if site_qc.get("segdup"):
        codes.append("segdup")
    if site_qc.get("umap_k100") is False:
        codes.append("low_mappability")
    m = site_qc.get("mapq") or {}
    if m.get("alt_reads", 0) >= 3 and m.get("ref_reads", 0) >= 3 and m.get("mean_ref") is not None and m.get("mean_alt") is not None \
            and m["mean_ref"] - m["mean_alt"] >= MAPQ_GAP:
        codes.append("alt_mapq_lower")
    return codes


def split_codes(codes):
    ex = [c for c in codes if c in EXCLUDING]
    return ex, [c for c in codes if c not in EXCLUDING]
