#!/usr/bin/env python
"""RNA IGV slice of one cell around its fusion breakpoints.

usage: sc_rna_fusion_slice.py <sorted.bam> <out.bam> <read_ids.txt> <chr:start-end,...>

Writes the reads of <sorted.bam> overlapping the regions (deduplicated) to
<out.bam>, coordinate sorted and indexed; reads named in read_ids.txt (Arriba
read_identifiers, with or without the "CELL|" prefix) get the tag ZF:i:1.
"""
import sys
import pysam

bam_in, bam_out, ids_file, regions = sys.argv[1:5]
ids = set()
with open(ids_file) as f:
    for line in f:
        line = line.strip()
        if line:
            ids.add(line)
            ids.add(line.split("|", 1)[-1])

src = pysam.AlignmentFile(bam_in, "rb")
tmp = bam_out + ".unsorted.bam"
dst = pysam.AlignmentFile(tmp, "wb", template=src)
seen = set()
for reg in regions.split(","):
    for r in src.fetch(region=reg):
        key = (r.query_name, r.flag, r.reference_id, r.reference_start)
        if key in seen:
            continue
        seen.add(key)
        name = r.query_name
        if name in ids or name.split("|", 1)[-1] in ids:
            r.set_tag("ZF", 1, value_type="i")
        dst.write(r)
dst.close()
src.close()
pysam.sort("-o", bam_out, tmp)
pysam.index(bam_out)
import os
os.remove(tmp)
