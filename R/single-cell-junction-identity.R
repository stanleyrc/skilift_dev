## Junction identity for single-cell gOS exports: one id per junction of the patient-level (re-collected) junction set,
## plus the per-cell supporting read counts. Matching of per-cell graph junctions to these ids and the event families
## are done by the srctools runner (inst/gos_single_cell/gos_sc_event_families.py), which reads the tables written here.

.sc_chrom_order = c(as.character(1:22), "X", "Y")

#' Pure helper: build the master junction table from breakpoint pairs.
#' @param bp data.table with chrom1,pos1,strand1,chrom2,pos2,strand2 (chrom with or without "chr") and merged_id
#' @return data.table (id, chrom1..strand2, merged_id); ends ordered by genome position; id as "5:5624784+|5:6302857+"
#' @export
sc_junction_master_from_bp = function(bp) {
  bp = data.table::copy(data.table::as.data.table(bp))
  bp[, `:=`(chrom1 = sub("^chr", "", as.character(chrom1)), chrom2 = sub("^chr", "", as.character(chrom2)))]
  o1 = match(bp$chrom1, .sc_chrom_order); o2 = match(bp$chrom2, .sc_chrom_order)
  o1[is.na(o1)] = 99L; o2[is.na(o2)] = 99L
  swap = o1 > o2 | (o1 == o2 & bp$pos1 > bp$pos2)
  out = bp[, .(chrom1 = ifelse(swap, chrom2, chrom1), pos1 = as.integer(ifelse(swap, pos2, pos1)),
               strand1 = as.character(ifelse(swap, strand2, strand1)),
               chrom2 = ifelse(swap, chrom1, chrom2), pos2 = as.integer(ifelse(swap, pos1, pos2)),
               strand2 = as.character(ifelse(swap, strand1, strand2)), merged_id)]
  out[, id := paste0(chrom1, ":", pos1, strand1, "|", chrom2, ":", pos2, strand2)]
  out[, `:=`(chrom1 = paste0("chr", chrom1), chrom2 = paste0("chr", chrom2))]
  data.table::setcolorder(out, c("id", "chrom1", "pos1", "strand1", "chrom2", "pos2", "strand2", "merged_id"))
  out[!duplicated(id)]
}

#' Pure helper: per-cell supporting reads for the master junctions.
#' @param counts data.table(pair, merged_id, max_reads) as in merged.counts.dt_*.rds
#' @param master output of sc_junction_master_from_bp
#' @return data.table(id, cell, reads); only cells with >= 1 read are listed (absent = 0 reads observed)
#' @export
sc_junction_support_from_counts = function(counts, master) {
  counts = data.table::as.data.table(counts)
  m = merge(master[, .(id, merged_id)], counts[, .(cell = pair, merged_id, reads = max_reads)], by = "merged_id")
  m = m[, .(reads = max(reads)), by = .(id, cell)]
  m[reads > 0]
}

#' Read a re-collected junction set (jj_*.rds + merged.counts.dt_*.rds) and write the master / support tables
#' @param jj_rds gGnome Junction object (jj_<date>_whole_genome*.rds); breakpoint strands are the gGnome junction strands
#' @param counts_rds merged.counts.dt_*.rds (columns pair, merged_id, max_reads = distinct supporting read names per cell)
#' @param out_dir directory for junction_master.tsv and junction_support.tsv
#' @param pair_map optional named character vector, counts `pair` -> gOS cell directory name
#' @export
sc_export_junction_tables = function(jj_rds, counts_rds, out_dir, pair_map = NULL) {
  jj = readRDS(jj_rds)
  grl = jj$grl
  dt = data.table::as.data.table(jj$dt)
  ix = which(lengths(grl) == 2)
  gl = grl[ix]
  first = unlist(gl)[seq(1, 2 * length(gl), 2)]; second = unlist(gl)[seq(2, 2 * length(gl), 2)]
  bp = data.table::data.table(chrom1 = as.character(GenomicRanges::seqnames(first)), pos1 = GenomicRanges::start(first),
                              strand1 = as.character(GenomicRanges::strand(first)),
                              chrom2 = as.character(GenomicRanges::seqnames(second)), pos2 = GenomicRanges::start(second),
                              strand2 = as.character(GenomicRanges::strand(second)),
                              merged_id = names(gl))
  ## jj$dt$merged_id is "JUNC_MERGED_<k>" and the grl element name is <k> (merged.ix)
  if (!is.null(dt$merged_id) && !is.null(dt$merged.ix)) {
    bp[, merged_id := dt$merged_id[match(as.integer(merged_id), dt$merged.ix)]]
  }
  master = sc_junction_master_from_bp(bp)
  counts = readRDS(counts_rds)
  if (!is.null(pair_map)) counts[, pair := ifelse(pair %in% names(pair_map), pair_map[pair], pair)]
  support = sc_junction_support_from_counts(counts, master)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  data.table::fwrite(master, file.path(out_dir, "junction_master.tsv"), sep = "\t")
  data.table::fwrite(support, file.path(out_dir, "junction_support.tsv"), sep = "\t")
  invisible(list(master = master, support = support))
}
