## Artefact filters of the splicing events table: a normal (non-tumour) cell panel and a low-complexity / simple-repeat
## flag at the new splice sites. Both feed the "recommended" tier (sc_splice_events, sc_ir_events).

SC_REPEAT_CLASSES <- c("Simple_repeat", "Low_complexity", "Satellite", "rRNA", "snRNA", "scRNA", "srpRNA", "tRNA")

#' Low-complexity / simple-repeat track (hg38)
#'
#' RepeatMasker classes Simple_repeat, Low_complexity, Satellite and the small-RNA families, plus the UCSC simpleRepeat (TRF)
#' track, from the lab's hg38 UCSC tables.
#'
#' @param rmsk,simple UCSC rmsk.txt.gz and simpleRepeat.txt.gz of hg38
#' @return GRanges (chromosome without "chr"), reduced
#' @export
#' @author Stanley Clarke
sc_repeat_track <- function(rmsk = "/gpfs/commons/groups/imielinski_lab/DB/hg38/ucsc/rmsk.txt.gz",
                            simple = "/gpfs/commons/groups/imielinski_lab/DB/hg38/ucsc/simpleRepeat.txt.gz") {
    r <- data.table::fread(rmsk, header = FALSE, select = c(6, 7, 8, 12), col.names = c("chr", "s", "e", "class"), showProgress = FALSE)
    r <- r[class %in% SC_REPEAT_CLASSES & !grepl("_", chr)]
    s <- data.table::fread(simple, header = FALSE, select = c(2, 3, 4), col.names = c("chr", "s", "e"), showProgress = FALSE)[!grepl("_", chr)]
    x <- data.table::rbindlist(list(r[, .(chr, s, e)], s))
    GenomicRanges::reduce(GenomicRanges::GRanges(sub("^chr", "", x$chr), IRanges::IRanges(x$s + 1L, x$e)))
}

#' Splice sites of junctions that sit in a low-complexity / simple repeat
#'
#' Only the sites that are not annotated splice sites are tested (an annotated site in a repeat is not a mapping suspect).
#'
#' @param ev data.table(chromosome, start, end) intron coordinates
#' @param track sc_repeat_track
#' @param ann_starts,ann_ends character "chromosome:position" of annotated intron starts / ends
#' @param pad bp either side of the site
#' @return logical per row
#' @export
sc_splice_repeat_flag <- function(ev, track, ann_starts, ann_ends, pad = 10L) {
    if (is.null(track) || !nrow(ev)) return(rep(FALSE, nrow(ev)))
    novel_l <- !paste0(ev$chromosome, ":", ev$start) %in% ann_starts
    novel_r <- !paste0(ev$chromosome, ":", ev$end) %in% ann_ends
    hit <- function(pos, sel) {
        out <- rep(FALSE, length(pos))
        if (any(sel)) out[sel] <- GenomicRanges::countOverlaps(GenomicRanges::GRanges(ev$chromosome[sel], IRanges::IRanges(pmax(1L, pos[sel] - pad), pos[sel] + pad)), track) > 0
        out
    }
    hit(ev$start, novel_l) | hit(ev$end, novel_r)
}

#' Normal-cell panel of junction usage
#'
#' @param per list patient -> per-cell junction table (rna_id, chromosome, start, end, count)
#' @param normal list patient -> rna ids of normal / non-tumour cells
#' @param cl junction table with cluster_id
#' @return data.table(cluster_id, start, end, reads, cells, cluster_reads, psi) pooled over every patient's normal cells
#' @export
sc_splice_normal_panel <- function(per, normal, cl) {
    nj <- data.table::rbindlist(lapply(names(per), function(p) {
        ids <- normal[[p]]
        if (!length(ids)) return(NULL)
        per[[p]][tolower(rna_id) %in% tolower(ids), .(rna_id = paste(p, rna_id), chromosome, start, end, count)]
    }))
    if (is.null(nj) || !nrow(nj)) return(data.table::data.table(cluster_id = character(0), start = integer(0), end = integer(0), reads = integer(0),
                                                               cells = integer(0), cluster_reads = integer(0), psi = numeric(0)))
    nj <- merge(nj, unique(cl[, .(chromosome, start, end, cluster_id)]), by = c("chromosome", "start", "end"))
    p <- nj[, .(reads = sum(count), cells = data.table::uniqueN(rna_id)), by = .(cluster_id, start, end)]
    p[, cluster_reads := sum(reads), by = cluster_id]
    p[, psi := reads / cluster_reads]
    p[]
}

## "seen in normal cells": the junction is used in the normal panel (>= 3 reads in >= 2 cells) at >= half the tumour PSI
sc_splice_normal_flag <- function(psi, normal_psi, normal_reads, normal_cells, min_reads = 3L, min_cells = 2L, ratio = 0.5) {
    !is.na(normal_reads) & normal_reads >= min_reads & !is.na(normal_cells) & normal_cells >= min_cells & !is.na(normal_psi) & normal_psi >= ratio * psi
}
