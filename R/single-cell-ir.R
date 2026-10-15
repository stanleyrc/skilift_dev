## Intron retention (IR) of single-cell RNA for the gOS splicing events table.
##
## IRFinder-style, per cell and annotated intron: reads that continue across the exon / intron boundary on either side of
## the intron (unspliced; `ret_l`, `ret_r`), reads inside the clean intron body (`body`) and reads spliced exactly out of the
## intron (`spl`, CIGAR N gaps). IR ratio = retained / (retained + spliced), retained = mean of the two boundaries.
## Only "clean" introns are scored (boundaries and body not covered by an exon of any other transcript, so an alternative
## exon is not read as retention). Targets: all introns of cancer genes plus introns with enough spliced-read coverage in a
## sample of cells. Counting is one BAM pass per cell (sc_ir_count_cell), summarised per patient by sc_ir_events.

SC_IR_FLANK <- 6L        ## bases either side of the boundary a read must cover to count as unspliced across it
SC_IR_CLEAN_FLANK <- 12L ## intron-side bases at each end that no other exon may touch
SC_IR_TRIM <- 20L        ## body = intron without 20 bp at each end

#' Target introns for intron retention
#'
#' Unique annotated introns (hg38 GTF) that are clean (no exon of any transcript within `SC_IR_CLEAN_FLANK` bp of either end,
#' >= 50 clean body bases) and are in a cancer gene, or have `min_sample_reads` spliced reads in >= `min_sample_cells` of the
#' sampled cells.
#'
#' @param gtf GTF (hg38, chr-prefixed)
#' @param cancer_genes gene symbols; all their clean introns are targets
#' @param cell_dirs named character: folders holding junctions.bed of the sampled cells (NULL: cancer genes only)
#' @param min_sample_reads,min_sample_cells coverage thresholds for non-cancer-gene introns
#' @param cores parallel reads
#' @return list(introns = data.table(idx, chromosome, start, end, strand, gene, cancer, sample_reads, sample_cells),
#'   pieces = data.table(idx, chromosome, start, end) clean body pieces)
#' @export
#' @author Stanley Clarke
sc_ir_targets <- function(gtf, cancer_genes, cell_dirs = NULL, min_sample_reads = 50L, min_sample_cells = 5L, cores = 4) {
    ex <- sc_gtf_exons(gtf)
    introns <- sc_gtf_introns(ex)
    ann <- unique(introns[, .(chromosome, start, end, strand, gene)], by = c("chromosome", "start", "end"))
    ann <- ann[end - start + 1L >= 2L * SC_IR_TRIM + 50L]
    exu <- GenomicRanges::reduce(GenomicRanges::GRanges(ex$chromosome, IRanges::IRanges(ex$start, ex$end)))
    left <- GenomicRanges::GRanges(ann$chromosome, IRanges::IRanges(ann$start, ann$start + SC_IR_CLEAN_FLANK - 1L))
    right <- GenomicRanges::GRanges(ann$chromosome, IRanges::IRanges(ann$end - SC_IR_CLEAN_FLANK + 1L, ann$end))
    ann[, clean := !(GenomicRanges::countOverlaps(left, exu) > 0 | GenomicRanges::countOverlaps(right, exu) > 0)]
    ann <- ann[clean == TRUE][, clean := NULL]
    ann[, idx := .I]
    body <- GenomicRanges::GRanges(ann$chromosome, IRanges::IRanges(ann$start + SC_IR_TRIM, ann$end - SC_IR_TRIM))
    body$idx <- ann$idx
    pcs <- GenomicRanges::setdiff(GenomicRanges::GRanges(ann$chromosome, IRanges::IRanges(ann$start + SC_IR_TRIM, ann$end - SC_IR_TRIM)), exu)
    ## setdiff merges overlapping introns: re-assign pieces to introns by overlap
    ov <- GenomicRanges::findOverlaps(pcs, body)
    pieces <- data.table::data.table(idx = body$idx[S4Vectors::subjectHits(ov)],
                                     chromosome = as.character(GenomicRanges::seqnames(pcs))[S4Vectors::queryHits(ov)],
                                     start = pmax(GenomicRanges::start(pcs)[S4Vectors::queryHits(ov)], GenomicRanges::start(body)[S4Vectors::subjectHits(ov)]),
                                     end = pmin(GenomicRanges::end(pcs)[S4Vectors::queryHits(ov)], GenomicRanges::end(body)[S4Vectors::subjectHits(ov)]))
    pieces <- pieces[end >= start]
    keep_len <- pieces[, .(clean_len = sum(end - start + 1L)), by = idx][clean_len >= 50L]
    ann <- ann[idx %in% keep_len$idx]
    ann[keep_len, on = "idx", clean_len := i.clean_len]
    ann[, cancer := toupper(gene) %in% toupper(cancer_genes)]
    ann[, `:=`(sample_reads = 0L, sample_cells = 0L)]
    if (length(cell_dirs)) {
        j <- sc_read_junctions(cell_dirs, cores = cores)
        s <- j[, .(sample_reads = sum(count), sample_cells = .N), by = .(chromosome, start, end)]
        ann[s, on = .(chromosome, start, end), `:=`(sample_reads = i.sample_reads, sample_cells = i.sample_cells)]
    }
    ann <- ann[cancer | (sample_reads >= min_sample_reads & sample_cells >= min_sample_cells)]
    ann[, idx_old := idx]
    ann[, idx := .I]
    pieces <- merge(pieces, ann[, .(idx_old, idx)], by.x = "idx", by.y = "idx_old")
    pieces[, idx := idx.y][, idx.y := NULL]
    ann[, idx_old := NULL]
    list(introns = ann[], pieces = pieces[, .(idx, chromosome, start, end)])
}

#' Intron retention counts of one cell
#'
#' @param bam coordinate-sorted, indexed STAR BAM (chr-prefixed)
#' @param tg sc_ir_targets
#' @param min_mapq 255 = STAR unique alignments only
#' @return data.table(idx, ret_l, ret_r, body, spl) of the targets with any read, attr "n_reads" (alignments examined)
#' @export
#' @author Stanley Clarke
sc_ir_count_cell <- function(bam, tg, min_mapq = 255L) {
    it <- tg$introns[paste0("chr", chromosome) %in% GenomeInfoDb::seqlevels(Rsamtools::BamFile(bam))]
    tg$pieces <- tg$pieces[idx %in% it$idx]
    reg <- GenomicRanges::reduce(GenomicRanges::GRanges(paste0("chr", it$chromosome), IRanges::IRanges(it$start - SC_IR_FLANK, it$end + SC_IR_FLANK)),
                                 min.gapwidth = 1000L)
    param <- Rsamtools::ScanBamParam(which = reg, mapqFilter = min_mapq,
                                     flag = Rsamtools::scanBamFlag(isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE, isUnmappedQuery = FALSE),
                                     what = c("qname", "flag"))
    ga <- GenomicAlignments::readGAlignments(bam, param = param, use.names = TRUE)
    ## a read spanning two merged query regions is returned twice
    ga <- ga[!duplicated(paste(names(ga), S4Vectors::mcols(ga)$flag, GenomicRanges::start(ga)))]
    out <- data.table::data.table(idx = integer(0), ret_l = integer(0), ret_r = integer(0), body = integer(0), spl = integer(0))
    if (!length(ga)) { attr(out, "n_reads") <- 0L; return(out) }
    chr <- sub("^chr", "", as.character(GenomicRanges::seqnames(ga)))
    cg <- GenomicAlignments::cigar(ga)
    mb <- GenomicAlignments::cigarRangesAlongReferenceSpace(cg, pos = GenomicRanges::start(ga), ops = c("M", "=", "X"))
    nb <- lengths(mb)
    blk <- GenomicRanges::GRanges(rep(chr, nb), unlist(mb, use.names = FALSE))
    ## unspliced across the boundary: one aligned block covers `flank` bases on each side
    lq <- GenomicRanges::GRanges(it$chromosome, IRanges::IRanges(it$start - SC_IR_FLANK, it$start + SC_IR_FLANK - 1L))
    rq <- GenomicRanges::GRanges(it$chromosome, IRanges::IRanges(it$end - SC_IR_FLANK + 1L, it$end + SC_IR_FLANK))
    ret_l <- GenomicRanges::countOverlaps(lq, blk, type = "within")
    ret_r <- GenomicRanges::countOverlaps(rq, blk, type = "within")
    pc <- tg$pieces
    pq <- GenomicRanges::GRanges(pc$chromosome, IRanges::IRanges(pc$start, pc$end))
    body <- as.integer(tapply(GenomicRanges::countOverlaps(pq, blk), factor(pc$idx, levels = it$idx), sum))
    body[is.na(body)] <- 0L
    ## spliced reads: N gaps = introns (first / last intronic base)
    nr <- GenomicAlignments::cigarRangesAlongReferenceSpace(cg, pos = GenomicRanges::start(ga), ops = "N")
    nn <- lengths(nr)
    spl <- integer(nrow(it))
    if (sum(nn)) {
        nrr <- unlist(nr, use.names = FALSE)
        g <- data.table::data.table(chromosome = rep(chr, nn), start = IRanges::start(nrr), end = IRanges::end(nrr))[, .N, by = .(chromosome, start, end)]
        m <- g[it[, .(chromosome, start, end, idx)], on = .(chromosome, start, end), nomatch = 0L]
        spl[match(m$idx, it$idx)] <- m$N
    }
    out <- data.table::data.table(idx = it$idx, ret_l = ret_l, ret_r = ret_r, body = body, spl = spl)[ret_l > 0 | ret_r > 0 | body > 0 | spl > 0]
    attr(out, "n_reads") <- length(ga)
    out[]
}

#' IR ratio statistics of one intron over cells
#'
#' @param ret_l,ret_r,spl per covering cell counts
#' @return list(k, n) per cell: retained (mean of the two boundaries) and retained + spliced
#' @export
sc_ir_kn <- function(ret_l, ret_r, spl) {
    k <- (ret_l + ret_r) / 2
    list(k = k, n = k + spl)
}

#' Intron retention events of a patient
#'
#' Same columns as sc_splice_events (event_type "intron_retention", type "intron_retention"): PSI = IR ratio (pooled over the
#' cells covering the intron, bootstrap CI over cells), cells using (>= 2 boundary reads, i.e. >= 1 retained read) of the cells
#' covering (>= `cover_reads` retained + spliced reads), clone enrichment (per-cell IR ratio, Wilcoxon), cancer gene, recurrence in
#' the other patients and the normal-cell panel, flags and the recommended tier. A retained intron must have boundary reads at both
#' ends in >= `min_using_cells` cells and body reads (so a mapping artefact at one end is not scored).
#'
#' @param patient_id patient id
#' @param counts data.table(rna_id, idx, ret_l, ret_r, body, spl) of the patient's cells (cells with no read at a target are absent)
#' @param cells rna ids of the patient's cells with RNA
#' @param tg sc_ir_targets
#' @param cl cluster junction table of sc_export_rna_splicing (to link an IR event to the sashimi cluster of its intron)
#' @param cancer cancer gene symbols (upper case); coding, tx_rank, ex, introns: as sc_splice_events
#' @param others list patient -> data.table(idx, ret, spl, cells_using) of the other patients (recurrence)
#' @param cg named character rna_id -> DNA clone
#' @param normal_counts counts (as `counts`) of the normal / non-tumour cells of every patient (the normal panel)
#' @return list(events, cells) like sc_splice_events (rows of intron retention only)
#' @export
#' @author Stanley Clarke
sc_ir_events <- function(patient_id, counts, cells, tg, cl = NULL, cancer = character(0), coding = character(0), others = list(), cg = NULL,
                         min_reads = 10L, min_using_cells = 3L, cover_reads = 3L, min_cover_cells = 10L, min_psi = 0.05, max_events = 500L,
                         common_patients = 3L, n_boot = 200L, seed = 1L, normal_counts = NULL) {
    set.seed(seed)
    it <- tg$introns
    cn <- data.table::copy(counts)[rna_id %in% cells]
    cn[, `:=`(k = (ret_l + ret_r) / 2)]
    cn[, n := k + spl]
    ## candidate introns: retained reads at both ends and in the body
    cov <- cn[n >= cover_reads, .(cells_covering = .N, cluster_reads = sum(n), reads = sum(k), spliced = sum(spl),
                                  cells_using = sum(ret_l + ret_r >= 2), end_l = sum(ret_l), end_r = sum(ret_r), body = sum(body)), by = idx]
    cov[, psi := reads / cluster_reads]
    cand <- cov[reads >= min_reads & cells_using >= min_using_cells & psi >= min_psi & end_l >= 3 & end_r >= 3 & body > 0]
    empty <- list(events = data.table::data.table(), cells = list())
    if (!nrow(cand)) return(empty)
    cand <- merge(cand, it, by = "idx")
    cand[, cancer_gene := !is.na(gene) & toupper(gene) %in% cancer]
    cand[, art_gene := sc_splice_artefact(gene)]
    cand[, prio := cancer_gene * 1e5 + (!art_gene) * 1e4 + pmin(reads, 9999)]
    data.table::setorder(cand, -prio)
    cand <- cand[seq_len(min(.N, max_events))]
    ## sashimi cluster of the intron (the annotated junction), when it passed the exporter's junction filters
    cand[, cluster_id := NA_character_]
    if (!is.null(cl)) cand[cl, on = .(chromosome, start, end), cluster_id := i.cluster_id]
    cand[, event_id := paste0("ev:ir:", chromosome, ":", start, "-", end)]
    ## recurrence in the other patients
    cand[, `:=`(n_other = length(others), n_other_with = 0L, max_other = NA_real_, other_patients = "")]
    if (length(others)) {
        ot <- data.table::rbindlist(lapply(names(others), function(q) others[[q]][, patient := q]), fill = TRUE)
        ot <- ot[idx %in% cand$idx & ret + spl >= 5]
        ot[, opsi := ret / (ret + spl)]
        r <- ot[, .(max_other = max(opsi), n_other_with = sum(ret >= 5 & cells_using >= 2), other_patients = paste(patient[ret >= 5 & cells_using >= 2], collapse = ",")), by = idx]
        cand[r, on = "idx", `:=`(max_other = i.max_other, n_other_with = i.n_other_with, other_patients = i.other_patients)]
    }
    cand[, dpsi_other := ifelse(!is.na(max_other), psi - max_other, NA_real_)]
    cand[, n_patients_with := n_other_with + 1L]
    cells_of <- split(cn[idx %in% cand$idx & n >= cover_reads], by = "idx", keep.by = TRUE)
    cell_data <- vector("list", nrow(cand)); names(cell_data) <- cand$event_id
    res <- vector("list", nrow(cand))
    for (i in seq_len(nrow(cand))) {
        e <- cand[i]
        cc <- cells_of[[as.character(e$idx)]]
        k <- as.integer(round(cc$k)); n <- as.integer(round(cc$n))
        cell_data[[i]] <- data.table::data.table(rna_id = cc$rna_id, k = k, n = n)
        ci <- sc_splice_psi_ci(cc$k, cc$n, n_boot)
        bc <- list(clone = NA_character_, clone_n = NA_integer_, clone_psi = NA_real_, rest_psi = NA_real_, clone_p = NA_real_)
        if (length(cg)) {
            g <- unname(cg[cc$rna_id]); psi_c <- cc$k / cc$n
            best <- NULL
            for (cl_i in unique(g[!is.na(g)])) {
                inn <- !is.na(g) & g == cl_i; out <- !is.na(g) & g != cl_i
                if (sum(inn) < 5 || sum(out) < 5) next
                pw <- tryCatch(stats::wilcox.test(psi_c[inn], psi_c[out], exact = FALSE)$p.value, error = function(err) NA_real_)
                if (is.finite(pw) && (is.null(best) || pw < best$clone_p))
                    best <- list(clone = cl_i, clone_n = sum(inn), clone_psi = sum(cc$k[inn]) / sum(cc$n[inn]), rest_psi = sum(cc$k[out]) / sum(cc$n[out]), clone_p = pw)
            }
            if (!is.null(best)) bc <- best
        }
        res[[i]] <- data.table::data.table(psi_lo = ci[1], psi_hi = ci[2], clone = bc$clone, clone_n = bc$clone_n, clone_psi = bc$clone_psi,
                                           rest_psi = bc$rest_psi, clone_p = bc$clone_p)
    }
    ev <- cbind(cand, data.table::rbindlist(res))
    ev[, clone_q := stats::p.adjust(clone_p, "BH")]
    ev[, `:=`(type = "intron_retention", event_type = "intron_retention", annotated = TRUE, n_skipped = NA_integer_)]
    ev[, cells_using_n := cells_using]
    ev[, cluster_reads := as.integer(round(cluster_reads))]
    ev[, reads := as.integer(round(reads))]
    ev[, `:=`(offset = NA_integer_, where = NA_character_, exon = NA_integer_, intron = NA_integer_, frameshift = NA, ref_start = start, ref_end = end,
              exons_skipped = NA_character_, frame = "unknown")]
    ev[, summary := sprintf("Intron retained in %s (chr%s:%s-%s, %s bp); %s of reads unspliced in %d of %d cells",
                            ifelse(is.na(gene), "an intergenic locus", gene), chromosome, format(start, big.mark = ",", trim = TRUE), format(end, big.mark = ",", trim = TRUE),
                            format(end - start + 1L, big.mark = ",", trim = TRUE), ifelse(psi > 0 & psi < 0.005, "<1%", paste0(round(100 * psi), "%")), cells_using, cells_covering)]
    ev[, `:=`(finding = NA_character_, finding_kind = NA_character_, finding_score = NA_real_, finding_severity = NA_character_)]
    ev[, flag_artefact := art_gene]
    ev[, flag_short := FALSE]
    ev[, flag_common := n_other_with >= common_patients & !is.na(max_other) & abs(psi - max_other) < 0.15]
    ev[, flag_few := cells_covering < min_cover_cells]
    ## normal-cell panel: IR ratio of the same intron in the normal cells of all patients
    ev[, `:=`(normal_reads = NA_integer_, normal_cells = NA_integer_, normal_psi = NA_real_)]
    if (!is.null(normal_counts) && nrow(normal_counts)) {
        nc <- normal_counts[idx %in% ev$idx][, `:=`(k = (ret_l + ret_r) / 2)][, n := k + spl][n >= 1]
        nn <- nc[, .(normal_reads = as.integer(round(sum(k))), normal_cells = sum(ret_l + ret_r >= 2), normal_psi = sum(k) / sum(n)), by = idx]
        ev[nn, on = "idx", `:=`(normal_reads = i.normal_reads, normal_cells = i.normal_cells, normal_psi = i.normal_psi)]
    }
    ev <- sc_splice_apply_filters(ev)
    ev[, differential := !is.na(dpsi_other) & dpsi_other >= 0.5]
    ev[, enriched := !is.na(clone_q) & clone_q < 0.1 & abs(clone_psi - rest_psi) >= 0.1]
    ev[, finding_notable := FALSE]
    ev[, reason := mapply(function(a, b, c) paste(c(if (a) "retained intron in cancer gene", if (b) "clone-enriched", if (c) "patient-specific"), collapse = ", "),
                          cancer_gene & psi >= 0.2, enriched, differential)]
    ev[, recommended := reason != "" & !flag_artefact & !flag_common & !flag_few & !flag_normal & !flag_repeat & cells_using >= min_using_cells & reads >= min_reads]
    ev[, patient := patient_id]
    ev[, intron_idx := idx]
    drop <- intersect(c("idx", "spliced", "end_l", "end_r", "body", "art_gene", "prio", "clean_len", "sample_reads", "sample_cells", "cancer"), names(ev))
    ev[, (drop) := NULL]
    list(events = ev[], cells = cell_data)
}
