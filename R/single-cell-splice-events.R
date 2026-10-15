## Event-level splicing for gOS: one row per spliced junction (event) with plain-language summary,
## explicit units, clone enrichment, recurrence, a "recommended" tier and links to nearby SNVs.
## Called by sc_export_rna_splicing; the pieces are exported for tests.

#' Event type of a junction (LeafCutter-style cluster member)
#'
#' @param type junction class of sc_junction_types
#' @param strand "+" / "-"
#' @return one of exon_skip, alt_5, alt_3, novel_junction, annotated_alt
#' @export
sc_splice_event_type <- function(type, strand) {
    ifelse(type == "exon_skip", "exon_skip",
    ifelse(type == "novel_acceptor", "alt_3",
    ifelse(type == "novel_donor", "alt_5",
    ifelse(type %in% c("novel", "novel_combination"), "novel_junction", "annotated_alt"))))
}

#' Plain-language one-line summary of a splice event
#'
#' @param gene,type,note,psi,n_using,frame_text see sc_export_rna_splicing; note is sc_splice_site_note
#' @param n_cover cells covering the locus
#' @return character(1)
#' @export
sc_splice_event_summary <- function(gene, type, note, psi, n_using, n_cover, chromosome, start, end) {
    g <- if (is.na(gene)) "an intergenic locus" else gene
    fs <- if (isTRUE(note$frameshift)) "frameshift" else if (identical(note$frameshift, FALSE)) "in frame" else NULL
    where <- switch(type,
        novel_acceptor = , novel_donor = {
            site <- if (type == "novel_acceptor") "acceptor" else "donor"
            if (identical(note$where, "distant") && is.finite(note$offset))
                sprintf("Novel %s %s kb from the annotated one in %s", site, format(round(note$offset / 1000, 1), nsmall = 1), g)
            else if (is.finite(note$offset) && identical(note$where, "intron") && is.finite(note$intron))
                sprintf("Novel %s %d bp into intron %d of %s", site, note$offset, note$intron, g)
            else if (is.finite(note$offset) && is.finite(note$exon))
                sprintf("Novel %s %d bp inside exon %d of %s", site, note$offset, note$exon, g)
            else sprintf("Novel %s in %s", site, g)
        },
        exon_skip = {
            ex <- note$exons
            what <- if (length(ex) == 1) sprintf("exon %d", ex) else if (length(ex) > 1) sprintf("exons %d-%d", min(ex), max(ex)) else "one or more exons"
            sprintf("Skips %s of %s", what, g)
        },
        novel_combination = sprintf("Novel pairing of annotated splice sites in %s", g),
        novel = sprintf("Novel junction in %s (chr%s:%s-%s)", g, chromosome, format(start, big.mark = ",", trim = TRUE), format(end, big.mark = ",", trim = TRUE)),
        sprintf("Annotated alternative junction in %s", g))
    paste0(where, "; used by ", if (is.finite(psi) && psi > 0 && psi < 0.005) "<1" else round(100 * psi), "% of reads in ", n_using, " of ", n_cover, " cells",
           if (length(fs)) paste0("; ", fs) else "")
}

## percentile bootstrap (over cells) CI of the pooled PSI sum(k) / sum(n)
sc_splice_psi_ci <- function(k, n, B = 200L) {
    m <- length(k)
    if (m < 2) return(c(NA_real_, NA_real_))
    idx <- matrix(sample.int(m, m * B, replace = TRUE), nrow = m)
    ps <- colSums(matrix(k[idx], nrow = m)) / pmax(1, colSums(matrix(n[idx], nrow = m)))
    unname(stats::quantile(ps, c(0.025, 0.975)))
}

#' Signed position of a genomic base relative to a splice site
#'
#' Intron coordinates (start, end = first / last intronic base). Returns the offset
#' in the transcript direction: positive = intronic (+1 is the first intron base of a donor,
#' for an acceptor the intron side is counted back from the last intron base, 1 = last),
#' <= 0 on the exon side (0 = the exon base next to the intron).
#'
#' @param q genomic position(s)
#' @param start,end intron coordinates of the junction
#' @param strand "+" / "-"
#' @param site "donor" or "acceptor"
#' @return integer offsets
#' @export
sc_splice_site_offset <- function(q, start, end, strand, site) {
    plus <- identical(strand, "+")
    if (site == "donor") (if (plus) q - start + 1L else end - q + 1L) else (if (plus) end - q + 1L else q - start + 1L)
}

#' Region class of an SNV relative to a splice site
#'
#' @param o offset of sc_splice_site_offset
#' @param site "donor" or "acceptor"
#' @return canonical (+-1-2 intronic), motif (extended motif: donor exon -3..-1, intron +3..+6;
#'   acceptor intron -3..-20 (polypyrimidine tract) and the first exon base), near (otherwise within the window)
#' @export
sc_splice_site_region <- function(o, site) {
    if (site == "donor")
        ifelse(o >= 1 & o <= 2, "canonical", ifelse((o >= 3 & o <= 6) | (o >= -2 & o <= 0), "motif", "near"))
    else
        ifelse(o >= 1 & o <= 2, "canonical", ifelse((o >= 3 & o <= 20) | o == 0, "motif", "near"))
}

## the splice sites of an intron (start, end) on a strand: genomic positions of donor and acceptor
sc_splice_sites <- function(start, end, strand) {
    if (identical(strand, "+")) c(donor = start, acceptor = end) else c(donor = end, acceptor = start)
}

#' SNV sites from a gOS snv_matrix.json
#'
#' @param file snv_matrix.json (format "compact": variants[], cells -> [[variant index (0-based), ref reads, alt reads, genotype]])
#' @return list(variants = data.table(index (0-based), id, chromosome, position, ref, alt, category, gene, consequence, protein, driver), cells = raw cells list)
#' @export
sc_read_snv_matrix <- function(file) {
    d <- jsonlite::fromJSON(file, simplifyVector = FALSE)
    vl <- d$variants
    chr1 <- function(f, d0 = NA_character_) vapply(vl, function(x) { y <- x[[f]]; if (is.null(y) || is.na(y)) d0 else as.character(y) }, "")
    id <- chr1("id")
    s <- data.table::tstrsplit(sub("^chr", "", id), "_", fixed = TRUE)
    v <- data.table::data.table(index = seq_along(vl) - 1L, id = id, chromosome = s[[1]], position = as.integer(s[[2]]), ref = s[[3]], alt = s[[4]],
                                category = chr1("category"), gene = chr1("gene"), consequence = chr1("consequence"), protein = chr1("protein"),
                                driver = vapply(vl, function(x) isTRUE(x$driver), TRUE))
    list(variants = v, cells = d$cells)
}

#' Genotype calls (ref / alt reads) of chosen cells and variants of sc_read_snv_matrix
#'
#' @param cells cells element of sc_read_snv_matrix
#' @param keep_cells cell ids
#' @param keep_index variant indexes (0-based)
#' @return data.table(cell, index, ref_reads, alt_reads)
#' @export
sc_snv_calls <- function(cells, keep_cells, keep_index) {
    data.table::rbindlist(lapply(intersect(names(cells), keep_cells), function(cid) {
        l <- cells[[cid]]
        if (!length(l)) return(NULL)
        idx <- vapply(l, function(e) as.integer(e[[1]]), 1L)
        sel <- which(idx %in% keep_index)
        if (!length(sel)) return(NULL)
        data.table::data.table(cell = cid, index = idx[sel], ref_reads = vapply(l[sel], function(e) as.integer(e[[2]]), 1L),
                               alt_reads = vapply(l[sel], function(e) as.integer(e[[3]]), 1L))
    }))
}

#' Genome-wide SNVs near the splice sites of events
#'
#' @param ev data.table(event_id, chromosome, strand, start, end, ref_start, ref_end) intron coordinates; ref_* = annotated intron the
#'   new site replaces (NA if none)
#' @param variants variants of sc_read_snv_matrix
#' @param window bp either side of a splice site
#' @param max_shift largest distance between the new and the annotated site for the segment between them to be searched
#' @return data.table(event_id, index, site, offset (signed, + intron / <= 0 exon), distance, region, in_shift)
#' @export
sc_splice_snv_hits <- function(ev, variants, window = 50L, max_shift = 300L) {
    out <- list()
    vk <- data.table::copy(variants)[, .(index, chromosome, position)]
    for (site_nm in c("donor", "acceptor")) {
        pos <- vapply(seq_len(nrow(ev)), function(i) sc_splice_sites(ev$start[i], ev$end[i], ev$strand[i])[[site_nm]], 1L)
        s <- data.table::data.table(event_id = ev$event_id, chromosome = ev$chromosome, strand = ev$strand, start = ev$start, end = ev$end,
                                    site_pos = as.integer(pos), lo = as.integer(pos) - window, hi = as.integer(pos) + window)
        h <- vk[s, on = .(chromosome, position >= lo, position <= hi), nomatch = 0L,
                .(event_id = i.event_id, index, position = x.position, strand = i.strand, start = i.start, end = i.end)]
        if (nrow(h)) {
            h[, offset := mapply(function(q, st, en, sd) sc_splice_site_offset(q, st, en, sd, site_nm), position, start, end, strand)]
            ## exon side: offsets <= 0 carry no intron sign; distance is |o| on the intron side, 1 - o on the exon side
            h[, `:=`(site = site_nm, distance = ifelse(offset >= 1L, offset, 1L - offset))]
            h[, region := sc_splice_site_region(offset, site_nm)]
            h[, `:=`(strand = NULL, start = NULL, end = NULL, position = NULL)]
            out[[site_nm]] <- h
        }
    }
    hits <- if (length(out)) data.table::rbindlist(out) else data.table::data.table(event_id = character(0), index = integer(0), offset = integer(0), site = character(0),
                                                                                  distance = integer(0), region = character(0))
    ## the segment between the new and the annotated site (retained intron / extended exon)
    sh <- ev[!is.na(ref_start) & (start != ref_start | end != ref_end) & abs(start - ref_start) + abs(end - ref_end) <= max_shift]
    if (nrow(sh)) {
        segs <- data.table::data.table(event_id = sh$event_id, chromosome = sh$chromosome,
                                       lo = ifelse(sh$start != sh$ref_start, pmin(sh$start, sh$ref_start), pmin(sh$end, sh$ref_end)),
                                       hi = ifelse(sh$start != sh$ref_start, pmax(sh$start, sh$ref_start), pmax(sh$end, sh$ref_end)))
        sg <- vk[segs, on = .(chromosome, position >= lo, position <= hi), nomatch = 0L, .(event_id = i.event_id, index, position = x.position)]
        if (nrow(sg)) {
            sg <- merge(sg, ev[, .(event_id, start, end, strand, ref_start, ref_end)], by = "event_id")
            ## offset to the nearest site of the new junction
            sg[, c("d_off", "d_site") := {
                dn <- mapply(function(q, st, en, sd) sc_splice_site_offset(q, st, en, sd, "donor"), position, start, end, strand)
                ac <- mapply(function(q, st, en, sd) sc_splice_site_offset(q, st, en, sd, "acceptor"), position, start, end, strand)
                dd <- ifelse(dn >= 1L, dn, 1L - dn); da <- ifelse(ac >= 1L, ac, 1L - ac)
                list(ifelse(dd <= da, dn, ac), ifelse(dd <= da, "donor", "acceptor"))
            }]
            sg <- sg[, .(event_id, index, offset = d_off, site = d_site)]
            sg[, `:=`(distance = ifelse(offset >= 1L, offset, 1L - offset))]
            sg[, region := ifelse(distance <= 2 & offset >= 1L, "canonical", "shifted_segment")]
            hits <- data.table::rbindlist(list(hits, sg), fill = TRUE)
        }
    }
    if (!nrow(hits)) return(hits)
    ## one row per (event, SNV): the nearest site / most specific region
    hits[, rk := match(region, c("canonical", "motif", "near", "shifted_segment"))]
    data.table::setorder(hits, event_id, index, rk, distance)
    hits <- hits[, .SD[1], by = .(event_id, index)][, rk := NULL]
    hits[]
}

#' Test the usage of a splice event against the DNA genotype of a SNV
#'
#' Cells with RNA and DNA: carriers (>= 1 alt read) vs non-carriers (ref reads only) of the SNV;
#' per-cell PSI (cells covering the locus) by Wilcoxon, cells using the junction (>= 1 read) by Fisher.
#'
#' @param k,n per covering cell: reads of the junction, reads of the cluster
#' @param carrier logical per covering cell: TRUE carrier, FALSE non-carrier, NA no DNA call
#' @param min_each minimum carrier and non-carrier cells
#' @return list(n_carriers, n_noncarriers, psi_carriers, psi_noncarriers, delta_psi, odds_ratio, p_wilcoxon, p_fisher, p) or NULL
#' @export
sc_splice_snv_test <- function(k, n, carrier, min_each = 2L) {
    ok <- !is.na(carrier) & n > 0
    k <- k[ok]; n <- n[ok]; carrier <- carrier[ok]
    nc <- sum(carrier); nn <- sum(!carrier)
    if (nc < min_each || nn < min_each) return(NULL)
    psi <- k / n
    pw <- tryCatch(stats::wilcox.test(psi[carrier], psi[!carrier], exact = FALSE)$p.value, error = function(e) NA_real_)
    tab <- matrix(c(sum(k[carrier] > 0), sum(k[carrier] == 0), sum(k[!carrier] > 0), sum(k[!carrier] == 0)), 2)
    ft <- tryCatch(stats::fisher.test(tab), error = function(e) NULL)
    pf <- if (is.null(ft)) NA_real_ else ft$p.value
    p <- suppressWarnings(min(2 * min(pw, pf, na.rm = TRUE), 1))
    list(n_carriers = nc, n_noncarriers = nn, psi_carriers = sum(k[carrier]) / sum(n[carrier]), psi_noncarriers = sum(k[!carrier]) / sum(n[!carrier]),
         delta_psi = sum(k[carrier]) / sum(n[carrier]) - sum(k[!carrier]) / sum(n[!carrier]),
         odds_ratio = if (is.null(ft)) NA_real_ else unname(ft$estimate), p_wilcoxon = pw, p_fisher = pf, p = if (is.finite(p)) p else NA_real_)
}

#' Event-level splicing table of a patient
#'
#' One event per junction of a LeafCutter-style cluster, with: event type, annotated vs novel, PSI (pooled reads of covering
#' cells) with a bootstrap-over-cells 95% CI, cells using the junction (>= 1 read) of the cells covering the locus (>= `cover_reads`
#' cluster reads), the clone with the most different per-cell PSI (Wilcoxon, BH q over events), frame effect, cancer gene,
#' recurrence in the other patients, flags (artefact classes) and the "recommended" tier.
#'
#' Recommended: enough coverage (>= `min_cover_cells` cells covering, >= `min_using_cells` using, >= `min_reads` reads), not flagged
#' (paralog / repeat family gene, a very short intron, or the same junction at a similar PSI in >= `common_patients` other patients), and
#' one of: novel junction in a cancer gene, a previous finding of the exporter (moderate or high), clone enrichment (q < 0.1),
#' patient-specific novel junction (PSI >= 0.5 above every other patient).
#'
#' @param patient_id patient id
#' @param x per-cell junction table (rna_id, chromosome, start, end, count, cluster_id, cl_total)
#' @param cl junction table of the clusters (cluster_id, chromosome, start, end, strand, gene, annotated, type, n_skipped)
#' @param allj junction reads per patient (patient, chromosome, start, end, reads, cells)
#' @param introns,ex,tx_rank,coding,cancer see sc_export_rna_splicing
#' @param cg named character rna_id -> DNA clone (normal cells removed)
#' @param findings findings of the patient (list)
#' @return list(events = data.table, cells = list event_id -> data.table(rna_id, k, n))
#' @export
sc_splice_events <- function(patient_id, x, cl, allj, introns, ex, tx_rank, coding, cancer, cg = NULL, findings = list(), min_reads = 10L,
                             min_using_cells = 3L, cover_reads = 3L, min_cover_cells = 10L, max_events = 2000L, common_patients = 3L,
                             short_intron = 30L, min_psi = 0.05, n_boot = 200L, seed = 1L) {
    set.seed(seed)
    ## cluster-level coverage
    cc <- unique(x[cl_total >= cover_reads, .(rna_id, cluster_id, cl_total)])
    cov <- cc[, .(cells_covering = .N, cluster_reads = sum(cl_total)), by = cluster_id]
    xc <- x[cl_total >= cover_reads]
    ju <- xc[, .(reads = sum(count), cells_using = data.table::uniqueN(rna_id)), by = .(cluster_id, start, end)]
    ju <- merge(ju, cov, by = "cluster_id")
    ju[, psi := reads / cluster_reads]
    info <- unique(cl[, .(cluster_id, chromosome, strand, start, end, gene, annotated, type, n_skipped)], by = c("cluster_id", "start", "end"))
    ju <- merge(ju, info, by = c("cluster_id", "start", "end"))
    ju[, cancer_gene := !is.na(gene) & toupper(gene) %in% cancer]
    ## previous findings of the exporter
    fkey <- data.table::rbindlist(lapply(findings, function(f) if (length(f$cluster_id)) data.table::data.table(
        cluster_id = f$cluster_id, start = as.integer(f$junction$start), end = as.integer(f$junction$end), finding = f$id, finding_kind = f$kind,
        finding_score = f$score, finding_severity = f$severity) else NULL))
    ju <- if (nrow(fkey)) merge(ju, fkey[order(-finding_score)][!duplicated(paste(cluster_id, start, end))], by = c("cluster_id", "start", "end"), all.x = TRUE) else
        cbind(ju, data.table::data.table(finding = NA_character_, finding_kind = NA_character_, finding_score = NA_real_, finding_severity = NA_character_))
    ## candidates: novel junctions, plus annotated alternatives in cancer genes or findings
    cand <- ju[reads >= min_reads & cells_using >= min_using_cells & psi >= min_psi &
               (!annotated | ((cancer_gene | !is.na(finding)) & psi <= 0.95))]
    if (!nrow(cand)) return(list(events = data.table::data.table(), cells = list()))
    cand[, art_gene := sc_splice_artefact(gene)]
    cand[, prio := (!is.na(finding)) * 1e6 + (cancer_gene & !annotated) * 1e5 + (!art_gene) * 1e4 + pmin(reads, 9999)]
    data.table::setorder(cand, -prio)
    cand <- cand[seq_len(min(.N, max_events))]
    cand[, event_id := paste0("ev:", cluster_id, ":", start, "-", end)]

    ## cross-patient recurrence of the junction
    pj <- merge(allj[, .(patient, chromosome, start, end, reads, cells)], unique(cl[, .(chromosome, start, end, cluster_id)]), by = c("chromosome", "start", "end"))
    pj[, ctot := sum(reads), by = .(patient, cluster_id)]
    cp <- unique(pj[patient != patient_id & ctot >= 30, .(patient, cluster_id, ctot)])
    rec <- merge(cand[, .(cluster_id, start, end)], cp, by = "cluster_id", allow.cartesian = TRUE)
    rec <- merge(rec, pj[, .(patient, cluster_id, start, end, jreads = reads, jcells = cells)], by = c("patient", "cluster_id", "start", "end"), all.x = TRUE)
    rec[is.na(jreads), `:=`(jreads = 0L, jcells = 0L)]
    rec[, opsi := jreads / ctot]
    rec <- rec[, .(n_other = .N, max_other = max(opsi), n_other_with = sum(jreads >= 5 & jcells >= 2),
                   other_patients = paste(patient[jreads >= 5 & jcells >= 2], collapse = ",")), by = .(cluster_id, start, end)]
    cand <- merge(cand, rec, by = c("cluster_id", "start", "end"), all.x = TRUE)
    cand[is.na(n_other), `:=`(n_other = 0L, n_other_with = 0L, other_patients = "")]
    cand[, dpsi_other := ifelse(n_other > 0, psi - max_other, NA_real_)]
    cand[, n_patients_with := n_other_with + 1L]
    data.table::setorder(cand, -prio)

    ## per-event cell vectors, CI, clone enrichment, site note
    cells_of <- split(cc, cc$cluster_id)
    rows_of <- split(xc[cand[, .(cluster_id, start, end)], on = .(cluster_id, start, end), nomatch = 0L, .(cluster_id, start, end, rna_id, count)],
                     by = c("cluster_id", "start", "end"), drop = TRUE, keep.by = FALSE)
    cell_data <- vector("list", nrow(cand)); names(cell_data) <- cand$event_id
    res <- vector("list", nrow(cand))
    for (i in seq_len(nrow(cand))) {
        e <- cand[i]
        cells <- cells_of[[e$cluster_id]]
        rr <- rows_of[[paste(e$cluster_id, e$start, e$end, sep = ".")]]
        k <- integer(nrow(cells))
        if (!is.null(rr)) { m <- match(rr$rna_id, cells$rna_id); k[m[!is.na(m)]] <- rr$count[!is.na(m)] }
        n <- cells$cl_total
        cell_data[[i]] <- data.table::data.table(rna_id = cells$rna_id, k = k, n = n)
        ci <- sc_splice_psi_ci(k, n, n_boot)
        ## clone with the most different per-cell PSI
        bc <- list(clone = NA_character_, clone_n = NA_integer_, clone_psi = NA_real_, rest_psi = NA_real_, clone_p = NA_real_)
        if (length(cg)) {
            g <- unname(cg[cells$rna_id]); psi_c <- k / n
            best <- NULL
            for (cl_i in unique(g[!is.na(g)])) {
                inn <- !is.na(g) & g == cl_i; out <- !is.na(g) & g != cl_i
                if (sum(inn) < 5 || sum(out) < 5) next
                pw <- tryCatch(stats::wilcox.test(psi_c[inn], psi_c[out], exact = FALSE)$p.value, error = function(e) NA_real_)
                if (is.finite(pw) && (is.null(best) || pw < best$clone_p))
                    best <- list(clone = cl_i, clone_n = sum(inn), clone_psi = sum(k[inn]) / sum(n[inn]), rest_psi = sum(k[out]) / sum(n[out]), clone_p = pw)
            }
            if (!is.null(best)) bc <- best
        }
        ## frame / site note
        note <- list(offset = NA_integer_, where = NA_character_, exon = NA_integer_, exons = NULL, frameshift = NA, intron = NA_integer_, ref_start = NA_integer_, ref_end = NA_integer_)
        if (e$type != "annotated" && e$type != "novel" && e$type != "novel_combination") {
            note <- tryCatch(sc_splice_site_note(e[, .(chromosome, start, end, strand, type)], introns, ex, tx_rank), error = function(err) note)
            if (!(e$gene %in% coding)) note$frameshift <- NA
        }
        res[[i]] <- data.table::data.table(
            cells_using_n = sum(k > 0), psi_lo = ci[1], psi_hi = ci[2], clone = bc$clone, clone_n = bc$clone_n, clone_psi = bc$clone_psi, rest_psi = bc$rest_psi,
            clone_p = bc$clone_p, offset = note$offset %||% NA_integer_, where = note$where %||% NA_character_, exon = note$exon %||% NA_integer_,
            intron = note$intron %||% NA_integer_, frameshift = if (is.na(note$frameshift)) NA else note$frameshift,
            ref_start = note$ref_start %||% NA_integer_, ref_end = note$ref_end %||% NA_integer_,
            exons_skipped = if (length(note$exons)) paste(note$exons, collapse = ",") else NA_character_, summary = sc_splice_event_summary(
                e$gene, e$type, note, e$psi, e$cells_using, e$cells_covering, e$chromosome, e$start, e$end))
    }
    ev <- cbind(cand, data.table::rbindlist(res))
    ev[, clone_q := stats::p.adjust(clone_p, "BH")]
    ev[, event_type := sc_splice_event_type(type, strand)]
    ev[, frame := ifelse(is.na(frameshift), "unknown", ifelse(frameshift, "frameshift", "in_frame"))]
    ## flags
    ev[, flag_artefact := art_gene]
    ev[, flag_short := (end - start + 1L) <= short_intron]
    ev[, flag_common := n_other_with >= common_patients & !is.na(max_other) & abs(psi - max_other) < 0.15]
    ev[, flag_few := cells_covering < min_cover_cells]
    ev[, flags := mapply(function(a, s, c, f) paste(c(if (a) "paralog/repeat family", if (s) "very short intron", if (c) "same junction at a similar PSI in other patients",
                                                       if (f) "few covering cells"), collapse = "; "), flag_artefact, flag_short, flag_common, flag_few)]
    ev[, differential := !is.na(dpsi_other) & dpsi_other >= 0.5 & !annotated]
    ev[, enriched := !is.na(clone_q) & clone_q < 0.1 & abs(clone_psi - rest_psi) >= 0.1]
    ev[, finding_notable := !is.na(finding_severity) & finding_severity %in% c("high", "moderate")]
    ev[, reason := mapply(function(a, b, c, d) paste(c(if (a) "novel in cancer gene", if (b) "clone-enriched", if (c) "patient-specific", if (d) "finding"), collapse = ", "),
                          !annotated & cancer_gene, enriched, differential, finding_notable)]
    ev[, recommended := reason != "" & !flag_artefact & !flag_short & !flag_common & !flag_few & cells_using >= min_using_cells & reads >= min_reads]
    ev[, patient := patient_id]
    list(events = ev[], cells = cell_data)
}

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

#' Link events to SNVs near their splice sites and test co-occurrence
#'
#' @param ev events of sc_splice_events
#' @param cells per-event per-cell data (sc_splice_events$cells)
#' @param snv_file snv_matrix.json of the patient
#' @param cell_map named character rna_id -> DNA cell id (NA for cells without DNA)
#' @param cg named character rna_id -> clone
#' @param window bp either side of a splice site
#' @param max_links per event kept
#' @return data.table of links (one row per event x SNV) with the test; q BH over tested pairs
#' @export
sc_splice_snv_links <- function(ev, cells, snv_file, cell_map, cg = NULL, window = 50L, max_links = 5L) {
    if (!nrow(ev) || !length(snv_file) || !file.exists(snv_file)) return(NULL)
    sn <- sc_read_snv_matrix(snv_file)
    sub <- ev[, .(event_id, chromosome, strand, start, end, ref_start = as.integer(ref_start), ref_end = as.integer(ref_end))]
    hits <- sc_splice_snv_hits(sub, sn$variants, window)
    if (!nrow(hits)) return(NULL)
    hits <- merge(hits, sn$variants, by = "index")
    ## genotype lookup
    calls <- sc_snv_calls(sn$cells, unique(stats::na.omit(unlist(cell_map))), unique(hits$index))
    if (!nrow(calls)) return(NULL)
    data.table::setkey(calls, index, cell)
    dna_of <- cell_map[!is.na(cell_map)]
    tests <- lapply(seq_len(nrow(hits)), function(i) {
        h <- hits[i]
        cd <- cells[[h$event_id]]
        dna <- unname(dna_of[cd$rna_id])
        g <- calls[.(h$index)]
        alt <- g$alt_reads[match(dna, g$cell)]; ref <- g$ref_reads[match(dna, g$cell)]
        carrier <- ifelse(is.na(dna) | is.na(alt), NA, ifelse(alt >= 1, TRUE, ifelse(ref >= 1, FALSE, NA)))
        tt <- sc_splice_snv_test(cd$k, cd$n, carrier)
        cars <- cd$rna_id[which(carrier %in% TRUE)]
        top_clone <- NA_character_; top_frac <- NA_real_
        if (length(cg) && length(cars)) { tb <- sort(table(cg[cars]), decreasing = TRUE); if (length(tb)) { top_clone <- names(tb)[1]; top_frac <- as.numeric(tb[1]) / length(cars) } }
        base <- list(n_carriers_all = sum(carrier %in% TRUE), n_noncarriers_all = sum(carrier %in% FALSE), carrier_clone = top_clone, carrier_clone_frac = top_frac,
                     carriers = cars[order(-cd$k[match(cars, cd$rna_id)])][seq_len(min(length(cars), 12))])
        if (is.null(tt)) tt <- list(n_carriers = base$n_carriers_all, n_noncarriers = base$n_noncarriers_all, psi_carriers = NA_real_, psi_noncarriers = NA_real_,
                                    delta_psi = NA_real_, odds_ratio = NA_real_, p_wilcoxon = NA_real_, p_fisher = NA_real_, p = NA_real_)
        c(tt, base[c("carrier_clone", "carrier_clone_frac")], list(carriers = list(base$carriers)))
    })
    tt <- data.table::rbindlist(lapply(tests, function(l) { l$carriers <- list(l$carriers[[1]]); data.table::as.data.table(l) }), fill = TRUE)
    out <- cbind(hits, tt)
    out[, q := NA_real_]
    out[!is.na(p), q := stats::p.adjust(p, "BH")]
    out[, canonical := region == "canonical"]
    out[, rk := (region == "canonical") * 1000 + (region == "motif") * 500 + ifelse(is.na(q), 0, 100 * (1 - pmin(q, 1))) - pmin(distance, 100) / 10]
    data.table::setorder(out, event_id, -rk)
    out <- out[, head(.SD, max_links), by = event_id]
    out[, rk := NULL]
    out[]
}

#' Events as a JSON-ready list
#'
#' @param evr list(events, cells) of sc_splice_events
#' @param links sc_splice_snv_links table or NULL
#' @param n_igv cells per event listed for the IGV view
#' @export
sc_splice_events_json <- function(evr, links = NULL, n_igv = 6L) {
    ev <- evr$events
    if (!nrow(ev)) return(list())
    num <- function(v) if (length(v) && is.finite(v)) as.numeric(v) else NULL
    lapply(seq_len(nrow(ev)), function(i) {
        e <- ev[i]
        cd <- evr$cells[[e$event_id]]
        lk <- if (!is.null(links)) links[event_id == e$event_id] else NULL
        list(id = e$event_id, cluster_id = e$cluster_id, gene = e$gene, chromosome = e$chromosome, strand = e$strand, start = e$start, end = e$end,
             type = e$type, event_type = e$event_type, annotated = e$annotated, summary = e$summary,
             reads = e$reads, cluster_reads = e$cluster_reads, psi = num(e$psi), psi_lo = num(e$psi_lo), psi_hi = num(e$psi_hi),
             cells_using = e$cells_using, cells_covering = e$cells_covering,
             clone = if (!is.na(e$clone)) list(clone = e$clone, n_cells = e$clone_n, psi = num(e$clone_psi), psi_rest = num(e$rest_psi), p = num(e$clone_p), q = num(e$clone_q)) else NULL,
             frame = e$frame, site = list(offset = if (is.na(e$offset)) NULL else e$offset, where = if (is.na(e$where)) NULL else e$where,
                                          exon = if (is.na(e$exon)) NULL else e$exon, intron = if (is.na(e$intron)) NULL else e$intron,
                                          exons_skipped = if (is.na(e$exons_skipped)) NULL else as.integer(strsplit(e$exons_skipped, ",")[[1]])),
             cancer_gene = e$cancer_gene, n_patients = e$n_patients_with, other_patients = I(if (nzchar(e$other_patients)) strsplit(e$other_patients, ",")[[1]] else character(0)),
             max_psi_other = num(e$max_other), dpsi_other = num(e$dpsi_other),
             flags = I(if (nzchar(e$flags)) strsplit(e$flags, "; ")[[1]] else character(0)), recommended = e$recommended,
             reasons = I(if (nzchar(e$reason)) strsplit(e$reason, ", ")[[1]] else character(0)),
             finding_id = if (is.na(e$finding)) NULL else e$finding, finding_score = num(e$finding_score),
             loci = list(list(chromosome = e$chromosome, position = e$start - 1L), list(chromosome = e$chromosome, position = e$end + 1L)),
             cells = I(cd[order(-k)][seq_len(min(.N, n_igv))][k > 0]$rna_id),
             snv_links = if (!is.null(lk) && nrow(lk)) lapply(seq_len(nrow(lk)), function(j) { l <- lk[j]
                 list(variant_id = l$id, variant_index = l$index, chromosome = l$chromosome, position = l$position, ref = l$ref, alt = l$alt, category = l$category,
                      gene = l$gene, consequence = l$consequence, protein = l$protein, site = l$site, offset = l$offset, distance = l$distance, region = l$region,
                      canonical = l$canonical, n_carriers = l$n_carriers, n_noncarriers = l$n_noncarriers, psi_carriers = num(l$psi_carriers),
                      psi_noncarriers = num(l$psi_noncarriers), delta_psi = num(l$delta_psi), odds_ratio = num(l$odds_ratio), p_wilcoxon = num(l$p_wilcoxon),
                      p_fisher = num(l$p_fisher), p = num(l$p), q = num(l$q), carrier_clone = if (is.na(l$carrier_clone)) NULL else l$carrier_clone,
                      carrier_clone_frac = num(l$carrier_clone_frac), carriers = I(l$carriers[[1]])) }) else list())
    })
}
