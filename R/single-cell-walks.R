#' Export ecDNA / amplicon walks of a single-cell patient for gOS
#'
#' Writes `walks.json` next to the patient's other single-cell files: every
#' walk of a gWalk (nodes in walk order with chromosome, start, end, strand
#' and graph copy number; junctions between consecutive nodes typed REF when
#' the next node continues the reference, ALT otherwise, the closing junction
#' for circular walks) together with the per-cell copy number of the walk
#' (`counts`: gw_id, pair, cn, amp), the amplicon summary (`coords`: span,
#' coordinates, cgc_genes, driver gene flags, ncells, total_cn) and the
#' curation flags of `summary` (ncells_filter, cn_filter, gene_label,
#' amp_id4). Only cells with cn >= `min_cn` are listed per walk.
#'
#' With `simplify` (default) consecutive nodes joined by a reference
#' adjacency are collapsed into one interval, as `gWalk$simplify()` does in
#' the blogs, so only the ALT junctions remain; the merged node's cn is the
#' width-weighted mean of its pieces. Internal pieces shorter than 1 kb
#' (templated insertions) are then dropped (listed in the `via` of the junction
#' that skips them) and same-strand neighbours within 1 kb are merged.
#'
#' @param walks gWalk or path to an rds holding one
#' @param counts data.table or rds path with columns gw_id, pair, cn, amp
#' @param coords optional data.table / rds path (amp_coords)
#' @param summary optional data.table / rds path (amplicon_summary_dt)
#' @param out_dir patient folder of the gOS dataset
#' @param patient patient id written into the file
#' @param min_cn copy number from which a cell counts as carrying the walk
#' @param simplify collapse reference-adjacent nodes, keeping only ALT junctions
#' @param cell_ids optional gOS cell ids of the patient; walk cell ids are renamed to them when they match up to underscores / case
#' @return path of the written json (invisibly)
#' @export
sc_export_walks <- function(walks, counts, coords = NULL, summary = NULL, out_dir, patient = NULL, min_cn = 1, cell_ids = NULL, simplify = TRUE) {
    load_rds <- function(x) if (is.character(x)) readRDS(path.expand(x)) else x
    gw <- load_rds(walks)
    if (!inherits(gw, "gWalk")) stop("walks must be a gWalk (got ", paste(class(gw), collapse = "/"), ")")
    counts <- data.table::as.data.table(load_rds(counts))
    counts[, gw_id := as.character(gw_id)]
    counts[, pair := as.character(pair)]
    ## cell ids in the walk tables may be spelled differently from the gOS cell ids
    ## (e.g. MGH302_MR2_pl3_10b vs MGH302_MR_2_pl3_10b): match on a key without underscores / case
    if (!is.null(cell_ids)) {
        key <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))
        lookup <- stats::setNames(as.character(cell_ids), key(cell_ids))
        hit <- lookup[key(counts$pair)]
        n_mapped <- sum(!is.na(hit))
        if (n_mapped) counts[!is.na(hit), pair := hit[!is.na(hit)]]
        message(sprintf("sc_export_walks: %d of %d walk cells matched to gOS cell ids", length(unique(counts$pair[!is.na(hit)])), length(unique(counts$pair))))
    }
    ## older patients' coords / summary tables are per cell (no gw_id): skip them
    by_walk <- function(x, what) {
        if (is.null(x)) return(NULL)
        x <- data.table::as.data.table(load_rds(x))
        if (!"gw_id" %in% names(x)) { message("sc_export_walks: ", what, " has no gw_id column, ignored"); return(NULL) }
        x[, gw_id := as.character(gw_id)]
    }
    coords <- by_walk(coords, "coords")
    if (!is.null(summary)) summary <- data.table::as.data.table(load_rds(summary))
    from_cells <- FALSE
    ## a per-cell summary (pair, amp_id4, amp_cn; no gw_id) names amplicons per cell:
    ## turn it into a per-walk summary by matching each label to the walk its cells carry
    if (!is.null(summary) && !"gw_id" %in% names(summary) && all(c("pair", "amp_id4", "amp_cn") %in% names(summary)))
        { summary <- walk_labels_from_cell_summary(summary, counts, min_cn); from_cells <- TRUE }
    summary <- by_walk(summary, "summary")
    dt <- data.table::as.data.table(gw$dt)
    if (!"gw_id" %in% names(dt)) dt[, gw_id := walk.id]
    dt[, gw_id := as.character(gw_id)]
    grl <- gw$grl
    cells <- sort(unique(counts$pair))
    known_flags <- c("EGFR", "PDGFRA", "MYCN", "CDK4", "MDM2", "MAP3K1", "RRAS2")
    walk_json <- lapply(seq_len(nrow(dt)), function(k) {
        gr <- grl[[k]]
        m <- as.data.frame(GenomicRanges::mcols(gr))
        nodes <- data.table::data.table(
            chromosome = sub("^chr", "", as.character(GenomicRanges::seqnames(gr))),
            start = GenomicRanges::start(gr), end = GenomicRanges::end(gr),
            strand = as.character(GenomicRanges::strand(gr)),
            cn = if ("cn" %in% names(m)) as.numeric(m$cn) else NA_real_)
        n_raw <- nrow(nodes)
        if (simplify && n_raw > 1) nodes <- simplify_walk_nodes(nodes, isTRUE(dt$circular[k]))
        n <- nrow(nodes)
        pairs <- if (n > 1) cbind(seq_len(n - 1), seq_len(n - 1) + 1) else matrix(integer(0), ncol = 2)
        ## circular: closing junction last -> first (a self-junction when the walk collapsed to one interval)
        if (isTRUE(dt$circular[k]) && n >= 1) pairs <- rbind(pairs, c(n, 1))
        junctions <- lapply(seq_len(nrow(pairs)), function(p) {
            a <- nodes[pairs[p, 1]]; b <- nodes[pairs[p, 2]]
            ref <- a$chromosome == b$chromosome && a$strand == b$strand &&
                ((a$strand != "-" && b$start == a$end + 1) || (a$strand == "-" && a$start == b$end + 1))
            j <- list(from = pairs[p, 1] - 1L, to = pairs[p, 2] - 1L, type = if (isTRUE(ref)) "REF" else "ALT")
            if ("via" %in% names(b) && nzchar(b$via)) j$via <- b$via
            j
        })
        gid <- dt$gw_id[k]
        cc <- counts[gw_id == gid]
        amp <- if ("amp" %in% names(cc) && nrow(cc)) cc$amp[!is.na(cc$amp) & cc$amp != ""][1] else NA_character_
        if (length(amp) == 0) amp <- NA_character_
        carriers <- cc[is.finite(cn) & cn >= min_cn]
        co <- if (!is.null(coords)) coords[gw_id == gid][1] else NULL
        su <- if (!is.null(summary)) summary[gw_id == gid][1] else NULL
        flags <- if (!is.null(co)) names(co)[vapply(names(co), function(nm) is.logical(co[[nm]]) && isTRUE(co[[nm]]), logical(1))] else character(0)
        cgc <- if (!is.null(co) && !is.na(co$cgc_genes)) trimws(strsplit(co$cgc_genes, ",")[[1]]) else character(0)
        list(
            id = gid, walk_id = dt$walk.id[k], name = dt$name[k],
            label = if (!is.na(amp)) amp else if (!is.null(su) && "gene_label" %in% names(su) && !is.na(su$gene_label)) su$gene_label else dt$name[k],
            circular = isTRUE(dt$circular[k]),
            span = if ("wid" %in% names(dt)) dt$wid[k] else sum(nodes$end - nodes$start + 1),
            n_nodes = n,
            n_nodes_raw = n_raw,
            coordinates = if (!is.null(co)) co$coordinates else NULL,
            genes = union(cgc, intersect(flags, known_flags)),
            driver_genes = intersect(flags, known_flags),
            ncells = nrow(carriers),
            total_cn = sum(carriers$cn),
            median_cn = if (nrow(carriers)) stats::median(carriers$cn) else 0,
            max_cn = if (nrow(carriers)) max(carriers$cn) else 0,
            curated = if (!is.null(su) && "cn_filter" %in% names(su)) isTRUE(su$cn_filter) else if (from_cells) FALSE else NULL,
            gene_label = if (!is.null(su) && "gene_label" %in% names(su)) su$gene_label else NULL,
            amp_id = if (!is.null(su) && "amp_id4" %in% names(su) && !is.na(su$amp_id4)) su$amp_id4 else NULL,
            nodes = nodes[, intersect(c("chromosome", "start", "end", "strand", "cn"), names(nodes)), with = FALSE],
            junctions = junctions,
            cells = as.list(stats::setNames(carriers$cn, carriers$pair)))
    })
    ## with labels matched from a per-cell summary, the unmatched walks get "<driver gene | Other> minor k"
    ## (by decreasing cell count) instead of their raw gWalk name
    if (from_cells) {
        raw <- which(vapply(walk_json, function(w) identical(w$label, w$name), logical(1)))
        fam <- vapply(walk_json[raw], function(w) if (length(w$driver_genes)) paste(w$driver_genes, collapse = "/") else "Other", character(1))
        ord <- order(fam, -vapply(walk_json[raw], function(w) w$ncells, numeric(1)))
        k <- stats::ave(seq_along(ord), fam[ord], FUN = seq_along)
        for (i in seq_along(ord)) walk_json[[raw[ord[i]]]]$label <- paste(fam[ord[i]], "minor", k[i])
    }
    out <- list(format = "gos-sc-walks/1", patient = patient, n_walks = length(walk_json),
                n_cells = length(cells), cells = cells, walks = walk_json)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    path <- file.path(out_dir, "walks.json")
    jsonlite::write_json(out, path, auto_unbox = TRUE, digits = NA, null = "null", na = "null")
    invisible(path)
}

## Per-walk labels from a per-cell amplicon summary (pair, amp_id4 e.g. "MGH285 EGFR 1 (81)", amp_cn).
## Each label goes to the walk whose carriers best match the label's cells (Jaccard >= 0.8) and whose
## per-cell copy numbers agree best with amp_cn (walks with identical carriers are told apart by cn);
## a walk takes at most one label. Returns gw_id, gene_label ("EGFR 1"), amp_id4, cn_filter = TRUE.
walk_labels_from_cell_summary <- function(summary, counts, min_cn = 1) {
    key <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))
    s <- data.table::data.table(cell = key(summary$pair), amp_id4 = as.character(summary$amp_id4), amp_cn = as.numeric(summary$amp_cn))
    s <- s[!is.na(amp_id4) & nzchar(amp_id4)]
    car <- counts[is.finite(cn) & cn >= min_cn, .(cell = key(pair), gw_id, cn)]
    cand <- data.table::rbindlist(lapply(unique(s$amp_id4), function(l) {
        sl <- s[amp_id4 == l]
        data.table::rbindlist(lapply(unique(car$gw_id), function(w) {
            cw <- car[gw_id == w]
            shared <- merge(sl, cw, by = "cell")
            data.table::data.table(amp_id4 = l, gw_id = w,
                                   jaccard = nrow(shared) / length(union(sl$cell, cw$cell)),
                                   cn_diff = if (nrow(shared)) mean(abs(shared$amp_cn - shared$cn)) else Inf)
        }))
    }))
    cand <- cand[jaccard >= 0.8][order(-round(jaccard, 2), cn_diff)]
    picked <- cand[0]
    for (i in seq_len(nrow(cand)))
        if (!cand$amp_id4[i] %in% picked$amp_id4 && !cand$gw_id[i] %in% picked$gw_id) picked <- rbind(picked, cand[i])
    message(sprintf("sc_export_walks: %d of %d per-cell amplicon labels matched to walks", nrow(picked), data.table::uniqueN(s$amp_id4)))
    out <- picked[, .(gw_id, amp_id4, cn_filter = TRUE,
                      gene_label = trimws(sub("\\s*\\(\\d+\\)$", "", sub("^\\S+\\s+", "", amp_id4))))]
    out
}

## TRUE where node b continues node a along the reference (same chromosome and strand, abutting)
## TRUE where node b continues node a along the reference (same chromosome and
## strand) after skipping at most `gap` bases
ref_adjacent <- function(a, b, gap = 0) {
    d <- ifelse(a$strand == "-", a$start - b$end - 1, b$start - a$end - 1)
    a$chromosome == b$chromosome & a$strand == b$strand & d >= 0 & d <= gap
}

## Merge runs of reference-adjacent nodes (within `gap`) into single intervals
## (width-weighted mean cn, the run keeps the `via` of its first node); a
## circular walk whose last node abuts its first is rotated so the run across
## the closing junction is merged too.
collapse_ref_runs <- function(nodes, circular = FALSE, gap = 0) {
    n <- nrow(nodes)
    if (n < 2) return(nodes)
    brk <- !ref_adjacent(nodes[-n], nodes[-1], gap)
    if (circular && any(brk) && isTRUE(ref_adjacent(nodes[n], nodes[1], gap))) {
        first <- which(brk)[1] + 1
        nodes <- nodes[c(first:n, seq_len(first - 1))]
        brk <- !ref_adjacent(nodes[-n], nodes[-1], gap)
    }
    nodes[, run := cumsum(c(TRUE, brk))]
    out <- nodes[, list(chromosome = chromosome[1], start = min(start), end = max(end), strand = strand[1],
                        cn = if (all(is.na(cn))) NA_real_ else stats::weighted.mean(cn, end - start + 1, na.rm = TRUE),
                        via = via[1]),
                 by = run]
    out[, run := NULL]
    out
}

## Simplify a walk for display, as gWalk$simplify() does in the blogs plus two
## steps for pieces too small to see at amplicon scale:
##  1. collapse reference-adjacent nodes, so only ALT junctions remain, and
##     same-strand neighbours within `merge_gap` (small deletions);
##  2. drop internal pieces shorter than `min_width` (templated insertions),
##     recording them in `via` of the junction that skips them;
##  3. merge again the neighbours that dropping brought together.
simplify_walk_nodes <- function(nodes, circular = FALSE, min_width = 1000, merge_gap = 1000) {
    nodes <- data.table::copy(nodes)
    nodes[, via := ""]
    nodes <- collapse_ref_runs(nodes, circular)
    if (merge_gap > 0) nodes <- collapse_ref_runs(nodes, circular, merge_gap)
    small <- (nodes$end - nodes$start + 1) < min_width
    if (min_width > 0 && any(small) && !all(small)) {
        desc <- sprintf("%s:%s-%s%s", nodes$chromosome, format(nodes$start, big.mark = ",", trim = TRUE),
                        format(nodes$end, big.mark = ",", trim = TRUE), nodes$strand)
        n <- nrow(nodes)
        pending <- character(0)
        ## walk around once (twice for circular, so pieces before the first kept node reach it)
        idx <- if (circular) c(seq_len(n), seq_len(n)) else seq_len(n)
        vias <- nodes$via
        seen <- logical(n)
        for (i in idx) {
            if (small[i]) { if (!seen[i]) pending <- c(pending, desc[i]); seen[i] <- TRUE; next }
            if (length(pending)) vias[i] <- paste(c(if (nzchar(vias[i])) vias[i], pending), collapse = "; ")
            pending <- character(0)
            seen[i] <- TRUE
            if (all(seen) && !length(pending)) break
        }
        nodes[, via := vias]
        nodes <- nodes[!small]
    }
    if (merge_gap > 0) nodes <- collapse_ref_runs(nodes, circular, merge_gap)
    nodes
}
