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
#' @param walks gWalk or path to an rds holding one
#' @param counts data.table or rds path with columns gw_id, pair, cn, amp
#' @param coords optional data.table / rds path (amp_coords)
#' @param summary optional data.table / rds path (amplicon_summary_dt)
#' @param out_dir patient folder of the gOS dataset
#' @param patient patient id written into the file
#' @param min_cn copy number from which a cell counts as carrying the walk
#' @return path of the written json (invisibly)
#' @export
sc_export_walks <- function(walks, counts, coords = NULL, summary = NULL, out_dir, patient = NULL, min_cn = 1) {
    load_rds <- function(x) if (is.character(x)) readRDS(path.expand(x)) else x
    gw <- load_rds(walks)
    if (!inherits(gw, "gWalk")) stop("walks must be a gWalk (got ", paste(class(gw), collapse = "/"), ")")
    counts <- data.table::as.data.table(load_rds(counts))
    counts[, gw_id := as.character(gw_id)]
    counts[, pair := as.character(pair)]
    coords <- if (!is.null(coords)) data.table::as.data.table(load_rds(coords)) else NULL
    if (!is.null(coords)) coords[, gw_id := as.character(gw_id)]
    summary <- if (!is.null(summary)) data.table::as.data.table(load_rds(summary)) else NULL
    if (!is.null(summary)) summary[, gw_id := as.character(gw_id)]
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
        n <- nrow(nodes)
        pairs <- if (n > 1) cbind(seq_len(n - 1), seq_len(n - 1) + 1) else matrix(integer(0), ncol = 2)
        if (isTRUE(dt$circular[k]) && n > 1) pairs <- rbind(pairs, c(n, 1))
        junctions <- lapply(seq_len(nrow(pairs)), function(p) {
            a <- nodes[pairs[p, 1]]; b <- nodes[pairs[p, 2]]
            ref <- a$chromosome == b$chromosome && a$strand == b$strand &&
                ((a$strand != "-" && b$start == a$end + 1) || (a$strand == "-" && a$start == b$end + 1))
            list(from = pairs[p, 1] - 1L, to = pairs[p, 2] - 1L, type = if (isTRUE(ref)) "REF" else "ALT")
        })
        gid <- dt$gw_id[k]
        cc <- counts[gw_id == gid]
        amp <- cc$amp[!is.na(cc$amp) & cc$amp != ""][1]
        carriers <- cc[is.finite(cn) & cn >= min_cn]
        co <- if (!is.null(coords)) coords[gw_id == gid][1] else NULL
        su <- if (!is.null(summary)) summary[gw_id == gid][1] else NULL
        flags <- if (!is.null(co)) names(co)[vapply(names(co), function(nm) is.logical(co[[nm]]) && isTRUE(co[[nm]]), logical(1))] else character(0)
        cgc <- if (!is.null(co) && !is.na(co$cgc_genes)) trimws(strsplit(co$cgc_genes, ",")[[1]]) else character(0)
        list(
            id = gid, walk_id = dt$walk.id[k], name = dt$name[k],
            label = if (!is.na(amp)) amp else dt$name[k],
            circular = isTRUE(dt$circular[k]),
            span = if ("wid" %in% names(dt)) dt$wid[k] else sum(nodes$end - nodes$start + 1),
            n_nodes = n,
            coordinates = if (!is.null(co)) co$coordinates else NULL,
            genes = union(cgc, intersect(flags, known_flags)),
            driver_genes = intersect(flags, known_flags),
            ncells = nrow(carriers),
            total_cn = sum(carriers$cn),
            median_cn = if (nrow(carriers)) stats::median(carriers$cn) else 0,
            max_cn = if (nrow(carriers)) max(carriers$cn) else 0,
            curated = if (!is.null(su) && "cn_filter" %in% names(su)) isTRUE(su$cn_filter) else NULL,
            gene_label = if (!is.null(su) && "gene_label" %in% names(su)) su$gene_label else NULL,
            amp_id = if (!is.null(su) && "amp_id4" %in% names(su) && !is.na(su$amp_id4)) su$amp_id4 else NULL,
            nodes = nodes,
            junctions = junctions,
            cells = as.list(stats::setNames(carriers$cn, carriers$pair)))
    })
    out <- list(format = "gos-sc-walks/1", patient = patient, n_walks = length(walk_json),
                n_cells = length(cells), cells = cells, walks = walk_json)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    path <- file.path(out_dir, "walks.json")
    jsonlite::write_json(out, path, auto_unbox = TRUE, digits = NA, null = "null", na = "null")
    invisible(path)
}
