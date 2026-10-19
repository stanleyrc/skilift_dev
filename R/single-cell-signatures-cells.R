#' Count SBS96 mutations on each cell's root-to-tip path
#'
#' Uses the exported mutation categories and MRCA anchors, including truncal,
#' subclonal and private sites. Read dropout does not remove ancestral sites;
#' ambiguous, unmapped and outside-tumor calls are excluded.
#' @param variants Variant records from snv_matrix.json.
#' @param tree An ape phylo with the same tip identifiers as the matrix.
#' @param cell_ids Cells to return; defaults to the tree tips.
#' @return An SBS96 by cell integer matrix.
#' @export
sc_cell_signature_counts <- function(variants, tree, cell_ids = tree$tip.label) {
    v <- data.table::as.data.table(variants)
    if (!all(c("id", "category", "anchor", "context") %in% names(v)))
        stop("variants require id, category, anchor and context")
    if (anyDuplicated(v$id)) stop("duplicate mutation identifiers")
    counts <- matrix(0L, nrow = length(SC_SBS96), ncol = length(cell_ids),
                     dimnames = list(SC_SBS96, cell_ids))
    keep <- which(v$category %in% c("truncal", "subclonal", "private") &
                  v$context %in% SC_SBS96 & !is.na(v$anchor) & nzchar(v$anchor))
    if (!length(keep)) stop("no mapped SBS96 mutations with root-to-tip anchors")
    parent <- stats::setNames(tree$edge[, 1], tree$edge[, 2])
    paths <- lapply(match(cell_ids, tree$tip.label), function(tip) {
        if (is.na(tip)) return(integer())
        path <- tip
        while (as.character(tip) %in% names(parent)) {
            tip <- unname(parent[as.character(tip)])
            path <- c(path, tip)
        }
        path
    })
    by_anchor <- split(keep, v$anchor[keep])
    for (anchor in names(by_anchor)) {
        tips <- strsplit(anchor, "|", fixed = TRUE)[[1]]
        if (!length(tips) || any(!tips %in% tree$tip.label)) next
        node <- if (length(tips) == 1) match(tips, tree$tip.label) else ape::getMRCA(tree, tips)
        if (is.null(node)) next
        members <- which(vapply(paths, function(path) node %in% path, logical(1)))
        if (!length(members)) next
        values <- as.integer(table(factor(v$context[by_anchor[[anchor]]], levels = SC_SBS96)))
        counts[, members] <- counts[, members, drop = FALSE] + values
    }
    counts
}

#' Export per-cell SBS signature fits on root-to-tip mutation sets
#'
#' Calls sc_fit_sbs_signatures (SigProfilerAssignment COSMIC SBS) once for all
#' eligible cells. Cells below min_mutations retain n and too_few but have no
#' exposure estimates. Writes only signatures_cells.json; existing files are
#' backed up before an atomic replacement.
#' @param patient_dir Patient folder containing tree.nwk and snv_matrix.json.
#' @param out_dir Output folder; can be a staging directory.
#' @param work_dir Private SigProfiler working directory.
#' @param min_mutations Minimum valid SBS96 mutations; defaults to 100.
#' @param counts Optional precomputed SBS96-by-cell counts matrix.
#' @param fit_fun Fitter, defaults to sc_fit_sbs_signatures.
#' @param quality_fun Fit-quality reader, defaults to sc_signature_fit_quality.
#' @param cfg sc_annotation_config().
#' @return Invisibly, the exported object.
#' @export
sc_export_cell_signatures <- function(patient_dir, out_dir = patient_dir,
                                     work_dir = file.path(out_dir, "cell_signature_work"),
                                     min_mutations = 100L, counts = NULL,
                                     fit_fun = sc_fit_sbs_signatures,
                                     quality_fun = sc_signature_fit_quality,
                                     cfg = sc_annotation_config()) {
    if (min_mutations < 100L) stop("min_mutations must be at least 100")
    input_md5 <- NULL
    if (is.null(counts)) {
        source_files <- file.path(patient_dir, c("snv_matrix.json", "tree.nwk"))
        input_md5 <- tools::md5sum(source_files)
        snv <- jsonlite::read_json(file.path(patient_dir, "snv_matrix.json"), simplifyVector = TRUE)
        tree <- ape::read.tree(file.path(patient_dir, "tree.nwk"))
        cell_ids <- unique(c(names(snv$cells), tree$tip.label))
        counts <- sc_cell_signature_counts(snv$variants, tree, cell_ids)
        if (!identical(input_md5, tools::md5sum(source_files))) stop("signature inputs changed while reading")
        names(input_md5) <- basename(names(input_md5))
    }
    counts <- as.matrix(counts)
    if (!all(SC_SBS96 %in% rownames(counts)) || is.null(colnames(counts)) || anyDuplicated(colnames(counts)))
        stop("counts need SBS96 rownames and unique cell column names")
    counts <- counts[SC_SBS96, , drop = FALSE]
    if (any(!is.finite(counts) | counts < 0 | counts != floor(counts))) stop("counts must be nonnegative integers")
    n <- colSums(counts)
    eligible <- names(n)[n >= min_mutations]
    sets <- lapply(eligible, function(id) rep(SC_SBS96, counts[, id]))
    names(sets) <- eligible
    fits <- if (length(sets)) fit_fun(sets, work_dir, min_mutations = min_mutations, cfg = cfg) else NULL
    fits <- if (!is.null(fits)) data.table::as.data.table(fits) else data.table::data.table()
    quality <- if (length(sets)) quality_fun(work_dir) else NULL
    cells <- lapply(names(n), function(id) {
        f <- if (nrow(fits)) fits[set == id & activity > 0] else NULL
        activities <- if (!is.null(f) && nrow(f)) stats::setNames(f$activity, as.character(f$signature)) else numeric()
        exposures <- if (length(activities) && sum(activities) > 0) as.list(activities / sum(activities)) else list()
        too_few <- n[[id]] < min_mutations
        cosine <- quality$sets[[id]]$stats$cosine
        list(cell_id = id, n = unname(n[[id]]), n_mutations = unname(n[[id]]), too_few = too_few,
             fit_status = if (too_few) "too_few" else if (!length(exposures)) "no_activity" else "ok",
             exposures = exposures, cosine = if (is.null(cosine)) NA_real_ else cosine)
    })
    all_signatures <- sort(unique(unlist(lapply(cells, function(x) names(x$exposures)))))
    out <- list(format = "gos-sc-cell-signatures/1", genome = "GRCh38", cosmic_version = "3.4",
        method = "SigProfilerAssignment cosmic_fit; root-to-tip truncal + subclonal + private mapped SBS96 sites",
        min_mutations = min_mutations, mutation_set = "mapped root-to-tip; ambiguous/unmapped/outside-tumor excluded",
        input_md5 = as.list(input_md5),
        signatures = as.list(all_signatures), cells = cells)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    target <- file.path(out_dir, "signatures_cells.json")
    if (file.exists(target)) {
        backup <- paste0(target, ".bak_", format(Sys.Date(), "%y%m%d"))
        if (file.exists(backup)) backup <- paste0(backup, "_", format(Sys.time(), "%H%M%S"))
        if (!file.copy(target, backup)) stop("could not back up existing signatures_cells.json")
    }
    temporary <- tempfile(".signatures_cells_", tmpdir = out_dir)
    on.exit(unlink(temporary), add = TRUE)
    jsonlite::write_json(out, temporary, auto_unbox = TRUE, digits = NA, na = "null", null = "null")
    if (!file.rename(temporary, target)) stop("could not install staged signatures_cells.json")
    invisible(out)
}
