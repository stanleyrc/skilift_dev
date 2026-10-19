#' @name sc_neftel_quadrant
#' @title sc_neftel_quadrant
#' @description
#'
#' Neftel 2019 two-dimensional cell-state coordinates from the four module
#' scores the single-cell Seurat objects carry (MES1, AC2, OPC3, NPC4).
#' D = max(OPC, NPC) - max(AC, MES) is the vertical axis (positive = OPC/NPC-like
#' half, negative = AC/MES-like half). The horizontal axis is the signed
#' within-half difference: log2(|NPC - OPC| + 1) for D > 0, log2(|MES - AC| + 1)
#' for D <= 0, positive towards NPC (upper half) and MES (lower half), negative
#' towards OPC and AC (the arrangement of Neftel Fig. 2: OPC top left, NPC top
#' right, AC bottom left, MES bottom right). The discrete state is the
#' arg-max module, which reproduces the Seurat `state` column exactly.
#'
#' @param MES,AC,OPC,NPC numeric vectors of module scores (same length)
#' @return data.table with D, x, state
#' @export
#' @author Stanley Clarke
sc_neftel_quadrant <- function(MES, AC, OPC, NPC) {
    stopifnot(length(MES) == length(AC), length(AC) == length(OPC), length(OPC) == length(NPC))
    D <- pmax(OPC, NPC) - pmax(AC, MES)
    top <- !is.na(D) & D > 0
    x <- ifelse(top,
                sign(NPC - OPC) * log2(abs(NPC - OPC) + 1),
                sign(MES - AC) * log2(abs(MES - AC) + 1))
    sc <- cbind(MES = MES, AC = AC, OPC = OPC, NPC = NPC)
    ok <- stats::complete.cases(sc)
    state <- rep(NA_character_, length(D))
    state[ok] <- colnames(sc)[max.col(sc[ok, , drop = FALSE], ties.method = "first")]
    data.table::data.table(D = D, x = x, state = state)
}

#' @name sc_export_state_scores
#' @title sc_export_state_scores
#' @description
#'
#' Continuous per-cell Neftel module scores and quadrant coordinates for the
#' paper figure 1e: <out_dir>/_cohort/rna/state_scores.json
#' { format, columns, cells: [{ patient, rna_id, cell_id, MES, AC, OPC, NPC, D, x, state }] }.
#' The scores are the Seurat metadata columns the discrete `state` of
#' _cohort/rna/cells.json was derived from (state == arg-max score; checked
#' and reported). Cells without scores (non-malignant) are not written.
#'
#' @param seurat_paths named character vector: patient id -> Seurat .rds path
#' @param out_dir gOS data directory (holds the patient folders and _cohort/)
#' @param score_cols names of the MES, AC, OPC, NPC columns in the Seurat metadata
#' @param cells_json _cohort/rna/cells.json used to link rna_id to the DNA cell_id (default: in out_dir)
#' @return the exported list (invisibly)
#' @export
#' @author Stanley Clarke
sc_export_state_scores <- function(seurat_paths, out_dir, score_cols = c(MES = "MES1", AC = "AC2", OPC = "OPC3", NPC = "NPC4"),
                                  cells_json = file.path(out_dir, "_cohort", "rna", "cells.json")) {
    links <- NULL
    if (file.exists(cells_json)) {
        cj <- jsonlite::read_json(cells_json, simplifyVector = TRUE)$cells
        links <- data.table::as.data.table(cj)[, c("patient", "rna_id", "cell_id", "state")]
        data.table::setnames(links, "state", "state_cells_json")
    }
    parts <- lapply(names(seurat_paths), function(p) {
        so <- readRDS(seurat_paths[[p]])
        meta <- methods::slot(so, "meta.data")
        miss <- setdiff(unname(score_cols), colnames(meta))
        if (length(miss)) {
            message(p, ": score columns missing (", paste(miss, collapse = ", "), "): not available")
            return(NULL)
        }
        dt <- data.table::data.table(patient = p, rna_id = rownames(meta))
        for (k in names(score_cols)) dt[[k]] <- as.numeric(meta[[score_cols[[k]]]])
        dt <- dt[stats::complete.cases(dt[, names(score_cols), with = FALSE])]
        q <- sc_neftel_quadrant(dt$MES, dt$AC, dt$OPC, dt$NPC)
        dt <- cbind(dt, q)
        dt[, state_seurat := as.character(meta[["state"]])[match(rna_id, rownames(meta))]]
        dt
    })
    cells <- data.table::rbindlist(parts, use.names = TRUE)
    if (!is.null(links)) {
        cells <- merge(cells, links, by = c("patient", "rna_id"), all.x = TRUE, sort = FALSE)
    } else {
        cells[, `:=`(cell_id = NA_character_, state_cells_json = NA_character_)]
    }
    agree <- cells[!is.na(state_seurat), mean(state == state_seurat)]
    message(sprintf("state_scores: %d cells, %d patients; arg-max score == Seurat state in %.4f of cells", nrow(cells), length(unique(cells$patient)), agree))
    out <- list(
        format = "gos-sc-state-scores/1",
        method = paste0("Seurat metadata module scores ", paste(sprintf("%s=%s", names(score_cols), score_cols), collapse = ", "),
                        "; state = arg-max score; D = max(OPC,NPC) - max(AC,MES); x = signed log2(|diff|+1) within half (+ NPC/MES, - OPC/AC)"),
        n_patients = length(unique(cells$patient)),
        state_equals_argmax_fraction = agree,
        columns = c("patient", "rna_id", "cell_id", "MES", "AC", "OPC", "NPC", "D", "x", "state"),
        cells = cells[, c("patient", "rna_id", "cell_id", "MES", "AC", "OPC", "NPC", "D", "x", "state")]
    )
    dir.create(file.path(out_dir, "_cohort", "rna"), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(out, file.path(out_dir, "_cohort", "rna", "state_scores.json"), auto_unbox = TRUE, digits = 6, na = "null", dataframe = "rows")
    invisible(out)
}
