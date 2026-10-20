## Joint junction balance copy number as an alternative per-cell CN source for gOS (2026-10-20).
## Writes data/<patient>/sc_cn_joint/{manifest.json, cells/<cell>/complex.json}. Each complex.json has the regular
## per-cell schema (settings / intervals / connections): the joint graph on the joint chromosomes, and the cell's regular
## JaBbA intervals elsewhere (every interval carries cn_source = "joint" | "regular"). Regular junctions joining a joint
## chromosome to another chromosome are dropped (their joint-side breakpoints are not on the joint segmentation) and
## counted in the manifest. The joint graphs come from the srctools joint-balance runner (feature/joint-balance, v3).

#' merge a joint genome JSON (joint chromosomes) with a regular genome JSON (all other chromosomes)
#' @param joint,regular lists parsed with jsonlite::fromJSON(simplifyVector = FALSE) (settings, intervals, connections)
#' @param joint_chrs chromosomes taken from the joint JSON
#' @return merged list; attr(, "dropped") = regular connections dropped because they touch a joint chromosome
sc_merge_joint_genome_json <- function(joint, regular, joint_chrs, description = "<h3>Joint junction balance</h3>") {
    joint_chrs <- as.character(joint_chrs)
    chr_of <- function(iv) sub("^chr", "", as.character(iv$chromosome))
    ji <- Filter(function(iv) chr_of(iv) %in% joint_chrs, joint$intervals %||% list())
    ri <- Filter(function(iv) !(chr_of(iv) %in% joint_chrs), regular$intervals %||% list())
    jmap <- stats::setNames(seq_along(ji), vapply(ji, function(iv) as.character(abs(iv$iid)), ""))
    rmap <- stats::setNames(length(ji) + seq_along(ri), vapply(ri, function(iv) as.character(abs(iv$iid)), ""))
    renum <- function(ivs, map, src) lapply(ivs, function(iv) {
        iv$iid <- unname(map[[as.character(abs(iv$iid))]]); iv$title <- as.character(iv$iid); iv$cn_source <- src; iv })
    intervals <- c(renum(ji, jmap, "joint"), renum(ri, rmap, "regular"))
    remap_end <- function(x, map) {
        if (is.null(x)) return(list(ok = TRUE, v = NULL))
        k <- as.character(abs(x))
        if (!(k %in% names(map))) return(list(ok = FALSE))
        list(ok = TRUE, v = sign(x) * unname(map[[k]]))
    }
    conns <- list(); dropped <- 0L
    take <- function(cs, map, src) {
        for (cc in cs %||% list()) {
            s <- remap_end(cc$source, map); k <- remap_end(cc$sink, map)
            if (!s$ok || !k$ok || (is.null(cc$source) && is.null(cc$sink))) { if (src == "regular" && xor(s$ok && !is.null(cc$source), k$ok && !is.null(cc$sink))) dropped <<- dropped + 1L; next }
            cc$source <- s$v; cc$sink <- k$v; cc$cn_source <- src
            conns[[length(conns) + 1L]] <<- cc
        }
    }
    take(joint$connections, jmap, "joint")
    take(regular$connections, rmap, "regular")
    conns <- lapply(seq_along(conns), function(i) { conns[[i]]$cid <- i; conns[[i]] })
    settings <- regular$settings %||% joint$settings %||% list()
    settings$description <- description
    out <- list(settings = settings, intervals = intervals, connections = conns)
    attr(out, "dropped") <- dropped
    out
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' export joint-balance per-cell CN for gOS
#' @param cells cell ids; joint_gg_files per-cell rds with $gg (gGraph, balanced) or a gGraph; regular_json_files per-cell
#'   regular complex.json; out_dir data/<patient>/sc_cn_joint; joint_chrs chromosomes solved jointly
#' @param verdict list(much_better, default_joint, criteria = list(list(key, label, pass, detail)), summary)
#' @param metrics,provenance free lists written into the manifest
#' @return manifest list (invisibly)
sc_export_joint_cn <- function(patient, cells, joint_gg_files, regular_json_files, out_dir, joint_chrs, verdict,
                               method = "joint junction balance", metrics = list(), provenance = list(), max.cn = 100,
                               cells_missing = list()) {
    stopifnot(length(cells) == length(joint_gg_files), length(cells) == length(regular_json_files))
    dir.create(file.path(out_dir, "cells"), recursive = TRUE, showWarnings = FALSE)
    done <- character(); missing <- cells_missing; dropped <- 0L
    for (i in seq_along(cells)) {
        cell <- cells[i]
        res <- tryCatch({
            x <- readRDS(joint_gg_files[i]); gg <- if (inherits(x, "gGraph")) x else x$gg
            gg <- gg[as.character(seqnames) %in% joint_chrs]
            tmp <- tempfile(fileext = ".json")
            gGnome::refresh(gg)$json(filename = tmp, verbose = FALSE, maxcn = max.cn)
            joint <- jsonlite::fromJSON(tmp, simplifyVector = FALSE); unlink(tmp)
            regular <- jsonlite::fromJSON(regular_json_files[i], simplifyVector = FALSE)
            m <- sc_merge_joint_genome_json(joint, regular, joint_chrs)
            d <- file.path(out_dir, "cells", cell); dir.create(d, recursive = TRUE, showWarnings = FALSE)
            jsonlite::write_json(m[c("settings", "intervals", "connections")], file.path(d, "complex.json"), auto_unbox = TRUE, digits = NA, null = "null")
            attr(m, "dropped")
        }, error = function(e) e)
        if (inherits(res, "error")) missing[[length(missing) + 1L]] <- list(cell = cell, reason = paste("export failed:", conditionMessage(res)))
        else { done <- c(done, cell); dropped <- dropped + res }
    }
    manifest <- list(
        format = "gos-sc-cn-joint/1", patient = patient, generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), method = method,
        chromosomes = I(as.character(joint_chrs)), cells = I(done), cells_missing = missing,
        outside_chromosomes = "regular per-cell JaBbA intervals (cn_source = regular)",
        dropped_regular_connections = dropped,
        default_source = if (isTRUE(verdict$default_joint)) "joint" else "regular",
        verdict = verdict, metrics = metrics, provenance = provenance)
    jsonlite::write_json(manifest, file.path(out_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
    invisible(manifest)
}
