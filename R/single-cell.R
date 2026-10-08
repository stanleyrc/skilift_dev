## Single-cell WGS ingestion for PGV and gOS
##
## Wrappers for lifting a single-cell patient (one "pair" per cell) into
##   (1) a PGV instance (Skilift object): per-cell coverage, total and allelic
##       genome graphs, and a patient-level phylogeny with SNV / junction
##       heatmaps; and
##   (2) a gOS single-cell dataset (see gOS_dev/SINGLE_CELL.md): per-cell case
##       folders, a patient folder with the tree and SNV matrix, Seurat RNA
##       export, and small per-cell BAM slices for IGV.
##
## Every function takes a `cells` data.table with one row per cell and a `pair`
## column; other inputs are named by column. None of it needs srctools.
##
## Example (how BWH70 was lifted):
##
##   pgv = Skilift$new(datafiles_json_path = ..., datadir = ..., settings = ...)
##   cells = all.dt[patient == "BWH70" & in_tree == TRUE]
##
##   ## PGV: coverage, total CN and allelic graphs (each calls pgv$add_plots)
##   lift_sc_coverage(pgv, cells, cov_col = "tumor_dryclean_cov_1kb",
##                    out_dir = rel2abs_dir, mask = egfr_mask, cores = 20)
##   lift_sc_genome_graphs(pgv, cells, gg_col = "jabba_wg_slack1e3", order = 10)
##   lift_sc_genome_graphs(pgv, cells, gg_col = "balanced_gg_rds", type = "allelic",
##                         order = 15, visible = FALSE)
##
##   ## PGV: phylogeny with SNV and junction heatmaps
##   jcn = sc_junction_matrix(cells, gg_col = "jabba_wg_slack1e3", cores = 20)
##   snvs = sc_snv_table(vaf = m$vaf, alt = m$alt_count, depth = m$depth)
##   lift_sc_phylogeny(pgv, "BWH70_phylogeny", tree = rooted_tree_file,
##                     normal_cells = normals, snvs = snvs, jcn = jcn)
##
##   ## gOS: dataset, RNA and BAM slices
##   build_gos_sc_dataset(cells[, .(pair, clone_id, region, ploidy)], "BWH70",
##                        tree = rooted_tree_file, out_dir = gos_dir,
##                        pgv_datadir = pgv$datadir, snvs = snvs)
##   export_gos_sc_rna(seurat_rds, file.path(gos_dir, "data", "BWH70"),
##                     datafiles = file.path(gos_dir, "datafiles.json"), patient = "BWH70")
##   lift_gos_sc_bams(cells, bam_col = "deduped_bam", gos_datadir = file.path(gos_dir, "data"),
##                    variants = unique(snvs$mutation), cores = 8)


`%||%` <- function(a, b) if (is.null(a)) b else a


## ------------------------------------------------------------------------
## Phylogeny helpers
## ------------------------------------------------------------------------

#' @name sc_tree_order
#' @title sc_tree_order
#' @description
#'
#' Tip labels of a phylogeny in plotted top-to-bottom order, matching the
#' default (ladderized) ggtree layout read from the highest y down, without
#' needing ggtree.
#'
#' @param tree path to a Newick file or an ape phylo object
#' @param convert_dash replace "-" with "_" in tip labels
#' @return character vector of tip labels
#' @export
#' @author Stanley Clarke
sc_tree_order <- function(tree, convert_dash = FALSE) {
    if (is.character(tree)) tree <- ape::read.tree(tree)
    tree <- ape::reorder.phylo(ape::ladderize(tree, right = FALSE), "cladewise")
    n_tips <- length(tree$tip.label)
    tips <- tree$edge[tree$edge[, 2] <= n_tips, 2]
    tip_order <- rev(tree$tip.label[tips])
    if (convert_dash) tip_order <- gsub("-", "_", tip_order, fixed = TRUE)
    return(tip_order)
}

#' @name sc_root_tree
#' @title sc_root_tree
#' @description
#'
#' Root a single-cell phylogeny (e.g. an unrooted CellPhy bestTree) on its
#' normal cells. When the normals are not monophyletic the tree is rooted on
#' the branch above the tumor cells' MRCA instead.
#'
#' @param tree path to a Newick file or an ape phylo object
#' @param outgroup tip labels of the normal cells
#' @param output_file optional Newick path to write
#' @return rooted phylo object
#' @export
#' @author Stanley Clarke
sc_root_tree <- function(tree, outgroup, output_file = NULL) {
    if (is.character(tree)) tree <- ape::read.tree(tree)
    outgroup <- intersect(outgroup, tree$tip.label)
    if (!length(outgroup)) stop("None of the outgroup cells are in the tree")
    if (ape::is.monophyletic(tree, outgroup, reroot = TRUE)) {
        rooted <- ape::root(tree, outgroup = outgroup, resolve.root = TRUE)
    } else {
        warning("Normal cells are not monophyletic; rooting above the tumor MRCA")
        tumor <- setdiff(tree$tip.label, outgroup)
        rooted <- ape::root(tree, outgroup = outgroup[1], resolve.root = TRUE)
        rooted <- ape::root(rooted, node = ape::getMRCA(rooted, tumor), resolve.root = TRUE)
    }
    if (!is.null(output_file)) {
        dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
        ape::write.tree(rooted, output_file)
    }
    rooted
}

#' @name sc_simplify_phylogeny
#' @title sc_simplify_phylogeny
#' @description
#'
#' Prepare a single-cell phylogeny for PGV: rotate nodes so tips render in
#' tree order (PGV draws the Newick vertically flipped relative to ggtree),
#' drop internal support labels, shrink the branch leading to the tumor MRCA
#' and optionally drop the normal cells.
#'
#' @param tree path to a Newick file or an ape phylo object
#' @param normal_cells tip labels (after convert_dash) of normal cells
#' @param tree_order desired top-to-bottom tip order; default sc_tree_order(tree)
#' @param shorten_tumor_branch scale the stem branch of the tumor clade
#' @param scale_tumor_branch factor applied to the tumor stem branch
#' @param convert_dash replace "-" with "_" in tip labels
#' @param remove_normal_cells drop normal_cells from the displayed tree
#' @param flip_vertical reverse tree_order to match PGV's orientation
#' @param output_file optional Newick path to write
#' @return the simplified phylo object (invisibly)
#' @export
#' @author Stanley Clarke
sc_simplify_phylogeny <- function(
    tree,
    normal_cells = character(0),
    tree_order = NULL,
    shorten_tumor_branch = TRUE,
    scale_tumor_branch = 0.1,
    convert_dash = TRUE,
    remove_normal_cells = TRUE,
    flip_vertical = TRUE,
    output_file = NULL) {
    if (is.character(tree)) tree <- ape::read.tree(tree)
    if (convert_dash) {
        tree$tip.label <- gsub("-", "_", tree$tip.label, fixed = TRUE)
    }
    if (is.null(tree_order)) tree_order <- sc_tree_order(tree)
    if (!setequal(tree_order, tree$tip.label) || length(tree_order) != length(tree$tip.label)) {
        stop("tree_order must contain exactly the tree's tip labels")
    }
    target_order <- if (flip_vertical) rev(tree_order) else tree_order
    tree <- ape::rotateConstr(tree, target_order)
    tree$node.label <- NULL

    tumor_tips <- which(!tree$tip.label %in% normal_cells)
    if (shorten_tumor_branch && length(normal_cells) && length(tumor_tips) > 1) {
        tumor_mrca <- ape::getMRCA(tree, tumor_tips)
        stem <- which(tree$edge[, 2] == tumor_mrca)
        tree$edge.length[stem] <- tree$edge.length[stem] * scale_tumor_branch
    }
    if (remove_normal_cells && length(normal_cells)) {
        tree <- ape::drop.tip(tree, intersect(normal_cells, tree$tip.label))
    }
    if (!is.null(output_file)) {
        dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
        ape::write.tree(tree, file = output_file)
        message("Phylogeny written to: ", output_file)
    }
    invisible(tree)
}


## ------------------------------------------------------------------------
## SNV and junction matrices
## ------------------------------------------------------------------------

#' @name sc_junction_matrix
#' @title sc_junction_matrix
#' @description
#'
#' Cells x junctions copy-number matrix from each cell's JaBbA gGraph. ALT
#' junctions are merged across cells with gGnome::ra.merge; a cell lacking a
#' junction gets 0. Column names are junction coordinates without "chr"
#' (e.g. "10:132180758-132180758+ <-> 11:50354829-50354829+").
#'
#' @param cells data.table with a pair column
#' @param gg_col column of gGraph rds paths (or gGraph objects)
#' @param cores number of cores
#' @param pad padding (bp) used when matching junctions across cells
#' @param min_cn minimum junction copy number kept per cell
#' @param seqnames_keep chromosomes to keep
#' @return numeric matrix, rows = cells, columns = junctions
#' @export
#' @author Stanley Clarke
sc_junction_matrix <- function(
    cells,
    gg_col,
    cores = 1,
    pad = 0,
    min_cn = 1,
    seqnames_keep = c(1:22, "X", "Y")) {
    cells <- data.table::as.data.table(cells)
    seqnames_keep <- c(seqnames_keep, paste0("chr", seqnames_keep))
    junc.lst <- parallel::mclapply(seq_len(nrow(cells)), function(i) {
        gg <- cells[[gg_col]][[i]]
        gg <- if (is.character(gg)) readRDS(gg) else gg$copy()
        gg <- gg[seqnames %in% seqnames_keep]
        juncs <- gg$junctions[type == "ALT"]
        if (!length(juncs)) return(NULL)
        grl <- juncs$grl
        grl <- grl[mcols(grl)$cn >= min_cn]
        if (!length(grl)) return(NULL)
        mcols(grl) <- mcols(grl)[, "cn", drop = FALSE]
        gUtils::gr.chr(grl)
    }, mc.cores = cores)
    names(junc.lst) <- cells$pair
    junc.lst <- junc.lst[!vapply(junc.lst, is.null, logical(1))]
    if (!length(junc.lst)) stop("No ALT junctions found in any cell")

    message("Merging junctions across ", length(junc.lst), " cells...")
    merged.jj <- gGnome::jJ(do.call(gGnome::ra.merge, c(junc.lst, pad = pad)))
    merged.dt <- merged.jj$dt
    merged.dt[, junc_coord := merged.jj$junc]
    cn_cols <- paste0("cn.", names(junc.lst))
    sub.dt <- unique(merged.dt[, c("junc_coord", cn_cols), with = FALSE])
    jcn <- as.matrix(sub.dt[, cn_cols, with = FALSE])
    jcn[is.na(jcn)] <- 0
    dimnames(jcn) <- list(gsub("chr", "", sub.dt$junc_coord), names(junc.lst))
    jcn <- jcn[!grepl("Un", rownames(jcn)), , drop = FALSE]

    ## cells with no junctions get a row of zeros
    out <- matrix(0, nrow = nrow(cells), ncol = nrow(jcn),
                  dimnames = list(cells$pair, rownames(jcn)))
    out[colnames(jcn), ] <- t(jcn)
    return(out)
}

#' @name sc_snv_table
#' @title sc_snv_table
#' @description
#'
#' Long SNV table (pair, mutation, vaf, ref.count.t, alt.count.t) from cells x
#' mutations matrices, e.g. the vaf / alt_count / depth matrices of a mapped
#' cellphy result. This is the snvs input of write_sc_phylogeny_json(),
#' lift_sc_phylogeny() and build_gos_sc_dataset().
#'
#' @param vaf cells x mutations matrix of VAFs
#' @param alt cells x mutations matrix of alt read counts
#' @param depth cells x mutations matrix of total depth
#' @param cells optional cell order (row names) to keep
#' @param mutations optional mutation order (column names) to keep
#' @return data.table
#' @export
#' @author Stanley Clarke
sc_snv_table <- function(vaf, alt, depth, cells = NULL, mutations = NULL) {
    if (is.null(cells)) cells <- rownames(vaf)
    if (is.null(mutations)) mutations <- colnames(vaf)
    pick <- function(m) m[cells, mutations, drop = FALSE]
    vaf <- pick(vaf)
    alt <- pick(alt)
    ref <- pick(depth) - alt
    data.table::data.table(
        pair = rep(cells, each = length(mutations)),
        mutation = rep(mutations, times = length(cells)),
        vaf = as.numeric(t(vaf)),
        ref.count.t = as.numeric(t(ref)),
        alt.count.t = as.numeric(t(alt))
    )
}

#' @name read_sc_mutations_json
#' @title read_sc_mutations_json
#' @description
#'
#' Read a PGV phylogeny mutations.sparse.json back into the long SNV table
#' written by write_sc_phylogeny_json().
#'
#' @param path path to mutations.sparse.json
#' @return data.table with pair, mutation, vaf, ref.count.t, alt.count.t
#' @export
#' @author Stanley Clarke
read_sc_mutations_json <- function(path) {
    sparse <- jsonlite::fromJSON(path)
    data.table::rbindlist(lapply(names(sparse$cells), function(cell) {
        x <- sparse$cells[[cell]]
        if (!length(x)) return(NULL)
        data.table::data.table(pair = cell, mutation = x$variantId, vaf = x$vaf,
                               ref.count.t = x$refCount, alt.count.t = x$altCount)
    }))
}


## ------------------------------------------------------------------------
## PGV
## ------------------------------------------------------------------------

## add_plots() only keeps columns already in pgv$plots, so add any missing ones
add_sc_plots <- function(pgv, plots, cores) {
    missing <- setdiff(c("order", "defaultChartType"), names(pgv$plots))
    if (length(missing)) {
        plots_all <- data.table::as.data.table(data.table::copy(pgv$plots))
        for (col in missing) plots_all[, (col) := if (col == "order") NA_real_ else NA_character_]
        pgv$plots <- plots_all
    }
    pgv$add_plots(plots, cores = cores)
}

#' @name write_sc_phylogeny_json
#' @title write_sc_phylogeny_json
#' @description
#'
#' Write the SNV and junction heatmap files used by the PGV phylogeny plot:
#' mutations.sparse.json (per-cell VAF / ref / alt at each variant) and
#' junctions.json (cells x junctions copy number). Both are validated by
#' reading them back before being moved into place.
#'
#' @param snvs data.table/data.frame with pair, mutation, vaf, ref.count.t, alt.count.t
#'   (mutation ids like chr1_18258813_G_A); see sc_snv_table()
#' @param jcn numeric matrix, rows = cells, columns = junction ids; see sc_junction_matrix()
#' @param output_dir patient folder in the PGV data directory
#' @param overwrite replace existing files
#' @return named vector of the written paths (invisibly)
#' @export
#' @author Stanley Clarke
write_sc_phylogeny_json <- function(snvs, jcn, output_dir, overwrite = FALSE) {
    output_names <- c("mutations.sparse.json", "junctions.json")
    output_paths <- file.path(output_dir, output_names)
    existing <- output_paths[file.exists(output_paths)]
    if (length(existing) && !overwrite) {
        stop("Refusing to overwrite existing output (use overwrite = TRUE): ",
             paste(existing, collapse = ", "))
    }

    json <- function(x, ...) {
        as.character(jsonlite::toJSON(x, digits = NA, na = "null", null = "null",
                                      json_verbatim = TRUE, ...))
    }
    numbers <- function(x) {
        out <- sprintf("%.17g", x)
        out[is.na(x)] <- "null"
        out
    }
    check_numbers <- function(x, label, maximum = Inf) {
        invalid <- !is.numeric(x) || any(is.nan(x)) ||
            any(!is.na(x) & (!is.finite(x) | x < 0 | x > maximum))
        if (invalid) stop("Invalid numeric values in ", label)
    }
    check_ids <- function(x, label, cells = FALSE) {
        if (!is.character(x) || anyNA(x) || any(!nzchar(trimws(x))) ||
            any(grepl("[[:cntrl:]]", x)) || (cells && any(grepl("[<>]", x)))) {
            stop("Invalid IDs in ", label)
        }
    }

    snvs <- as.data.frame(snvs, stringsAsFactors = FALSE)
    required <- c("pair", "mutation", "vaf", "ref.count.t", "alt.count.t")
    if (!all(required %in% names(snvs))) {
        stop("snvs requires columns: ", paste(required, collapse = ", "))
    }
    for (col in c("vaf", "ref.count.t", "alt.count.t")) {
        ## a mixed character/numeric matrix arrives as character
        if (!is.numeric(snvs[[col]])) snvs[[col]] <- as.numeric(as.character(snvs[[col]]))
    }
    cell <- as.character(snvs$pair)
    site <- as.character(snvs$mutation)
    check_ids(cell, "SNV pair", cells = TRUE)
    check_ids(site, "SNV mutation")
    if (any(!grepl("^[A-Za-z0-9][A-Za-z0-9.-]*_[0-9]+_[A-Za-z*-]+_[A-Za-z*-]+$", site))) {
        stop("Malformed SNV mutation ID (expected chrom_pos_ref_alt)")
    }
    check_numbers(snvs$vaf, "SNV vaf", 1)
    check_numbers(snvs$ref.count.t, "SNV ref.count.t")
    check_numbers(snvs$alt.count.t, "SNV alt.count.t")
    if (anyDuplicated(paste(cell, site, sep = "\r"))) {
        stop("Duplicate SNV cell/site observation")
    }

    if (!is.matrix(jcn) || !is.numeric(jcn)) stop("jcn must be a numeric matrix")
    check_ids(rownames(jcn), "JCN row names", cells = TRUE)
    check_ids(colnames(jcn), "JCN column names")
    if (anyDuplicated(rownames(jcn)) || anyDuplicated(colnames(jcn))) {
        stop("JCN matrix needs unique row and column names")
    }
    check_numbers(jcn, "JCN")

    cell_ids <- unique(cell)
    variant_ids <- unique(site)
    rows_by_cell <- split(seq_len(nrow(snvs)), factor(cell, levels = cell_ids))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    staged <- vapply(output_names, function(name) {
        tempfile(paste0(".", name, "-"), tmpdir = output_dir)
    }, character(1))
    on.exit(unlink(staged), add = TRUE)

    ## mutations.sparse.json, streamed one cell at a time
    encoded_sites <- vapply(variant_ids, json, character(1), auto_unbox = TRUE)
    encoded_cells <- vapply(cell_ids, json, character(1), auto_unbox = TRUE)
    site_index <- match(site, variant_ids)
    con <- file(staged[[1]], open = "wt", encoding = "UTF-8")
    tryCatch({
        writeLines(paste0('{"schemaVersion":1,"variants":', json(data.frame(id = variant_ids)),
                          ',"countFields":["refCount","altCount"],"cells":{'), con)
        for (i in seq_along(cell_ids)) {
            rows <- rows_by_cell[[i]]
            records <- paste0('{"variantId":', encoded_sites[site_index[rows]],
                              ',"vaf":', numbers(snvs$vaf[rows]),
                              ',"refCount":', numbers(snvs$ref.count.t[rows]),
                              ',"altCount":', numbers(snvs$alt.count.t[rows]), "}")
            writeLines(paste0(if (i > 1) "," else "", encoded_cells[[i]], ":[",
                              paste(records, collapse = ","), "]"), con)
        }
        writeLines("}}", con)
    }, finally = close(con))

    ## junctions.json
    values <- lapply(seq_len(nrow(jcn)), function(r) {
        structure(paste0("[", paste(numbers(jcn[r, ]), collapse = ","), "]"), class = "json")
    })
    writeLines(json(list(schemaVersion = jsonlite::unbox(1L), cellIds = rownames(jcn),
                         junctions = data.frame(id = colnames(jcn)), values = values)),
               staged[[2]], useBytes = TRUE)

    ## round trip
    same <- function(a, b) isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 0))
    decoded <- jsonlite::fromJSON(staged[[1]])
    if (!identical(names(decoded$cells), cell_ids) || !identical(decoded$variants$id, variant_ids)) {
        stop("SNV order or labels changed on round trip")
    }
    for (i in seq_along(cell_ids)) {
        got <- decoded$cells[[i]]
        rows <- rows_by_cell[[i]]
        if (!identical(got$variantId, site[rows]) || !same(got$vaf, snvs$vaf[rows]) ||
            !same(got$refCount, snvs$ref.count.t[rows]) || !same(got$altCount, snvs$alt.count.t[rows])) {
            stop("SNV values changed on round trip for ", cell_ids[[i]])
        }
    }
    decoded_jcn <- jsonlite::fromJSON(staged[[2]])
    if (!identical(decoded_jcn$cellIds, rownames(jcn)) ||
        !identical(decoded_jcn$junctions$id, colnames(jcn)) ||
        !same(decoded_jcn$values, jcn)) {
        stop("JCN values changed on round trip")
    }

    for (i in seq_along(output_paths)) {
        if (!file.rename(staged[[i]], output_paths[[i]])) stop("Could not write ", output_paths[[i]])
    }
    message(sprintf("Wrote %s (%d cells, %d variants) and %s (%d cells, %d junctions)",
                    output_paths[[1]], length(cell_ids), length(variant_ids),
                    output_paths[[2]], nrow(jcn), ncol(jcn)))
    invisible(setNames(output_paths, c("mutations", "junctions")))
}

#' @name sc_rel2abs_coverage
#' @title sc_rel2abs_coverage
#' @description
#'
#' Convert a single-cell dryclean coverage to absolute copy number with
#' rel2abs, optionally removing bins that overlap a mask.
#'
#' @param cov dryclean coverage GRanges or rds path
#' @param ploidy cell ploidy
#' @param purity purity (1 for single cells)
#' @param mask optional GRanges (or rds path) of bins to drop
#' @param field coverage column to convert
#' @param new_col column for the converted values
#' @return GRanges
#' @export
#' @author Stanley Clarke
sc_rel2abs_coverage <- function(
    cov,
    ploidy,
    purity = 1,
    mask = NULL,
    field = "foreground",
    new_col = "foregroundabs") {
    if (!inherits(cov, "GRanges")) cov <- readRDS(cov)
    if (!is.null(mask)) {
        if (is.character(mask)) mask <- readRDS(mask)
        cov <- gUtils::gr.nochr(cov)
        mask <- gUtils::gr.nochr(mask)
        cov <- cov[!gUtils::`%^%`(cov, mask)]
    }
    mcols(cov)[[new_col]] <- rel2abs(gr = cov, purity = purity, ploidy = ploidy, field = field)
    return(cov)
}

#' @name lift_sc_coverage
#' @title lift_sc_coverage
#' @description
#'
#' Add rel2abs coverage scatterplots for each cell to a PGV instance. Each
#' cell's coverage is converted with sc_rel2abs_coverage() using its ploidy,
#' saved to out_dir (reused if already there) and added with pgv$add_plots().
#'
#' @param pgv Skilift object
#' @param cells data.table with pair, the coverage column and the ploidy column
#' @param cov_col column of dryclean coverage rds paths
#' @param out_dir directory for the per-cell rel2abs rds files
#' @param ploidy_col column with cell ploidy
#' @param purity purity (1 for single cells)
#' @param mask optional GRanges (or rds path) of bins to drop
#' @param ref reference for new patients in PGV
#' @param order plot order
#' @param visible whether the plot is shown by default
#' @param title plot titles; default "<pair> Rel2abs coverage (Ploidy: x)"
#' @param overwrite recompute existing rel2abs files and overwrite PGV plots
#' @param cores number of cores
#' @return the data.table of plots added (invisibly)
#' @export
#' @author Stanley Clarke
lift_sc_coverage <- function(
    pgv,
    cells,
    cov_col,
    out_dir,
    ploidy_col = "ploidy",
    purity = 1,
    mask = NULL,
    ref = "hg38",
    order = 20,
    visible = FALSE,
    title = NULL,
    overwrite = FALSE,
    cores = 1) {
    cells <- data.table::as.data.table(cells)
    cells <- cells[!is.na(get(cov_col)) & file.exists(get(cov_col))]
    if (!nrow(cells)) stop("No cells with an existing ", cov_col)
    if (!is.null(mask) && is.character(mask)) mask <- readRDS(mask)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    rel2abs_paths <- file.path(out_dir, paste0(cells$pair, "_rel2abs_cov.rds"))
    todo <- which(overwrite | !file.exists(rel2abs_paths))
    parallel::mclapply(todo, function(i) {
        cov <- sc_rel2abs_coverage(cells[[cov_col]][i], ploidy = cells[[ploidy_col]][i],
                                   purity = purity, mask = mask)
        saveRDS(cov, rel2abs_paths[i])
        NULL
    }, mc.cores = cores, mc.preschedule = FALSE)
    missing <- !file.exists(rel2abs_paths)
    if (any(missing)) {
        warning("rel2abs failed for: ", paste(cells$pair[missing], collapse = ", "))
        cells <- cells[!missing]
        rel2abs_paths <- rel2abs_paths[!missing]
    }
    if (is.null(title)) {
        title <- paste0(cells$pair, " Rel2abs coverage (Ploidy: ", cells[[ploidy_col]], ")")
    }
    plots <- bw_temp(patient_id = cells$pair, x = rel2abs_paths, ref = ref, order = order,
                     field = "foregroundabs", title = title, chart_type = "scatterplot",
                     type = "scatterplot", visible = visible, overwrite = overwrite)
    add_sc_plots(pgv, plots, cores)
    invisible(plots)
}

#' @name lift_sc_genome_graphs
#' @title lift_sc_genome_graphs
#' @description
#'
#' Add each cell's total copy-number (type = "genome") or allelic
#' (type = "allelic") JaBbA graph to a PGV instance with pgv$add_plots().
#' Allelic graphs are also linked to the cell's genome plot via allelicSource.
#'
#' @param pgv Skilift object
#' @param cells data.table with pair, the graph column and optionally the ploidy column
#' @param gg_col column of gGraph rds paths
#' @param type "genome" or "allelic"
#' @param ploidy_col column with cell ploidy, added to titles; NULL to omit
#' @param ref reference for new patients in PGV
#' @param order plot order
#' @param visible whether the plot is shown by default
#' @param title plot titles; default "<pair> Total CN|Allelic; Ploidy: x"
#' @param max.cn maximum copy number shown
#' @param link_allelic for type = "allelic", set allelicSource on each cell's genome plot
#' @param overwrite overwrite existing json files
#' @param cores number of cores
#' @return the data.table of plots added (invisibly)
#' @export
#' @author Stanley Clarke
lift_sc_genome_graphs <- function(
    pgv,
    cells,
    gg_col,
    type = c("genome", "allelic"),
    ploidy_col = "ploidy",
    ref = "hg38",
    order = if (type == "genome") 10 else 15,
    visible = type == "genome",
    title = NULL,
    max.cn = 1000,
    link_allelic = TRUE,
    overwrite = FALSE,
    cores = 1) {
    type <- match.arg(type)
    cells <- data.table::as.data.table(cells)
    cells <- cells[!is.na(get(gg_col)) & file.exists(get(gg_col))]
    if (!nrow(cells)) stop("No cells with an existing ", gg_col)
    source <- if (type == "genome") "genome.json" else "allelic.json"
    if (is.null(title)) {
        title <- paste0(cells$pair, if (type == "genome") " Total CN" else " Allelic")
        if (!is.null(ploidy_col)) title <- paste0(title, "; Ploidy: ", cells[[ploidy_col]])
    }
    plots <- genome_temp(patient_id = cells$pair, x = cells[[gg_col]], ref = ref, order = order,
                         type = type, source = source, visible = visible, title = title,
                         max.cn = max.cn, overwrite = overwrite)
    add_sc_plots(pgv, plots, cores)
    if (type == "allelic" && link_allelic) {
        plots_all <- data.table::copy(pgv$plots)
        plots_all[patient.id %in% cells$pair & type == "genome" & source == "genome.json",
                  allelicSource := "allelic.json"]
        pgv$plots <- plots_all
        pgv$update_datafiles_json()
    }
    invisible(plots)
}

#' @name lift_sc_phylogeny
#' @title lift_sc_phylogeny
#' @description
#'
#' Add a single-cell phylogeny patient to a PGV instance: writes
#' phylogeny.newick (via sc_simplify_phylogeny), and when snvs and jcn are
#' given the SNV / junction heatmap files (via write_sc_phylogeny_json), then
#' adds the patient and its phylogeny plot and updates datafiles.json.
#' Leaves in the tree should match the cell pair names used for the cells'
#' own plots.
#'
#' @param pgv Skilift object
#' @param patient_id patient.id for the phylogeny entry (e.g. "BWH70_phylogeny")
#' @param tree path to a rooted Newick file or an ape phylo object
#' @param normal_cells normal cell tip labels (dropped from the display tree)
#' @param snvs optional long SNV table; see sc_snv_table()
#' @param jcn optional cells x junctions matrix; see sc_junction_matrix()
#' @param ref reference
#' @param title plot title
#' @param order plot order
#' @param overwrite replace existing heatmap json files
#' @param ... passed to sc_simplify_phylogeny()
#' @return path to the patient folder (invisibly)
#' @export
#' @author Stanley Clarke
lift_sc_phylogeny <- function(
    pgv,
    patient_id,
    tree,
    normal_cells = character(0),
    snvs = NULL,
    jcn = NULL,
    ref = "hg38",
    title = paste(patient_id, "Phylogeny"),
    order = 1,
    overwrite = FALSE,
    ...) {
    patient_dir <- file.path(pgv$datadir, patient_id)
    dir.create(patient_dir, recursive = TRUE, showWarnings = FALSE)
    sc_simplify_phylogeny(tree, normal_cells = normal_cells,
                          output_file = file.path(patient_dir, "phylogeny.newick"), ...)
    with_heatmap <- !is.null(snvs) && !is.null(jcn)
    if (with_heatmap) {
        write_sc_phylogeny_json(snvs, jcn, patient_dir, overwrite = overwrite)
    } else if (!is.null(snvs) || !is.null(jcn)) {
        warning("Both snvs and jcn are needed for the heatmap; writing the tree only")
    }

    if (!patient_id %in% pgv$metadata$patient.id) {
        new_meta <- data.table::data.table(patient.id = patient_id, ref = ref,
                                           description = list(setNames(list(), character(0))))
        pgv$metadata <- rbind(pgv$metadata, new_meta, fill = TRUE)
    }
    plot <- data.table::data.table(patient.id = patient_id, type = "phylogeny", visible = TRUE,
                                   source = "phylogeny.newick", title = title, order = order)
    if (with_heatmap) {
        plot[, `:=`(heatmap.mutationSource = "mutations.sparse.json",
                    heatmap.mutationFormat = "sparse",
                    heatmap.junctionSource = "junctions.json")]
    }
    pgv$plots <- rbind(pgv$plots[!(patient.id == patient_id & type == "phylogeny")], plot, fill = TRUE)
    pgv$update_datafiles_json()
    invisible(patient_dir)
}


## ------------------------------------------------------------------------
## gOS single-cell datasets
## ------------------------------------------------------------------------

write_gos_json <- function(x, path) {
    jsonlite::write_json(x, path, auto_unbox = TRUE, digits = NA, null = "null", na = "null")
}

## row of a data.table -> named list without NA / empty values
gos_record <- function(row) {
    rec <- as.list(row)
    rec[vapply(rec, function(v) length(v) == 1 && !is.na(v) && nzchar(as.character(v)), logical(1))]
}

## gOS per-cell mutations.json (bulk mutations format) from SNV observations
gos_mutations_json <- function(obs) {
    settings <- list(y_axis = list(title = "copy number", visible = TRUE))
    if (!nrow(obs)) return(list(settings = settings, intervals = list()))   # e.g. normals with no called sites
    parts <- data.table::tstrsplit(obs$mutation, "_", fixed = TRUE)
    pos <- as.integer(parts[[2]])
    intervals <- data.table::data.table(
        chromosome = sub("^chr", "", parts[[1]]),
        startPoint = pos,
        endPoint = pos + 1L,
        iid = seq_len(nrow(obs)),
        title = seq_len(nrow(obs)),
        type = "interval",
        y = 2,
        annotation = sprintf("Type: SNV; Genomic_variant: %s>%s; VAF: %.3f; Alt_count: %s; Ref_count: %s; ",
                             parts[[3]], parts[[4]], obs$vaf, obs$alt.count.t, obs$ref.count.t)
    )
    list(settings = settings, intervals = intervals)
}

#' @name sc_node_anchors
#' @title sc_node_anchors
#' @description
#' Two tips whose most recent common ancestor is each node (one tip for a
#' leaf), so a node can be found again in any copy of the same topology,
#' whatever its node numbering.
#'
#' @param tree phylo object the node numbers refer to
#' @param nodes integer node numbers (NA allowed)
#' @return character vector "tipA|tipB", "tip" or NA
#' @export
#' @author Stanley Clarke
sc_node_anchors <- function(tree, nodes) {
    n <- ape::Ntip(tree)
    first_tip <- function(k) {
        while (k > n) k <- tree$edge[tree$edge[, 1] == k, 2][1]
        tree$tip.label[k]
    }
    anchor <- function(k) {
        if (is.na(k)) return(NA_character_)
        if (k <= n) return(tree$tip.label[k])
        kids <- tree$edge[tree$edge[, 1] == k, 2]
        paste(first_tip(kids[1]), first_tip(kids[length(kids)]), sep = "|")
    }
    lookup <- unique(stats::na.omit(nodes))
    anchors <- vapply(lookup, anchor, character(1))
    unname(anchors[match(nodes, lookup)])
}

#' @name build_gos_sc_dataset
#' @title build_gos_sc_dataset
#' @description
#'
#' Build (or update) a gOS single-cell dataset for one patient, following
#' gOS SINGLE_CELL.md. Each cell gets a case folder with metadata.json,
#' mutations.json (called SNVs) and, unless already written by
#' lift_gos_sc_cells(), complex.json / allelic.json / coverage.arrow linked (or
#' copied) from the cell's PGV folder. The patient
#' folder gets tree.nwk, metadata.json and snv_matrix.json (reads at every
#' site, called or not, for the SNV heatmap). The patient's entries in the
#' dataset's datafiles.json are replaced; other patients are kept.
#'
#' @param cells data.table with pair and any per-cell attributes to show in gOS
#'   (e.g. clone_id, region, ploidy, state, Phase); every non-NA column is copied
#'   into the cell's datafiles.json entry
#' @param patient_id patient id
#' @param out_dir gOS dataset folder (holding datafiles.json and data/)
#' @param tree optional rooted Newick path or phylo object; leaves = cell pairs.
#'   When given, cells are ordered by the tree.
#' @param pgv_datadir optional PGV data directory with per-cell genome.json, allelic.json,
#'   coverage.arrow to link in; leave NULL when the cell files were written with
#'   lift_gos_sc_cells()
#' @param snvs optional long SNV table (pair, mutation, vaf, ref.count.t, alt.count.t)
#' @param variant_info optional per-site table (mutation plus any fields, e.g. category,
#'   mapped_node, cellphy_input) written onto each variant of snv_matrix.json
#' @param attributes named list of fields added to every entry (e.g. tumor_type, disease)
#' @param link symlink PGV files (TRUE) or copy them (FALSE)
#' @return data.table of the cell records (invisibly)
#' @export
#' @author Stanley Clarke
build_gos_sc_dataset <- function(
    cells,
    patient_id,
    out_dir,
    tree = NULL,
    pgv_datadir = NULL,
    snvs = NULL,
    variant_info = NULL,
    attributes = list(tumor_type = "GBM", disease = "Glioblastoma", primary_site = "brain"),
    link = TRUE) {
    cells <- data.table::copy(data.table::as.data.table(cells))
    data_dir <- file.path(out_dir, "data")
    patient_dir <- file.path(data_dir, patient_id)
    dir.create(patient_dir, recursive = TRUE, showWarnings = FALSE)

    if (!is.null(tree)) {
        if (is.character(tree)) tree <- ape::read.tree(tree)
        ape::write.tree(tree, file.path(patient_dir, "tree.nwk"))
        leaves <- tree$tip.label
        missing <- setdiff(leaves, cells$pair)
        if (length(missing)) cells <- rbind(cells, data.table::data.table(pair = missing), fill = TRUE)
        cells <- cells[match(leaves, pair)]
    }
    if (!is.null(snvs)) {
        snvs <- data.table::copy(data.table::as.data.table(snvs))
        data.table::setkey(snvs, pair)
    }
    pgv_files <- c(complex.json = "genome.json", allelic.json = "allelic.json",
                   coverage.arrow = "coverage.arrow")

    records <- lapply(seq_len(nrow(cells)), function(i) {
        cell <- cells$pair[i]
        rec <- c(list(pair = cell, entry_type = "cell", patient_id = patient_id),
                 attributes, gos_record(cells[i, !"pair"]))
        cell_dir <- file.path(data_dir, cell)
        dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)

        if (!is.null(pgv_datadir)) {
            src <- file.path(pgv_datadir, cell, pgv_files)
            if (all(file.exists(src))) {
                dst <- file.path(cell_dir, names(pgv_files))
                unlink(dst)
                if (link) file.symlink(normalizePath(src), dst) else file.copy(src, dst)
            }
        }
        complex <- file.path(cell_dir, "complex.json")
        if (file.exists(complex)) {
            connections <- jsonlite::fromJSON(complex)$connections
            rec$junction_count <- if (length(connections)) sum(connections$type == "ALT") else 0L
        }
        if (!is.null(snvs) && cell %in% snvs$pair) {
            called <- snvs[.(cell)][alt.count.t > 0]
            rec$snv_count <- nrow(called)
            write_gos_json(gos_mutations_json(called), file.path(cell_dir, "mutations.json"))
        }
        rec$summary <- paste("Single cell", rec$clone_id, rec$region, sep = " · ")
        write_gos_json(list(c(rec, cov_slope = 1, cov_intercept = 0)),
                       file.path(cell_dir, "metadata.json"))
        rec
    })

    if (!is.null(snvs)) {
        ## compact matrix: cell -> [[variant index (0-based), ref, alt], ...]
        count <- function(x) ifelse(is.na(x), "null", format(x, scientific = FALSE, trim = TRUE))
        variant_ids <- unique(snvs$mutation)
        snvs[, variant := match(mutation, variant_ids) - 1L]
        cell_json <- snvs[, .(json = paste0(
            jsonlite::toJSON(pair[1], auto_unbox = TRUE), ":[",
            paste0("[", variant, ",", count(ref.count.t), ",", count(alt.count.t), "]", collapse = ","), "]")),
            by = pair]
        variants <- data.table::data.table(id = variant_ids)
        if (!is.null(variant_info)) {
            info <- data.table::as.data.table(variant_info)
            variants <- cbind(variants, info[match(variant_ids, info$mutation), !"mutation"])
            ## node numbers only mean something for this phylo object, so also
            ## write two tips whose MRCA is the node ("tipA|tipB"; one tip for a leaf)
            if (!is.null(tree) && "node" %in% names(variants)) {
                variants[, anchor := sc_node_anchors(tree, node)]
            }
        }
        writeLines(paste0('{"schemaVersion":1,"format":"compact","variants":',
                          jsonlite::toJSON(variants, na = "null", digits = NA),
                          ',"cells":{', paste(cell_json$json, collapse = ","), "}}"),
                   file.path(patient_dir, "snv_matrix.json"))
    }

    clones <- table(vapply(records, function(r) as.character(r$clone_id %||% "NA"), character(1)))
    patient <- c(list(pair = patient_id, entry_type = "patient", patient_id = patient_id),
                 attributes,
                 list(cell_count = length(records),
                      summary = paste0("Single-cell WGS patient\nCells: ", length(records),
                                       "\nClones: ", paste0(names(clones), " (", clones, " cells)", collapse = ", "),
                                       if (!is.null(tree)) "\nTree: included")))
    write_gos_json(list(patient), file.path(patient_dir, "metadata.json"))

    datafiles <- file.path(out_dir, "datafiles.json")
    entries <- if (file.exists(datafiles)) jsonlite::fromJSON(datafiles, simplifyVector = FALSE) else list()
    entries <- Filter(function(e) !identical(e$patient_id, patient_id) && !identical(e$pair, patient_id), entries)
    write_gos_json(c(entries, list(patient), records), datafiles)
    message(sprintf("gOS dataset %s: %d cells (%d with CN, %d with SNVs)", patient_id, length(records),
                    sum(vapply(records, function(r) !is.null(r$junction_count), logical(1))),
                    sum(vapply(records, function(r) !is.null(r$snv_count), logical(1)))))
    invisible(data.table::rbindlist(records, fill = TRUE))
}

#' @name lift_gos_sc_cells
#' @title lift_gos_sc_cells
#' @description
#'
#' Write each cell's gOS case files straight from its source objects:
#' complex.json (total CN JaBbA graph), allelic.json (allelic graph) and
#' coverage.arrow (dryclean coverage converted with sc_rel2abs_coverage(),
#' plotted as foregroundabs). Uses the same writers as pgv$add_plots().
#'
#' @param cells data.table with pair, the source columns and the ploidy column
#' @param gos_datadir gOS data folder (cell folders are created inside it)
#' @param total_gg_col column of total CN gGraph rds paths (NULL to skip)
#' @param allelic_gg_col column of allelic gGraph rds paths (NULL to skip)
#' @param cov_col column of dryclean coverage rds paths (NULL to skip)
#' @param ploidy_col column with cell ploidy, used for rel2abs
#' @param purity purity for rel2abs (1 for single cells)
#' @param mask optional GRanges (or rds path) of bins to drop before rel2abs
#' @param rel2abs_dir optional folder to keep the rel2abs coverage rds files
#' @param settings PGV-style settings.json with the reference seqlengths
#' @param ref reference name in settings
#' @param max.cn maximum copy number in the graph jsons
#' @param overwrite rewrite files that already exist
#' @param cores number of cells processed in parallel
#' @return data.table with pair and which files were written (invisibly)
#' @export
#' @author Stanley Clarke
lift_gos_sc_cells <- function(
    cells,
    gos_datadir,
    total_gg_col = NULL,
    allelic_gg_col = NULL,
    cov_col = NULL,
    ploidy_col = "ploidy",
    purity = 1,
    mask = NULL,
    rel2abs_dir = NULL,
    settings = internal_settings_path,
    ref = "hg38",
    max.cn = 1000,
    overwrite = FALSE,
    cores = 1) {
    cells <- data.table::as.data.table(cells)
    if (!is.null(mask) && is.character(mask)) mask <- readRDS(mask)
    if (!is.null(rel2abs_dir)) dir.create(rel2abs_dir, recursive = TRUE, showWarnings = FALSE)
    has <- function(col, i) !is.null(col) && !is.na(cells[[col]][i]) && file.exists(cells[[col]][i])
    done <- function(path) file.exists(path) && !overwrite
    attempt <- function(label, cell, expr) {
        tryCatch({ expr; TRUE }, error = function(e) {
            message(label, " failed for ", cell, ": ", conditionMessage(e))
            FALSE
        })
    }
    res <- parallel::mclapply(seq_len(nrow(cells)), function(i) {
        cell <- cells$pair[i]
        cell_dir <- file.path(gos_datadir, cell)
        dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)
        out <- list(pair = cell, complex = NA, allelic = NA, coverage = NA)
        if (has(total_gg_col, i)) {
            out$complex <- done(file.path(cell_dir, "complex.json")) || attempt("complex.json", cell,
                create_ggraph_json(genome_temp(patient_id = cell, x = cells[[total_gg_col]][i], ref = ref,
                                               source = "complex.json", max.cn = max.cn, overwrite = TRUE),
                                   gos_datadir, settings))
        }
        if (has(allelic_gg_col, i)) {
            out$allelic <- done(file.path(cell_dir, "allelic.json")) || attempt("allelic.json", cell,
                create_allelic_json(genome_temp(patient_id = cell, x = cells[[allelic_gg_col]][i], ref = ref,
                                                type = "allelic", source = "allelic.json", max.cn = max.cn,
                                                overwrite = TRUE),
                                    gos_datadir, settings))
        }
        if (has(cov_col, i)) {
            out$coverage <- done(file.path(cell_dir, "coverage.arrow")) || attempt("coverage.arrow", cell, {
                cov <- sc_rel2abs_coverage(cells[[cov_col]][i], ploidy = cells[[ploidy_col]][i],
                                           purity = purity, mask = mask)
                if (!is.null(rel2abs_dir)) saveRDS(cov, file.path(rel2abs_dir, paste0(cell, "_rel2abs_cov.rds")))
                create_scatterplot_arrow(arrow_temp(patient_id = cell, x = list(cov), ref = ref,
                                                    field = "foregroundabs", source = "coverage.arrow",
                                                    overwrite = TRUE),
                                         gos_datadir, settings)
            })
        }
        data.table::as.data.table(out)
    }, mc.cores = cores, mc.preschedule = FALSE)
    failed <- !vapply(res, data.table::is.data.table, logical(1))
    if (any(failed)) warning("Worker error for: ", paste(cells$pair[failed], collapse = ", "))
    res <- data.table::rbindlist(res[!failed])
    message(sprintf("gOS cell files: %d complex.json, %d allelic.json, %d coverage.arrow for %d cells",
                    sum(res$complex %in% TRUE), sum(res$allelic %in% TRUE), sum(res$coverage %in% TRUE), nrow(res)))
    invisible(res)
}

#' @name sc_seurat_metadata
#' @title sc_seurat_metadata
#' @description
#'
#' Per-cell metadata columns of a Seurat object as a data.table (pair = cell
#' barcode), e.g. to join onto the cells table given to build_gos_sc_dataset().
#'
#' @param seurat Seurat object or rds path
#' @param cols metadata columns to keep (missing ones are skipped)
#' @return data.table
#' @export
#' @author Stanley Clarke
sc_seurat_metadata <- function(
    seurat,
    cols = c("Clone_Annotation", "wEGFR_Species_Clone", "state", "Region",
             "Region_Annotation", "Phase", "Cell_Type", "seurat_clusters")) {
    if (is.character(seurat)) seurat <- readRDS(seurat)
    meta <- seurat[[]]
    cols <- intersect(cols, colnames(meta))
    out <- data.table::data.table(pair = rownames(meta))
    for (col in cols) out[[col]] <- as.character(meta[[col]])
    out
}

#' @name export_gos_sc_rna
#' @title export_gos_sc_rna
#' @description
#'
#' Export a patient's Seurat object to the gOS rna/ layout (format
#' gos-sc-rna/1, same as gOS services/sc-analysis/r/export_seurat.R): the
#' log-normalized data layer as a gene-major CSR matrix
#' (matrix.indptr.i32 / matrix.indices.i32 / matrix.data.f32), genes.tsv,
#' cells.json (cell metadata and UMAP) and manifest.json. RNA cells link to gOS
#' cells by pair, or by an rna_id on the cell's datafiles.json entry.
#'
#' @param seurat Seurat object or rds path
#' @param out_dir gOS patient folder (rna/ is written inside it)
#' @param datafiles optional gOS datafiles.json used to link RNA barcodes to cells
#' @param patient patient id used to filter datafiles.json entries
#' @param assay assay to export
#' @param reduction 2D embedding to include as umap_1 / umap_2
#' @param cluster_col cluster column
#' @param cell_type_col cell type column
#' @param meta_cols further metadata columns for cells.json
#' @param embeddings optional named list of extra 2D embeddings (matrices with
#'   barcodes as row names, e.g. from sc_recompute_umap()); each is written to
#'   cells.json as <name>_1 / <name>_2, NA for cells not in it
#' @return path to the rna/ folder (invisibly)
#' @export
#' @author Stanley Clarke
export_gos_sc_rna <- function(
    seurat,
    out_dir,
    datafiles = NULL,
    patient = NULL,
    assay = "RNA",
    reduction = "umap",
    cluster_col = "seurat_clusters",
    cell_type_col = "Cell_Type",
    meta_cols = character(0),
    embeddings = list()) {
    if (!requireNamespace("SeuratObject", quietly = TRUE) || !requireNamespace("Matrix", quietly = TRUE)) {
        stop("export_gos_sc_rna needs the SeuratObject and Matrix packages")
    }
    source_rds <- if (is.character(seurat)) normalizePath(seurat) else NA_character_
    if (is.character(seurat)) seurat <- readRDS(seurat)
    if (inherits(seurat[[assay]], "Assay5")) {
        seurat[[assay]] <- SeuratObject::JoinLayers(seurat[[assay]])
        expr <- SeuratObject::LayerData(seurat, assay = assay, layer = "data")
    } else {
        expr <- SeuratObject::GetAssayData(seurat, assay = assay, slot = "data")
    }
    expr <- methods::as(expr, "CsparseMatrix")                     # genes x cells
    by_gene <- methods::as(Matrix::t(expr), "CsparseMatrix")       # column-compressed by gene
    if (length(by_gene@x) > .Machine$integer.max) stop("too many non-zeros for int32 indptr")

    rna_dir <- file.path(out_dir, "rna")
    dir.create(rna_dir, recursive = TRUE, showWarnings = FALSE)
    write_bin <- function(x, name, int) {
        con <- file(file.path(rna_dir, name), "wb")
        on.exit(close(con))
        writeBin(if (int) as.integer(x) else as.numeric(x), con, size = 4, endian = "little")
    }
    write_bin(by_gene@p, "matrix.indptr.i32", TRUE)
    write_bin(by_gene@i, "matrix.indices.i32", TRUE)
    write_bin(by_gene@x, "matrix.data.f32", FALSE)
    writeLines(rownames(expr), file.path(rna_dir, "genes.tsv"))

    barcodes <- colnames(expr)
    cell_id <- rep(NA_character_, length(barcodes))
    if (!is.null(datafiles)) {
        for (r in jsonlite::fromJSON(datafiles, simplifyVector = FALSE)) {
            if (!identical(tolower(r$entry_type %||% ""), "cell")) next
            if (!is.null(patient) && !identical(as.character(r$patient_id), patient)) next
            id <- as.character(r$pair %||% r$case_id %||% r$id)
            cell_id[barcodes == as.character(r$rna_id %||% id)] <- id
        }
    }

    meta <- seurat[[]]
    cells <- data.frame(rna_id = barcodes, cell_id = cell_id, stringsAsFactors = FALSE)
    absent <- setdiff(meta_cols, colnames(meta))
    if (length(absent)) warning("metadata columns not found: ", paste(absent, collapse = ", "))
    for (col in unique(c(cluster_col, cell_type_col, "nCount_RNA", "nFeature_RNA", "percent.mt", meta_cols))) {
        if (!col %in% colnames(meta)) next
        v <- meta[barcodes, col]
        cells[[gsub("[^A-Za-z0-9_]", "_", col)]] <- if (is.factor(v)) as.character(v) else v
    }
    if (reduction %in% SeuratObject::Reductions(seurat)) {
        emb <- SeuratObject::Embeddings(seurat, reduction = reduction)[barcodes, 1:2, drop = FALSE]
        cells$umap_1 <- emb[, 1]
        cells$umap_2 <- emb[, 2]
    }
    for (name in names(embeddings)) {
        emb <- embeddings[[name]]
        k <- match(barcodes, rownames(emb))
        cells[[paste0(name, "_1")]] <- emb[k, 1]
        cells[[paste0(name, "_2")]] <- emb[k, 2]
    }
    jsonlite::write_json(list(cells = cells), file.path(rna_dir, "cells.json"),
                         auto_unbox = TRUE, na = "null", digits = NA, dataframe = "rows")
    manifest <- list(
        format = "gos-sc-rna/1",
        n_cells = ncol(expr),
        n_genes = nrow(expr),
        n_matched = sum(!is.na(cell_id)),
        normalization = paste0("Seurat ", assay, " data layer"),
        source_rds = source_rds,
        assay = assay,
        seurat_version = as.character(utils::packageVersion("SeuratObject")),
        created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    )
    jsonlite::write_json(manifest, file.path(rna_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
    message(sprintf("wrote %d cells x %d genes (%d linked to gOS cells) to %s",
                    ncol(expr), nrow(expr), sum(!is.na(cell_id)), rna_dir))
    invisible(rna_dir)
}

#' @name sc_recompute_umap
#' @title sc_recompute_umap
#' @description
#'
#' Recompute a UMAP on a subset of a Seurat object's cells (e.g. only the cells
#' that also have DNA), repeating the standard workflow on the subset:
#' variable features, scaling, PCA, then UMAP with the given settings (the
#' defaults match the GBM objects' own RunUMAP call).
#'
#' @param seurat Seurat object or rds path
#' @param cells barcodes to keep
#' @param dims PCs used for the UMAP
#' @param n.neighbors,min.dist,metric,seed.use passed to RunUMAP
#' @param nfeatures variable features
#' @return matrix (barcodes x 2)
#' @export
#' @author Stanley Clarke
sc_recompute_umap <- function(
    seurat,
    cells,
    dims = 1:15,
    n.neighbors = 30,
    min.dist = 0.3,
    metric = "cosine",
    seed.use = 42,
    nfeatures = 2000) {
    if (!requireNamespace("Seurat", quietly = TRUE)) stop("sc_recompute_umap needs the Seurat package")
    if (is.character(seurat)) seurat <- readRDS(seurat)
    cells <- intersect(cells, colnames(seurat))
    if (length(cells) < max(dims) + 2) stop("Too few cells (", length(cells), ") for a UMAP on ", max(dims), " PCs")
    sub <- subset(seurat, cells = cells)
    sub <- Seurat::FindVariableFeatures(sub, nfeatures = nfeatures, verbose = FALSE)
    sub <- Seurat::ScaleData(sub, verbose = FALSE)
    sub <- Seurat::RunPCA(sub, npcs = min(50, length(cells) - 1), verbose = FALSE)
    sub <- Seurat::RunUMAP(sub, dims = dims, n.neighbors = min(n.neighbors, length(cells) - 1),
                           min.dist = min.dist, metric = metric, seed.use = seed.use, verbose = FALSE)
    SeuratObject::Embeddings(sub, reduction = "umap")[, 1:2, drop = FALSE]
}

#' @name sc_igv_regions
#' @title sc_igv_regions
#' @description
#'
#' Regions to keep in a cell's IGV BAM slice: every SNV site plus the
#' breakpoints of the cell's ALT junctions (from its gOS complex.json /
#' PGV genome.json), each padded by pad bp.
#'
#' @param variants SNV ids like chr1_18258813_G_A
#' @param genome_json optional path to the cell's complex.json / genome.json
#' @param pad padding in bp
#' @return data.table with chrom, start, end (BED, sorted)
#' @export
#' @author Stanley Clarke
sc_igv_regions <- function(variants, genome_json = NULL, pad = 200, junction_pad = 1000) {
    parts <- data.table::tstrsplit(variants, "_", fixed = TRUE)
    pos <- as.integer(parts[[2]])
    regions <- data.table::data.table(chrom = parts[[1]], pos = pos, pad = pad)
    if (!is.null(genome_json) && file.exists(genome_json)) {
        g <- jsonlite::fromJSON(genome_json)
        alt <- g$connections
        if (length(alt)) alt <- alt[alt$type == "ALT" & !is.na(alt$source) & !is.na(alt$sink), ]
        if (length(alt) && nrow(alt)) {
            iv <- data.table::as.data.table(g$intervals)[, .(iid, chromosome, startPoint, endPoint)]
            ## a positive source leaves from the interval's end, a negative sink enters at its end
            ends <- rbind(data.table::data.table(iid = abs(alt$source), at_end = alt$source > 0),
                          data.table::data.table(iid = abs(alt$sink), at_end = alt$sink < 0))
            ends <- merge(ends, iv, by = "iid")
            ## junction breakpoints get a wider window so split / discordant reads are visible
            regions <- rbind(regions, ends[, .(
                chrom = paste0("chr", sub("^chr", "", chromosome)),
                pos = ifelse(at_end, endPoint, startPoint),
                pad = junction_pad)])
        }
    }
    regions[, `:=`(start = pmax(0L, as.integer(pos) - pad), end = as.integer(pos) + pad)]
    regions[order(chrom, start), .(chrom, start, end)]
}

#' @name lift_gos_sc_bams
#' @title lift_gos_sc_bams
#' @description
#'
#' Cut a small BAM per cell for IGV in gOS: reads overlapping sc_igv_regions()
#' (SNV sites plus the cell's junction breakpoints) are written to
#' <gos_datadir>/<pair>/reads.bam and indexed. Uses samtools.
#'
#' @param cells data.table with pair and the BAM column
#' @param bam_col column of per-cell BAM/CRAM paths
#' @param gos_datadir gOS data folder holding the cell folders
#' @param variants SNV ids like chr1_18258813_G_A
#' @param pad padding in bp around SNV sites
#' @param junction_pad padding in bp around junction (fusion / SV) breakpoints
#' @param bed_dir folder for the per-cell BED files (default: a temp dir)
#' @param reference reference fasta, needed when the inputs are CRAMs
#' @param samtools samtools binary
#' @param cores number of cells processed in parallel
#' @return data.table with pair, bam and ok (invisibly)
#' @export
#' @author Stanley Clarke
lift_gos_sc_bams <- function(
    cells,
    bam_col,
    gos_datadir,
    variants,
    pad = 200,
    junction_pad = 1000,
    bed_dir = tempfile("igv_regions_"),
    reference = NULL,
    samtools = "samtools",
    cores = 1) {
    cells <- data.table::as.data.table(cells)
    cells <- cells[!is.na(get(bam_col)) & file.exists(get(bam_col)) & dir.exists(file.path(gos_datadir, pair))]
    dir.create(bed_dir, recursive = TRUE, showWarnings = FALSE)
    ok <- parallel::mclapply(seq_len(nrow(cells)), function(i) {
        cell <- cells$pair[i]
        cell_dir <- file.path(gos_datadir, cell)
        bed <- file.path(bed_dir, paste0(cell, ".bed"))
        data.table::fwrite(sc_igv_regions(variants, file.path(cell_dir, "complex.json"), pad = pad, junction_pad = junction_pad),
                           bed, sep = "\t", col.names = FALSE)
        out <- file.path(cell_dir, "reads.bam")
        tmp <- paste0(out, ".tmp")
        ref_args <- if (!is.null(reference)) c("-T", shQuote(reference)) else NULL
        status <- system2(samtools, c("view", "-b", "-M", ref_args, "-L", shQuote(bed), "-o", shQuote(tmp),
                                      shQuote(cells[[bam_col]][i])))
        if (status != 0 || !file.rename(tmp, out)) return(FALSE)
        if (system2(samtools, c("index", shQuote(out))) != 0) return(FALSE)
        Sys.chmod(c(out, paste0(out, ".bai")), mode = "0664", use_umask = FALSE)
        TRUE
    }, mc.cores = cores, mc.preschedule = FALSE)
    result <- data.table::data.table(pair = cells$pair, bam = file.path(gos_datadir, cells$pair, "reads.bam"),
                                     ok = vapply(ok, isTRUE, logical(1)))
    if (any(!result$ok)) warning("BAM slicing failed for: ", paste(result[ok == FALSE]$pair, collapse = ", "))
    invisible(result)
}
