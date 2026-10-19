# Paper Fig 3A/B data: amplicon classes (per-cell amplified gene clusters of the JaBbA graphs) and graph-level
# "amplicon structures" (cells of a class clustered by the Jaccard distance of their ALT-junction sets), plus the
# Fig 4B ecDNA SNV set (pileup sites inside an amplicon with no alt reads in normal cells).
# Inputs are the per-patient tables of the amplicon copy-number database (gene clusters with mean CN and PASS filter,
# per-cell junction table `jj.dt`); no patient data lives here.

#' Graph states of one gene cluster: complete-linkage clusters of cells by the Jaccard distance of their junction sets
#'
#' Port of srctools `assign_graph_states`: cells sharing (nearly) the same junctions form one state; the tree is
#' cut at `max_difference`.
#' @param x data.table with columns pair (cell) and coord (junction id)
#' @param max_difference Jaccard distance cut (0.25 in the paper)
#' @return data.table(pair, graph_state)
#' @export
sc_graph_states <- function(x, max_difference = 0.25) {
  x <- unique(data.table::as.data.table(x)[!is.na(coord), .(pair, coord)])
  cells <- unique(x$pair)
  if (length(cells) == 0) return(data.table::data.table(pair = character(), graph_state = integer()))
  if (length(cells) == 1) return(data.table::data.table(pair = cells, graph_state = 1L))
  junctions <- unique(x$coord)
  m <- Matrix::sparseMatrix(i = match(x$pair, cells), j = match(x$coord, junctions), x = 1L,
                            dims = c(length(cells), length(junctions)), dimnames = list(cells, junctions))
  inter <- as.matrix(Matrix::tcrossprod(m))
  nj <- Matrix::rowSums(m)
  d <- 1 - inter / (outer(nj, nj, "+") - inter)
  diag(d) <- 0
  hc <- stats::hclust(stats::as.dist(d), method = "complete")
  data.table::data.table(pair = cells, graph_state = as.integer(stats::cutree(hc, h = max_difference)))
}

#' Graph states of every (patient, gene cluster) in a junction table
#'
#' Port of srctools `assign_graph_states_to_jj_dt`. With `use_junction_universe`, a cell is compared on the junctions
#' seen in any cell of the cluster (dropout-tolerant), keeping junctions with CN > 0, the "No_ALT_junctions"
#' placeholder, and every junction of cells whose cluster junctions all have CN 0.
#' @param jj data.table: patient, pair, coord, cn, genes_in_cluster, alt_junctions
#' @return data.table(patient, genes_in_cluster, pair, graph_state)
#' @export
sc_graph_states_jj <- function(jj, max_difference = 0.25, use_junction_universe = TRUE) {
  jj <- data.table::copy(data.table::as.data.table(jj))
  in_cluster <- function(d) d[!is.na(genes_in_cluster) & genes_in_cluster != "Not_in_subgraph"]
  if (!use_junction_universe) {
    sel <- in_cluster(jj[(cn > 0 | alt_junctions %in% "No_ALT_junctions")])
    return(sel[, sc_graph_states(.SD, max_difference), by = .(patient, genes_in_cluster)])
  }
  jj[is.na(cn), cn := 0]
  jj[, any_pos := any(cn > 0), by = .(pair, genes_in_cluster)]
  jj[, keep := (cn > 0) | (alt_junctions %in% "No_ALT_junctions") | !any_pos]
  kept <- jj[keep == TRUE]
  cc <- unique(in_cluster(kept)[, .(patient, pair, genes_in_cluster)])
  uni <- in_cluster(kept)[, .(universe = list(unique(coord))), by = .(patient, genes_in_cluster)]
  ex <- merge(kept[!is.na(coord), .(patient, pair, coord)], cc, by = c("patient", "pair"), allow.cartesian = TRUE)
  if (!nrow(uni)) return(data.table::data.table(patient = character(), genes_in_cluster = character(), pair = character(), graph_state = integer()))
  uni[, sc_graph_states(ex[patient == .BY$patient & genes_in_cluster == .BY$genes_in_cluster & coord %in% universe[[1]]], max_difference),
      by = .(patient, genes_in_cluster)]
}

#' Junction table prepared as in the paper's Fig 3A: junctions spanning >= min_span (or the no-ALT placeholder),
#' and cells of a cluster left without any such junction kept as one "No_ALT_junctions" state
#' @export
sc_amplicon_structure_input <- function(jj, patient, min_span = 1e5) {
  jj <- data.table::copy(data.table::as.data.table(jj))
  pg <- unique(jj[!is.na(genes_in_cluster) & genes_in_cluster != "Not_in_subgraph", .(pair, genes_in_cluster)])
  pg[, pgc := paste0(pair, "_", genes_in_cluster)]
  keep <- unique(jj[span >= min_span | alt_junctions %in% "No_ALT_junctions", coord])
  jj <- jj[coord %in% keep]
  jj[, pgc := paste0(pair, "_", genes_in_cluster)]
  add <- pg[!pgc %in% jj$pgc]
  if (nrow(add)) {
    add[, `:=`(patient = patient, coord = "No_ALT_junctions", type = "ALT")]
    jj <- rbind(jj, add, fill = TRUE)
  }
  jj[, patient := patient]
  jj[, pgc := NULL]
  jj[]
}

sc_cluster_class <- function(g) vapply(strsplit(as.character(g), ","), function(z) paste(sort(trimws(z)), collapse = "+"), "")

#' Amplicon classes and structures of one patient
#'
#' @param clusters per-cell gene clusters: pair, genes_in_cluster (comma-separated), mean_cn, filter
#' @param jj per-cell junction table (see sc_graph_states_jj)
#' @param min_cells classes need >= this many PASS cells
#' @return list(classes = data.table(cls, n_cells, n_structures), cells = data.table(cls, cell, cn, structure), dropped)
#' @export
sc_amplicon_structures <- function(clusters, jj, patient, min_cells = 5, max_difference = 0.25, min_span = 1e5) {
  gc <- data.table::as.data.table(clusters)
  gc <- gc[filter %in% "PASS", .(cell = pair, cls = sc_cluster_class(genes_in_cluster), cn = as.numeric(mean_cn))]
  n <- gc[, .(n_cells = data.table::uniqueN(cell)), by = cls]
  keep <- n[n_cells >= min_cells]
  dropped <- n[n_cells < min_cells]
  gc <- gc[cls %in% keep$cls]
  st <- if (is.null(jj) || !nrow(jj)) NULL else sc_graph_states_jj(sc_amplicon_structure_input(jj, patient, min_span), max_difference, TRUE)
  if (!is.null(st) && nrow(st)) {
    st[, cls := sc_cluster_class(genes_in_cluster)]
    ns <- st[, .(n_structures = data.table::uniqueN(graph_state)), by = cls]
    cells <- merge(gc, st[, .(cls, cell = pair, structure = graph_state)], by = c("cls", "cell"), all.x = TRUE)
  } else {
    ns <- data.table::data.table(cls = character(), n_structures = integer())
    cells <- gc[, structure := NA_integer_]
  }
  classes <- merge(keep, ns, by = "cls", all.x = TRUE)
  list(classes = classes[], cells = cells[order(cls, cell)], dropped = dropped[])
}

#' Write _cohort/amplicon_structures.json (gos-sc-amplicon-structures/1) for the gOS Paper figures 3A/B
#'
#' @param amp_db data.frame of the amplicon copy-number database: patient, gene_clusters_amp_thresh4 (rds path),
#'   gene_cluster_sets_thresh4_lst (rds path of a list with jj.dt)
#' @param out_dir gOS data root (the file goes to <out_dir>/_cohort/)
#' @param n_tumour optional named vector: tumour cells with a copy-number profile per patient (denominator of 3B)
#' @export
sc_export_amplicon_structures <- function(amp_db, out_dir, patients = NULL, n_tumour = NULL, min_cells = 5, max_difference = 0.25, min_span = 1e5) {
  db <- data.table::as.data.table(amp_db)
  if (!is.null(patients)) db <- db[patient %in% patients]
  ok <- function(f) length(f) == 1 && !is.na(f) && nzchar(f) && file.exists(f)
  out <- list()
  na <- list()
  for (P in unique(db$patient)) {
    r <- db[patient == P][1]
    if (!ok(r$gene_clusters_amp_thresh4)) {
      na[[P]] <- "no gene-cluster table"
      next
    }
    jj <- if (ok(r$gene_cluster_sets_thresh4_lst)) data.table::as.data.table(readRDS(r$gene_cluster_sets_thresh4_lst)$jj.dt) else NULL
    a <- sc_amplicon_structures(readRDS(r$gene_clusters_amp_thresh4), jj, P, min_cells, max_difference, min_span)
    classes <- lapply(seq_len(nrow(a$classes)), function(i) {
      k <- a$classes$cls[i]
      cc <- a$cells[cls == k]
      list(cls = k, genes = strsplit(k, "+", fixed = TRUE)[[1]], n_cells = a$classes$n_cells[i],
           n_structures = if (is.na(a$classes$n_structures[i])) NULL else a$classes$n_structures[i],
           cells = cc$cell, cn = round(cc$cn, 2), structure = cc$structure)
    })
    out[[P]] <- list(n_tumour_cells = if (!is.null(n_tumour) && P %in% names(n_tumour)) unname(n_tumour[[P]]) else NULL,
                     junction_table = !is.null(jj), classes = classes,
                     dropped = lapply(seq_len(nrow(a$dropped)), function(i) list(cls = a$dropped$cls[i], n_cells = a$dropped$n_cells[i])))
  }
  res <- list(format = "gos-sc-amplicon-structures/1",
              method = paste("classes: per-cell amplified gene clusters of the JaBbA graphs (nodes > 4x ploidy), PASS = cluster mean CN >= 10,",
                             sprintf(">= %d cells; structures: cells of a class clustered by Jaccard distance of ALT-junction sets", min_cells),
                             sprintf("(junction span >= %g bp, junction universe of the class, complete linkage cut at %g)", min_span, max_difference)),
              min_cells = min_cells, max_difference = max_difference, min_span = min_span,
              patients = out, not_available = na)
  dir.create(file.path(out_dir, "_cohort"), showWarnings = FALSE, recursive = TRUE)
  f <- file.path(out_dir, "_cohort", "amplicon_structures.json")
  jsonlite::write_json(res, f, auto_unbox = TRUE, null = "null", na = "null", digits = NA)
  invisible(res)
}

#' ecDNA SNV sites: pileup sites with alt reads in >= min_cells cells and no alt read in the normal cells at
#' >= min_normal_depth normal reads (germline excluded)
#'
#' @param pileup data.table: cell, site ("<pos>_<ref>_<alt>" or with chr prefix), ref_count, alt_count
#' @param normal data.table: pos, normal_alt, normal_depth
#' @return sites (ordered by position) passing the filter
#' @export
sc_ecdna_snv_sites <- function(pileup, normal, min_cells = 2, min_normal_depth = 5) {
  p <- data.table::as.data.table(pileup)
  car <- p[alt_count > 0, .(n = data.table::uniqueN(cell)), by = site][n >= min_cells, site]
  pos <- as.integer(sub("_.*", "", car))
  nm <- data.table::as.data.table(normal)
  ok <- pos %in% nm[normal_alt == 0 & normal_depth >= min_normal_depth, pos]
  s <- car[ok]
  s[order(as.integer(sub("_.*", "", s)))]
}

#' Write <out_dir>/<patient>/ecdna_snvs.json (gos-sc-ecdna-snvs/1): the filtered sites and per-cell ref / alt reads
#' @export
sc_export_ecdna_snvs <- function(patient, pileup, normal, out_dir, chromosome, region = NULL, min_cells = 2, min_normal_depth = 5, min_reads = 2, source = "") {
  p <- data.table::as.data.table(pileup)
  s <- sc_ecdna_snv_sites(p, normal, min_cells, min_normal_depth)
  d <- p[site %in% s & ref_count + alt_count >= min_reads]
  nm <- data.table::as.data.table(normal)
  sites <- lapply(s, function(z) {
    parts <- strsplit(z, "_", fixed = TRUE)[[1]]
    x <- d[site == z]
    nn <- nm[pos == as.integer(parts[1])][1]
    list(id = paste(chromosome, z, sep = "_"), chromosome = chromosome, pos = as.integer(parts[1]), ref = parts[2], alt = parts[3],
         normal_alt = nn$normal_alt, normal_depth = nn$normal_depth,
         cells = x$cell, ref_count = as.integer(x$ref_count), alt_count = as.integer(x$alt_count))
  })
  res <- list(format = "gos-sc-ecdna-snvs/1", patient = patient, chromosome = chromosome, region = region,
              filter = sprintf("alt reads in >= %d cells; 0 alt reads in normal cells at >= %d normal reads; cells with >= %d reads", min_cells, min_normal_depth, min_reads),
              source = source, n_sites = length(s), sites = sites)
  dir.create(file.path(out_dir, patient), showWarnings = FALSE, recursive = TRUE)
  jsonlite::write_json(res, file.path(out_dir, patient, "ecdna_snvs.json"), auto_unbox = TRUE, null = "null", na = "null", digits = NA)
  invisible(res)
}
