#' @name sc_export_rna_fusions
#' @title sc_export_rna_fusions
#' @description
#'
#' RNA fusions of one patient for gOS: the per-cell Arriba calls (STAR chimeric
#' alignment of each cell, fusions.tsv) with confidence in `keep_confidence`,
#' aggregated per fusion (gene pair + breakpoints) with the read support of every
#' cell, the gOS cell id of cells that also have DNA, and whether the patient's
#' DNA filtered events hold the same gene pair. Most single-cell calls are
#' one-cell artefacts, so a fusion is exported when it matches a DNA fusion, is
#' an Arriba known / recurrent fusion (tags), is seen in `min_cells` cells, or is
#' a high-confidence call of one cell with `single_min_reads` reads. With `slices`, every carrier
#' cell gets an IGV slice of its RNA reads within `pad` bp of its fusion
#' breakpoints (rna/reads/<rna_id>.bam; supporting reads tagged ZF:i:1).
#' Writes data/<P>/rna/fusions.json (gos-sc-rna-fusions/1).
#'
#' @param patient patient id
#' @param cell_dirs named character: rna_id -> per-cell folder (fusions.tsv, Aligned.sorted.bam)
#' @param out_dir patient folder of the gOS dataset (data/<P>)
#' @param cell_map named character rna_id -> gOS cell id (NA when the cell has no DNA)
#' @param dna_events optional path to the patient's filtered.events.json
#' @param keep_confidence Arriba confidence levels kept
#' @param slices write per-cell RNA BAM slices
#' @param min_cells fusions seen in at least this many cells are kept
#' @param single_min_reads a fusion of one cell is kept when high confidence with at least this many reads
#' @param max_fusions cap on exported fusions (DNA-matched ones always kept)
#' @param pad bp around each breakpoint in the slices
#' @param cores parallel cells
#' @return the fusions data.table (invisibly)
#' @export
#' @author Stanley Clarke
sc_export_rna_fusions <- function(patient, cell_dirs, out_dir, cell_map = NULL, dna_events = NULL,
                                  keep_confidence = c("high", "medium"), min_cells = 2, single_min_reads = 10,
                                  max_fusions = 2000, slices = TRUE, pad = 300, cores = 4) {
    files <- file.path(cell_dirs, "fusions.tsv")
    keep <- file.exists(files)
    calls <- data.table::rbindlist(parallel::mclapply(which(keep), function(i) {
        x <- data.table::fread(files[i], sep = "\t", quote = "", colClasses = "character", showProgress = FALSE)
        if (!nrow(x)) return(NULL)
        data.table::setnames(x, sub("^#", "", names(x)))
        x[, rna_id := names(cell_dirs)[i]]
        x
    }, mc.cores = cores), fill = TRUE)
    rna_dir <- file.path(out_dir, "rna")
    dir.create(rna_dir, showWarnings = FALSE, recursive = TRUE)
    empty <- list(format = "gos-sc-rna-fusions/1", patient = patient, n_cells_rna = sum(keep), fusions = list())
    if (!nrow(calls)) {
        jsonlite::write_json(empty, file.path(rna_dir, "fusions.json"), auto_unbox = TRUE)
        return(invisible(data.table::data.table()))
    }
    calls <- calls[confidence %in% keep_confidence]
    num <- function(v) suppressWarnings(as.integer(v))
    calls[, `:=`(split1 = num(split_reads1), split2 = num(split_reads2), discordant = num(discordant_mates))]
    calls[, fusion_id := paste0(gene1, "::", gene2, "|", breakpoint1, "|", breakpoint2)]
    calls[, cell_id := if (is.null(cell_map)) NA_character_ else unname(cell_map[rna_id])]
    rank <- c(high = 3L, medium = 2L, low = 1L)

    ## DNA fusions of the patient (gene pairs, either order)
    dna_pairs <- character(0)
    if (!is.null(dna_events) && file.exists(dna_events)) {
        ev <- data.table::as.data.table(jsonlite::fromJSON(dna_events))
        if ("vartype" %in% names(ev)) dna_pairs <- unique(ev[vartype %in% c("fusion", "outframe_fusion")]$gene)
    }
    pair_key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "::")
    dna_keys <- if (length(dna_pairs)) {
        g <- data.table::tstrsplit(dna_pairs, "::", fixed = TRUE)
        stats::setNames(dna_pairs, pair_key(g[[1]], g[[2]]))
    } else character(0)
    clean_gene <- function(g) sub("\\(.*$", "", sub(",.*$", "", g))   ## intergenic "A(123),B(456)" -> A

    fus <- calls[, .(gene1 = gene1[1], gene2 = gene2[1], breakpoint1 = breakpoint1[1], breakpoint2 = breakpoint2[1],
                     strand1 = `strand1(gene/fusion)`[1], strand2 = `strand2(gene/fusion)`[1],
                     site1 = site1[1], site2 = site2[1], type = type[1], reading_frame = reading_frame[1],
                     confidence = names(rank)[match(max(rank[confidence]), rank)],
                     n_cells = data.table::uniqueN(rna_id), n_cells_dna = data.table::uniqueN(cell_id[!is.na(cell_id)]),
                     split_reads = sum(split1 + split2, na.rm = TRUE), discordant_mates = sum(discordant, na.rm = TRUE)),
                 by = .(id = fusion_id)]
    fus[, dna_gene := unname(dna_keys[pair_key(clean_gene(gene1), clean_gene(gene2))])]
    known <- unique(calls[!is.na(tags) & tags != "." & grepl("Mitelman|known|COSMIC|CCLE", tags, ignore.case = TRUE)]$fusion_id)
    fus[, known := id %in% known]
    n_all <- nrow(fus)
    fus <- fus[!is.na(dna_gene) | known | n_cells >= min_cells |
               (confidence == "high" & split_reads + discordant_mates >= single_min_reads)]
    data.table::setorder(fus, -n_cells, -split_reads)
    fus <- fus[!is.na(dna_gene) | seq_len(.N) <= max_fusions]
    message(sprintf("sc_export_rna_fusions %s: %d of %d fusions kept (%d DNA-matched, %d known)", patient, nrow(fus), n_all,
                    sum(!is.na(fus$dna_gene)), sum(fus$known)))
    calls <- calls[fusion_id %in% fus$id]

    ## per-cell RNA slices around the breakpoints, supporting reads tagged ZF:i:1
    if (slices) {
        reads_dir <- file.path(rna_dir, "reads")
        dir.create(reads_dir, showWarnings = FALSE)
        by_cell <- split(calls, calls$rna_id)
        script <- system.file("extdata", "scripts", "sc_rna_fusion_slice.py", package = "Skilift")
        if (!nzchar(script)) script <- file.path(getNamespaceInfo("Skilift", "path"), "inst", "extdata", "scripts", "sc_rna_fusion_slice.py")
        py <- Sys.getenv("SC_PYTHON", "/gpfs/commons/groups/imielinski_lab/Software/miniforge3/envs/mskilab_ne1/bin/python")
        ok <- parallel::mclapply(names(by_cell), function(r) {
            bam <- file.path(cell_dirs[[r]], "Aligned.sorted.bam")
            if (!file.exists(paste0(bam, ".bai"))) return(FALSE)
            x <- by_cell[[r]]
            regions <- unique(unlist(lapply(c(x$breakpoint1, x$breakpoint2), function(b) {
                p <- strsplit(b, ":", fixed = TRUE)[[1]]
                pos <- as.integer(p[2])
                paste0(p[1], ":", max(1, pos - pad), "-", pos + pad)
            })))
            ids <- unique(unlist(strsplit(x$read_identifiers[!is.na(x$read_identifiers) & x$read_identifiers != "."], ",", fixed = TRUE)))
            idf <- tempfile(fileext = ".txt"); on.exit(unlink(idf))
            writeLines(ids, idf)
            out <- file.path(reads_dir, paste0(r, ".bam"))
            rc <- system2(py, c(shQuote(script), shQuote(bam), shQuote(out), shQuote(idf), shQuote(paste(regions, collapse = ","))),
                          stdout = FALSE, stderr = FALSE)
            rc == 0 && file.exists(paste0(out, ".bai"))
        }, mc.cores = cores, mc.preschedule = FALSE)
        sliced <- names(by_cell)[unlist(ok) %in% TRUE]
        message(sprintf("sc_export_rna_fusions %s: RNA slices for %d of %d carrier cells", patient, length(sliced), length(by_cell)))
    } else sliced <- character(0)

    cells_of <- split(calls, calls$fusion_id)
    out <- lapply(seq_len(nrow(fus)), function(i) {
        f <- fus[i]
        cc <- cells_of[[f$id]][order(-(split1 + split2 + discordant))]
        c(as.list(f[, !"dna_gene"]),
          list(dna_match = if (is.na(f$dna_gene)) NULL else list(kind = "gene_pair", event_gene = f$dna_gene),
               cells = lapply(seq_len(nrow(cc)), function(k) list(
                   rna_id = cc$rna_id[k], cell_id = if (is.na(cc$cell_id[k])) NULL else cc$cell_id[k],
                   split1 = cc$split1[k], split2 = cc$split2[k], discordant = cc$discordant[k], confidence = cc$confidence[k],
                   bam = if (cc$rna_id[k] %in% sliced) paste0("rna/reads/", cc$rna_id[k], ".bam") else NULL))))
    })
    res <- c(empty[c("format", "patient", "n_cells_rna")], list(fusions = out))
    jsonlite::write_json(res, file.path(rna_dir, "fusions.json"), auto_unbox = TRUE, null = "null", na = "null", digits = NA)
    invisible(fus)
}

## ------------------------------------------------------------------ splicing

#' @name sc_read_junctions
#' @title sc_read_junctions
#' @description
#'
#' Per-cell splice junctions from regtools `junctions extract` BED files:
#' intron coordinates (1-based first and last intronic base, chromosome
#' without "chr") and read counts.
#'
#' @param cell_dirs named character rna_id -> folder holding junctions.bed
#' @param cores parallel reads
#' @return data.table(rna_id, chromosome, start, end, count)
#' @export
#' @author Stanley Clarke
sc_read_junctions <- function(cell_dirs, cores = 4) {
    files <- file.path(cell_dirs, "junctions.bed")
    data.table::rbindlist(parallel::mclapply(which(file.exists(files)), function(i) {
        ## blockSizes "77,33" must stay text (fread would read it as the decimal 77.33)
        x <- data.table::fread(files[i], header = FALSE, select = c(1:3, 5, 11), showProgress = FALSE,
                               colClasses = list(character = 11), col.names = c("chr", "s", "e", "count", "blocks"))
        if (!nrow(x)) return(NULL)
        b <- data.table::tstrsplit(x$blocks, ",", fixed = TRUE, type.convert = TRUE)
        x[, .(rna_id = names(cell_dirs)[i], chromosome = sub("^chr", "", chr),
              start = s + b[[1]] + 1L, end = e - b[[2]], count = as.integer(count))]
    }, mc.cores = cores, mc.preschedule = FALSE))
}

## introns of the GTF transcripts: chromosome, start, end, strand, gene, transcript, exon numbers
sc_gtf_introns <- function(gtf) {
    ex <- data.table::as.data.table(rtracklayer::import(gtf, feature.type = "exon"))
    ex <- ex[, .(chromosome = sub("^chr", "", as.character(seqnames)), start, end, strand = as.character(strand),
                 gene = gene_name, transcript = transcript_id, exon = as.integer(exon_number))]
    data.table::setorder(ex, transcript, start)
    ex[, .(chromosome = chromosome[-.N], start = end[-.N] + 1L, end = start[-1] - 1L, strand = strand[1], gene = gene[1],
           exon_left = exon[-.N], exon_right = exon[-1]), by = transcript][start <= end]
}

## LeafCutter-style clusters: junctions of one chromosome/strand sharing a donor or acceptor
sc_junction_clusters <- function(j) {
    j <- data.table::copy(j)
    j[, jid := .I]
    parent <- seq_len(nrow(j))
    find <- function(x) { while (parent[x] != x) { parent[x] <<- parent[parent[x]]; x <- parent[x] }; x }
    for (key in list(c("chromosome", "strand", "start"), c("chromosome", "strand", "end"))) {
        grp <- j[, .(ids = list(jid)), by = key]$ids
        for (g in grp) if (length(g) > 1) for (k in g[-1]) { a <- find(g[1]); b <- find(k); if (a != b) parent[b] <- a }
    }
    j[, cluster := vapply(jid, find, 1L)]
    j[, n_in_cluster := .N, by = cluster]
    j[n_in_cluster > 1][, `:=`(jid = NULL, n_in_cluster = NULL)][]
}

## canonical-transcript exon junctions of known GBM splice variants (alt vs reference junction)
SC_SPLICE_VARIANTS <- list(
    list(id = "EGFRvIII", gene = "EGFR", transcript = "ENST00000275493", alt = c(1, 8), ref = c(1, 2),
         description = "EGFR exon 1 → exon 8 (Δ exons 2–7, EGFRvIII)"),
    list(id = "EGFRvII", gene = "EGFR", transcript = "ENST00000275493", alt = c(13, 16), ref = c(13, 14),
         description = "EGFR exon 13 → exon 16 (Δ exons 14–15, EGFRvII)"),
    list(id = "MET_ex14_skip", gene = "MET", transcript = "ENST00000397752", alt = c(13, 15), ref = c(13, 14),
         description = "MET exon 13 → exon 15 (exon 14 skipping)"),
    list(id = "PDGFRA_d8_9", gene = "PDGFRA", transcript = "ENST00000257290", alt = c(7, 10), ref = c(7, 8),
         description = "PDGFRA exon 7 → exon 10 (Δ exons 8–9)"))

sc_variant_junctions <- function(gtf_exons) {
    rbindlist_safe <- function(x) data.table::rbindlist(Filter(Negate(is.null), x))
    rbindlist_safe(lapply(SC_SPLICE_VARIANTS, function(v) {
        e <- gtf_exons[sub("\\..*$", "", transcript) == v$transcript]
        if (!nrow(e)) return(NULL)
        jn <- function(a, b) {
            l <- e[exon == a]; r <- e[exon == b]
            if (!nrow(l) || !nrow(r)) return(c(NA, NA))
            if (l$strand[1] == "-") c(r$end[1] + 1L, l$start[1] - 1L) else c(l$end[1] + 1L, r$start[1] - 1L)
        }
        a <- jn(v$alt[1], v$alt[2]); r <- jn(v$ref[1], v$ref[2])
        data.table::data.table(id = v$id, gene = v$gene, description = v$description, chromosome = e$chromosome[1],
                               alt_start = a[1], alt_end = a[2], ref_start = r[1], ref_end = r[2])
    }))
}

#' @name sc_export_rna_splicing
#' @title sc_export_rna_splicing
#' @description
#'
#' Splice junction usage of the single-cell RNA for gOS. Junctions of every
#' cell (regtools) are pooled per patient; junctions with at least `min_reads`
#' reads in `min_cells` cells (cohort) are grouped into LeafCutter-style
#' clusters (shared donor or acceptor, same strand). Strand and gene come from
#' the GTF introns (annotated) or, for novel junctions, the GTF gene they fall
#' in. Writes:
#'  - data/_cohort/rna/splicing.json: clusters whose junction usage differs
#'    between patients (chi-square on patient x junction pseudobulk counts,
#'    BH q, max delta-PSI between patients with >= `min_patient_reads`),
#'  - data/<P>/rna/splicing.json per patient: per-cell counts of the clusters
#'    most variable between the patient's cells plus the cohort clusters, and
#'    per-cell alt / reference reads of known GBM splice variants (EGFRvIII,
#'    EGFRvII, MET exon 14 skipping, PDGFRA delta 8-9).
#'
#' @param cell_dirs list patient -> named character (rna_id -> per-cell folder with junctions.bed)
#' @param data_dir data folder of the gOS dataset
#' @param gtf GTF of the alignment reference
#' @param cell_maps list patient -> named character rna_id -> gOS cell id
#' @param min_reads,min_cells cohort junction filter
#' @param min_patient_reads patient cluster reads for delta-PSI
#' @param n_cohort,n_patient clusters exported
#' @param cores parallel reads
#' @return list(cohort = cohort clusters, patients = per-patient cluster ids) invisibly
#' @export
#' @author Stanley Clarke
sc_export_rna_splicing <- function(cell_dirs, data_dir, gtf, cell_maps = list(), min_reads = 30, min_cells = 10,
                                   min_patient_reads = 30, n_cohort = 300, n_patient = 150, cores = 8) {
    message("sc_export_rna_splicing: GTF introns")
    ex <- data.table::as.data.table(rtracklayer::import(gtf, feature.type = "exon"))
    ex <- ex[, .(chromosome = sub("^chr", "", as.character(seqnames)), start, end, strand = as.character(strand),
                 gene = gene_name, transcript = transcript_id, exon = as.integer(exon_number))]
    data.table::setorder(ex, transcript, start)
    introns <- ex[, .(chromosome = chromosome[-.N], start = end[-.N] + 1L, end = start[-1] - 1L, strand = strand[1], gene = gene[1]),
                  by = transcript][start <= end]
    ann <- unique(introns[, .(chromosome, start, end, strand, gene)], by = c("chromosome", "start", "end"))
    genes <- ex[, .(start = min(start), end = max(end), strand = strand[1], chromosome = chromosome[1]), by = gene]
    variants <- sc_variant_junctions(ex)

    ## per patient: per-cell junction counts
    per <- list()
    for (p in names(cell_dirs)) {
        message("sc_export_rna_splicing: reading junctions of ", p, " (", length(cell_dirs[[p]]), " cells)")
        per[[p]] <- sc_read_junctions(cell_dirs[[p]], cores = cores)
    }
    allj <- data.table::rbindlist(lapply(names(per), function(p) per[[p]][, .(reads = sum(count), cells = .N), by = .(chromosome, start, end)][, patient := p]))
    tot <- allj[, .(reads = sum(reads), cells = sum(cells)), by = .(chromosome, start, end)][reads >= min_reads & cells >= min_cells]
    ## strand / gene: annotated intron, else the gene the junction lies in (unique strand), else drop
    tot <- merge(tot, ann, by = c("chromosome", "start", "end"), all.x = TRUE)
    tot[, annotated := !is.na(strand)]
    nov <- tot[is.na(strand)]
    if (nrow(nov)) {
        g <- genes[, .(gchr = chromosome, gs = start, ge = end, gstrand = strand, gname = gene)]
        hit <- g[nov, on = .(gchr = chromosome, gs <= start, ge >= end), .(chromosome = i.chromosome, start = i.start, end = i.end, gstrand, gname),
                 allow.cartesian = TRUE, nomatch = 0L]
        hit <- hit[, .(strand = if (data.table::uniqueN(gstrand) == 1) gstrand[1] else NA_character_, gene = gname[1]), by = .(chromosome, start, end)]
        tot[hit, on = .(chromosome, start, end), `:=`(strand = i.strand, gene = i.gene)]
    }
    tot <- tot[!is.na(strand)]
    cl <- sc_junction_clusters(tot[, .(chromosome, strand, start, end, gene, annotated)])
    cl[, gene := { g <- gene[!is.na(gene)]; if (length(g)) names(sort(table(g), decreasing = TRUE))[1] else NA_character_ }, by = cluster]
    cl[, cluster_id := paste0("clu_", chromosome, "_", min(start), "_", max(end), "_", strand), by = cluster]
    data.table::setorder(cl, cluster_id, start, end)
    message("sc_export_rna_splicing: ", data.table::uniqueN(cl$cluster_id), " clusters from ", nrow(cl), " junctions")

    ## cohort: patient x junction pseudobulk per cluster
    pb <- merge(allj, cl[, .(chromosome, start, end, cluster_id)], by = c("chromosome", "start", "end"))
    patients <- names(per)
    stats <- pb[, {
        m <- data.table::dcast(.SD, patient ~ start + end, value.var = "reads", fill = 0, fun.aggregate = sum)
        mat <- as.matrix(m[, -1, with = FALSE]); rownames(mat) <- m$patient
        mat <- mat[rowSums(mat) >= min_patient_reads, , drop = FALSE]
        if (nrow(mat) < 2 || ncol(mat) < 2) list(p = NA_real_, max_dpsi = NA_real_, n_patients = nrow(mat)) else {
            psi <- mat / rowSums(mat)
            list(p = suppressWarnings(stats::chisq.test(mat)$p.value), max_dpsi = max(apply(psi, 2, function(v) diff(range(v)))),
                 n_patients = nrow(mat))
        }
    }, by = cluster_id, .SDcols = c("patient", "start", "end", "reads")]
    stats <- stats[!is.na(p)]
    stats[, q := stats::p.adjust(p, "BH")]
    top <- stats[q < 0.05][order(-max_dpsi, q)][seq_len(min(.N, n_cohort))]
    cohort_clusters <- lapply(top$cluster_id, function(id) {
        jj <- cl[cluster_id == id]
        u <- pb[cluster_id == id]
        s <- top[cluster_id == id]
        list(id = id, gene = jj$gene[1], chromosome = jj$chromosome[1], strand = jj$strand[1],
             start = min(jj$start), end = max(jj$end), p = s$p, q = s$q, max_dpsi = s$max_dpsi, n_patients = s$n_patients,
             junctions = lapply(seq_len(nrow(jj)), function(k) list(start = jj$start[k], end = jj$end[k], annotated = jj$annotated[k])),
             usage = stats::setNames(lapply(patients, function(p) {
                 cnt <- vapply(seq_len(nrow(jj)), function(k) as.numeric(sum(u[patient == p & start == jj$start[k] & end == jj$end[k]]$reads)), 0)
                 list(counts = cnt, total = sum(cnt), n_cells = sum(per[[p]][chromosome == jj$chromosome[1] & start %in% jj$start & end %in% jj$end, data.table::uniqueN(rna_id)]))
             }), patients))
    })
    cdir <- file.path(data_dir, "_cohort", "rna")
    dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(list(format = "gos-sc-splicing-cohort/1", patients = patients, n_clusters_tested = nrow(stats), clusters = cohort_clusters),
                         file.path(cdir, "splicing.json"), auto_unbox = TRUE, digits = NA, null = "null", na = "null")

    ## per patient: clusters most variable between cells (+ cohort clusters), per-cell counts; known variants
    pat_ids <- list()
    for (p in patients) {
        x <- merge(per[[p]], cl[, .(chromosome, start, end, cluster_id)], by = c("chromosome", "start", "end"))
        x[, cl_total := sum(count), by = .(rna_id, cluster_id)]
        ## variability: SD across cells (>= 5 reads in the cluster) of the cluster's top-junction PSI
        var <- x[cl_total >= 5, {
            topj <- .SD[, .(r = sum(count)), by = .(start, end)][which.max(r)]
            v <- .SD[start == topj$start & end == topj$end, .(psi = count / cl_total), by = rna_id]
            ## cells without the top junction have PSI 0
            n <- data.table::uniqueN(rna_id)
            psi <- c(v$psi, rep(0, n - nrow(v)))
            list(n_cells = n, sd = if (n >= 10) stats::sd(psi) else NA_real_)
        }, by = cluster_id][!is.na(sd)]
        ids <- union(intersect(top$cluster_id, unique(x$cluster_id)), var[order(-sd)][seq_len(min(.N, n_patient))]$cluster_id)
        pat_ids[[p]] <- ids
        clusters <- lapply(ids, function(id) {
            jj <- cl[cluster_id == id]
            xx <- x[cluster_id == id]
            w <- data.table::dcast(xx, rna_id ~ start + end, value.var = "count", fill = 0, fun.aggregate = sum)
            keyj <- paste(jj$start, jj$end, sep = "_")
            for (k in setdiff(keyj, names(w))) w[, (k) := 0L]
            m <- as.matrix(w[, keyj, with = FALSE])
            list(id = id, gene = jj$gene[1], chromosome = jj$chromosome[1], strand = jj$strand[1],
                 junctions = lapply(seq_len(nrow(jj)), function(k) list(start = jj$start[k], end = jj$end[k], annotated = jj$annotated[k])),
                 cells = stats::setNames(lapply(seq_len(nrow(m)), function(i) unname(as.integer(m[i, ]))), w$rna_id))
        })
        vv <- lapply(seq_len(nrow(variants)), function(i) {
            v <- variants[i]
            alt <- per[[p]][chromosome == v$chromosome & start == v$alt_start & end == v$alt_end, .(alt = sum(count)), by = rna_id]
            ref <- per[[p]][chromosome == v$chromosome & start == v$ref_start & end == v$ref_end, .(ref = sum(count)), by = rna_id]
            ar <- merge(alt, ref, by = "rna_id", all = TRUE)
            ar[is.na(alt), alt := 0L][is.na(ref), ref := 0L]
            list(id = v$id, gene = v$gene, description = v$description,
                 alt_junction = list(chromosome = v$chromosome, start = v$alt_start, end = v$alt_end),
                 ref_junction = list(chromosome = v$chromosome, start = v$ref_start, end = v$ref_end),
                 n_cells_alt = sum(ar$alt > 0),
                 cells = stats::setNames(lapply(seq_len(nrow(ar)), function(k) c(ar$alt[k], ar$ref[k])), ar$rna_id))
        })
        cm <- cell_maps[[p]]
        pdir <- file.path(data_dir, p, "rna")
        dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
        jsonlite::write_json(list(format = "gos-sc-splicing/1", patient = p, n_cells = data.table::uniqueN(per[[p]]$rna_id),
                                  variants = vv, clusters = clusters,
                                  cell_map = if (length(cm)) as.list(cm) else structure(list(), names = character(0))),
                             file.path(pdir, "splicing.json"), auto_unbox = TRUE, digits = NA, null = "null", na = "null")
        message("sc_export_rna_splicing: ", p, ": ", length(ids), " clusters, variants: ",
                paste(vapply(vv, function(v) paste0(v$id, " ", v$n_cells_alt), ""), collapse = ", "))
    }
    invisible(list(cohort = top, patients = pat_ids))
}
