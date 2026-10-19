## Transcript-level gene models for the gOS gene track (collapsed / squished / expanded isoform modes).

## transcript ids (no version) of a GTF carrying `tag` on their transcript line (empty when absent or unreadable)
sc_gtf_tagged <- function(gtf, tag) {
    if (is.null(gtf) || !nzchar(gtf) || !file.exists(gtf)) return(character(0))
    x <- tryCatch(system2("grep", c("-F", shQuote(sprintf('tag "%s"', tag)), shQuote(gtf)), stdout = TRUE, stderr = FALSE),
                  error = function(e) character(0))
    x <- grep("\ttranscript\t", x, value = TRUE)
    unique(sub("\\..*$", "", sub('.*transcript_id "([^"]+)".*', "\\1", x)))
}

sc_json_str <- function(x) paste0('"', gsub('(["\\\\])', "\\\\\\1", ifelse(is.na(x), "", x)), '"')

#' @name sc_gene_models_table
#' @title sc_gene_models_table
#' @description
#'
#' Gene / transcript / exon tables of a GENCODE GRanges (as read from
#' cfg$gencode: a `type` column with gene, transcript, exon, CDS and
#' stop_codon rows). The CDS span of a transcript includes its stop codon
#' (GENCODE CDS rows exclude it). One transcript per gene is canonical: the
#' MANE Select transcript, else the Ensembl canonical one (ids from
#' `mane` / `ensembl_canonical`, or the GTF tags when present), else the
#' best-supported one (appris principal, CCDS, basic, longest CDS, longest
#' exons). Chromosomes lose any "chr" prefix.
#'
#' @param gr GENCODE GRanges
#' @param mane transcript ids (no version) tagged MANE_Select
#' @param ensembl_canonical transcript ids (no version) tagged Ensembl_canonical
#' @param drop_types transcript types to leave out (canonical transcripts are always kept)
#' @return list(genes, transcripts, exons) data.tables
#' @export
#' @author Stanley Clarke
sc_gene_models_table <- function(gr, mane = character(0), ensembl_canonical = character(0), drop_types = character(0)) {
    d <- data.table::as.data.table(gr)
    d <- d[type %in% c("gene", "transcript", "exon", "CDS", "stop_codon")]
    d[, chromosome := sub("^chr", "", as.character(seqnames))]
    d[, strand := as.character(strand)]
    nover <- function(x) sub("\\..*$", "", x)
    if (!"tag" %in% names(d)) d[, tag := NA_character_]
    if (!"transcript_name" %in% names(d)) d[, transcript_name := NA_character_]
    if (!"transcript_type" %in% names(d)) d[, transcript_type := NA_character_]
    tag <- function(x) { x <- as.character(x); x[is.na(x)] <- ""; x }

    genes <- d[type == "gene", .(gene_id = nover(gene_id), name = gene_name, gene_type = as.character(gene_type),
                                 chromosome, strand, start = as.integer(start), end = as.integer(end))]
    genes <- unique(genes, by = "gene_id")
    tx <- d[type == "transcript", .(transcript_id = nover(transcript_id), gene_id = nover(gene_id), name = transcript_name,
                                    tx_type = as.character(transcript_type), tag = tag(tag),
                                    start = as.integer(start), end = as.integer(end))]
    tx <- unique(tx, by = "transcript_id")
    ex <- d[type == "exon", .(transcript_id = nover(transcript_id), start = as.integer(start), end = as.integer(end))]
    ex <- unique(ex)
    data.table::setorder(ex, transcript_id, start)
    cds <- d[type %in% c("CDS", "stop_codon"), .(cds_start = min(as.integer(start)), cds_end = max(as.integer(end))),
             by = .(transcript_id = nover(transcript_id))]
    tx <- merge(tx, cds, by = "transcript_id", all.x = TRUE)
    elen <- ex[, .(exon_len = sum(end - start + 1L), n_exons = .N), by = transcript_id]
    tx <- merge(tx, elen, by = "transcript_id", all.x = TRUE)
    clen <- d[type == "CDS", .(cds_len = sum(as.integer(end) - as.integer(start) + 1L)), by = .(transcript_id = nover(transcript_id))]
    tx <- merge(tx, clen, by = "transcript_id", all.x = TRUE)
    tx[is.na(cds_len), cds_len := 0L]
    tx[is.na(exon_len), exon_len := 0L]

    tx[, mane := transcript_id %in% mane | grepl("(^|,)MANE_Select(,|$)", tag)]
    tx[, ens_canonical := transcript_id %in% ensembl_canonical | grepl("(^|,)Ensembl_canonical(,|$)", tag)]
    tx[, rank_score := 1e10 * mane + 1e9 * ens_canonical + 1e8 * grepl("appris_principal", tag) + 1e7 * grepl("CCDS", tag) +
             1e6 * grepl("(^|,)basic(,|$)", tag) + pmin(cds_len, 999999) / 10 + pmin(exon_len, 999999) / 1e7]
    data.table::setorder(tx, gene_id, -rank_score, transcript_id)
    tx[, canonical := !duplicated(gene_id)]
    tx[, canonical_basis := ifelse(!canonical, NA_character_, ifelse(mane, "MANE_Select", ifelse(ens_canonical, "Ensembl_canonical", "inferred")))]
    if (length(drop_types)) tx <- tx[canonical | !(tx_type %in% drop_types)]
    tx <- tx[transcript_id %in% ex$transcript_id]
    ex <- ex[transcript_id %in% tx$transcript_id]
    genes <- genes[gene_id %in% tx$gene_id]
    tx[, `:=`(mane = NULL, ens_canonical = NULL, rank_score = NULL)]
    list(genes = genes[], transcripts = tx[], exons = ex[])
}

## JSON text of one chromosome's models (genes ordered by start; canonical transcript first, then coding, then by length)
sc_gene_models_json <- function(m, chr, source = "") {
    g <- m$genes[m$genes$chromosome == chr]
    data.table::setorder(g, start, end)
    tx <- m$transcripts[gene_id %in% g$gene_id]
    ex <- m$exons[transcript_id %in% tx$transcript_id]
    exs <- ex[, .(exons = paste0("[", paste0("[", start, ",", end, "]", collapse = ","), "]")), by = transcript_id]
    tx <- merge(tx, exs, by = "transcript_id")
    data.table::setorder(tx, gene_id, -canonical, -cds_len, -exon_len, transcript_id)
    tx[, json := paste0('{"id":', sc_json_str(transcript_id), ',"name":', sc_json_str(name), ',"type":', sc_json_str(tx_type),
                        ',"canonical":', ifelse(canonical, 1L, 0L),
                        ifelse(canonical, paste0(',"basis":', sc_json_str(canonical_basis)), ""),
                        ',"cds":', ifelse(is.na(cds_start), "null", paste0("[", cds_start, ",", cds_end, "]")),
                        ',"exons":', exons, "}")]
    txs <- tx[, .(transcripts = paste0("[", paste(json, collapse = ","), "]"), n = .N), by = gene_id]
    g <- merge(g, txs, by = "gene_id")
    data.table::setorder(g, start, end)
    gj <- paste0('{"id":', sc_json_str(g$gene_id), ',"name":', sc_json_str(g$name), ',"type":', sc_json_str(g$gene_type),
                 ',"strand":', sc_json_str(g$strand), ',"start":', g$start, ',"end":', g$end, ',"transcripts":', g$transcripts, "}")
    paste0('{"version":1,"chromosome":', sc_json_str(chr), ',"source":', sc_json_str(source),
           ',"n_genes":', nrow(g), ',"n_transcripts":', nrow(tx), ',"genes":[\n', paste(gj, collapse = ",\n"), "\n]}\n")
}

#' @name sc_export_gene_models
#' @title sc_export_gene_models
#' @description
#'
#' Writes transcript-level gene models for the gOS gene track, one JSON per
#' chromosome (out_dir/chr7.json): genes [{ id, name, type, strand, start,
#' end, transcripts: [{ id, name, type, canonical (0/1; basis MANE_Select /
#' Ensembl_canonical / inferred), cds [start, end] (with stop codon) or null,
#' exons [[start, end], ...] }] }], 1-based closed coordinates, plus
#' out_dir/index.json (chromosomes, counts, bytes, dropped transcript types).
#' A chromosome whose file exceeds max_mb is rewritten without its
#' non-canonical `drop_types` transcripts (retained intron, NMD, ...), and
#' the index says so. All gene types are kept (lncRNA included).
#'
#' @param out_dir output folder
#' @param cfg sc_annotation_config() (gencode: GENCODE rds; canonical_gtf: GTF with MANE / Ensembl_canonical tags)
#' @param gr optional GENCODE GRanges (default readRDS(cfg$gencode))
#' @param chromosomes chromosomes to write (default 1-22, X, Y, M present in the data)
#' @param max_mb per-file size cap
#' @param drop_types transcript types dropped from chromosomes over the cap
#' @return index data.table (invisibly)
#' @export
#' @author Stanley Clarke
sc_export_gene_models <- function(out_dir, cfg = sc_annotation_config(), gr = NULL, chromosomes = NULL, max_mb = 4,
                                  drop_types = c("retained_intron", "nonsense_mediated_decay", "non_stop_decay", "processed_transcript")) {
    if (is.null(gr)) gr <- readRDS(cfg$gencode)
    mane <- sc_gtf_tagged(cfg$canonical_gtf, "MANE_Select")
    ens <- sc_gtf_tagged(cfg$canonical_gtf, "Ensembl_canonical")
    source <- paste0(basename(as.character(cfg$gencode)), if (length(mane) || length(ens)) paste0("; canonical tags from ", basename(cfg$canonical_gtf)))
    m <- sc_gene_models_table(gr, mane = mane, ensembl_canonical = ens)
    if (is.null(chromosomes)) chromosomes <- intersect(c(1:22, "X", "Y", "M"), unique(m$genes$chromosome))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    small <- NULL
    idx <- data.table::rbindlist(lapply(chromosomes, function(chr) {
        txt <- sc_gene_models_json(m, chr, source)
        dropped <- ""
        if (nchar(txt, type = "bytes") > max_mb * 1e6) {
            if (is.null(small)) small <<- sc_gene_models_table(gr, mane = mane, ensembl_canonical = ens, drop_types = drop_types)
            txt <- sc_gene_models_json(small, chr, source)
            dropped <- paste(drop_types, collapse = ",")
        }
        f <- file.path(out_dir, paste0("chr", chr, ".json"))
        writeLines(txt, f, sep = "", useBytes = TRUE)
        mm <- if (nzchar(dropped)) small else m
        g <- mm$genes[chromosome == chr]
        data.table::data.table(chromosome = chr, file = basename(f), bytes = file.size(f), n_genes = nrow(g),
                               n_coding = sum(g$gene_type == "protein_coding"), n_lncrna = sum(g$gene_type == "lncRNA"),
                               n_transcripts = sum(mm$transcripts$gene_id %in% g$gene_id), dropped_types = dropped)
    }))
    canon <- m$transcripts[canonical == TRUE, .N, by = canonical_basis]
    meta <- list(version = 1L, source = source, coordinates = "1-based closed, hg38, no chr prefix",
                 canonical = stats::setNames(as.list(canon$N), canon$canonical_basis), chromosomes = idx)
    jsonlite::write_json(meta, file.path(out_dir, "index.json"), auto_unbox = TRUE, pretty = TRUE, dataframe = "rows")
    invisible(idx)
}
