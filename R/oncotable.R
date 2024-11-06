#' @title process_gencode
#' @description
#'
#' Helper script to process gencode parameter
#'
#' @param gencode path to gencode file. Gencode file must be either rds or some format accepted by rtracklayer::import (e.g. GTF) with just a single entry for each gene (so gencode entries for each gene are collapse to a single range). The input could be .gtf or .rds with GRanges object, or a GRanges object i.e. resulting from importing the (appropriate) GENCODE .gtf via rtracklayer, note: this input is only used in CNA to gene mapping.
#' @return gencode_gr GRanges
#' @author Marcin Imielinski
process_gencode = function(gencode = NULL){
  if (is.null(gencode))
    stop('gencode file must be provided')
  else if (is.character(gencode))
  {
    if (grepl('.rds$', gencode))
      gencode = readRDS(gencode)
    else
      gencode = rtracklayer::import(gencode)
  }
  return(gencode)
}

#' @title collect_gene_fusions
#' @description
#' Collects gene fusion data from a specified file and processes it.
#'
#' @param fusions Path to the fusions.rds file.
#' @param pge GRanges object with gencode annotations.
#' @param verbose Logical flag to indicate if messages should be printed.
#' @return A data.table containing processed gene fusion information.
collect_gene_fusions <- function(fusions, pge, verbose = TRUE) {
  if (is.null(fusions) || !file.exists(fusions)) {
    if (verbose) message('Fusions file is missing or does not exist.')
    return(data.table(type = NA, source = 'fusions'))
  }
  
  if (verbose) message('pulling fusions')
  fus <- readRDS(fusions)$meta
  
  if (nrow(fus) == 0) {
    if (verbose) message('No fusions found in the file.')
    return(data.table(type = NA, source = 'fusions'))
  }
  
  non_silent_fusions <- fus[silent == FALSE, ]
  unique_fusions <- non_silent_fusions[!duplicated(genes), ]
  unique_fusions[, vartype := ifelse(in.frame == TRUE, 'fusion', 'outframe_fusion')]
  
  split_genes <- function(genes) unlist(strsplit(genes, ','))
  gene_lengths <- function(genes) sapply(strsplit(genes, ','), length)
  
  fus <- unique_fusions[, .(
    gene = split_genes(genes),
    vartype = rep(vartype, gene_lengths(genes)),
    fusion_genes = rep(genes, gene_lengths(genes))
  )][, `:=`(track = 'variants', type = vartype, source = 'fusions')]
  
  get_gene_coords <- function(genes) {
    coords <- lapply(genes, function(gene) {
      gene_ranges <- pge[mcols(pge)$gene_name == gene]
      paste0(seqnames(gene_ranges), ":", start(gene_ranges), "-", end(gene_ranges))
    })
    paste(unlist(coords), collapse = ",")
  }
  
  fus[, fusion_gene_coords := unlist(lapply(strsplit(fusion_genes, ','), get_gene_coords))]
  
  return(fus)
}

#' @name oncotable
#' @title oncotable
#' @description
#'
#' @param annotated_bcf Path to annotated.bcf file
#' @param fusions Path to fusion.rds file
#' @param jabba_rds Path to jabba.simple.rds file
#' @param complex Path to complex.rds file
#' @param signature_counts Path to signature_counts.txt file
#' @param karyograph Optional path to the karyograph.rds file
#' @param gencode_gr GRanges object with gencode annotations 
#' @param amp.thresh SCNA amplification threshold to call an amp as a function of ploidy (4)
#' @param del.thresh SCNA deletion threshold for (het) del as a function of ploidy (by default cn = 1 will be called del, but this allows additoinal regions in high ploidy tumors to be considered het dels)
#' @param verbose logical flag 
#' @export
oncotable = function(
  pair,
  annotated_bcf = NULL,
  fusions = NULL,
  jabba_rds = NULL,
  complex = NULL,
  signature_counts = NULL,
  gencode,
  verbose = TRUE,
  amp.thresh = 4,
  filter = 'PASS',
  del.thresh = 0.5,
  karyograph = NULL
) {
  out <- data.table()

  if ('type' %in% names(mcols(gencode))) {
    pge <- gencode %Q% (type == 'gene' & gene_type == 'protein_coding')
  } else {
    pge <- gencode %Q% (gene_type == 'protein_coding')
  }

  ## collect gene fusions
  out <- rbind(out, collect_gene_fusions(
    fusions,
    pge,
    verbose
  ), fill = TRUE, use.names = TRUE)

  ## collect complex events
  out <- rbind(
    out,
    collect_complex_events(complex, verbose
  ), fill = TRUE, use.names = TRUE)

  ## collect copy number / jabba
  out <- rbind(
    out,
    collect_copy_number_jabba(jabba_rds, pge, amp.thresh, del.thresh, verbose, karyograph),
    fill = TRUE,
    use.names = TRUE
  )

  # collect signatures
  out <- rbind(out, collect_signatures(signature_counts, verbose), fill = TRUE, use.names = TRUE)

  ## collect gene mutations
  if (!is.null(annotated_bcf) && file.exists(annotated_bcf)) {
    if (verbose) message('pulling annotated_bcf using FILTER=', filter)
    local_bcftools_path <- Sys.which("bcftools")
    local_bcftools_path <- ifelse(local_bcftools_path == "", stop("bcftools not found in the system PATH. Please install or module load bcftools."), local_bcftools_path)
    message("bcftools found at: ", local_bcftools_path)
    bcf <- skitools::grok_bcf(annotated_bcf, label = "sample", long = TRUE, filter = filter, bpath = local_bcftools_path)
    if (verbose) message(length(bcf), ' variants pass filter')
    genome.size <- sum(seqlengths(bcf), na.rm = TRUE) / 1e6
    if (is.na(genome.size)) genome.size <- sum(seqlengths(gG(jabba = jabba_rds)), na.rm = TRUE) / 1e6
    nmut <- data.table(as.character(seqnames(bcf)), start(bcf), end(bcf), bcf$REF, bcf$ALT) %>% unique %>% nrow
    mut.density <- data.table(value = c(nmut, nmut / genome.size), type = c('count', 'density'), track = 'tmb', source = 'annotated_bcf')
    out <- rbind(out, mut.density, fill = TRUE, use.names = TRUE)
    keepeff <- c('trunc', 'cnadel', 'cnadup', 'complexsv', 'splice', 'inframe_indel', 'fusion', 'missense', 'promoter', 'regulatory', 'mir')
    bcf <- bcf[bcf$short %in% keepeff]
    if (verbose) message(length(bcf), ' variants pass keepeff')
    vars <- NULL
    if (length(bcf)) {
      bcf$variant.g <- paste0(seqnames(bcf), ':', start(bcf), '-', end(bcf), ' ', bcf$REF, '>', bcf$ALT)
      vars <- gr2dt(bcf)[, .(gene, vartype, variant.g, variant.p, distance, annotation, type = short, track = 'variants', source = 'annotated_bcf')]
      setkey(vars, variant.g)
      vars <- vars[, .SD[1], by = variant.g]
    }
    out <- rbind(out, vars, fill = TRUE, use.names = TRUE)
  } else {
    out <- rbind(out, data.table(type = NA, source = 'annotated_bcf'), fill = TRUE, use.names = TRUE)
  }

  out$id = pair

  if (verbose) message('done processing sample')
  return(out)
}

#' @title collect_signatures
#' @description
#' Collects signature data from a specified file and processes it.
#'
#' @param signature_counts Path to the signature_counts.txt file.
#' @param verbose Logical flag to indicate if messages should be printed.
#' @return A data.table containing processed signature information.
collect_signatures <- function(signature_counts, verbose = TRUE) {
  #TODO: update this function to use sigprofiler
  if (!is.null(signature_counts) && file.exists(signature_counts)) {
    if (verbose) message('pulling signature_counts')
    sig <- fread(signature_counts)
    sig <- sig[, .(value = num_events, type = Signature, etiology = Etiology, frac = frac.events, track = 'signature', source = 'signature_counts')]
    return(sig)
  } else {
    return(data.table(type = NA, source = 'signature_counts'))
  }
}
#' @description
#' Collects complex events from a specified file and processes them.
#'
#' @param complex Path to the complex.rds file.
#' @param verbose Logical flag to indicate if messages should be printed.
#' @return A data.table containing processed complex event information.
collect_complex_events <- function(complex, verbose = TRUE) {
  if (is.null(complex) || !file.exists(complex)) {
    if (verbose) message('Complex events file is missing or does not exist.')
    return(data.table(type = NA, source = 'complex'))
  }
  
  if (verbose) message('pulling complex events')
  sv <- readRDS(complex)$meta$events
  
  if (nrow(sv) == 0) {
    if (verbose) message('No complex events found in the file.')
    return(data.table(type = NA, source = 'complex'))
  }
  
  sv_summary <- sv[, .(value = .N), by = type]
  simple_sv_types <- c('del', 'dup', 'invdup', 'tra', 'inv')
  sv_summary[, track := ifelse(type %in% simple_sv_types, 'simple sv', 'complex sv')]
  sv_summary[, source := 'complex']
  
  return(sv_summary)
}

#' @title collect_copy_number_jabba
#' @description
#' Collects copy number and jabba data from a specified file and processes it.
#'
#' @param jabba_rds Path to the jabba.simple.rds file.
#' @param pge GRanges object with gencode annotations.
#' @param amp.thresh SCNA amplification threshold.
#' @param del.thresh SCNA deletion threshold.
#' @param verbose Logical flag to indicate if messages should be printed.
#' @param karyograph Optional path to the karyograph.rds file.
#' @return A data.table containing processed copy number and jabba information.
collect_copy_number_jabba <- function(
  jabba_rds,
  pge,
  amp.thresh,
  del.thresh,
  verbose = TRUE,
  karyograph = NULL
) {
  if (is.null(jabba_rds) || !file.exists(jabba_rds)) {
    if (verbose) message('Jabba RDS file is missing or does not exist.')
    return(data.table(type = NA, source = 'jabba_rds'))
  }
  
  if (verbose) message('pulling jabba_rds to get SCNA and purity / ploidy')
  jab <- readRDS(jabba_rds)
  jabpurity <- if (is.null(jab$meta$purity)) { jab$meta$purity } else { jab$purity }
  jabploidy <- if (is.null(jab$meta$ploidy)) { jab$meta$ploidy } else { jab$ploidy }
  result <- data.table(
    value = c(jabpurity, jabploidy),
    type = c('purity', 'ploidy'),
    track = 'pp'
  )
  
  # get the ncn data from jabba
  nseg <- NULL

  # Don't fail out if karyograph isn't found
  # Just provide nseg as NULL, and get_gene_ampdels_from_jabba
  # assumes ncn = 2 (pretty safe assumption)
  if (!is.null(karyograph) && file.exists(karyograph)) {
    nseg <- readRDS(karyograph)$segstats[, c("ncn")]
  }
  scna <- skitools::get_gene_ampdels_from_jabba(
    jab,
    amp.thresh = amp.thresh,
    del.thresh = del.thresh,
    pge = pge,
    nseg = nseg
  )
  
  if (nrow(scna)) {
    scna[, track := 'variants'][, source := 'jabba_rds'][, vartype := 'scna']
    result <- rbind(
      result,
      scna[, .(value = min_cn, type, track, gene = gene_name)],
      fill = TRUE,
      use.names = TRUE
    )
  }
  
  return(result)
}
