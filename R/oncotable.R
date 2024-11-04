#' @title process_gencode
#' @description
#'
#' Helper script to process gencode parameter
#'
#' @param gencode path to gencode file. Gencode file must be either rds or some format accepted by rtracklayer::import (e.g. GTF)
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

#' @name oncotable
#' @title oncotable
#' @description
#'
#' THIS ONCOTABLE IS ADAPTED FOR SKILIFT
#' in that it will report only the putatively impacted transcript based on SnpEff outputs to remove duplicated SNV entries
#'
#' Takes as input (keyed) "tumors" (aka pairs) table which a metadata table with specific
#' columns pointing to paths corresponding to one or more of the following pipeline outputs:
#'
#' $annotated_bcf  Path to annotated.bcf file that is the primary output of SnpEff module from which TMB and basic mutation
#' descriptions are extracted along with their basic categories (these will comprising the core of the oncoplot are computed)
#' 
#' $fusions  Path to fusion.rds file that is the primary output of the Fusions modjle, from which protein coding fusions will
#' be computed for
#' 
#' $jabba_rds  Path to jabba.simple.rds output representing primary output of JaBbA module from which SCNA and overall
#' junction burden are computed
#' 
#' $complex    Path to complex.rds gGnome cached object that is the primary output of Events module, from which simple
#' and complex event burdens are computed
#' 
#' $signature_counts Path to signature_counts.txt that is the primary output of Signatures module from which SNV signature
#' counts are computed
#' 
#' The function then outputs a melted data.table of "interesting" features that can be saved and/or immediately output
#' into oncoprint.  This data.table will at the very least have fields $id $type (event type), $track, and  $source
#' populated in addition to a few other data type specific columns.
#'
#' The $source column is the name of the column of tumors from which that data was extracted, and track is a grouping
#' variable that allows separation of the various data types. 
#'
#' All the paths above can be NA or non existent, in which case a dummy row is inserted into the table so that downstream
#' applications know that data is missing for that sample. 
#'
#' @param tumors keyed data.table i.e. keyed by unique tumor id with specific columns corresponding to  paths to pipeline outputs(see description)
#' @param gencode path to gencode with just a single entry for each gene (so gencode entries for each gene are collapse to a single range). The input could be .gtf or .rds with GRanges object, or a GRanges object i.e. resulting from importing the (appropriate) GENCODE .gtf via rtracklayer, note: this input is only used in CNA to gene mapping. If nothing is provided then 'http://mskilab.com/fishHook/hg19/gencode.v19.genes.gtf' is used by default.
#' @param amp.thresh SCNA amplification threshold to call an amp as a function of ploidy (4)
#' @param del.thresh SCNA deletion threshold for (het) del as a function of ploidy (by default cn = 1 will be called del, but this allows additoinal regions in high ploidy tumors to be considered het dels)
#' @param mc.cores number of cores for multithreading
#' @param verbose logical flag 
#' @author Marcin Imielinski
#' @export
oncotable = function(
  tumors,
  gencode = '~/DB/GENCODE/gencode.v29lift37.annotation.nochr.rds',
  verbose = TRUE,
  amp.thresh = 4,
  filter = 'PASS',
  del.thresh = 0.5,
  mc.cores = 1
) {
  gencode = process_gencode(gencode)

  if ('type' %in% names(mcols(gencode))){
      # This is a bit hacky. The hg38 object does not contain the "type" column so we check if it is there and only use it when it is present
      pge = gencode %Q% (type  == 'gene' & gene_type == 'protein_coding')
  } else {
      pge = gencode %Q% (gene_type == 'protein_coding')
  }

  .oncotable = function(
    dat,
    x = dat[[key(dat)]][1],
    pge,
    verbose = TRUE,
    amp.thresh = 2,
    del.thresh = 0.5,
    filter = 'PASS'
  ) {
    out = data.table()

    ## collect gene fusions
    if (!is.null(dat$fusions) && file.exists(dat[x, fusions]))
    {
      if (verbose)
        message('pulling $fusions for ', x)
      fus = readRDS(dat[x, fusions])$meta
      if (nrow(fus))
      {
        fus = fus[silent == FALSE, ][!duplicated(genes), ]
        fus[, vartype := ifelse(in.frame == TRUE, 'fusion', 'outframe_fusion')] # annotate out of frame fusions
        fus = fus[, .(
            gene = strsplit(genes, ',') %>% unlist,
            vartype = rep(vartype, sapply(strsplit(genes, ','), length)),
            fusion_genes = rep(genes, sapply(strsplit(genes, ','), length))
        )][, id := x][, track := 'variants'][, type := vartype][, source := 'fusions']
        # get coordinates for fusion genes
        fus[, fusion_gene_coords := unlist(lapply(strsplit(fusion_genes, ','), function(genes) {
          coords <- lapply(genes, function(gene) {
            gene_ranges <- pge[mcols(pge)$gene_name == gene]
            paste0(seqnames(gene_ranges), ":", start(gene_ranges), "-", end(gene_ranges))
          })
          paste(unlist(coords), collapse = ",")
        }))]
        out = rbind(out, fus, fill = TRUE, use.names = TRUE)
      }
    } 
    else ## signal missing result
      out = rbind(out, data.table(id = x, type = NA, source = 'fusions'), fill = TRUE, use.names = TRUE)

    ## collect complex events
    if (!is.null(dat$complex) && file.exists(dat[x, complex]))
    {
      if (verbose)
        message('pulling $complex events for ', x)
      sv = readRDS(dat[x, complex])$meta$events
      if (nrow(sv))
      {
        sv = sv[, .(value = .N), by = type][, id := x][, track := ifelse(type %in% c('del', 'dup', 'invdup', 'tra', 'inv'), 'simple sv', 'complex sv')][, source := 'complex']
        out = rbind(out, sv, fill = TRUE, use.names = TRUE)
      }
    }
    else
      out = rbind(out, data.table(id = x, type = NA, source = 'complex'), fill = TRUE, use.names = TRUE)

    ## collect copy number / jabba
    if (!is.null(dat$jabba_rds) && file.exists(dat[x, jabba_rds]))
    {
      if (verbose)
        message('pulling $jabba_rds to get SCNA and purity / ploidy for ', x)
      jab = readRDS(dat[x, jabba_rds])
      jabpurity = NULL
      jabploidy = NULL
      jabpurity = if (is.null(jab$meta$purity)) { jab$meta$purity } else { jab$purity }
      jabploidy = if (is.null(jab$meta$ploidy)) { jab$meta$ploidy } else { jab$ploidy }
      out = rbind(out,
                  data.table(id = x, value = c(jabpurity, jabploidy), type = c('purity', 'ploidy'), track = 'pp'),
                  fill = TRUE, use.names = TRUE)

      # get the ncn data from jabba
      nseg = NULL
      
      ## Don't fail out if karyograph isn't found 
      ## Just provide nseg is NULL, and get_gene_ampdels_from_jabba
      ## assumes ncn = 2 (pretty safe assumption)
      if (!is.null(dat$karyograph) && file.exists(dat[x, karyograph])) {
        nseg = readRDS(dat$karyograph)$segstats[,c("ncn")]
      }

      scna = skitools::get_gene_ampdels_from_jabba(jab, amp.thresh = amp.thresh,
                                     del.thresh = del.thresh, pge = pge, nseg = nseg)

      if (nrow(scna))
      {
        scna[, track := 'variants'][, source := 'jabba_rds'][, vartype := 'scna']
        out = rbind(out,
                    scna[, .(id = x, value = min_cn, type, track, gene = gene_name)],
                    fill = TRUE, use.names = TRUE)
      }
    } else {
      out = rbind(out, data.table(id = x, type = NA, source = 'jabba_rds'), fill = TRUE, use.names = TRUE)
    }

    ## collect signatures
    if (!is.null(dat$signature_counts) && file.exists(dat[x, signature_counts]))
    {
      if (verbose)
        message('pulling $signature_counts for ', x)
      sig = fread(dat[x, signature_counts])
      sig = sig[, .(id = x, value = num_events, type = Signature, etiology = Etiology, frac = frac.events, track = 'signature', source = 'signature_counts')]
      out = rbind(out, sig, fill = TRUE, use.names = TRUE)
    }
    else
      out = rbind(out, data.table(id = x, type = NA, source = 'signature_counts'), fill = TRUE, use.names = TRUE)

    ## collect gene mutations
    if (!is.null(dat$annotated_bcf) && file.exists(dat[x, annotated_bcf]))
    {
      if (verbose)
        message('pulling $annotated_bcf for ', x, ' using FILTER=', filter)
      local_bcftools_path <- Sys.which("bcftools")
      local_bcftools_path <- ifelse(local_bcftools_path == "", stop("bcftools not found in the system PATH. Please install or moudule load bcftools."), local_bcftools_path)
      message("bcftools found at: ", local_bcftools_path)
      bcf = skitools::grok_bcf(dat[x, annotated_bcf], label = x, long = TRUE, filter = filter, bpath=local_bcftools_path)
      if (verbose)
        message(length(bcf), ' variants pass filter')
      genome.size = sum(seqlengths(bcf), na.rm = TRUE)/1e6
      if (is.na(genome.size)) ## something went wrong with vcf
        genome.size = sum(seqlengths(gG(jabba = dat[x, jabba_rds])), na.rm = TRUE)/1e6
      nmut = data.table(as.character(seqnames(bcf)), start(bcf), end(bcf), bcf$REF, bcf$ALT) %>% unique %>% nrow
      mut.density = data.table(id = x, value = c(nmut, nmut/genome.size), type = c('count', 'density'),  track = 'tmb', source = 'annotated_bcf')
      out = rbind(out, mut.density, fill = TRUE, use.names = TRUE)
      keepeff = c('trunc', 'cnadel', 'cnadup', 'complexsv', 'splice', 'inframe_indel', 'fusion', 'missense', 'promoter', 'regulatory','mir')
      bcf = bcf[bcf$short %in% keepeff]
      if (verbose)
        message(length(bcf), ' variants pass keepeff')
      vars = NULL
      if (length(bcf))
      {
        bcf$variant.g = paste0(seqnames(bcf), ':', start(bcf), '-', end(bcf), ' ', bcf$REF, '>', bcf$ALT)
        vars = gr2dt(bcf)[, .(id = x, gene, vartype, variant.g, variant.p, distance, annotation, type = short, track = 'variants', source = 'annotated_bcf')]
        setkey(vars, variant.g)
        vars = vars[, .SD[1], by = variant.g]
      }
      out = rbind(out, vars, fill = TRUE, use.names = TRUE)
    }
    else
      out = rbind(out, data.table(id = x, type = NA, source = 'annotated_bcf'), fill = TRUE, use.names = TRUE)

    if (verbose)
      message('done ', x)

    return(out)
  }

  if (is.null(key(tumors)))
  {
    if (is.null(tumors$id))
      stop('Input tumors table must be keyed or have column $id')
    else
      setkey(tumors, id)
  }

  out = mclapply(tumors[[key(tumors)]], .oncotable,
                 dat = tumors, pge = pge, amp.thresh = amp.thresh, filter = filter, del.thresh = del.thresh, verbose = verbose, mc.cores = mc.cores)
  out = rbindlist(out, fill = TRUE, use.names = TRUE)

  setnames(out, 'id', key(tumors))
  return(out)
}

