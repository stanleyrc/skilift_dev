#' @description
#' Creates a multiplicity data table from SNV copy number data, handling both somatic and germline cases.
#'
#' @param snv_cn Path to SNV copy number RDS file or data.table object
#' @param is_germline Logical indicating if mutations are germline (default: FALSE)
#' @param field Column name to use for copy number values (default: "total_copies")
#'
#' @return data.table containing processed mutation data
#' @export
create_multiplicity <- function(snv_cn, oncokb_snv=NULL, is_germline = FALSE, field = "altered_copies") {
  mutations.gr <- tryCatch({
    if (!grepl("\\.rds$", snv_cn)) {
      message("Expected .rds ending for mutations. Attempting to read anyway: ", snv_cn)
    }
    readRDS(snv_cn)
  }, error = function(e) {
    message(paste0("Input was not .rds; failed with '", e$message, "'\nAssuming input is .maf"))
    return(fread(snv_cn) %>% dt2gr())
  }, finally = {
    message("Finished attempting to load input.")
  })

  if (!inherits(mutations.gr, "GRanges")) {
    mutations.gr.tmp = tryCatch(
      {
        dt2gr(mutations.gr)
      }, error = function(e) tryCatch(
      {
        as(mutations.gr, "GRanges")
      }, error = function(e) NULL)
    )
    if (is.null(mutations.gr.tmp)) {
      stop("input snv_cn must be a maf/snv file coercible into GRanges")
    }
    mutations.gr = mutations.gr.tmp
    rm(mutations.gr.tmp)
  }

  if (is.null(mutations.gr)) {
    stop("Failed to assign a valid value to mutations.gr.")
  } else {
    message("Successfully loaded input snv_cn.")
  }

  
  mcols(mutations.gr)$snpeff_annotation = mutations.gr$annotation
  annotationsplit = strsplit(mcols(mutations.gr)$snpeff_annotation, "&")
  annotationsplit = gGnome::dunlist(annotationsplit)
  annotationsplit[, ix := seq_len(.N), by = listid]
  annotationsplit[, num := .N, by = listid]

  ## Normalize everything to just the 1st variant
  ## type that appears if we get
  ## splice&intron_variant nonsense.
  mcols(mutations.gr)$snpeff_annotation = annotationsplit[ix == 1]$V1
  rm("annotationsplit")

  mut_ann <- ""

  # Annotation fields that should always be present
  # We're using 
  annotation_fields <- list(
    # Variant_Classification = "Type"
    # Let's just use the snpeff annotations here because
    # it should still be possible to populate with just SnpEff/Multiplicity.
    snpeff_annotation = "Type",
    Gene = "Gene",
    # HGVSc = "Variant",
    variant.c = "Variant",
    # HGVSp = "Protein_variant",
    variant.p = "Protein_variant",
    variant.g = "Genomic_variant",
    vaf = "VAF",
    alt = "Alt_count",
    ref = "Ref_count",
    normal.alt = "Normal_alt_count",
    normal.ref = "Normal_ref_count",
    FILTER = "Filter"
  )
  
  mutations.dt <- gr2dt(mutations.gr)
  if (!is.null(oncokb_snv)) {
    message("oncokb_snv provided, processing input")
    is_path_character = is.character(oncokb_snv)
    is_length_one = NROW(oncokb_snv) == 1
    is_file_exists = is_path_character && is_length_one && file.exists(oncokb_snv)
    is_rds = is_file_exists && grepl("rds$", oncokb_snv)
    is_txt = is_file_exists && grepl("maf$|(c|t)sv$|txt$", oncokb_snv)
    if (is_rds) {
      oncokb_snv = readRDS(oncokb_snv)
    } else if (is_txt) {
      oncokb_snv = fread(oncokb_snv)
    } else {
      stop("Provided oncokb_snv file extension not supported")
    }

    is_data_frame = inherits(oncokb_snv, "data.frame")
    if (is_data_frame) {
      oncokb_snv_tmp = tryCatch(
        {
          dt2gr(oncokb_snv)
        }, error = function(e) tryCatch(
        {
          as(oncokb_snv, "GRanges")
        }, error = function(e) NULL)
      )
      if (is.null(oncokb_snv_tmp)) {
        stop("oncokb_snv must be coercible to GRanges")
      }
      oncokb_snv = oncokb_snv_tmp
    }
    if (!inherits(oncokb_snv, "GRanges")) {
      stop("final oncokb_snv not a GRanges object")
    }

    mutations.gr.annotated = merge_oncokb_multiplicity(
      oncokb_snv,
      mutations.gr,
      overwrite = TRUE,
      other.cols.keep = c("snpeff_annotation")
    )
    mutations.gr.annotated$gene = mutations.gr.annotated$Hugo_Symbol
    mutations.dt = gr2dt(mutations.gr.annotated)

    ## Process mutations
    setnames(mutations.dt, old = "VAF", new = "vaf", skip_absent = TRUE)
    mutations.dt <- mutations.dt[!is.na(get(field)), ]
    mutations.dt[start == end, end := end + 1]
    mutations.dt[, vaf := round(vaf, 3)] ## round for frontend legibility
    mutations.dt[, ONCOGENIC := fcase(
      is.na(ONCOGENIC), "",
      grepl("Unknown", ONCOGENIC), "", ## necessitated by frontend implementation
      default = ONCOGENIC
    )]
    mutations.dt[, MUTATION_EFFECT := fcase(
      is.na(MUTATION_EFFECT), "",
      grepl("Unknown", MUTATION_EFFECT), "", ## extraneous string
      default = MUTATION_EFFECT
    )]
    mutations.dt[, HIGHEST_LEVEL := fcase(
      HIGHEST_LEVEL == "" | is.na(HIGHEST_LEVEL), "",
      default = gsub("LEVEL_", "", HIGHEST_LEVEL) ## extraneous string
    )]

    mutations.dt <- mutations.dt[FILTER == "PASS"] #### TEMPORARY before implementation of fast coverage

    if ("strand" %in% colnames(mutations.dt)) {
      mutations.dt[, strand := NULL]
    }

    ## Add OncoKB specific annotation fields
    annotation_fields = c(
      annotation_fields,
      list(
        ONCOGENIC = "Oncogenicity",
        MUTATION_EFFECT = "Effect",
        HIGHEST_LEVEL = "Level"
      )
    )

    # Converting to the OncoKB HGVSc and p variants.
    # for internal consistency
    names(annotation_fields)[
      names(annotation_fields) == "variant.c"
    ] = "HGVSc"

    names(annotation_fields)[
      names(annotation_fields) == "variant.p"
    ] = "HGVSp"

    
  } 

  if (!any(class(mutations.dt) == "data.table")) {
    stop("Input must be a data.table.")
  }

  ## parse mut_ann variable
  for (col in names(annotation_fields)) {
    if (col %in% colnames(mutations.dt)) {
      mut_ann <- paste0(mut_ann, annotation_fields[[col]], ": ", mutations.dt[[col]], "; ")
    }
  }

  mutations.dt[, annotation := mut_ann]

  return(mutations.dt)
}

#' @title Convert Multiplicity Data to Intervals
#' @description
#' Converts multiplicity data into a list of intervals and settings for visualization
#'
#' @param multiplicity data.table containing mutation data with seqnames, start, and end columns
#' @param field Column name to use for y-axis values (default: "total_copies")
#' @param settings Path to settings JSON file
#' @param node_metadata Additional columns to include in node data (default: NULL)
#' @param reference_name Reference genome name (default: "hg19")
#'
#' @return List containing settings and intervals for visualization
#' @export
multiplicity_to_intervals <- function(
    multiplicity,
    field = "altered_copies",
    settings = Skilift:::default_settings_path,
    node_metadata = NULL,
    reference_name = "hg19",
    cohort_type
) {
    # Load chromosome lengths from settings
    settings_data <- jsonlite::fromJSON(settings)
    chrom_lengths <- as.data.table(
        settings_data$coordinates$sets[[reference_name]])[
        , .(chromosome, startPoint, endPoint)
    ]
    setnames(chrom_lengths, c("seqnames", "start", "end"))

  # Ensure chromosome naming consistency
  if (nrow(chrom_lengths[grepl("chr", seqnames), ]) > 0) {
    chrom_lengths[
      !grepl("chr", seqnames),
      seqnames := paste0("chr", seqnames)
    ]
  }

  # Set y-values from specified field
  multiplicity[, y_value := get(field)]

    # Convert to GRanges
    gr = dt2gr(multiplicity[order(seqnames, start), ]) %>% sortSeqlevels()
    if (nrow(chrom_lengths[grepl("chr", seqnames), ]) > 0) {
        GenomeInfoDb::seqlevelsStyle(gr) = "UCSC"
    } else {
        GenomeInfoDb::seqlevelsStyle(gr) = "NCBI"
    }

    if (is.null(cohort_type)) stop("Cohort type is missing")

    if (cohort_type == "heme") {
        hemedb = readRDS(Skilift:::HEMEDB)
        # HEMEDB = "/gpfs/data/imielinskilab/projects/Clinical_NYU/db/master_heme_database.20250128_095937.790322.rds"
        # FIXME: HARDCODED PATH!
        gencode = "/gpfs/data/imielinskilab/DB/GENCODE/gencode.v19.annotation.gtf.nochr.rds"
        gencode <- Skilift:::process_gencode(gencode)
        genes = gencode[gencode$type == "gene"]
        gr_heme_genes = genes[na.omit(match(hemedb$GENE, genes$gene_name))]
        message("Filtering multiplicity to heme relevant genes")
        is_heme = (gr %^% gr_heme_genes) & (gr$gene %in% gr_heme_genes$gene_name)
        gr_heme = gr[is_heme]
        gr_other = gr[!is_heme]
        remaining = 1e4 - NROW(gr_heme)
        if (remaining > 0) {
            otherix = 1:NROW(gr_other)
            set.seed(42)
            sampled_otherix = sample(otherix, size = remaining, replace = FALSE)
            gr_subsampled_other = gr_other[sampled_otherix]
            gr = c(gr_heme, gr_subsampled_other)
        } else {
            gr = gr_heme
        }
        
    }
    
    

  # Validate ranges
  if (any(gr@seqinfo@seqlengths >
    chrom_lengths[seqnames %in% names(seqlengths(gr))]$end)) {
    stop(paste("Ranges exceed chromosome lengths in", reference_name))
  }

  # Create graph object and convert to data.table
  jab <- gG(nodes = gr)
  node_cols <- c("snode.id", "y_value", "annotation") ### node_metadata)
  node_dt <- gr2dt(jab$nodes$gr[, node_cols])

  # Create final node data
  intervals <- node_dt[, .(
    chromosome = seqnames,
    startPoint = start,
    endPoint = end,
    iid = snode.id,
    title = snode.id,
    type = "interval",
    y = y_value,
    annotation = node_dt$annotation
  )]

  # Construct return list
  list(
    settings = list(
      y_axis = list(
        title = "copy number",
        visible = TRUE
      )
    ),
    intervals = intervals,
    connections = data.table()
  )
}


#' @name lift_multiplicity
#' @title lift_multiplicity
#' @description
#' Create multiplicity JSON files for all samples in a cohort
#'
#' @param cohort Cohort object containing sample information
#' @param is_germline Logical indicating if mutations are germline (default: FALSE)
#' @param output_data_dir Base directory for output files
#' @param cores Number of cores for parallel processing (default: 1)
#' @return None
#' @export
lift_multiplicity <- function(
    cohort,
    output_data_dir,
    is_germline = FALSE,
    node_metadata = c("gene", "feature_type", "annotation", "REF", "ALT", "variant.c", "variant.p", "vaf", "transcript_type", "impact", "rank"),
    field = "altered_copies",
    cores = 1
) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }
    
    if (!dir.exists(output_data_dir)) {
        dir.create(output_data_dir, recursive = TRUE)
    }
    
    # Determine which column to use based on is_germline
    snv_cn_col <- if(is_germline) "germline_multiplicity" else "multiplicity"

    oncokb_snv_col <- if(is_germline) NULL else "oncokb_snv"
    
    # Validate required column exists
    if (!snv_cn_col %in% names(cohort$inputs)) {
        stop(sprintf("Missing required column in cohort: %s", snv_cn_col))
    }
    
    # Get reference name from cohort
    reference_name <- cohort$reference_name
    if (is.null(reference_name)) {
        stop("Reference name not found in cohort object")
    }
    
    # Process each sample in parallel
    mclapply(seq_len(nrow(cohort$inputs)), function(i) {
        row <- cohort$inputs[i,]
        pair_dir <- file.path(output_data_dir, row$pair)
        
        if (!dir.exists(pair_dir)) {
            dir.create(pair_dir, recursive = TRUE)
        }
        
        # Determine output filename based on is_germline
        out_file <- file.path(
            pair_dir,
            if(is_germline) "germline_mutations.json" else "mutations.json"
        )
        
        snv_cn_path = row[[snv_cn_col]]
        oncokb_snv_path = NULL
        if (!is.null(oncokb_snv_col)) oncokb_snv_path = row[[oncokb_snv_col]]
        futile.logger::flog.threshold("ERROR")
        tryCatchLog({
            # Create multiplicity data.table
            mult_dt <- create_multiplicity(
                snv_cn = snv_cn_path,
                oncokb_snv = oncokb_snv_path,
                is_germline = is_germline,
                field = field
            )
            
            # Convert to intervals
            intervals_list <- multiplicity_to_intervals(
                multiplicity = mult_dt,
                reference_name = reference_name,
                node_metadata = node_metadata,
                field = field,
                cohort_type = cohort$type
            )
            
            # Write to JSON
            jsonlite::write_json(
                intervals_list,
                out_file,
                pretty = TRUE,
                auto_unbox = TRUE,
                digits = 4
            )
            
        }, error = function(e) {
            print(sprintf("Error processing %s: %s", row$pair, e$message))
            NULL
        })
    }, mc.cores = cores, mc.preschedule = TRUE)
    
    invisible(NULL)
}
