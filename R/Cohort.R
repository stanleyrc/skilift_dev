#' @import R6
#' @import data.table
#' @export
Cohort <- R6Class("Cohort",
  public = list(
    #' @field inputs data.table containing cohort information
    inputs = NULL,
    
    #' @field reference_name character string specifying genome reference, defaults to hg19
    reference_name = "hg19",

    #' @field cohort_cols_to_x_cols mapping of cohort columns to possible input columns
    cohort_cols_to_x_cols = NULL,

    #' @description
    #' Initialize a new Cohort object
    #' @param x Either a data.table or path to pipeline output directory
    #' @param reference_name character string specifying genome reference
    #' @param col_mapping Optional list mapping cohort columns to possible input columns
    initialize = function(x, reference_name = "hg19", col_mapping = NULL) {
      self$reference_name <- reference_name
      
      default_col_mapping <- list(
        pair = c("pair", "patient_id", "pair_id", "sample"),
        tumor_type = c("tumor_type"),
        disease = c("disease"),
        primary_site = c("primary_site"),
        inferred_sex = c("inferred_sex"),
        structural_variants = c("structural_variants", "gridss_somatic", "gridss_sv", "svaba_sv", "sv", "svs"),
        tumor_coverage = c("tumor_coverage", "dryclean_tumor", "tumor_dryclean_cov"),
        somatic_snvs = c("somatic_snvs", "sage_somatic_vcf", "strelka_somatic_vcf", "strelka2_somatic_vcf", "somatic_snv", "snv_vcf", "somatic_snv_vcf"),
        germline_snvs = c("germline_snvs", "sage_germline_vcf", "germline_snv", "germline_snv_vcf"),
        het_pileups = c("het_pileups", "hets", "sites_txt", "hets_sites"),
        somatic_snv_cn = c("somatic_snv_cn", "multiplicity"),
        germline_snv_cn = c("germline_snv_cn", "germline_multiplicity"),
        somatic_variant_annotations = c("somatic_variant_annotations", "annotated_bcf"),
        germline_variant_annotations = c("germline_variant_annotations", "annotated_vcf_germline"),
        oncokb_snv = c("oncokb_snv", "oncokb_maf", "maf"),
        oncokb_cna = c("oncokb_cna", "cna"),
        oncokb_fusions = c("oncokb_fusions", "oncokb_fusion", "fusion_maf"),
        jabba_gg = c("jabba_gg", "jabba_simple", "jabba_rds", "jabba_simple_gg"),
        karyograph = c("karyograph"),
        balanced_jabba_gg = c("balanced_jabba_gg", "non_integer_balance", "balanced_gg"),
        events = c("events", "complex"),
        fusions = c("fusions"),
        allelic_jabba_gg = c("allelic_jabba_gg", "lp_phased_balance", "allelic_gg"),
        activities_sbs_signatures = c("activities_sbs_signatures", "sbs_activities"),
        matrix_sbs_signatures = c("matrix_sbs_signatures", "sbs_matrix"),
        decomposed_sbs_signatures = c("decomposed_sbs_signatures", "sbs_decomposed"),
        activities_indel_signatures = c("activities_indel_signatures", "indel_activities"),
        matrix_indel_signatures = c("matrix_indel_signatures", "indel_matrix"),
        decomposed_indel_signatures = c("decomposed_indel_signatures", "indel_decomposed"),
        hrdetect = c("hrdetect", "hrd"),
        oncotable = c("oncotable"),
        estimate_library_complexity = c("estimate_library_complexity", "library_complexity_metrics", "est_lib_complex"),
        alignment_summary_metrics = c("alignment_summary_metrics", "alignment_metrics"),
        insert_size_metrics = c("insert_size_metrics", "insert_metrics"),
        wgs_metrics = c("wgs_metrics", "wgs_stats")
      )
      
      # Merge user-provided mapping with default mapping
      if (!is.null(col_mapping)) {
        for (col_name in names(col_mapping)) {
          if (col_name %in% names(default_col_mapping)) {
            # col_mapping goes first to override default_mapping
            default_col_mapping[[col_name]] <- unique(c(col_mapping[[col_name]], default_col_mapping[[col_name]]))
          } else {
            default_col_mapping[[col_name]] <- col_mapping[[col_name]]
          }
        }
      }
      
      self$cohort_cols_to_x_cols <- default_col_mapping
      
      if (is.character(x) && length(x) == 1) {
        self$inputs <- private$construct_from_path(x)
      } else if (is.data.table(x)) {
        self$inputs <- private$construct_from_datatable(x)
      } else {
        stop("Input must be either a path (character) or data.table")
      }
    },

    #' Print method for cohort object
    print = function(...) {
      cat("Cohort methods:\n")
      cat(format(self, ...), sep = "\n")
      cat("\n")
      base::print(self$inputs)
    }
    ,

    #' @description
    #' Validate inputs data.table for missing values and files
    #' @return NULL if all inputs are valid, or a data.table with details about missing data
    validate_inputs = function() {
      if (is.null(self$inputs) || nrow(self$inputs) == 0) {
        stop("No inputs data available to validate")
      }
      
      # Initialize results data.table
      missing_data <- data.table(
        pair = character(),
        column = character(),
        reason = character()
      )
      
      # Define metadata fields that should only check for NA/NULL
      metadata_fields <- c("tumor_type", "disease", "primary_site", "inferred_sex")

      # Check each column for missing values
      for (col in names(self$inputs)) {
        # Skip pair column as it's our identifier
        if (col == "pair") next
        
        # Check for NULL or NA values
        na_pairs <- self$inputs[is.na(get(col)) | is.null(get(col)), pair]
        if (length(na_pairs) > 0) {
          missing_data <- rbindlist(list(
            missing_data,
            data.table(
              pair = na_pairs,
              column = col,
              reason = "NULL or NA value"
            )
          ))
        }
        
        # For non-metadata columns that should contain file paths, check if files exist
        if (col %in% names(self$cohort_cols_to_x_cols) && !(col %in% metadata_fields)) {
          file_paths <- self$inputs[!is.na(get(col)), get(col)]
          for (path in file_paths) {
            if (!file.exists(path)) {
              missing_data <- rbindlist(list(
                missing_data,
                data.table(
                  pair = self$inputs[get(col) == path, pair],
                  column = col,
                  reason = "File does not exist"
                )
              ))
            }
          }
        }
      }
      
      # Return results
      if (nrow(missing_data) == 0) {
        message("All inputs are valid - no missing values or files found")
        return(NULL)
      } else {
        # Sort by pair and column for better readability
        setorder(missing_data, pair, column)
        return(missing_data)
      }
    }
  ),
  
  private = list(
    construct_from_path = function(pipeline_outdir) {
      if (!dir.exists(pipeline_outdir)) {
        stop("Pipeline directory does not exist: ", pipeline_outdir)
      }
      
      # Get sample metadata first - this contains our patient IDs
      sample_metadata <- private$get_pipeline_samples_metadata(pipeline_outdir)

      if (is.null(sample_metadata)) {
        stop("Could not get sample metadata - this is required for patient IDs")
      }
      
      # Initialize data.table with all pairs from metadata
      ## outputs <- data.table(pair = sample_metadata$pair)

      
      # Get all file paths recursively
      pipeline_output_paths <- list.files(pipeline_outdir, recursive = TRUE, full.names = TRUE)
      
      outputs <- parse_pipeline_paths(
        pipeline_output_paths,
        initial_dt = sample_metadata
      )
      
      if (nrow(outputs) == 0) {
        warning("No data could be extracted from pipeline directory")
      }
      
      return(outputs)
    },
    
    get_pipeline_samples_metadata = function(pipeline_outdir) {
      report_path <- file.path(pipeline_outdir, "pipeline_info/pipeline_report.txt")
      if (!file.exists(report_path)) {
        warning("Pipeline report not found: ", report_path)
        return(NULL)
      }
      nflogs = list.files(pipeline_outdir, pattern = "\\.nextflow\\.log.*", all.files = TRUE, full.names = TRUE)
      parsed_meta = lapply(nflogs, private$read_nf_log)
      ## transpose the list
      parsed_meta = do.call(Map, c(f = c, parsed_meta))
      samplesheets = unique(parsed_meta$samplesheet)
      samplesheets = samplesheets[nzchar(samplesheets)]
      launchdir = unique(parsed_meta$launchdir)
      launchdir = launchdir[nzchar(launchdir)]

      unique_patients = Reduce(union, lapply(samplesheets, function(x) fread(x)$patient))
      samplesheet = data.table::rbindlist(lapply(samplesheets, fread), fill = TRUE)

      metadata <- data.table(
        pair = samplesheet$patient,
        tumor_type = samplesheet$tumor_type,
        disease = samplesheet$disease,
        primary_site = samplesheet$primary_site,
        inferred_sex = samplesheet$sex
      )
      
      # remove duplicates
      metadata <- unique(metadata)
      return(metadata)
      
      browser()

      # Read pipeline report
      report_lines <- readLines(report_path)

      # Extract launch directory and samplesheet path
      launch_dir <- grep("launchDir:", report_lines, value = TRUE)
      samplesheet_path <- grep("input:", report_lines, value = TRUE)
      
      if (length(launch_dir) == 0 || length(samplesheet_path) == 0) {
        warning("Could not find launch directory or samplesheet path in pipeline report")
        return(NULL)
      }

      # Clean up paths
      browser()
      launch_dir <- gsub(".*launchDir: ", "", launch_dir)
      samplesheet_filename <- gsub(".*input: ", "", samplesheet_path)
      samplesheet_filename <- gsub("^\\./", "", samplesheet_filename)
      if(grepl("^/", samplesheet_filename)) { #TRUE if the samplesheet is already a full path
          samplesheet_path <- samplesheet_filename
      } else {
          samplesheet_path <- file.path(launch_dir, gsub("^\\./", "", samplesheet_filename))
      }
      
      if (!file.exists(samplesheet_path)) {
        warning("Samplesheet not found: ", samplesheet_path)
        return(NULL)
      }
      
      # Read samplesheet and extract metadata
      samplesheet <- fread(samplesheet_path)
      metadata <- data.table(
        pair = samplesheet$patient,
        tumor_type = samplesheet$tumor_type,
        disease = samplesheet$disease,
        primary_site = samplesheet$primary_site,
        inferred_sex = samplesheet$sex
      )
      
      # remove duplicates
      metadata <- unique(metadata)

      return(metadata)
    },
    
    # Determine column based on file path
    construct_from_datatable = function(dt) {
      if (!is.data.table(dt)) {
        stop("Input must be a data.table")
      }
      
      result_dt <- data.table()
      
      for (cohort_col in names(self$cohort_cols_to_x_cols)) {
        possible_cols <- self$cohort_cols_to_x_cols[[cohort_col]]
        
        found_col <- NULL
        # First try exact matches
        for (col in possible_cols) {
          if (col %in% names(dt)) {
            found_col <- col
            break
          }
        }
        
        # If no exact match found, try prefix matches
        if (is.null(found_col)) {
          for (col in possible_cols) {
            for (name in names(dt)) {
              if (startsWith(name, col)) {
                found_col <- name
                break
              }
            }
            if (!is.null(found_col)) break
          }
        }
        
        if (!is.null(found_col)) {
          result_dt[, (cohort_col) := dt[[found_col]]]
        } else {
          warning(sprintf("No matching column found for '%s'. Expected one of: %s", 
                        cohort_col, paste(possible_cols, collapse = ", ")))
        }
      }
      
      if (nrow(result_dt) == 0) {
        warning("No data could be extracted from input data.table")
      }
      
      if ("pair" %in% names(result_dt) && inherits(result_dt, "data.table")) {
        data.table::setkey(result_dt, pair)
      }
      
      return(result_dt)
    },
    
    filter_inputs = function() {
      stop("Method not implemented yet")
    }
  )
)

#' Expected nf-casereport patterns
#' 
#' Encoding the expected path variables
#' 
#' This variable contains the paths to be expected from
#' the nf-casereport directory structure. 
#' You can output the default template to modify it using
#' base::dput(Skilift::nf_path_patterns)
#' @export
nf_path_patterns <- list(
  balanced_jabba_gg = "non_integer_balance/.*/non_integer.balanced.gg.rds$",
  tumor_coverage = "dryclean_tumor/.*/drycleaned.cov.rds$",
  het_pileups = "(hetpileups|amber)/.*/sites.txt$",
  jabba_gg = "jabba/.*/jabba.simple.gg.rds$",
  events = "events/.*/complex.rds$",
  fusions = "fusions/.*/fusions.rds$",
  structural_variants = c("gridss.*/.*/.*high_confidence_somatic.vcf.bgz$", "tumor_only_junction_filter/.*/.*somatic.filtered.sv.rds$", "gridss.*/.*.gridss.filtered.vcf.gz$"),
  karyograph = "jabba/.*/karyograph.rds$",
  allelic_jabba_gg = "lp_phased_balance/.*/lp_phased.balanced.gg.rds$",
  somatic_snvs = c("sage/somatic/tumor_only_filter/.*/.*.sage.pass_filtered.tumoronly.vcf.gz$", "sage/somatic/.*/.*sage.somatic.vcf.gz$"),
  somatic_variant_annotations = "snpeff/somatic/.*/.*ann.bcf$",
  somatic_snv_cn = "snv_multiplicity3/.*/.*est_snv_cn_somatic.rds",
  activities_sbs_signatures = "signatures/sigprofilerassignment/somatic/.*/sbs_results/Assignment_Solution/Activities/sbs_Assignment_Solution_Activities.txt",
  matrix_sbs_signatures = "signatures/sigprofilerassignment/somatic/.*/SBS/sigmat_results.SBS96.all",
  decomposed_sbs_signatures = "signatures/sigprofilerassignment/somatic/.*/sbs_results/Assignment_Solution/Activities/Decomposed_MutationType_Probabilities*.txt",
  activities_indel_signatures = "signatures/sigprofilerassignment/somatic/.*/indel_results/Assignment_Solution/Activities/indel_Assignment_Solution_Activities.txt",
  matrix_indel_signatures = "signatures/sigprofilerassignment/somatic/.*/ID/sigmat_results.ID83.all",
  decomposed_indel_signatures = "signatures/sigprofilerassignment/somatic/.*/indel_results/.*/Decomposed_MutationType_Probabilities.txt",
  hrdetect = "hrdetect/.*/hrdetect_results.rds",
  estimate_library_complexity = "qc_reports/gatk/.*/.*metrics",
  alignment_summary_metrics = "qc_reports/picard/.*/.*alignment_summary_metrics",
  insert_size_metrics = "qc_reports/picard/.*/.*insert_size_metrics",
  wgs_metrics = "qc_reports/picard/.*/.*coverage_metrics"
)

#' nf-casereport Parser
#' 
#' Parses outputs of the nf-casereport directory structure
#' 
#' This is a function that extracts paths from a 
#' pipeline. The default structure is that the
#' path patterns expected follow that of nf-casereport run.
#' The columns it extracts are defined by
#' variable Skilift::nf_column_map.
#' You can output the default template to modify it using
#' base::dput(Skilift::nf_path_patterns)
#' @export
parse_pipeline_paths = function(
  paths,
  initial_dt = data.table(pair = character(0)),
  path_patterns = Skilift::nf_path_patterns,
  id_name = "pair"
) {
  if (length(paths) == 1 && dir.exists(paths)) {
    paths = list.files(
      paths, 
      recursive = TRUE,
      full.names = TRUE
    )
  }
  
  # Get the list of pairs we're looking for
  pairs_to_match <- unique(initial_dt[[id_name]])
  # Map of regex patterns to column names
  for (col_name in names(path_patterns)) {
    patterns = path_patterns[[col_name]]
    for (pattern in patterns) {
      present_paths = grep(pattern, paths, value = TRUE)
      if (!length(present_paths) > 0) next
      
      # Create mapping between paths and pairs
      dt <- data.table()
      for (pair in pairs_to_match) {
        # Find paths that contain this pair name
        pair_paths <- grep(paste0("/", pair, "/"), present_paths, value = TRUE)
        if (length(pair_paths) > 0) {
          dt <- rbindlist(list(dt, data.table(
            pair = pair,
            path = pair_paths
          )))
        }
      }
      
      if (nrow(dt) > 0) {
        setnames(dt, "path", col_name)
        initial_dt = merge.data.table(
          initial_dt,
          dt,
          by = "pair",
          all.x = TRUE,
          all.y = FALSE,  # Only keep pairs that were in initial_dt
          allow.cartesian = FALSE
        )
        break
      }
    }
  }
  
  return(initial_dt)
}


#' deparse1 copy
#'
#' Robustly deparse object
#'
deparse1 = function(expr, collapse = " ", width.cutoff = 500L, ...) {
    paste(deparse(expr, width.cutoff, ...), collapse = collapse)
}

#' Subset Cohort object
#'
#' Overloads subset operator for Cohort
#'
#' @export 
'[.Cohort' = function(obj, i = NULL, j = NULL, with = TRUE, ...) {
  expri = deparse1(substitute(i))
  is_i_given = any(expri != "NULL")
  vector_of_column_names = deparse1(substitute(j))
  is_j_given = any(vector_of_column_names != "NULL")
  tbl = data.table::copy(obj$inputs)
  tblj = tbl
  if (is_j_given) {
    selectj = unique(c("pair", vector_of_column_names))
    tblj = base::subset(tblj, select = selectj) 
  }
  colmap = as.list(names(tblj))
  names(colmap) = names(tblj)
  tbli = tblj
  if (is_i_given) {
    tbli = tblj[i,,with = with]
  }
  # obj_out$inputs = tbli
  obj_out = Skilift::Cohort$new(x = data.table(), reference_name = obj$reference_name)
  obj_out$inputs = tbli[]
  invisible(obj_out$inputs[])
  return(obj_out)
}


#' Parse .nextflow.log
#'
#' Parse nextflow log files for launchdir and input options
#' to get samplesheet
#'
#' @author Kevin Hadi 
read_nf_log = function(nflog_path, max_lines = 1000000L) {
    f = file(nflog_path, open = "r")
    on.exit({close(f)})
    line = readLines(f, n = 1L)
    is_line_at_options = grepl("Input/output options", line)
    is_line_at_launchdir = grepl("launchDir", line)
    are_all_lines_parsed = is_line_at_options && is_line_at_launchdir
    counter = 1
    ansi_pattern = '\\\033\\[((?:\\d|;)*)([a-zA-Z])'
    input_samplesheet = ""
    launchdir = ""
    while (!are_all_lines_parsed && counter <= max_lines) {
        line = readLines(f, n = 1L)
        if (!identical(is_line_at_options, TRUE)) {
            is_line_at_options = length(line) && grepl("Input/output options", line)
            if (identical(is_line_at_options, TRUE)) {
                line = readLines(f, n = 1L)
                input_samplesheet = gsub(ansi_pattern, "", line)
                input_samplesheet = gsub("input[[:space:]]+:[[:space:]]+", "", trimws(input_samplesheet))
            }
        }
        if (!identical(is_line_at_launchdir, TRUE)) {
            is_line_at_launchdir = length(line) && grepl("launchDir", line)
            if (identical(is_line_at_launchdir, TRUE)) {
                launchdir_line = gsub(ansi_pattern, "", line)
                launchdir = gsub("launchDir[[:space:]]+:[[:space:]]+", "", trimws(launchdir_line))
            }
        } 
        are_all_lines_parsed = is_line_at_options && is_line_at_launchdir
        counter = counter + 1
    }
    if (counter == max_lines && !are_all_lines_parsed) {
        stop("Could not find Input/output options!")
    }
    return(
        list(
            samplesheet = input_samplesheet,
            launchdir = launchdir
        )
    )
}


## Assign new methods
Cohort$private_methods[["read_nf_log"]] = read_nf_log
