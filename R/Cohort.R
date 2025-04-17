default_settings_path <- system.file("extdata", "test_data", "settings.json", package = "Skilift")

#' @export
Cohort <- R6Class("Cohort",
  public = list(
    #' @field inputs data.table containing cohort information
    inputs = NULL,

    #' @field reference_name character string specifying genome reference, defaults to hg19
    reference_name = "hg19",

    #' @field cohort_cols_to_x_cols mapping of cohort columns to possible input columns
    cohort_cols_to_x_cols = NULL,

    #' @field unified list of cohort column names and nf-casereport map
    unified_map = NULL,

    #' @field unified list of cohort column names and nf-casereport map
    path_patterns = NULL,

    #' @field enum("paired", "heme", "tumor_only")
    type = NULL,

    #' @field character
    nextflow_results_path = NULL,

    #' @description
    #' Initialize a new Cohort object
    #' @param x Either a data.table or path to pipeline output directory
    #' @param reference_name character string specifying genome reference
    #' @param col_mapping Optional list mapping cohort columns to possible input columns
    initialize = function(x,
                          reference_name = "hg19",
                          settings = Skilift:::default_settings_path,
                          col_mapping = NULL,
                          path_patterns = Skilift::nf_path_patterns,
                          cohort_type = "paired") {
      self$reference_name <- reference_name

      default_cohort_types <- c("paired", "heme", "tumor_only")
      if (!cohort_type %in% default_cohort_types) {
        stop("cohort_type must be one of: ", paste(default_cohort_types, collapse = ", "))
      }
      self$type <- cohort_type

      # Merge user-provided mapping with default mapping
      default_col_mapping <- Skilift::default_col_mapping
      if (!is.null(col_mapping)) {
        for (col_name in names(col_mapping)) {
          col_mapping
          if (col_name %in% names(default_col_mapping)) {
            # col_mapping goes first to override default_mapping
            default_col_mapping[[col_name]] <- unique(
              c(
                col_mapping[[col_name]],
                default_col_mapping[[col_name]]
              )
            )
          } else {
            default_col_mapping[[col_name]] <- col_mapping[[col_name]]
          }
        }
      }

      self$path_patterns <- path_patterns

      self$cohort_cols_to_x_cols <- default_col_mapping


      if (is.character(x) && length(x) == 1) {
        self$inputs <- private$construct_from_path(x)[]
        self$nextflow_results_path <- x
      } else if (is.data.table(x)) {
        self$inputs <- private$construct_from_datatable(x)[]
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
    },
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

      # Define metadata fields and configuration parameter fields that should only check for NA/NULL
      metadata_fields <- c("tumor_type", "disease", "primary_site", "inferred_sex", "metadata_is_visible")
      config_fields <- Skilift:::config_parameter_names # Use the existing config parameter names

      # Combine metadata and config fields
      na_check_only_fields <- c(metadata_fields, config_fields)

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

        # For non-metadata/non-config columns that should contain file paths, check if files exist
        if (col %in% names(self$cohort_cols_to_x_cols) && !(col %in% na_check_only_fields)) {
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
      sample_metadata$pair_original = sample_metadata$pair

      is_castable = (
          all(!is.na(sample_metadata$sample_type))
      )
      is_paired = is_castable && all(c("normal", "tumor") %in% sample_metadata$sample_type)


      id_to_parse = "pair"
      if (is_castable) {
          sample_metadata = Skilift::dcastski(sample_metadata, id_columns = c("pair", "tumor_type", "disease", "primary_site", "inferred_sex", "pair_original"), type_columns = "sample_type", cast_columns = c("sample", "bam"), sep = "_")
      }

      if (is_paired) {
          sample_metadata$realpair = paste(sample_metadata$tumor_sample, "_vs_", sample_metadata$normal_sample, sep = "")
          id_to_parse = c("realpair", "tumor_sample", "normal_sample")
      }

      is_null_tumor_bam = is.null(sample_metadata$tumor_bam)
      is_all_na_tumor_bam = !is_null_tumor_bam && all(is.na(sample_metadata$tumor_bam))
      is_null_normal_bam = is.null(sample_metadata$normal_bam)
      is_all_na_normal_bam = !is_null_normal_bam && all(is.na(sample_metadata$normal_bam))
      if (is_all_na_tumor_bam) {
          sample_metadata$tumor_bam = NULL
      }
      if (is_all_na_normal_bam) {
          sample_metadata$normal_bam = NULL
      }
      
      

      if (is.null(sample_metadata)) {
        stop("Could not get sample metadata - this is required for patient IDs")
      }

      # Add default configuration parameters to sample_metadata
      for (cohort_col in names(self$cohort_cols_to_x_cols)) {
        possible_cols <- self$cohort_cols_to_x_cols[[cohort_col]]
        default_value <- attr(possible_cols, "default")

        if (!is.null(default_value)) {
          if (is.list(default_value) || length(default_value) > 1) {
            sample_metadata[, (cohort_col) := list(list(default_value))]
          } else {
            sample_metadata[, (cohort_col) := default_value]
          }
        }
      }

      # Get all file paths recursively
      pipeline_output_paths <- list.files(pipeline_outdir, recursive = TRUE, full.names = TRUE)

      outputs_lst = mclapply(id_to_parse, function(id_name) {
          sample_metadata$pair = sample_metadata[[id_name]]
          parse_pipeline_paths(
              pipeline_output_paths,
              initial_dt = sample_metadata,
              path_patterns = self$path_patterns,
              id_name = id_name
          )
      }, mc.cores = 2)

      nm = names(sample_metadata)

      remove_cols = nm[!nm %in% c("pair_original")]

      for (i in seq_len(NROW(outputs_lst))) {
          nm_outputs = names(outputs_lst[[i]])
          outputs_lst[[i]] = base::subset(outputs_lst[[i]], select = !nm_outputs %in% remove_cols)
      }

      outputs = outputs_lst[[1]]
      if (length(outputs_lst) > 1) outputs = Reduce(function(x,y) merge(x,y, by = "pair_original", all = TRUE), outputs_lst)

      # outputs = Reduce(function(x,y) merge(x,y, by = "pair_original", all = TRUE), outputs_lst)
      outputs = merge(sample_metadata, outputs, by = "pair_original", all = TRUE)

      ## outputs <- parse_pipeline_paths(
      ##     pipeline_output_paths,
      ##     initial_dt = sample_metadata,
      ##     path_patterns = self$path_patterns
      ## )

      if (nrow(outputs) == 0) {
        warning("No data could be extracted from pipeline directory")
      }

      return(outputs)
    },
    get_pipeline_samples_metadata = function(pipeline_outdir) {
      report_path <- file.path(pipeline_outdir, "pipeline_info/pipeline_report.txt")
      launch_dir = character(0)
      samplesheet_path = character(0)
      if (!file.exists(report_path)) {
        message("Pipeline report not found: ", report_path)
        message("Attempting to read from nextflow logs")
        nflogs = list.files(pipeline_outdir, pattern = "\\.nextflow.*\\.log.*", all.files = TRUE, full.names = TRUE)
        if (NROW(nflogs) == 0) {
          pipeline_outdir_upone = dirname(normalizePath(pipeline_outdir))
          nflogs = list.files(
            pipeline_outdir_upone, 
            pattern = "\\.nextflow.*\\.log.*", 
            all.files = TRUE, 
            full.names = TRUE
          )
        }
        parsed_meta = lapply(nflogs, private$read_nf_log)
        parsed_meta = do.call(Map, c(f = c, parsed_meta))
        launch_dir = unique(parsed_meta$launchdir)
        launch_dir = launch_dir[nzchar(launch_dir)]
        samplesheet_path = unique(parsed_meta$samplesheet)
        samplesheet_path = samplesheet_path[nzchar(samplesheet_path)]
        if (NROW(launch_dir) > 1) stop("More than one launch_dir found - something's wrong")
        if (NROW(samplesheet_path) > 1) stop("More than one samplesheet_path found in pipeline directory")
      } else {
        # Read pipeline report
        report_lines <- readLines(report_path)

        # Extract launch directory and samplesheet path
        launch_dir <- grep("launchDir:", report_lines, value = TRUE)
        samplesheet_path <- grep("input:", report_lines, value = TRUE)

      }

      #### Edge case: multiple runs (with different samplesheets) point to the same results
      #### this needs to be handled a bit more carefully as the nflogs have a cap
      #### solution: coerce unique results directory per run (with common results directory)
      ####  merge the samplesheets with the same parent
      
      # nflogs = list.files(pipeline_outdir, pattern = "\\.nextflow.*\\.log.*", all.files = TRUE, full.names = TRUE)
      # parsed_meta = lapply(nflogs, private$read_nf_log)
      # ## transpose the list
      # parsed_meta = do.call(Map, c(f = c, parsed_meta))
      # samplesheets = unique(parsed_meta$samplesheet)
      # samplesheets = samplesheets[nzchar(samplesheets)]
      # launchdir = unique(parsed_meta$launchdir)
      # launchdir = launchdir[nzchar(launchdir)]
      #
      # # unique_patients = Reduce(union, lapply(samplesheets, function(x) fread(x)$patient))
      #
      # samplesheets = samplesheets[file.exists(samplesheets)]
      # samplesheet = data.table::rbindlist(lapply(samplesheets, fread), fill = TRUE)
      #
      # metadata <- data.table(
      #   pair = samplesheet$patient,
      #   sample = samplesheet$sample,
      #   status = samplesheet$status,
      #   sex = samplesheet$sex,
      #   tumor_type = samplesheet$tumor_type,
      #   disease = samplesheet$disease,
      #   primary_site = samplesheet$primary_site
      # )
      #
      # # remove duplicates
      # metadata <- unique(metadata)
      #
      # return(metadata)

      

      if (length(launch_dir) == 0 || length(samplesheet_path) == 0) {
        warning("Could not find launch directory or samplesheet path in pipeline report")
        return(NULL)
      }


      # Clean up paths
      launch_dir <- gsub(".*launchDir: ", "", launch_dir)
      samplesheet_filename <- gsub(".*input: ", "", samplesheet_path)
      samplesheet_filename <- gsub("^\\./", "", samplesheet_filename)
      if (grepl("^/", samplesheet_filename)) { # TRUE if the samplesheet is already a full path
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
      metacols = c("patient", "sample", "tumor_type", "status", "disease", "primary_site", "sex", "bam")
      metavars = base::mget(
        metacols, 
        as.environment(as.list(samplesheet)),
        ifnotfound = rep_len(
          list(rep_len(NA_character_, NROW(samplesheet))),
          NROW(metacols)
        )
      )
      metadata <- data.table(
        pair = metavars$patient,
        sample = metavars$sample,
        tumor_type = metavars$tumor_type,
        status = metavars$status,
        disease = metavars$disease,
        primary_site = metavars$primary_site,
        inferred_sex = metavars$sex,
        bam = metavars$bam        
      )
      metadata$sample_type = ifelse(metadata$status == 0, "normal", ifelse(metadata$status == 1, "tumor", NA_character_))

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
      remaining_dt_cols <- names(dt) # Track available columns

      # First pass: exact matches only
      for (cohort_col in names(self$cohort_cols_to_x_cols)) {
        possible_cols <- self$cohort_cols_to_x_cols[[cohort_col]]

        # Get default value if it exists
        default_value <- attr(possible_cols, "default")

        # Try exact matches first
        found_col <- intersect(possible_cols, remaining_dt_cols)[1]

        if (!is.null(found_col)) {
          result_dt[, (cohort_col) := dt[[found_col]]]
          # Remove matched column from candidates
          remaining_dt_cols <- setdiff(remaining_dt_cols, found_col)
        }
      }

      # Second pass: prefix matches for columns that weren't exactly matched
      for (cohort_col in names(self$cohort_cols_to_x_cols)) {
        if (!is.null(result_dt[[cohort_col]]) && !all(is.na(result_dt[[cohort_col]]))) next # Skip if already matched

        possible_cols <- self$cohort_cols_to_x_cols[[cohort_col]]
        default_value <- attr(possible_cols, "default")

        found_col <- NULL
        # Try prefix matches
        for (col in possible_cols) {
          for (name in remaining_dt_cols) {
            if (startsWith(name, col)) {
              found_col <- name
              break
            }
          }
          if (!is.null(found_col)) break
        }

        if (!is.null(found_col)) {
          result_dt[, (cohort_col) := dt[[found_col]]]
          # Remove matched column from candidates
          remaining_dt_cols <- setdiff(remaining_dt_cols, found_col)
        } else if (!is.null(default_value) && nrow(dt) > 0) {
          # Only add default values if the input data.table is not empty
          if (is.list(default_value) || length(default_value) > 1) {
            result_dt[, (cohort_col) := list(replicate(nrow(dt), list(default_value), simplify = FALSE))]
          } else {
            result_dt[, (cohort_col) := default_value]
          }
        } else if (nrow(dt) > 0 && !cohort_col %in% Skilift:::config_parameter_names) {
          warning(sprintf(
            "No matching column found for '%s'. Expected one of: %s",
            cohort_col, paste(possible_cols, collapse = ", ")
          ))
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
  ),
  active = list()
)

Cohort$private_methods[["length"]] = function() {
  return(base::NROW(self$inputs))
}

Cohort$private_methods[["key"]] = function() {
  return(data.table::key(self$inputs))
}

Cohort$active[["length"]] = function() {
  return(private$length())
}

Cohort$active[["getkeys"]] = function() {
  input_key = private$key()
  return(self$inputs[[input_key]])
}

Cohort$public_methods[["setkey"]] = function(..., verbose = getOption("datatable.verbose"), physical = TRUE) {
  Skilift::setkey.Cohort(x = self, ... = ..., verbose = verbose, physical = physical)
  return(self)
}

Cohort$active[["key"]] = function() {
  input_key = private$key()
  return(input_key)
}

#' Cohort length
#' 
#' Get Cohort length (use inputs)
#' @export 
#' @export length.Cohort
length.Cohort = function(x) {
  return(x$length)
}

#' @export 
setkey <- function(x, ...) {
  UseMethod("setkey")
}

#' @export 
setkey.default = function(x, ..., verbose = getOption("datatable.verbose"), physical = TRUE) {
  return(
    data.table::setkey(x, ..., verbose = getOption("datatable.verbose"), physical = TRUE)
  )
}

#' Cohort setkey
#' 
#' setkey on Cohort inputs
#' @export 
#' @export setkey.Cohort
setkey.Cohort = function(x, ..., verbose = getOption("datatable.verbose"), physical = TRUE) {
  return(
    data.table::setkey(x$inputs, ..., verbose = verbose, physical = physical)
  )
}

#' @export 
key <- function(x, ...) {
  UseMethod("key")
}

#' @export 
key.default = function(x) {
  return(data.table::key(x))
}

#' Cohort key
#' 
#' Get key from Cohort inputs
#' @export 
key.Cohort = function(x) {
  return(
    data.table::key(x$inputs)
  )
}

# Cohort$private_fields = list()
# Cohort$private_fields[[".inputs"]] = data.table::data.table()


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
  purple_pp_range = "purple/.*purple.purity.range.tsv",
  purple_pp_bestFit = "purple/.*purple.purity.tsv",
  jabba_gg = "jabba/.*/jabba.simple.gg.rds$",
  events = "events/.*/complex.rds$",
  fusions = "fusions/.*/fusions.rds$",
  fragcounter_normal = "fragcounter_normal/.*/.*cov.rds$",
  fragcounter_tumor = "fragcounter_tumor/.*/.*cov.rds$",
  segments_cbs = "cbs/.*/.*(?<!n)seg.rds$",
  structural_variants = c("gridss.*/.*/.*high_confidence_somatic.vcf.bgz$", "tumor_only_junction_filter/.*/.*somatic.filtered.sv.rds$"),
  structural_variants_unfiltered = "gridss.*/.*.gridss.filtered.vcf.gz$",
  karyograph = "jabba/.*/karyograph.rds$",
  allelic_jabba_gg = "lp_phased_balance/.*/lp_phased.balanced.gg.rds$",
  somatic_snvs = c("sage/somatic/tumor_only_filter/.*/.*.sage.pass_filtered.tumoronly.vcf.gz$", "sage/somatic/.*/.*.sage.pass_filtered.vcf.gz$"),
  somatic_snvs_unfiltered = c("sage/somatic/.*/.*sage.somatic.vcf.gz$"),
  somatic_variant_annotations = "snpeff/somatic/.*/.*ann.bcf$",
  multiplicity = c("snv_multiplicity/.*/.*est_snv_cn_somatic.rds", "snv_multiplicity3/.*/.*est_snv_cn_somatic.rds"),
  germline_multiplicity = c("snv_multiplicity/.*/.*est_snv_cn_germline.rds", "snv_multiplicity3/.*/.*est_snv_cn_germline.rds"), ### TO DO FIX ME
  hetsnps_multiplicity = c("snv_multiplicity/.*/.*est_snv_cn_hets.rds", "snv_multiplicity3/.*/.*est_snv_cn_hetsnps.rds"), ### TO DO FIX ME
  ## discard for heme/tumor-only
  activities_sbs_signatures = "signatures/sigprofilerassignment/somatic/.*/sbs_results/Assignment_Solution/Activities/sbs_Assignment_Solution_Activities.txt",
  matrix_sbs_signatures = "signatures/sigprofilerassignment/somatic/.*/SBS/sigmat_results.SBS96.all",
  decomposed_sbs_signatures = "signatures/sigprofilerassignment/somatic/.*/sbs_results/Assignment_Solution/Activities/Decomposed_MutationType_Probabilities*.txt",
  activities_indel_signatures = "signatures/sigprofilerassignment/somatic/.*/indel_results/Assignment_Solution/Activities/indel_Assignment_Solution_Activities.txt",
  matrix_indel_signatures = "signatures/sigprofilerassignment/somatic/.*/ID/sigmat_results.ID83.all",
  decomposed_indel_signatures = "signatures/sigprofilerassignment/somatic/.*/indel_results/.*/Decomposed_MutationType_Probabilities.txt",
  hrdetect = "hrdetect/.*/hrdetect_results.rds",
  ## ^^ discard for heme/tumor-only
  estimate_library_complexity = "qc_reports/gatk/.*/.*metrics",
  alignment_summary_metrics = "qc_reports/picard/.*/.*alignment_summary_metrics",
  insert_size_metrics = "qc_reports/picard/.*/.*insert_size_metrics",
  wgs_metrics = "qc_reports/picard/.*/.*coverage_metrics",
  msisensor_pro = "msisensorpro/.*/.*msisensor_pro_results.tsv", ## TODO FILL ME
  oncokb_snv = "oncokb/.*/merged_oncokb.maf",
  oncokb_cna = "oncokb/.*/merged_oncokb_cna.tsv",
  oncokb_fusions = "oncokb/.*/merged_oncokb_fusions.tsv"
)

#' Default column mappings
#'
#' Defines naming convention cohort inputs table
#'
#' This variable contains the possible regex patterns
#' from an input data.table with columns corresponding to processed data.
#' You can output the default template to modify it using
#' base::dput(Skilift::default_col_mapping)
#' @export
default_col_mapping <- list(
  pair = c("pair", "patient_id", "pair_id", "sample"),
  tumor_type = c("tumor_type", "status"),
  disease = c("disease"),
  primary_site = c("primary_site"),
  inferred_sex = c("inferred_sex"),
  tumor_bam = c("tumor_bam"),
  normal_bam = c("normal_bam"),
  structural_variants = c("structural_variants", "gridss_somatic", "gridss_sv", "svaba_sv", "sv", "svs", "vcf"),
  structural_variants_unfiltered = "structural_variants_unfiltered",
  fragcounter_normal = c("fragcounter_normal"),
  fragcounter_tumor = c("fragcounter_tumor"),
  segments_cbs = c("cbs_seg_rds", "seg_rds", "cbs_seg"),
  tumor_coverage = c("tumor_coverage", "dryclean_tumor", "tumor_dryclean_cov", "dryclean_cov"),
  somatic_snvs = c("somatic_snvs", "sage_somatic_vcf", "strelka_somatic_vcf", "strelka2_somatic_vcf", "somatic_snv", "snv_vcf", "somatic_snv_vcf", "snv_somatic_vcf"),
  somatic_snvs_unfiltered = "somatic_snvs_unfiltered",
  germline_snvs = c("germline_snvs", "sage_germline_vcf", "germline_snv", "germline_snv_vcf"),
  het_pileups = c("het_pileups", "hets", "sites_txt", "hets_sites"),
  multiplicity = c("multiplicity", "somatic_snv_cn", "snv_multiplicity"),
  germline_multiplicity = c("germline_multiplicity", "germline_snv_cn"),
  hetsnps_multiplicity = c("hetsnps_multiplicity", "hets_snv_cn"),
  somatic_variant_annotations = c("somatic_variant_annotations", "annotated_bcf", "variant_somatic_ann"),
  germline_variant_annotations = c("germline_variant_annotations", "annotated_vcf_germline"),
  oncokb_snv = c("oncokb_snv", "oncokb_maf", "maf"),
  oncokb_cna = c("oncokb_cna", "cna"),
  oncokb_fusions = c("oncokb_fusions", "oncokb_fusion", "fusion_maf"),
  jabba_gg = c("jabba_gg", "jabba_simple", "jabba_rds", "jabba_simple_gg"),
  karyograph = c("karyograph"),
  balanced_jabba_gg = c("balanced_jabba_gg", "non_integer_balance", "balanced_gg", "ni_balanced_gg"),
  events = c("events", "complex"),
  fusions = c("fusions"),
  allelic_jabba_gg = c("allelic_jabba_gg", "lp_phased_balance", "allelic_gg", "lp_balanced_gg"),
  activities_sbs_signatures = c("activities_sbs_signatures", "sbs_activities"),
  matrix_sbs_signatures = c("matrix_sbs_signatures", "sbs_matrix"),
  decomposed_sbs_signatures = c("decomposed_sbs_signatures", "sbs_decomposed"),
  activities_indel_signatures = c("activities_indel_signatures", "indel_activities"),
  matrix_indel_signatures = c("matrix_indel_signatures", "indel_matrix"),
  decomposed_indel_signatures = c("decomposed_indel_signatures", "indel_decomposed"),
  hrdetect = c("hrdetect", "hrd"),
  onenesstwoness = c("onenesstwoness","oneness_twoness"),
  oncotable = c("oncotable"),
  estimate_library_complexity = c("estimate_library_complexity", "library_complexity_metrics", "est_lib_complex"),
  alignment_summary_metrics = c("alignment_summary_metrics", "alignment_metrics"),
  insert_size_metrics = c("insert_size_metrics", "insert_metrics"),
  wgs_metrics = c("wgs_metrics", "wgs_stats"),
  tumor_wgs_metrics = c("tumor_wgs_metrics", "tumor_wgs_stats"),
  normal_wgs_metrics = c("normal_wgs_metrics", "normal_wgs_stats"),
  purple_pp_range = c("purple_pp_range", "purple_range"),
  purple_pp_bestFit = c("purple_pp_bestFit", "purple_bestFit", "purple_solution"),
  msisensorpro = c("msisensor_pro", "msisensor_pro_results", "msisensor_results", "msisensorpro"),
  # Configuration parameters with default values
  metadata_is_visible = structure(c("metadata_is_visible"), default = TRUE),
  copy_number_graph_max_cn = structure(c("copy_number_graph_max_cn"), default = 100),
  copy_number_graph_annotations = structure(c("copy_number_graph_annotations"), default = list(c("bfb", "chromoplexy", "chromothripsis", "del", "dm", "cpxdm", "dup", "pyrgo", "rigma", "simple", "tic", "tyfonas"))),
  multiplicity_node_metadata = structure(c("multiplicity_node_metadata"), default = c("gene", "feature_type", "annotation", "REF", "ALT", "variant.c", "variant.p", "vaf", "transcript_type", "impact", "rank")),
  multiplicity_field = structure(c("multiplicity_field"), default = "altered_copies"),
  denoised_coverage_apply_mask = structure(c("denoised_coverage_apply_mask"), default = TRUE),
  denoised_coverage_field = structure(c("denoised_coverage_field"), default = "foreground"),
  denoised_coverage_color_field = structure(c("denoised_coverage_color_field"), default = NULL),
  denoised_coverage_bin_width = structure(c("denoised_coverage_bin_width"), default = 1e3L), ## More headaches with testing setting to NA than worth changing
  hetsnps_field = structure(c("hetsnps_field"), default = "count"),
  hetsnps_color_field = structure(c("hetsnps_color_field"), default = "col"),
  hetsnps_bin_width = structure(c("hetsnps_bin_width"), default = NA),
  hetsnps_mask = structure(c("hetsnps_mask"), default = system.file("extdata", "data", "maskA_re.rds", package = "Skilift")),
  hetsnps_subsample_size = structure(c("hetsnps_subsample_size"), default = 100000),
  hetsnps_min_normal_freq = structure(c("hetsnps_min_normal_freq"), default = 0.2),
  hetsnps_max_normal_freq = structure(c("hetsnps_max_normal_freq"), default = 0.8),
  segment_width_distribution_annotations = structure(c("segment_width_distribution_annotations"), default = NULL)
)

#' Configuration parameter names
#'
#' List of configuration parameter column names
#'
#' This variable contains the names of columns that are used
#' for configuration parameters in the Cohort inputs data.table
#' @export
config_parameter_names <- c(
  "copy_number_graph_max_cn",
  "copy_number_graph_annotations",
  "multiplicity_node_metadata",
  "multiplicity_field",
  "denoised_coverage_apply_mask",
  "denoised_coverage_field",
  "denoised_coverage_color_field",
  "denoised_coverage_bin_width",
  "hetsnps_field",
  "hetsnps_color_field",
  "hetsnps_bin_width",
  "hetsnps_mask",
  "hetsnps_subsample_size",
  "hetsnps_min_normal_freq",
  "hetsnps_max_normal_freq",
  "segment_width_distribution_annotations"
)

#' Create unified column and nf casereports map
#'
#' Parse column mapping and nf-casereport path mapping
#'
#' Takes in a list of named column regex patterns and
#' a list of named nf-casereport path regex patterns and
#' coerces them into a unified object for later use.
#' To be added as an active binding.
#' Not exported.
unify_colmap_nfmap <- function(col_mapping, path_patterns) {
  if (!inherits(col_mapping, "list") || !inherits(path_patterns, "list")) {
    stop("col_mapping and path_patterns objects must be lists")
  }
  expected_map_fields <- c(
    "column_name_regex",
    "nf_path_regex"
  )

  unified_map <- list()

  all_col_names <- union(
    names(col_mapping),
    names(path_patterns)
  )

  for (col in all_col_names) {
    unified_map[[col]] <- list(
      column_name_regex = default_col_mapping[[col]],
      nf_path_regex = nf_path_patterns[[col]]
    )
  }
  return(unified_map)
}

Cohort$private_methods[["unify_colmap_nfmap"]] <- unify_colmap_nfmap


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
parse_pipeline_paths <- function(
    paths,
    initial_dt = data.table(pair = character(0)),
    path_patterns = Skilift::nf_path_patterns,
    id_name = "pair") {
  if (length(paths) == 1 && dir.exists(paths)) {
    paths <- list.files(
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
      present_paths <- grep(pattern, paths, value = TRUE, perl = TRUE)
      if (!length(present_paths) > 0) next

      # Create mapping between paths and pairs
      dt <- data.table()
      for (pair in pairs_to_match) {
        # Find paths that contain this pair name
        # FIXME: hack for snpeff and other paths that have "snpeff/somatic/PAIR-lane_X/..." directory structure
        # Seems not to be localized to just snpeff
        # balanced_jabba_gg, cbs, jabba_gg, somatic_variant_annotations, events, allelic_jabba_gg
        pair_paths <- grep(paste0("/", pair, "(-lane_.*)?", "/"), present_paths, value = TRUE, perl = TRUE)
        # if (col_name %in% c("somatic_variant_annotations", "jabba_gg")) {
        #   pair_paths <- grep(paste0("/", pair, "(-lane.*)?", "/"), present_paths, value = TRUE, perl = TRUE)
        # } else {
        #   pair_paths <- grep(paste0("/", pair, "/"), present_paths, value = TRUE, perl = TRUE)
        # }
        
        if (length(pair_paths) > 0) {
          dt <- rbindlist(list(dt, data.table(
            pair = pair,
            path = pair_paths
          )))
        }
      }

      if (nrow(dt) > 0) {
        setnames(dt, "path", col_name)
        initial_dt <- merge.data.table(
          initial_dt,
          dt,
          by = "pair",
          all.x = TRUE,
          all.y = FALSE, # Only keep pairs that were in initial_dt
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
deparse1 <- function(expr, collapse = " ", width.cutoff = 500L, ...) {
  paste(deparse(expr, width.cutoff, ...), collapse = collapse)
}

cohort_attributes <- c(
  "inputs",
  "reference_name",
  "type",
  "nextflow_results_path"
)


#' Subset Cohort object
#'
#' Overloads subset operator for Cohort
#'
#' @export
"[.Cohort" <- function(obj, i = NULL, j = NULL, with = TRUE, ...) {
  expri <- deparse1(substitute(i))
  is_i_given <- any(expri != "NULL")
  vector_of_column_names <- deparse1(substitute(j))
  is_j_given <- any(vector_of_column_names != "NULL")
  tbl <- data.table::copy(obj$inputs)
  tblj <- tbl
  if (is_j_given) {
    if (any(grepl("\"|\'", vector_of_column_names))) vector_of_column_names <- j
    selectj <- unique(c("pair", vector_of_column_names))
    tblj <- base::subset(tblj, select = selectj)
  }
  colmap <- as.list(names(tblj))
  names(colmap) <- names(tblj)
  tbli <- tblj[]
  if (is_i_given) {
      expr = substitute(i)
      if (is.call(expr)) {
          tbli <- tblj[eval(expr), , with = with][]
      } else {
          tbli <- tblj[expr, , with = with][]
      }
      
  }
  # obj_out$inputs = tbli
  suppressWarnings({
    obj_out <- Skilift::Cohort$new(
      x = data.table()
    )
  })
  attributes_to_copy <- Skilift:::cohort_attributes
  attributes_to_copy <- attributes_to_copy[!attributes_to_copy %in% "inputs"]
  for (attribute in attributes_to_copy) {
    obj_out[[attribute]] <- obj[[attribute]]
  }
  obj_out$inputs <- tbli[]
  invisible(obj_out$inputs[])
  return(obj_out)
}


#' Parse .nextflow.log
#'
#' Parse nextflow log files for launchdir and input options
#' to get samplesheet
#'
#' @author Kevin Hadi
read_nf_log <- function(nflog_path, max_lines = 1000000L) {
  f <- file(nflog_path, open = "r")
  on.exit({
    close(f)
  })
  line <- readLines(f, n = 1L)
  is_line_at_options <- grepl("Input/output options", line)
  is_line_at_launchdir <- grepl("launchDir", line)
  are_all_lines_parsed <- is_line_at_options && is_line_at_launchdir
  counter <- 1
  ansi_pattern <- "\\\033\\[((?:\\d|;)*)([a-zA-Z])"
  input_samplesheet <- ""
  launchdir <- ""
  while (!are_all_lines_parsed && counter <= max_lines) {
    line <- readLines(f, n = 1L)
    if (!identical(is_line_at_options, TRUE)) {
      is_line_at_options <- length(line) && grepl("Input/output options", line)
      if (identical(is_line_at_options, TRUE)) {
        line <- readLines(f, n = 1L)
        input_samplesheet <- gsub(ansi_pattern, "", line)
        input_samplesheet <- gsub("input[[:space:]]+:[[:space:]]+", "", trimws(input_samplesheet))
      }
    }
    if (!identical(is_line_at_launchdir, TRUE)) {
      is_line_at_launchdir <- length(line) && grepl("launchDir", line)
      if (identical(is_line_at_launchdir, TRUE)) {
        launchdir_line <- gsub(ansi_pattern, "", line)
        launchdir <- gsub("launchDir[[:space:]]+:[[:space:]]+", "", trimws(launchdir_line))
      }
    }
    are_all_lines_parsed <- is_line_at_options && is_line_at_launchdir
    counter <- counter + 1
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
Cohort$private_methods[["read_nf_log"]] <- read_nf_log


## Already - we need to deal with backwards compatibility
## It makes sense to refer to Cohort$type and in the arguments to functions, specify "cohort_type"
## But there are cached cohorts where we've already instantiated the cohort_type attribute.
## Skilift::refresh() can restore this, but we would need to make sure that the right attributes are
## propagated forward.
refresh_attributes = list(
  list("type", "cohort_type") ## attributes that could be present, but need to be updated
)

#' Refresh Cohort object
#'
#' Reinstantiate Cohort object
#'
#' @export
refresh_cohort <- function(cohort) {
  former_inputs = data.table::copy(cohort$inputs) # Create a deep copy of the inputs
  former_inputs_colnames = names(former_inputs)
  obj_out <- Skilift::Cohort$new(
    x = former_inputs 
  )
  for (attribute in setdiff(Skilift:::cohort_attributes, "inputs")) {
    obj_out[[attribute]] <- cohort[[attribute]]
  }
  for (attribute_lst in Skilift:::refresh_attributes) {
    if ( !is.null(cohort[[ attribute_lst[[2]] ]]) ) {
      obj_out[[ attribute_lst[[1]] ]] <- cohort[[ attribute_lst[[2]] ]]
    }
  }
  for (remaining_col in former_inputs_colnames[! former_inputs_colnames %in% names(obj_out$inputs)]) {
    obj_out$inputs[[remaining_col]] = former_inputs[[remaining_col]]
  }
  return(obj_out)
}


#' pairify_cohort_inputs
#' 
#' Create paired cohort inputs if tumor and normal provided
#' for paired outputs
#' 
#' @return paired table from cohort inputs
#' @export
pairify_cohort_inputs = function(cohort, tumor_status = 1L, normal_status = 0L, tumor_normal_columns, sep_cast = "___", keep_remaining = FALSE, ...) {
  nextflow_results_path = NULL
  if (R6::is.R6(cohort) && !is.null(cohort$inputs)) {
    return_inputs = inputs = data.table::copy(cohort$inputs)
    nextflow_results_path = cohort$nextflow_results_path
  } else if (inherits(cohort, "data.table")) {
    return_inputs = inputs = data.table::copy(cohort)
  } else {
    stop("Cohort object or inputs data.table must be provided")
  }

  default_tumor_type_cols = Skilift::default_col_mapping[["tumor_type"]]

  tumor_type_vec = character(NROW(inputs))
  is_tumor_type_found = FALSE
  for (i in seq_along(default_tumor_type_cols)) {
    tumor_type_col = default_tumor_type_cols[i]
    tumor_type_vals = inputs[[tumor_type_col]]
    is_character = is.character(tumor_type_vals)
    is_numeric = is.integer(tumor_type_vals) || is.numeric(tumor_type_vals)
    is_0or1 = is_numeric && all(tumor_type_vals %in% c(0, 1))  
    if (is_character) {
      tumor_type_vec[
        grepl("tumor", tumor_type_vals, ignore.case = TRUE)
      ] = "tumor"
      tumor_type_vec[
        grepl("normal", tumor_type_vals, ignore.case = TRUE)
      ] = "normal"
    } else if (is_0or1) {
      tumor_type_vec[
        tumor_type_vals == tumor_status
      ] = "tumor"
      tumor_type_vec[
        tumor_type_vals == normal_status
      ] = "normal"
    }
    if (all(tumor_type_vec %in% c("tumor", "normal"))) is_tumor_type_found = TRUE
    if (is_tumor_type_found) {
      break
    }
  }

  # if (!is_tumor_type_found) stop("Tumor type was not found!")
  inputs$tumor_type = tumor_type_vec

	# is_status_in_cohort = !is.null(inputs$status)

	is_nextflow_results_path_present = !is.null(nextflow_results_path)
	# is_tumor_bam_in_cohort = !is.null(inputs$tumor_bam)
	# is_normal_bam_in_cohort = !is.null(inputs$normal_bam)

	# is_unpaired = is_status_in_cohort && is_nextflow_results_path_present
	# is_paired = is_tumor_bam_in_cohort && is_normal_bam_in_cohort

  is_unpaired = is_tumor_type_found ## TODO: test if this is sufficient.

  tumor_normal_columns = tumor_normal_columns[tumor_normal_columns %in% names(inputs)]

	if (is_unpaired && NROW(tumor_normal_columns) > 0) {
		# inputs$tumor_type = ifelse(inputs$status == tumor_status, "tumor", "normal")
		return_inputs = Skilift::dcastski(
			inputs,
			id_columns = "pair",
			type_columns = "tumor_type",
			cast_columns = tumor_normal_columns,
			sep_cast = sep_cast,
			keep_remaining = keep_remaining,
      ...
		)
	}
	return(return_inputs)
}

#' Cast table
#' 
#' Cast using base R.
#' 
#' Using base R for robustness and flexibility.
#' 
#' @export
dcastski = function(
	tbl, 
	id_columns, 
	type_columns, 
	cast_columns, 
	keep_remaining = FALSE, 
	sep_cast = "___", 
	prefix_type = TRUE,
	drop = FALSE,
	use_regex = FALSE
) {
  if (identical(use_regex, TRUE)) {
	.NotYetImplemented()
  }
  if (base::anyDuplicated(names(tbl)) > 0) {
    stop("Duplicated names present in table! Dedup first: names(tbl) = base::make.unique(names(tbl))")
  }
  remaining_cols = character(0)
  columns_to_process = unique(c(id_columns, type_columns))
  if (identical(keep_remaining, TRUE)) {
    remaining_cols = names(tbl)[!names(tbl) %in% c(columns_to_process, cast_columns)]
    columns_to_process = c(columns_to_process, remaining_cols)
  }
  columns_to_process = c(columns_to_process, cast_columns)
  ## A very unfortunate circumstance.
  ## Have to use unique.data.frame due to lists present in cohort inputs
  ## The data.frame base R method is slow, but works on list elements in columns.
  ## If tbl is a data.table, skeleton will remain a data.table
  skeleton = unique.data.frame(
    base::subset(tbl, select = names(tbl) %in% c(id_columns, remaining_cols))
  )
  reduced_tbl = base::subset(tbl, select = names(tbl) %in% columns_to_process)
  reduced_tbl$types = reduced_tbl[[type_columns[1]]]
  types = unique(reduced_tbl$types)
  if (is.factor(reduced_tbl$types) && !identical(drop, TRUE)) {
    types = base::levels(reduced_tbl$types)
  }
  if (length(type_columns) > 1) {
    .NotYetImplemented()
    lst = as.list(base::subset(reduced_tbl, select = names(reduced_tbl) %in% type_columns))
    types = do.call(interaction, lst)
    skeleton$types = types
  }
  for (type in types) {
    merge_type = reduced_tbl[reduced_tbl$types == type,]
    merge_type = base::subset(merge_type, select = names(merge_type) %in% c(id_columns, cast_columns))
    colnms = names(merge_type)
    names_to_change = colnms[colnms %in% cast_columns]
	if (prefix_type) {
		names_to_change = paste(type, names_to_change, sep = sep_cast)
	} else {
		names_to_change = paste(names_to_change, type, sep = sep_cast)
	}
    
    colnms[colnms %in% cast_columns] = names_to_change
    names(merge_type) = colnms
    skeleton = merge(skeleton, merge_type, by = id_columns, all.x = TRUE)
  }
  return(skeleton)
}


#' Melt table
#' 
#' Melt using base R.
#' 
#' Using base R for robustness and flexibility, not speed necessarily.
#' Since this will be used primarily for sample/pairs table manipulation.
#' 
#' @export
meltski = function(
    tbl,
    id.vars,
    measure.vars,
    is_measure_regex = FALSE,
    var.name = "variable",
    value.name = "value",
    group_regex = "",
    group_select = "\\1",
    replacement_regex = "",
    replacement = "",
    drop = FALSE,
    keep_first_variable_col = FALSE
) {
  if (missing(id.vars)) {
    id.vars = names(tbl)
    id.vars = id.vars[!id.vars %in% measure.vars]
  }
  skel = NULL
  all_names = names(tbl)
  collected_var_names = character(0)
  if (is_measure_regex) {
    message("!!WARNING!!: Regex using measure.vars does not guarantee ordered outputs.")
    message("!!WARNING!!: You must ensure that input columns are ordered with respect to the regex that you provide")
    message("!!WARNING!!: Otherwise, you should check that your output value columns are not mixed.")
    message("!!WARNING!!: Attempting to help by sorting column names first.")
    all_names = sort(all_names)
  }
  for (var_i in seq_along(measure.vars)) {
    var = measure.vars[var_i]
    if (is_measure_regex) var = grep(var, all_names, value = TRUE)
    if (is.list(var)) var = unlist(var)
    tbl_to_rbind = base::subset(tbl, select = c(id.vars, var))
    if (drop) {
      is_any_na = base::complete.cases(
        base::subset(tbl_to_rbind, select = var)
      )
      tbl_to_rbind = base::subset(tbl_to_rbind, is_any_na)
    }
    nm = names(tbl_to_rbind)
    for (i in seq_along(var)) {
      len_value_name = NROW(value.name)
      if (!(len_value_name > 1 && len_value_name == length(var))) {
        value_col_name = paste(value.name, "_", i, sep = "")
      } else {
        value_col_name = value.name[i]
      }
      nm[nm %in% var[i]] = value_col_name
    }
    names(tbl_to_rbind) = nm
    for (i in seq_along(var)) {
      var.col = paste(var.name, "_", i, sep = "")
      tbl_to_rbind[[var.col]] = var[i]
      collected_var_names = c(collected_var_names, var.col)
    }    
    skel = rbind(skel, tbl_to_rbind)
  }
  collected_var_names = unique(collected_var_names)
  if (nzchar(group_regex)) {
    for (variable in collected_var_names) {
      skel[[variable]] = gsub(group_regex, group_select, skel[[variable]])
    }
  }
  if (nzchar(replacement_regex)) {
    for (variable in collected_var_names) {
      skel[[variable]] = gsub(replacement_regex, replacement, skel[[variable]])
    }    
  }
  remaining_var_names = collected_var_names[-1]
  remaining_names = names(skel)[!names(skel) %in% remaining_var_names]
  if (identical(keep_first_variable_col, TRUE) && NROW(remaining_names) > 0) {
    skel = base::subset(skel, select = remaining_names)
    names(skel)[names(skel) == collected_var_names[1]] = gsub("_[0-9]+$", "", collected_var_names[1])
  }
  return(skel)
}


#' Merge Cohort objects
#'
#' Combines two or more Cohort objects into a single Cohort
#'
#' @param x Cohort object to merge
#' @param y Cohort object to merge
#' @param ... Arbitrary remaining Cohort objects to merge
#' @param warn_duplicates Logical indicating whether to warn about duplicate pairs (default TRUE)
#' @param rename_duplicates Logical indicating whether to rename duplicate pairs instead of overwriting (default FALSE)
#' @return A new Cohort object containing all data from input Cohorts
#' @export
merge.Cohort <- function(x, y, ..., warn_duplicates = TRUE, rename_duplicates = FALSE) {
  cohorts_dots <- list(...)
  is_x_missing = missing(x)
  is_y_missing = missing(y)
  is_x_or_y_missing = is_x_missing || is_y_missing
  number_cohorts_dots = NROW(cohorts_dots)
  length_cohorts = sum(as.integer(!is_x_missing), as.integer(!is_y_missing), number_cohorts_dots)
  # Validate inputs are all Cohorts
  if (length_cohorts < 2) {
    stop("At least two Cohort objects must be provided")
  }
  if (is_x_or_y_missing) {
    if (is_x_missing) print("x argument not provided!")
    if (is_y_missing) print("y argument not provided!")
    stop("Both x and y must be provided as Cohort objects")
  }
  cohorts = c(list(x), list(y), cohorts_dots)
  for (cohort in cohorts) {
    if (!inherits(cohort, "Cohort")) {
      stop("All arguments must be Cohort objects")
    }
  }

  # Get first cohort to initialize merged result
  merged_dt <- data.table::copy(cohorts[[1]]$inputs)

  # Merge remaining cohorts
  for (i in 2:length(cohorts)) {
    current_dt <- data.table::copy(cohorts[[i]]$inputs)

    # Check for duplicate pairs
    duplicate_pairs <- intersect(merged_dt$pair, current_dt$pair)
    if (length(duplicate_pairs) > 0) {
      if (warn_duplicates) {
        warning(sprintf(
          "Found %d duplicate pair(s): %s",
          length(duplicate_pairs),
          paste(duplicate_pairs, collapse = ", ")
        ))
      }

      if (rename_duplicates) {
        # Add incrementing suffix to duplicate pairs
        for (dup_pair in duplicate_pairs) {
          suffix <- 1
          while (paste0(dup_pair, "_", suffix) %in% merged_dt$pair) {
            suffix <- suffix + 1
          }
          current_dt[pair == dup_pair, pair := paste0(pair, "_", suffix)]
        }
      } else {
        # Remove duplicates from current_dt that would overwrite existing pairs
        current_dt <- current_dt[!pair %in% duplicate_pairs]
      }
    }

    # Merge inputs data.tables
    merged_dt <- merge(
      merged_dt, 
      current_dt,
      by="pair", 
      all.x = TRUE, ## Need these two on to do outer join, which seems to be expected behavior
      all.y = TRUE
    )
  }

  # Create new Cohort with merged data
  result <- Cohort$new(merged_dt, cohort_type=cohorts[[1]]$type)

  return(result)
}
