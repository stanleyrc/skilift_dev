#' @import R6
#' @import data.table
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
        jabba_gg = c("jabba_gg", "jabba_simple", "jabba_rds", "jabba_simple_gg"),
        karyograph = c("karyograph"),
        balanced_jabba_gg = c("balanced_jabba_gg", "non_integer_balance", "balanced_gg"),
        events = c("events", "complex"),
        fusions = c("fusions"),
        allelic_jabba_gg = c("allelic_jabba_gg", "lp_phased_balance", "allelic_gg"),
        activities_sbs_signatures = c("activities_sbs_signatures", "sbs_activities"),
        matrix_sbs_signatures = c("matrix_sbs_signatures", "sbs_matrix"),
        activities_indel_signatures = c("activities_indel_signatures", "indel_activities"),
        matrix_indel_signatures = c("matrix_indel_signatures", "indel_matrix"),
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
      outputs <- data.table(pair = sample_metadata$pair)
      
      # Get all file paths recursively
      pipeline_output_paths <- list.files(pipeline_outdir, recursive = TRUE, full.names = TRUE)
      
      # Process each path and build the data.table
      for (path in pipeline_output_paths) {
        column <- private$extract_pipeline_outpath_to_column(path)
        if (!is.null(column)) {
          # Extract the pair ID from the path components
          path_parts <- strsplit(path, "/")[[1]]
          # Find which pair this file belongs to by matching the path against each pair ID
          for (curr_pair in outputs$pair) {
            if (any(grepl(curr_pair, path_parts))) {
              outputs[pair == curr_pair, (column) := path]
              break
            }
          }
        }
      }
      
      # Merge with rest of sample metadata
      outputs <- merge(outputs, sample_metadata, by = "pair", all.x = TRUE)
      
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
      launch_dir <- gsub(".*launchDir: ", "", launch_dir)
      samplesheet_path <- gsub(".*input: ", "", samplesheet_path)
      samplesheet_path <- file.path(launch_dir, gsub("^\\./", "", samplesheet_path))
      
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
    extract_pipeline_outpath_to_column = function(path) {
      # Map of regex patterns to column names
      path_patterns <- list(
        "balanced_jabba_gg" = "non_integer_balance/.*/balanced.gg.rds$",
        "tumor_coverage" = "dryclean_tumor/.*/drycleaned.cov.rds$",
        "het_pileups" = "hetpileups/.*/sites.txt$",
        "jabba_gg" = "jabba/.*/jabba.simple.gg.rds$",
        "events" = "events/.*/complex.rds$",
        "fusions" = "fusions/.*/fusions.rds$",
        "structural_variants" = "gridss_somatic/.*/.*high_confidence_somatic.vcf.bgz$",
        "karyograph" = "jabba/.*/karyograph.rds$",
        "allelic_jabba_gg" = "lp_phased_balance/.*/balanced.gg.rds$",
        "somatic_snvs" = "sage/somatic/.*/.*sage.somatic.vcf.gz$",
        "somatic_variant_annotations" = "snpeff/somatic/.*/.*ann.bcf$",
        "somatic_snv_cn" = "snv_multiplicity3/.*/.*est_snv_cn_somatic.rds",
        "activities_sbs_signatures" = "signatures/sigprofilerassignment/somatic/.*/sbs_results/.*/Assignment_Solution_Activities.txt",
        "activities_indel_signatures" = "signatures/sigprofilerassignment/somatic/.*/indel_results/.*/Assignment_Solution_Activities.txt",
        "matrix_sbs_signatures" = "signatures/sigprofilerassignment/somatic/.*/SBS/sigmat_results.SBS96.all",
        "matrix_indel_signatures" = "signatures/sigprofilerassignment/somatic/.*/ID/sigmat_results.ID28.all",
        "hrdetect" = "hrdetect/.*/hrdetect_results.rds"
      )
      
      for (col_name in names(path_patterns)) {
        if (grepl(path_patterns[[col_name]], path)) {
          return(col_name)
        }
      }
      
      return(NULL)
    },

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
      
      return(result_dt)
    },
    
    filter_inputs = function() {
      stop("Method not implemented yet")
    }
  )
)
