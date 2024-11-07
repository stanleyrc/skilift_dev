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
        oncotable = c("oncotable")
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
    }
  ),
  
  private = list(
    construct_from_path = function(path) {
      stop("Method not implemented yet")
    },
    
    construct_from_datatable = function(dt) {
      if (!is.data.table(dt)) {
        stop("Input must be a data.table")
      }
      
      result_dt <- data.table()
      
      # browser()
      for (cohort_col in names(self$cohort_cols_to_x_cols)) {
        possible_cols <- self$cohort_cols_to_x_cols[[cohort_col]]
        
        found_col <- NULL
        for (col in possible_cols) {
          for (name in names(dt)) {
            if (startsWith(name, col)) {
              found_col <- name
              break
            }
          }
          if (!is.null(found_col)) break
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
    
    validate_inputs = function() {
      stop("Method not implemented yet")
    },
    
    filter_inputs = function() {
      stop("Method not implemented yet")
    }
  )
)
