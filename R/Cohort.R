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
      
      # Set up default column mapping
      default_mapping <- list(
        pair = c("pair", "sample", "sample_id", "pair_id"),
        tumor = c("tumor", "tumor_id", "tumor_sample"),
        normal = c("normal", "normal_id", "normal_sample"),
        annotated_bcf = c("annotated_bcf", "bcf", "vcf"),
        fusions = c("fusions", "fusion_file"),
        jabba_simple = c("jabba_simple", "jabba", "jabba_file"),
        karyograph = c("karyograph", "kg", "kg_file"),
        events = c("events", "complex", "complex_events"),
        signature_counts = c("signature_counts", "signatures"),
        oncokb_maf = c("oncokb_maf", "maf"),
        oncokb_cna = c("oncokb_cna", "cna")
      )
      
      # Merge user-provided mapping with default mapping
      if (!is.null(col_mapping)) {
        for (col_name in names(col_mapping)) {
          if (col_name %in% names(default_mapping)) {
            default_mapping[[col_name]] <- unique(c(default_mapping[[col_name]], col_mapping[[col_name]]))
          } else {
            default_mapping[[col_name]] <- col_mapping[[col_name]]
          }
        }
      }
      
      self$cohort_cols_to_x_cols <- default_mapping
      
        self$inputs <- private$construct_from_path(x)
      } else if (is.data.table(x)) {
        self$inputs <- private$construct_from_datatable(x)
      } else {
        stop("Input must be either a path (character) or data.table")
      }
    },
  ),
  
  private = list(
    construct_from_path = function(path) {
      if (!is.data.table(dt)) {
        stop("Input must be a data.table")
      }
      
      # Create new data.table to store results
      result_dt <- data.table()
      
      # For each desired column in the cohort
      for (cohort_col in names(self$cohort_cols_to_x_cols)) {
        # Get possible column names in input
        possible_cols <- self$cohort_cols_to_x_cols[[cohort_col]]
        
        # Find first matching column in input data
        found_col <- NULL
        for (col in possible_cols) {
          if (col %in% names(dt)) {
            found_col <- col
            break
          }
        }
        
        if (!is.null(found_col)) {
          # Add column to result
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
    
    construct_from_datatable = function(dt) {
      stop("Method not implemented yet")
    },
    
    validate_inputs = function() {
      stop("Method not implemented yet")
    },
    
    filter_inputs = function() {
      stop("Method not implemented yet")
    }
  )
)
