#' @import R6
#' @import data.table
NULL

#' CaseReportsData Class
#'
#' @description
#' A class to handle case reports data
#'
#' @export
CaseReportsData <- R6::R6Class(
  "CaseReportsData",
  
  public = list(
    #' @field case_reports data.table containing case reports data
    case_reports = NULL,

    #' @description
    #' Initialize a new CaseReportsData object
    #' @param path Character string specifying path to case reports data directory
    #' @return A new CaseReportsData object
    initialize = function(case_reports_data_path) {
      private$construct_from_path(case_reports_data_path)
    }
  ),

  private = list(
    get_case_dirs = function(path) {
      if (!dir.exists(path)) {
        stop("Case reports directory does not exist: ", path)
      }
      list.dirs(path, full.names = TRUE, recursive = FALSE)
    },

    get_case_id = function(case_dir) {
      basename(case_dir)
    },

    find_case_files = function(case_dir, file_patterns = "\\.(json|arrow)$") {
      files <- list.files(case_dir, pattern = file_patterns, 
                          full.names = TRUE, recursive = FALSE)
      names(files) <- tools::file_path_sans_ext(basename(files))
      files
    },

    create_case_row = function(case_dir) {
      case_id <- self$get_case_id(case_dir)
      files <- self$find_case_files(case_dir)
      
      # Create a named list with NA as default value
      expected_files <- c(
        "metadata",
        "coverage",
        "sage_qc",
        "filtered.events",
        "hetsnps",
        "allelic",
        "complex",
        "ppfit",
        "mutation_catalog",
        "id_mutation_catalog",
        "sbs_decomposed_prob",
        "id_decomposed_prob"
      )
      row_data <- setNames(rep(NA_character_, length(expected_files)), expected_files)
      
      # Update with found files
      for (name in names(files)) {
        if (name %in% names(row_data)) {
          row_data[[name]] <- files[[name]]
        } else {
          warning(sprintf("Unexpected file found in case %s: %s", case_id, name))
        }
      }
      
      # Add case_id as pair
      row_data$pair <- case_id
      
      # Convert to data.table row
      as.data.table(row_data)
    },

    construct_from_path = function(path) {
      case_dirs <- private$get_case_dirs(path)
      
      # Create data.table from all case directories
      self$case_reports <- rbindlist(
        lapply(case_dirs, private$create_case_row),
        fill = TRUE
      )
      
      # Check for missing files
      for (col in setdiff(names(self$case_reports), "pair")) {
        missing <- is.na(self$case_reports[[col]])
        if (any(missing)) {
          warning(sprintf("Missing %s file in cases: %s", 
                         col, 
                         paste(self$case_reports[missing]$pair, collapse = ", ")))
        }
      }
    },

    validate_case_reports_files = function() {
      # Will be implemented in next step
    },

    filter_cases = function() {
      # Will be implemented in next step
    }
  )
)
