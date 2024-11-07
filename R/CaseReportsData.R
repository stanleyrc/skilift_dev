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
    initialize = function(path) {
      private$construct_from_path(path)
    }
  ),

  private = list(
    construct_from_path = function(path) {
      # Will be implemented in next step
    },

    validate_case_reports_files = function() {
      # Will be implemented in next step
    },

    filter_cases = function() {
      # Will be implemented in next step
    }
  )
)
