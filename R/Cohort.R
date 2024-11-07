#' @import R6
#' @import data.table
Cohort <- R6Class("Cohort",
  public = list(
    #' @field inputs data.table containing cohort information
    inputs = NULL,
    
    #' @field reference_name character string specifying genome reference, defaults to hg19
    reference_name = "hg19",

    #' @description
    #' Initialize a new Cohort object
    #' @param x Either a data.table or path to pipeline output directory
    #' @param reference_name character string specifying genome reference
    initialize = function(x, reference_name = "hg19") {
      self$reference_name <- reference_name
      
      if (is.character(x) && length(x) == 1) {
        # x is a path
        self$inputs <- private$construct_from_path(x)
      } else if (is.data.table(x)) {
        # x is a data.table
        self$inputs <- private$construct_from_datatable(x)
      } else {
        stop("Input must be either a path (character) or data.table")
      }
    }
  ),
  
  private = list(
    # Placeholder methods to be implemented later
    construct_from_path = function(path) {
      stop("Method not implemented yet")
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
