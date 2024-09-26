library(R6)
library(data.table)

PairsAdapter <- R6Class("PairsAdapter",
  public = list(
    pairslike = NULL,
    column_map = NULL,
    
    initialize = function(pairslike, column_map) {
      self$pairslike <- pairslike
      self$column_map <- column_map

      default_column_map = c(
        "patient_id", "tumor_type", "disease", "primary_site", "inferred_sex",
        "non_integer_balance", "dryclean_tumor", "events", "fusions", "gridss_somatic",
        "karyograph", "lp_phased_balance", "dryclean_normal", "fragcounter_normal",
        "fragcounter_tumor", "het_pileups_wgs", "hrdetect", "jabba_simple",
        "indel_activities", "sbs_activities", "indel_matrix", "sbs_matrix",
        "somatic_snv_cn", "germline_snv_cn", "normal_sample", "tumor_sample",
        "tumor_bam", "normal_bam", "gridss_somatic_all", "sage_somatic_vcf",
        "sage_germline_vcf", "annotated_vcf", "annotated_bcf", "annotated_vcf_germline",
        "annotated_bcf_germline", "cbs_cov", "cbs_nseg", "cbs_seg"
      )
      
      # Merge user-provided column map with default
      self$column_map <- if (is.null(column_map)) {
        default_column_map
      } else {
        column_map
      }
    },
    
    handlePairsRds = function() {
      # Placeholder for handling RDS input
      # This function should read the RDS file and return a data.table
      return(data.table::data.table())
    },
    
    handlePairsTable = function() {
      # Placeholder for handling data.table input
      # This function should process the data.table and return it
      return(self$pairslike)
    },
    
    handlePipelineResultsPath = function() {
      # Placeholder for handling path to results directory
      # This function should read files from the directory and return a data.table
      return(data.table::data.table())
    },
    
    getPatients = function() {
      if (is.character(self$pairslike)) {
        return(self$handlePipelineResultsPath())
      } else if (inherits(self$pairslike, "data.table")) {
        return(self$handlePairsTable())
      } else if (file.exists(self$pairslike)) {
        return(self$handlePairsRds())
      } else {
        stop("Unsupported pairslike input")
      }
    }
  )
)
