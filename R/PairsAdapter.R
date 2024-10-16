library(R6)
library(data.table)

PairsAdapter <- R6Class("PairsAdapter",
  public = list(
    pairslike = NULL,
    column_map = NULL,
    
    initialize = function(pairslike, column_map) {
      self$pairslike <- pairslike
      self$column_map <- column_map

      # Default column map based on spec.md
      default_column_map = list(
        patient_id = "patient_id",
        tumor_type = "tumor_type",
        disease = "disease",
        primary_site = "primary_site",
        inferred_sex = "inferred_sex",
        non_integer_balance = "non_integer_balance",
        dryclean_tumor = "dryclean_tumor",
        events = "events",
        fusions = "fusions",
        gridss_somatic = "gridss_somatic",
        karyograph = "karyograph",
        lp_phased_balance = "lp_phased_balance",
        dryclean_normal = "dryclean_normal",
        fragcounter_normal = "fragcounter_normal",
        fragcounter_tumor = "fragcounter_tumor",
        het_pileups_wgs = "het_pileups_wgs",
        hrdetect = "hrdetect",
        jabba_simple = "jabba_simple",
        indel_activities = "indel_activities",
        sbs_activities = "sbs_activities",
        indel_matrix = "indel_matrix",
        sbs_matrix = "sbs_matrix",
        somatic_snv_cn = "somatic_snv_cn",
        germline_snv_cn = "germline_snv_cn",
        normal_sample = "normal_sample",
        tumor_sample = "tumor_sample",
        tumor_bam = "tumor_bam",
        normal_bam = "normal_bam",
        gridss_somatic_all = "gridss_somatic_all",
        sage_somatic_vcf = "sage_somatic_vcf",
        sage_germline_vcf = "sage_germline_vcf",
        annotated_vcf = "annotated_vcf",
        annotated_bcf = "annotated_bcf",
        annotated_vcf_germline = "annotated_vcf_germline",
        annotated_bcf_germline = "annotated_bcf_germline",
        cbs_cov = "cbs_cov",
        cbs_nseg = "cbs_nseg",
        cbs_seg = "cbs_seg"
      )
      
      # Merge user-provided column map with default
      self$column_map <- modifyList(default_column_map, column_map)
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
