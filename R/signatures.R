#' @name read_decomposed_mutationtype_probabilities
#' @title read_decomposed_mutationtype_probabilities
#' @description
#' Reads and formats decomposed probabilities from sigprofiler outputs for a single sample
#'
#' @param decomposed_probs_path Path to sigprofiler decomposed probabilities file
#' @param is_indel Boolean to indicate if the input matrix is for indels
#' @return data.frame containing formatted decomposed probabilities
#' @export
#' @author Johnathan Rafailov, Shihab Dider
read_decomposed_mutationtype_probabilities <- function(
    decomposed_probs_path,
    is_indel = FALSE
) {
    decomposed_probs <- fread(decomposed_probs_path)
    colnames(decomposed_probs)[1] <- "samples"
    decomposed_probs[, samples := gsub("_somatic", "", samples)]
    decomposed_probs_long <- melt(decomposed_probs,
        measure.vars = c(3:ncol(decomposed_probs)),
        variable.name = "signature",
        value.name = "p"
    )
    if (is_indel) {
        result <- data.frame(
            signature = decomposed_probs_long[, signature],
            insdel = decomposed_probs_long[[2]],
            p = decomposed_probs_long[, p]
        )
    } else {
        result <- data.frame(
            signature = decomposed_probs_long[, signature],
            tnc = decomposed_probs_long[[2]],
            p = decomposed_probs_long[, p]
        )
    }
    return(result)
}

#' @name create_mutations_catalog
#' @title create_mutations_catalog
#' @description
#' Creates a mutations catalog data frame from sigprofiler matrix generator output
#'
#' @param sigmat_path Path to either the SBS or ID matrix
#' @param is_indel Boolean to indicate if the matrix is for indels
#' @return list containing:
#'   - data: data.frame with columns:
#'     - id: row identifier
#'     - insdel/tnc: mutation type (insdel for indels, tnc for SBS)
#'     - mutations: count of mutations
#' @export
#' @author Shihab Dider, Sukanya Panja
create_mutations_catalog <- function(sigmat_path, is_indel = FALSE) {
    # Read and process input matrix
    sig_matrix <- fread(sigmat_path)
    sig_matrix_dt <- as.data.frame(sig_matrix)
    
    # Create the data frame structure
    result <- if (is_indel) {
        data.frame(
            id = 1:nrow(sig_matrix_dt),
            insdel = sig_matrix_dt[, 1],
            mutations = sig_matrix_dt[, 2]
        )
    } else {
        data.frame(
            id = 1:nrow(sig_matrix_dt),
            tnc = sig_matrix_dt[, 1],
            mutations = sig_matrix_dt[, 2]
        )
    }
    
    # Return as a list with 'data' element
    return(list(data = result))
}

#' @name lift_signatures
#' @title lift_signatures
#' @description
#' Create signature files for all samples in a cohort
#'
#' @param cohort Cohort object containing sample information
#' @param output_data_dir Base directory for output files
#' @param cores Number of cores for parallel processing (default: 1)
#' @return None
#' @export
lift_signatures <- function(cohort, output_data_dir, cores = 1) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }
    
    if (!dir.exists(output_data_dir)) {
        dir.create(output_data_dir, recursive = TRUE)
    }
    
    # Define processing configurations
    signature_configs <- list(
        sbs = list(
            catalog = list(
                input_col = "matrix_sbs_signatures",
                output_file = "mutation_catalog.json",
                is_indel = FALSE,
                processor = create_mutations_catalog
            ),
            decomposed = list(
                input_col = "decomposed_sbs_signatures",
                output_file = "sbs_decomposed_prob.json",
                is_indel = FALSE,
                processor = read_decomposed_mutationtype_probabilities
            )
        ),
        indel = list(
            catalog = list(
                input_col = "matrix_indel_signatures",
                output_file = "id_mutation_catalog.json",
                is_indel = TRUE,
                processor = create_mutations_catalog
            ),
            decomposed = list(
                input_col = "decomposed_indel_signatures",
                output_file = "id_decomposed_prob.json",
                is_indel = TRUE,
                processor = read_decomposed_mutationtype_probabilities
            )
        )
    )

    # Validate required columns exist
    required_cols <- unique(sapply(signature_configs, function(type) {
        sapply(type, function(config) config$input_col)
    }))
    missing_cols <- required_cols[!required_cols %in% names(cohort$inputs)]
    if (length(missing_cols) > 0) {
        stop("Missing required columns in cohort: ", paste(missing_cols, collapse = ", "))
    }
    
    # Process each sample in parallel
    mclapply(seq_len(nrow(cohort$inputs)), function(i) {
        row <- cohort$inputs[i,]
        pair_dir <- file.path(output_data_dir, row$pair)
        
        if (!dir.exists(pair_dir)) {
            dir.create(pair_dir, recursive = TRUE)
        }
        
        futile.logger::flog.threshold("ERROR")
        tryCatchLog({
            # Process all signature types
            for (sig_type in names(signature_configs)) {
                for (process_type in names(signature_configs[[sig_type]])) {
                    config <- signature_configs[[sig_type]][[process_type]]
                    input_file <- row[[config$input_col]]
                    
                    if (file.exists(input_file)) {
                        result <- config$processor(
                            input_file,
                            is_indel = config$is_indel
                        )
                        jsonlite::write_json(
                            result,
                            file.path(pair_dir, config$output_file),
                            pretty = TRUE
                        )
                    } else {
                        warning(sprintf("Missing input file for %s: %s", row$pair, input_file))
                    }
                }
            }
        }, error = function(e) {
            print(sprintf("Error processing %s: %s", row$pair, e$message))
            NULL
        })
    }, mc.cores = cores, mc.preschedule = TRUE)
    
    invisible(NULL)
}
