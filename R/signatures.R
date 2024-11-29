#' @name sigprofiler_decomposed_probs_json
#' @title sigprofiler_decomposed_probs_json
#' @description
#'
#' function to create decomposed probabilities json for case reports as adapted to sigprofiler outputs
#'
#' @param probs path to sigprofiler decomposed probabilities
#' @param is_indel Boolean to indicate if the input matrix is for indels
#' @param data_dir path to data directory of which each subdirectory contains jsons per patient
#' @param pairs vector of samples for which to create jsons, by default = NULL, indicating process all patients
#' @param cores number of cores to parallelize with
#' @return data.table or NULL
#' @export
#' @author Johnathan Rafailov
sigprofiler_decomposed_probs_json <- function(probs, is_indel, data_dir, pairs, cores = 1) {
    fread(probs) -> decomposed.probs
    ### lets change first column to something more predictable
    colnames(decomposed.probs)[1] <- "samples"
    ### seems like _somatic is systematically ended to the end of all sample names in the beginning
    decomposed.probs[, samples := gsub("_somatic", "", samples)]
    if (!is.null(pairs)) {
        decomposed.probs <- decomposed.probs[samples %in% pairs, ]
    }
    decomposed.probs.per.sample <- melt(decomposed.probs,
        measure.vars = c(3:ncol(decomposed.probs)),
        variable.name = "signature",
        value.name = "p"
    )
    samples <- unique(decomposed.probs.per.sample$samples)
    mclapply(samples, function(pair) {
        pair_data <- decomposed.probs.per.sample[samples == pair]
        if (is_indel) {
            samp_data <- data.frame(
                signature = pair_data[, signature],
                insdel = pair_data[[2]],
                p = pair_data[, p]
            )
            catalog_file_name <- "id_decomposed_prob.json"
        } else {
            samp_data <- data.frame(
                signature = pair_data[, signature],
                tnc = pair_data[[2]],
                p = pair_data[, p]
            )
            catalog_file_name <- "sbs_decomposed_prob.json"
        }
        write_json(samp_data, paste0(data_dir, "/", pair, "/", catalog_file_name))
    }, mc.cores = cores)
}

#' @name create_mutations_catalog_json
#' @title create_mutations_catalog_json
#' @description
#'
#' Builds the mutations catalog for case reports using the output of sigprofiler matrix generator
#'
#' @param sig_matrix_path Path to either the SBS or ID matrix
#' @param is_indel Boolean to indicate if the matrix is for indels
#' @param output_dir directory to write the jsons, will create directory for each patient id
#' @return
#' @export
#' @author Shihab Dider, Sukanya Panja
create_mutations_catalog_json <- function(
    sig_matrix_path,
    is_indel,
    output_dir) {
    sig_matrix <- fread(sig_matrix_path)
    sig_matrix_dt <- as.data.frame(sig_matrix)
    samples <- colnames(sig_matrix_dt)
    for (i in 2:(ncol(sig_matrix_dt))) { # skip first column as it has the tnc
        pair <- samples[i]

        if (is_indel) {
            samp_data <- data.frame(
                id = 1:nrow(sig_matrix_dt),
                insdel = sig_matrix_dt[, 1],
                mutations = sig_matrix_dt[, i]
            )
            catalog_file_name <- "id_mutation_catalog.json"
        } else {
            samp_data <- data.frame(
                id = 1:nrow(sig_matrix_dt),
                tnc = sig_matrix_dt[, 1],
                mutations = sig_matrix_dt[, i]
            )
            catalog_file_name <- "mutation_catalog.json"
        }
        pair_data <- list(pair = pair, data = samp_data)
        system(paste("mkdir -p", paste0(output_dir, "/", pair)))
        output_path <- paste0(output_dir, "/", catalog_file_name)
        write_json(
            pair_data,
            output_path,
            pretty = TRUE,
            auto_unbox = TRUE
        )
    }
}


#' @name create_distributions
#' @title create_distributions
#' @description
#'
#' Creates the background distributions for aggregated metadata KPIs across all patients
#'
#' @param case_reports_datadir directory containing all patient directories
#' @param filter_patients list() - subset of samples to filter on to create the distributions
#' @param cores number of cores to use for parallel processing
#' @param write_to_json boolean flag to indicate whether to write the files to JSON or not
#' @param common_dir path to the common directory in which to write the distribution JSONs (required if write_to_json is TRUE)
#' @param overwrite boolean flag to indicate whether to overwrite existing files (required if write_to_json is TRUE)
#' @return List of data tables containing distributions if write_to_json is FALSE, otherwise NULL
#' @export
create_distributions <- function(
    case_reports_datadir,
    filter_patients = NULL,
    cores = 1,
    write_to_json = FALSE,
    common_dir = NULL,
    overwrite = FALSE,
    haveSignatures = TRUE) {
    case_reports_data_folder <- paste0(case_reports_datadir, "/")

    files <- list.files(
        case_reports_data_folder,
        pattern = "metadata.json",
        recursive = TRUE,
        full.names = TRUE
    )

    patient_ids <- basename(dirname(files))

    meta.dt <- data.table(
        meta_json = files,
        patient_id = patient_ids
    )

    meta.dt <- meta.dt[file.exists(meta_json), ]

    if (!is.null(filter_patients)) {
        meta.dt <- meta.dt[patient_id %in% filter_patients, ]
    }

    jsons.lst <- mclapply(1:nrow(meta.dt), function(i) {
        meta.sub.dt <- meta.dt[i, ]
        read_meta_data_json(meta_json = meta.sub.dt$meta_json, patient_id = meta.sub.dt$patient_id)
    }, mc.cores = cores)
    jsons.dt <- rbindlist(jsons.lst, fill = TRUE)

    subset_and_setnames <- function(dt, cols, new_names) {
        if (all(cols %in% names(dt))) {
            return(dt[, ..cols] %>% setnames(., new_names))
        }
        return(NULL)
    }

    # add new metadata attributes here
    column_map <- list(
        snv.dt = list(
            cols = c("pair", "snv_count", "tumor_type"),
            new_names = c("pair", "value", "tumor_type")
        ),
        sv.dt = list(
            cols = c("pair", "sv_count", "tumor_type"),
            new_names = c("pair", "value", "tumor_type")
        ),
        loh.dt = list(
            cols = c("pair", "tumor_type", "loh_fraction", "loh_seglen", "loh_total_genome"),
            new_names = c("pair", "tumor_type", "value", "LOH_seg_len", "genome_width")
        ),
        ploidy.dt = list(
            cols = c("pair", "tumor_type", "ploidy", "purity"),
            new_names = c("pair", "tumor_type", "value", "purity")
        ),
        purity.dt = list(
            cols = c("pair", "tumor_type", "ploidy", "purity"),
            new_names = c("pair", "tumor_type", "ploidy", "value")
        ),
        cov_var.dt = list(
            cols = c("pair", "tumor_type", "dlrs"),
            new_names = c("pair", "tumor_type_mod", "value")
        ),
        tmb.dt = list(
            cols = c("pair", "tmb", "tumor_type"),
            new_names = c("pair", "value", "tumor_type")
        )
    )

    results <- lapply(names(column_map), function(name) {
        subset_and_setnames(jsons.dt, column_map[[name]]$cols, column_map[[name]]$new_names)
    })

    names(results) <- names(column_map)

    # signatures
    if (haveSignatures) {
        sigs.lst <- convert_signature_meta_json(jsons.dt = jsons.dt, cores = cores)

        sigs_names <- c(
            "deconstruct_sigs.dt",
            "sigprofilerassignment_indels.dt",
            "sigprofilerassignment_indels_counts.dt",
            "sigprofilerassignment_sbs_counts.dt",
            "sigprofilerassignment_sbs.dt"
        )

        sigs_results <- lapply(sigs_names, function(name) {
            if (name %in% names(sigs.lst)) {
                return(sigs.lst[[name]])
            }
            return(NULL)
        })

        json.lst <- c(results, sigs_results)

        names(json.lst) <- c(
            "snvCount",
            "svCount",
            "lohFraction",
            "ploidy",
            "purity",
            "coverageVariance",
            "tmb",
            "sbs",
            "sigprofiler_indel_fraction",
            "sigprofiler_indel_count",
            "sigprofiler_sbs_count",
            "sigprofiler_sbs_fraction"
        )
    } else {
        json.lst <- results
        names(json.lst) <- c(
            "snvCount",
            "svCount",
            "lohFraction",
            "ploidy",
            "purity",
            "coverageVariance",
            "tmb"
        )
    }

    if (write_to_json) {
        if (is.null(common_dir)) {
            stop("common_dir must be provided if write_to_json is TRUE")
        }
        write_distributions_to_json(json.lst, common_dir, cores, overwrite, haveSignatures)
        return(NULL)
    } else {
        return(json.lst)
    }
}

#' @name write_signature_jsons
#' @title write_signature_jsons
#' @description
#' Writes the different distribution JSONs for signatures from the list outputted by convert_signature_meta_json
#'
#' @param signatures_list list outputted from convert_sign
#' @param common_folder location of the signatures folder
#' @param cores cores for generating the signature distribution files
#' @return NULL
#' @export
#' @author Shihab Dider, Stanley Clarke
write_signature_jsons <- function(signatures_list, common_folder, cores = 1) {
    write_signature_json <- function(sig.dt, folder_path, sig) {
        sig.dt[, sig := sig]
        write_json(
            sig.dt,
            paste0(common_folder, folder_path, sig, ".json"),
            pretty = TRUE
        )
    }

    process_signatures <- function(sig_data, folder_path) {
        sig_add <- names(sig_data) %>% grep("pair|tumor_type", ., invert = TRUE, value = TRUE)
        empty.lst <- mclapply(sig_add, function(sig) {
            cols_keep <- c("pair", "tumor_type", sig)
            sig.dt <- sig_data[, ..cols_keep] %>% setnames(., c("pair", "tumor_type", "value"))
            write_signature_json(sig.dt, folder_path, sig)
            return(NULL)
        }, mc.cores = cores)
    }

    if ("sbs" %in% names(signatures_list)) {
        process_signatures(
            signatures_list[["sbs"]],
            "signatures/deconstructsigs_sbs_fraction/"
        )
    }
    if ("sigprofiler_indel_fraction" %in% names(signatures_list)) {
        process_signatures(
            signatures_list[["sigprofiler_indel_fraction"]],
            "signatures/sigprofiler_indel_fraction/"
        )
    }
    if ("sigprofiler_indel_count" %in% names(signatures_list)) {
        process_signatures(
            signatures_list[["sigprofiler_indel_count"]],
            "signatures/sigprofiler_indel_count/"
        )
    }
    if ("sigprofiler_sbs_count" %in% names(signatures_list)) {
        process_signatures(
            signatures_list[["sigprofiler_sbs_count"]],
            "signatures/sigprofiler_sbs_count/"
        )
    }
    if ("sigprofiler_sbs_fraction" %in% names(signatures_list)) {
        process_signatures(
            signatures_list[["sigprofiler_sbs_fraction"]],
            "signatures/sigprofiler_sbs_fraction/"
        )
        process_signatures(
            signatures_list[["sigprofiler_sbs_fraction"]],
            "signatures/sbs/"
        )
    }
}

#' @name write_distributions_to_json
#' @title write_distributions_to_json
#' @description
#'
#' Writes the distributions data table to JSON files
#'
#' @param distributions List of data tables containing distributions
#' @param common_dir path to the common directory in which to write the distribution JSONs
#' @param cores number of cores to use for parallel processing
#' @param overwrite boolean flag to indicate whether to overwrite existing files
#' @return NULL
#' @export
write_distributions_to_json <- function(
    distributions,
    common_dir,
    cores = 1,
    overwrite = FALSE,
    haveSignatures = TRUE) {
    common_folder <- paste0(common_dir, "/")

    if (!file.exists(common_folder)) {
        stop("common_folder does not exist. Make the directory first")
    }

    signature_folder_paths <- paste0(
        common_folder,
        c(
            "signatures/sbs/",
            "signatures/deconstructsigs_sbs_fraction/",
            "signatures/insertionDeletion/",
            "signatures/sigprofiler_indel_fraction/",
            "signatures/sigprofiler_indel_count/",
            "signatures/sigprofiler_sbs_count/",
            "signatures/sigprofiler_sbs_fraction/"
        )
    )

    if (!all(file.exists(signature_folder_paths))) {
        empty.lst <- mclapply(signature_folder_paths, function(folder_path) {
            cmd <- paste0("mkdir -p ", folder_path)
            print(paste0("Making directory ", folder_path))
            system(cmd)
            return(NULL)
        }, mc.cores = cores)
    }

    write_json_file <- function(data, file_path) {
        if (overwrite || !file.exists(file_path)) {
            write_json(data, file_path, pretty = TRUE)
        } else {
            message(paste0("File ", file_path, " already exists. Skipping."))
        }
    }

    message(paste0("writing jsons to ", common_folder))
    write_json_file(distributions$snvCount, paste0(common_folder, "/snvCount.json"))
    write_json_file(distributions$svCount, paste0(common_folder, "/svCount.json"))
    write_json_file(distributions$lohFraction, paste0(common_folder, "/lohFraction.json"))
    write_json_file(distributions$ploidy, paste0(common_folder, "/ploidy.json"))
    write_json_file(distributions$purity, paste0(common_folder, "/purity.json"))
    write_json_file(distributions$coverageVariance, paste0(common_folder, "/coverageVariance.json"))
    write_json_file(distributions$tmb, paste0(common_folder, "/tmb.json"))

    if (haveSignatures) {
        write_signature_jsons(
            signatures_list = distributions,
            common_folder = common_folder,
            cores = cores
        )
    }
}

#' @name read_meta_data_json
#' @title read_meta_data_json
#' @description
#' Reads in meta data jsons and gets signatures that are lists into the correct list format in the data.table. Used in create_distributions
#'
#' @param meta_json path to metadata.json for case reports
#' @param patient_id sample name for metadata_json
#'
#' @return data.table of the metadata.json for each sample
#' @export
#' @author Stanley Clarke
read_meta_data_json <- function(meta_json, patient_id) {
    json.dt <- as.data.table(jsonlite::read_json(meta_json, simplifyVector = TRUE)) ## have to change from this to work with signatures
    keep_cols <- names(json.dt) %>% grep("signatures|deletionInsertion|sigprofiler|deconstructsigs", ., value = TRUE, invert = TRUE)
    json.dt <- json.dt[, ..keep_cols]
    json.lst <- jsonlite::read_json(meta_json, simplifyVector = FALSE)
    if (!is.null(json.lst[[1]]$signatures)) {
        json.dt$signatures <- list(json.lst[[1]]$signatures)
    }
    if (!is.null(json.lst[[1]]$deletionInsertion)) {
        json.dt$deletionInsertion <- list(json.lst[[1]]$deletionInsertion)
    }
    if (!is.null(json.lst[[1]]$deconstructsigs_sbs_fraction)) {
        json.dt$deconstructsigs_sbs_fraction <- list(json.lst[[1]]$deconstructsigs_sbs_fraction)
    }
    if (!is.null(json.lst[[1]]$sigprofiler_indel_fraction)) {
        json.dt$sigprofiler_indel_fraction <- list(json.lst[[1]]$sigprofiler_indel_fraction)
    }
    if (!is.null(json.lst[[1]]$sigprofiler_indel_count)) {
        json.dt$sigprofiler_indel_count <- list(json.lst[[1]]$sigprofiler_indel_count)
    }
    if (!is.null(json.lst[[1]]$sigprofiler_sbs_fraction)) {
        json.dt$sigprofiler_sbs_fraction <- list(json.lst[[1]]$sigprofiler_sbs_fraction)
    }
    if (!is.null(json.lst[[1]]$sigprofiler_sbs_count)) {
        json.dt$sigprofiler_sbs_count <- list(json.lst[[1]]$sigprofiler_sbs_count)
    }
    return(json.dt)
}

#' @name convert_signature_meta_json
#' @title convert_signature_meta_json
#' @description
#' converts signatures from the list objects in the meta data to data.tables to write jsons for distributions. Used in create_distributions
#'
#' @param jsons.dt json dt after using read_meta_data_json
#' @param cores cores for generating each signature
#' @return list object with all of the signature distribution data.tables
#' @export
#' @author Stanley Clarke
convert_signature_meta_json <- function(jsons.dt, cores = 1) {
    ## signatures-sbs deconstruct sigs ## in column signatures but also in deconstructsigs_sbs_fraction
    if ("signatures" %in% names(jsons.dt)) {
        jsons.dt2 <- copy(jsons.dt)
        jsons.dt2[, length_signatures := sapply(signatures, function(x) length(unlist(x)))]
        jsons.dt2 <- jsons.dt2[length_signatures != 0, ]
        sigs.lst <- mclapply(1:nrow(jsons.dt2), function(x) {
            jsons.sub.dt2 <- jsons.dt2[x, ]
            sigs.dt <- jsons.sub.dt2$signatures[[1]] %>% as.data.table()
            sigs.dt[, pair := jsons.sub.dt2$pair]
            sigs.dt[, tumor_type := jsons.sub.dt2$tumor_type]
            return(sigs.dt)
        }, mc.cores = cores)
        deconstruct_sigs.dt <- rbindlist(sigs.lst, fill = TRUE)
    }
    ## signatures-indels sigprofiler assignment## in column deletionInsertion and sigprofiler_indel_fraction
    if ("deletionInsertion" %in% names(jsons.dt)) {
        jsons.dt2 <- copy(jsons.dt)
        jsons.dt2[, length_signatures := sapply(deletionInsertion, function(x) length(unlist(x)))]
        jsons.dt2 <- jsons.dt2[length_signatures != 0, ]
        sigs.lst <- mclapply(1:nrow(jsons.dt2), function(x) {
            jsons.sub.dt2 <- jsons.dt2[x, ]
            sigs.dt <- jsons.sub.dt2$deletionInsertion[[1]] %>% as.data.table()
            sigs.dt[, pair := jsons.sub.dt2$pair]
            sigs.dt[, tumor_type := jsons.sub.dt2$tumor_type]
            return(sigs.dt)
        }, mc.cores = cores)
        sigprofilerassignment_indels.dt <- rbindlist(sigs.lst, fill = TRUE)
    }
    ## signatures-indel counts sigprofiler assignment
    if ("sigprofiler_indel_count" %in% names(jsons.dt)) {
        jsons.dt2 <- copy(jsons.dt)
        jsons.dt2[, length_signatures := sapply(sigprofiler_indel_count, function(x) length(unlist(x)))]
        jsons.dt2 <- jsons.dt2[length_signatures != 0, ]
        sigs.lst <- mclapply(1:nrow(jsons.dt2), function(x) {
            jsons.sub.dt2 <- jsons.dt2[x, ]
            sigs.dt <- jsons.sub.dt2$sigprofiler_indel_count[[1]] %>% as.data.table()
            sigs.dt[, pair := jsons.sub.dt2$pair]
            sigs.dt[, tumor_type := jsons.sub.dt2$tumor_type]
            return(sigs.dt)
        }, mc.cores = cores)
        sigprofilerassignment_indels_counts.dt <- rbindlist(sigs.lst, fill = TRUE)
    } else {
        sigprofilerassignment_indels_counts.dt <- NULL
    }
    ## signatures-sigprofiler sbs fraction
    if ("sigprofiler_sbs_fraction" %in% names(jsons.dt)) {
        jsons.dt2 <- copy(jsons.dt)
        jsons.dt2[, length_signatures := sapply(sigprofiler_sbs_fraction, function(x) length(unlist(x)))]
        jsons.dt2 <- jsons.dt2[length_signatures != 0, ]
        sigs.lst <- mclapply(1:nrow(jsons.dt2), function(x) {
            jsons.sub.dt2 <- jsons.dt2[x, ]
            sigs.dt <- jsons.sub.dt2$sigprofiler_sbs_fraction[[1]] %>% as.data.table()
            sigs.dt[, pair := jsons.sub.dt2$pair]
            sigs.dt[, tumor_type := jsons.sub.dt2$tumor_type]
            return(sigs.dt)
        }, mc.cores = cores)
        sigprofilerassignment_sbs.dt <- rbindlist(sigs.lst, fill = TRUE)
    } else {
        sigprofilerassignment_sbs.dt <- NULL
    }

    ## signatures-sigprofiler sbs count
    if ("sigprofiler_sbs_count" %in% names(jsons.dt)) {
        jsons.dt2 <- copy(jsons.dt)
        jsons.dt2[, length_signatures := sapply(sigprofiler_sbs_count, function(x) length(unlist(x)))]
        jsons.dt2 <- jsons.dt2[length_signatures != 0, ]
        sigs.lst <- mclapply(1:nrow(jsons.dt2), function(x) {
            jsons.sub.dt2 <- jsons.dt2[x, ]
            sigs.dt <- jsons.sub.dt2$sigprofiler_sbs_count[[1]] %>% as.data.table()
            sigs.dt[, pair := jsons.sub.dt2$pair]
            sigs.dt[, tumor_type := jsons.sub.dt2$tumor_type]
            return(sigs.dt)
        }, mc.cores = cores)
        sigprofilerassignment_sbs_count.dt <- rbindlist(sigs.lst, fill = TRUE)
    } else {
        sigprofilerassignment_sbs_count.dt <- NULL
    }
    ## create list object to return
    list_return <- list(deconstruct_sigs.dt, sigprofilerassignment_indels.dt, sigprofilerassignment_indels_counts.dt, sigprofilerassignment_sbs.dt, sigprofilerassignment_sbs_count.dt)
    names(list_return) <- c("deconstruct_sigs.dt", "sigprofilerassignment_indels.dt", "sigprofilerassignment_indels_counts.dt", "sigprofilerassignment_sbs.dt", "sigprofilerassignment_sbs_counts.dt")
    return(list_return)
}
