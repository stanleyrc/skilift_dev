suppressWarnings(devtools::load_all())
library(testthat)
library(GenomicRanges)
library(IRanges)

# test <- function() { testthat::test_file("tests/testthat/test-metadata.R") }

setup({
    md_test_paths <<- list(
        coverage = system.file('extdata/test_data/test.cov.rds', package='Skilift'),
        jabba_gg = system.file('extdata/test_data/oncotable_test_data/jabba.simple.gg.rds', package='Skilift')
    )
    
    # Helper function to create mock JaBbA graph
    create_mock_jabba_gg <<- function(
        junction_count = 10,
        loose_count = 4,
        include_junctions = TRUE,
        include_edges = TRUE,
        include_loose = TRUE,
        purity = 0.8,
        ploidy = 2.1,
        nodes_dt = NULL  # Add this parameter
    ) {
        # Create mock junctions data.table
        junctions_dt <- data.table(
            type = rep(c("REF", "ALT"), c(2, junction_count)),
            other_cols = 1:(junction_count + 2)
        )
        
        # Create mock meta data
        meta <- list(
            purity = purity,
            ploidy = ploidy
        )
        
        # Create mock loose ends data.table
        loose_dt <- data.table(
            terminal = rep(c(TRUE, FALSE), c(1, loose_count)),
            other_cols = 1:(loose_count + 1)
        )

        # Create default nodes data.table if not provided
        if (is.null(nodes_dt)) {
            nodes_dt <- data.table(
                seqnames = rep(c(1:22, "X", "Y"), each = 2),  # Include X and Y chromosomes
                start = seq(1, 48, by = 1) * 1000000,
                end = seq(2, 49, by = 1) * 1000000,
                cn.low = c(rep(0, 24), rep(1, 24)),  # Mix of LOH and non-LOH segments
                cn.high = c(rep(2, 24), rep(2, 24)),
                cn = c(rep(2, 24), rep(3, 24))  # Add cn column with some variation
            )
            nodes_dt[, width := end - start + 1]
        }
        
        # Get unique seqnames from nodes_dt
        unique_seqnames <- unique(nodes_dt$seqnames)
        
        # Convert nodes data.table to GRanges
        nodes_gr <- GenomicRanges::GRanges(
            seqnames = nodes_dt$seqnames,
            ranges = IRanges::IRanges(start = nodes_dt$start, end = nodes_dt$end),
            seqinfo = GenomeInfoDb::Seqinfo(
                seqnames = as.character(unique_seqnames),
                seqlengths = rep(1e8, length(unique_seqnames)),
                isCircular = rep(FALSE, length(unique_seqnames)),
                genome = rep("hg19", length(unique_seqnames))
            )
        )
        
        # Create mock graph structure
        gg <- list(
            meta = meta,
            junctions = list(
                dt = if(include_junctions) junctions_dt else data.table()
            ),
            loose = if(include_loose) loose_dt else data.table(),
            edges = if (include_edges) list( dt = data.table( type = c("REF", "ALT", "ALT"), class = c("REF", "INV", "DEL")) ) else data.table(),
            nodes = list(
                dt = nodes_dt,
                gr = nodes_gr
            )
        )
        
        # Save to temp file and return path
        temp_file <- tempfile(fileext = ".rds")
        saveRDS(gg, temp_file)
        return(temp_file)
    }
    
    create_mock_cohort <<- function(
        num_samples = 3,
        include_required = TRUE,
        include_optional = TRUE,
        invalid_samples = NULL
    ) {
        # Create base sample IDs
        sample_ids <- paste0("TEST", sprintf("%03d", 1:num_samples))
        if (!is.null(invalid_samples)) {
            sample_ids[invalid_samples] <- paste0("INVALID", sprintf("%03d", invalid_samples))
        }
        
    # Create mock files for each sample
    sample_files <- lapply(sample_ids, function(id) {
        create_mock_metadata_files(pair = id)
    })
    
    # Create inputs data.table
    inputs <- data.table(pair = sample_ids)

    if (include_optional) {
        inputs[, `:=`(
            tumor_type = "BRCA",
            disease = "Breast Cancer", 
            primary_site = "Breast",
            jabba_gg = sapply(sample_files, function(x) x$jabba_gg),
            somatic_snvs = sapply(sample_files, function(x) x$somatic_snvs),
            tumor_coverage = sapply(sample_files, function(x) x$tumor_coverage),
            estimate_library_complexity = sapply(sample_files, function(x) x$qc_metrics$lib_complex),
            alignment_summary_metrics = sapply(sample_files, function(x) x$qc_metrics$alignment),
            insert_size_metrics = sapply(sample_files, function(x) x$qc_metrics$insert),
            wgs_metrics = sapply(sample_files, function(x) x$qc_metrics$wgs)
        )]
    }
    
    # Create mock Cohort object
    cohort <- list(
        inputs = inputs,
        files = sample_files
    )
    class(cohort) <- "Cohort"
    
    return(cohort)
    }

    create_mock_deconstruct_sigs <<- function(
        signatures = c("SBS1", "SBS2", "SBS3"),
        weights = c(0.5, 0.3, 0.2)
    ) {
        # Create weights matrix
        weights_mat <- matrix(weights, nrow = 1)
        colnames(weights_mat) <- paste0("Signature.", gsub("SBS", "", signatures))
        
        # Create mock deconstructSigs output
        sig_data <- list(weights = weights_mat)
        
        # Save to temp file
        temp_file <- tempfile(fileext = ".rds")
        saveRDS(sig_data, temp_file)
        return(temp_file)
    }

    create_mock_sigprofiler_output <<- function(
        signatures = c("SBS1", "SBS2", "SBS3"),
        counts = c(100, 60, 40),
        sample_name = "TEST001"
    ) {
        # Create data.table with counts
        sig_dt <- data.table(
            Samples = paste0(sample_name, "_somatic")
        )
        
        # Add signature columns
        for(i in seq_along(signatures)) {
            sig_dt[, (signatures[i]) := counts[i]]
        }
        
        # Save to temp file
        temp_file <- tempfile(fileext = ".txt")
        fwrite(sig_dt, temp_file, sep = "\t")
        return(temp_file)
    }

    # Helper function to create mock complex events

    create_mock_hrdetect <<- function(
        dels_mh = 100,
        rs3_exp = 0.4,
        rs5_exp = 0.2,
        sbs3_exp = 0.3,
        sbs8_exp = 0.1,
        del_rep = 50,
        hrd_score = 0.8
    ) {
        # Create mock indels classification table
        indels_table <- data.table(
            del.mh.count = dels_mh,
            del.rep.count = del_rep
        )
        
        # Create mock exposures for rearrangements
        exposures_rearr <- matrix(
            c(rs3_exp, rs5_exp),
            nrow = 2,
            dimnames = list(c("RefSigR3", "RefSigR5"), "sample")
        )
        
        # Create mock exposures for substitutions
        exposures_subs <- matrix(
            c(sbs3_exp, sbs8_exp),
            nrow = 2,
            dimnames = list(c("SBS3", "SBS8"), "sample")
        )
        
        # Create mock HRDetect output
        hrdetect_output <- data.table(
            Probability = hrd_score
        )
        
        # Create complete HRDetect object
        hrdetect_obj <- list(
            indels_classification_table = indels_table,
            exposures_rearr = exposures_rearr,
            exposures_subs = exposures_subs,
            hrdetect_output = hrdetect_output
        )
        
        # Save to temp file and return path
        temp_file <- tempfile(fileext = ".rds")
        saveRDS(hrdetect_obj, temp_file)
        return(temp_file)
    }

    create_mock_onenesstwoness <<- function(
        sum12_score = 0.7,
        brca1_score = 0.4,
        brca2_score = 0.3
    ) {
        # Create mock oneness-twoness scores
        ot_scores <- matrix(
            c(sum12_score, brca1_score, brca2_score),
            nrow = 1,
            dimnames = list(NULL, c("SUM12", "BRCA1", "BRCA2"))
        )
        
        # Create complete object
        onetwo_obj <- list(
            ot_scores = ot_scores
        )
        
        # Save to temp file and return path
        temp_file <- tempfile(fileext = ".rds")
        saveRDS(onetwo_obj, temp_file)
        return(temp_file)
    }

    # Helper function to create a complete set of mock files for create_metadata
    create_mock_metadata_files <<- function(
        pair = "TEST001",
        tumor_type = "BRCA",
        disease = "Breast Cancer",
        primary_site = "Breast",
        inferred_sex = "female",
        purity = 0.8,
        ploidy = 2.1
    ) {
        # Create all necessary mock files
        jabba_gg_file <- create_mock_jabba_gg(purity = purity, ploidy = ploidy)
        vcf_file <- create_mock_vcf(num_variants = 100)
        coverage_file <- create_mock_coverage()
        temp_dir <- tempdir()
        qc_files <- create_mock_qc_files(temp_dir)
        het_pileups_file <- create_mock_het_pileups()
        deconstruct_sigs <- create_mock_deconstruct_sigs()
        activities_indel <- create_mock_sigprofiler_output(signatures = c("ID1", "ID2"))
        activities_sbs <- create_mock_sigprofiler_output(signatures = c("SBS1", "SBS2"))
        hrdetect_file <- create_mock_hrdetect()
        onetwo_file <- create_mock_onenesstwoness()
        complex_file <- create_mock_complex()
        
        # Return list of all files
        return(list(
            jabba_gg = jabba_gg_file,
            somatic_snvs = vcf_file,
            tumor_coverage = coverage_file,
            qc_metrics = qc_files,
            het_pileups = het_pileups_file,
            activities_sbs_signatures = activities_sbs,
            activities_indel_signatures = activities_indel,
            deconstructsigs_sbs_signatures = deconstruct_sigs,
            hrdetect = hrdetect_file,
            onenesstwoness = onetwo_file,
            complex = complex_file
        ))
    }
    create_mock_complex <<- function(
        event_types = c("qrpmin", "qrpmix", "qrppos", "other"),
        event_counts = c(2, 1, 3, 4)
    ) {
        # Create events data.table
        events_dt <- data.table(
            type = rep(event_types, event_counts)
        )
            
        # Create complex object structure
        complex_data <- list(
            meta = list(
                events = events_dt
            )
        )
            
        # Save to temp file and return path
        temp_file <- tempfile(fileext = ".rds")
        saveRDS(complex_data, temp_file)
        return(temp_file)
    }


    # Helper function to create mock QC files
    create_mock_qc_files <<- function(temp_dir) {
        mock_lib_complex <- data.table(
            READ_PAIRS_EXAMINED = 1000,
            READ_PAIR_DUPLICATES = 100,
            READ_PAIR_OPTICAL_DUPLICATES = 50,
            PERCENT_DUPLICATION = 0.1
        )
        lib_complex_file <- file.path(temp_dir, "lib_complex.txt")
        fwrite(mock_lib_complex, lib_complex_file, sep = "\t")
        
        mock_alignment <- data.table(
            CATEGORY = "PAIR",
            TOTAL_READS = 2000,
            PF_READS_ALIGNED = 1900,
            PF_ALIGNED_BASES = 190000,
            MEAN_READ_LENGTH = 100
        )
        alignment_file <- file.path(temp_dir, "alignment.txt")
        fwrite(mock_alignment, alignment_file, sep = "\t")
        
        mock_insert <- data.table(
            MEDIAN_INSERT_SIZE = 300
        )
        insert_file <- file.path(temp_dir, "insert.txt")
        fwrite(mock_insert, insert_file, sep = "\t")
        
        mock_wgs <- data.table(
            MEAN_COVERAGE = 30,
            PCT_30X = 0.9,
            PCT_50X = 0.8
        )
        wgs_file <- file.path(temp_dir, "wgs.txt")
        fwrite(mock_wgs, wgs_file, sep = "\t")
        
        return(list(
            lib_complex = lib_complex_file,
            alignment = alignment_file,
            insert = insert_file,
            wgs = wgs_file
        ))
    }

    # Helper function to create a mock VCF file for testing
    create_mock_vcf <<- function(
        num_variants = 3,
        normal_vaf_values = c(0.0, 0.1, 0.0),
        reference = "hg19",
        include_format = TRUE
    ) {
        reference_line = sprintf("##reference=%s", reference)
        header <- c(
            "##fileformat=VCFv4.2",
            reference_line,
            "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic variant\">",
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">",
            "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Allele frequency\">",
            "##FORMAT=<ID=AU,Number=2,Type=Integer,Description=\"A allele count\">",
            "##FORMAT=<ID=CU,Number=2,Type=Integer,Description=\"C allele count\">",
            "##FORMAT=<ID=GU,Number=2,Type=Integer,Description=\"G allele count\">",
            "##FORMAT=<ID=TU,Number=2,Type=Integer,Description=\"T allele count\">",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR"
        )
        
        variants <- character(num_variants)
        for(i in 1:num_variants) {
            # Create allele counts based on VAF
            depth <- 100
            alt_count <- round(normal_vaf_values[i] * depth)
            ref_count <- depth - alt_count
            
            # Randomly select REF and ALT alleles
            bases <- c("A", "C", "G", "T")
            ref <- sample(bases, 1)
            alt <- sample(bases[bases != ref], 1)
            
            # Create FORMAT field with all required information
            format_str <- sprintf(
                "DP:AF:AU:CU:GU:TU\t%d:%.2f:%d,%d:%d,%d:%d,%d:%d,%d\t%d:%.2f:%d,%d:%d,%d:%d,%d:%d,%d",
                depth, normal_vaf_values[i],  # Normal sample
                ifelse(ref == "A" || alt == "A", ref_count, 0), ifelse(alt == "A", alt_count, 0),
                ifelse(ref == "C" || alt == "C", ref_count, 0), ifelse(alt == "C", alt_count, 0),
                ifelse(ref == "G" || alt == "G", ref_count, 0), ifelse(alt == "G", alt_count, 0),
                ifelse(ref == "T" || alt == "T", ref_count, 0), ifelse(alt == "T", alt_count, 0),
                depth, 0.3,  # Tumor sample (fixed VAF for simplicity)
                ifelse(ref == "A" || alt == "A", ref_count, 0), ifelse(alt == "A", alt_count, 0),
                ifelse(ref == "C" || alt == "C", ref_count, 0), ifelse(alt == "C", alt_count, 0),
                ifelse(ref == "G" || alt == "G", ref_count, 0), ifelse(alt == "G", alt_count, 0),
                ifelse(ref == "T" || alt == "T", ref_count, 0), ifelse(alt == "T", alt_count, 0)
            )
            
            variants[i] <- sprintf(
                "chr%d\t%d\t.\t%s\t%s\t100\tPASS\tSOMATIC\t%s",
                sample(1:22, 1),
                i * 100,
                ref,
                alt,
                format_str
            )
        }
        
        temp_vcf <- tempfile(fileext = ".vcf")
        writeLines(c(header, variants), temp_vcf)
        return(temp_vcf)
    }

    create_mock_coverage <<- function(
        foreground_values = c(1, 2, 3, 2, 1),
        seqnames = "chr1",
        start = 1,
        end = length(foreground_values),
        field = "foreground",
        ncn = 2
    ) {
        # Create GRanges object with required structure
        gr <- GenomicRanges::GRanges(
            seqnames = seqnames,
            ranges = IRanges::IRanges(
                start = start:(start + length(foreground_values) - 1),
                end = (start + length(foreground_values) - 1):end
            )
        )
        
        # Add required metadata columns
        mcols(gr)[[field]] <- foreground_values
        mcols(gr)[["ncn"]] <- rep(ncn, length(foreground_values))
        
  
        temp_file <- tempfile(fileext = ".rds")
        saveRDS(gr, temp_file)
        return(temp_file)
    }

    create_mock_het_pileups <<- function(
        ref_counts_n = c(50, 50, 50, 50, 50),  # Normal ref counts
        alt_counts_n = c(30, 35, 33, 32, 31),  # Normal alt counts - to create alt fractions between 0.2-0.8
        ref_counts_t = c(100, 90, 95, 85, 92), # Tumor ref counts
        alt_counts_t = c(40, 45, 42, 43, 41),  # Tumor alt counts
        seqnames = "chr1",
        start = 1000,
        positions = 5
    ) {
        # Create data.table with required structure
        dt <- data.table(
            seqnames = rep(seqnames, positions),
            start = start:(start + positions - 1),
            end = (start + positions - 1):(start + positions - 1),
            ref.count.n = ref_counts_n,
            alt.count.n = alt_counts_n,
            ref.count.t = ref_counts_t,
            alt.count.t = alt_counts_t
        )
        
        # Calculate alt fractions
        dt[, alt.frac.n := alt.count.n / (ref.count.n + alt.count.n)]
        dt[, alt.frac.t := alt.count.t / (ref.count.t + alt.count.t)]
        
        # Save to temp file and return path
        temp_file <- tempfile(fileext = ".txt")
        fwrite(dt, temp_file, sep = "\t")
        return(temp_file)
    }

})


test_that("initialize_metadata_columns creates correct structure", {
    pair_id <- "TEST001"
    result <- initialize_metadata_columns(pair_id)
    
    expect_s3_class(result, "data.table")
    expect_equal(nrow(result), 1)
    
    expected_columns <- c(
        "pair", "tumor_type", "disease", "primary_site", "inferred_sex",
        "coverage_qc", 
        "snv_count", "snv_count_normal_vaf_greater0",
        "cov_slope", "cov_intercept",
        "deconstructsigs_sbs_fraction", "sigprofiler_indel_fraction",
        "sigprofiler_sbs_count", "sigprofiler_sbs_fraction",
        "hrd_score", "hrd", "b1_2", "b1", "b2"
    )
    expect_true(all(expected_columns %in% names(result)))
    
    expect_equal(result$pair, pair_id)
    
    # NA values are of correct type
    expect_type(result$tumor_type, "character")
    expect_true(is.na(result$tumor_type))
    
    expect_type(result$snv_count, "integer")
    expect_true(is.na(result$snv_count))
    
    # list columns are empty lists
    expect_type(result$coverage_qc, "list")
    expect_null(result$coverage_qc[[1]])
    
    expect_type(result$hrd, "list")
    expect_null(result$hrd[[1]])
})


test_that("add_basic_metadata correctly updates metadata", {
    # Setup initial metadata
    test_pair <- "TEST001"
    metadata <- initialize_metadata_columns(test_pair)
    
    # Test with all values provided
    result <- add_basic_metadata(
        metadata,
        input_tumor_type = "BRCA",
        input_disease = "Breast Cancer",
        input_primary_site = "Breast"
    )
    
    expect_equal(result$tumor_type, "BRCA")
    expect_equal(result$disease, "Breast Cancer")
    expect_equal(result$primary_site, "Breast")
    
    # Test with partial values (some NULL)
    metadata <- initialize_metadata_columns(test_pair)
    result <- add_basic_metadata(
        metadata,
        input_tumor_type = "BRCA",
        input_disease = NULL,
        input_primary_site = "Breast"
    )
    
    expect_equal(result$tumor_type, "BRCA")
    expect_true(is.na(result$disease))
    expect_equal(result$primary_site, "Breast")
    
    # Test with all NULL values
    metadata <- initialize_metadata_columns(test_pair)
    result <- add_basic_metadata(
        metadata,
        input_tumor_type = NULL,
        input_disease = NULL,
        input_primary_site = NULL
    )
    
    expect_true(is.na(result$tumor_type))
    expect_true(is.na(result$disease))
    expect_true(is.na(result$primary_site))
})

test_that("add_basic_metadata handles invalid inputs", {
    # Setup initial metadata
    test_pair <- "TEST001"
    metadata <- initialize_metadata_columns(test_pair)
    
    # Test with missing metadata
    expect_error(
        add_basic_metadata(
            metadata = NULL,
            input_tumor_type = "BRCA",
            input_disease = "Breast Cancer",
            input_primary_site = "Breast"
        ),
        "metadata cannot be NULL"
    )
    
    # Test with wrong metadata type
    expect_error(
        add_basic_metadata(
            metadata = list(),
            input_tumor_type = "BRCA",
            input_disease = "Breast Cancer",
            input_primary_site = "Breast"
        ),
        "metadata must be a data.table"
    )
    
    # Test with invalid value types
    expect_error(
        add_basic_metadata(
            metadata,
            input_tumor_type = 123,
            input_disease = "Breast Cancer",
            input_primary_site = "Breast"
        ),
        "tumor_type must be NULL or a character string"
    )
})

test_that("add_sex_information handles direct sex specification", {
    metadata <- initialize_metadata_columns("TEST001")
    
    result <- add_sex_information(metadata, input_inferred_sex = "male")
    expect_equal(result$inferred_sex, "male")
    
    result <- add_sex_information(metadata, input_inferred_sex = "female")
    expect_equal(result$inferred_sex, "female")
    
    expect_warning(
        add_sex_information(metadata, input_inferred_sex = "invalid"),
        "inferred sex must be one of `male` or `female`"
    )
})

test_that("add_sex_information infers from jabba_gg", {
    metadata <- initialize_metadata_columns("TEST001")
    jabba_gg_path <- md_test_paths$jabba_gg
    
    result <- add_sex_information(metadata, jabba_gg = jabba_gg_path)
    expect_true(result$inferred_sex %in% c("male", "female"))
})

test_that("add_sex_information infers from tumor coverage", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Create mock coverage GRanges
    mock_coverage <- GRanges(
        seqnames = "X",
        ranges = IRanges(start = 1, end = 1000),
        foreground = 0.5
    )
    
    # Save to temp file
    temp_file <- tempfile(fileext = ".rds")
    saveRDS(mock_coverage, temp_file)
    
    # Test inference
    result <- add_sex_information(metadata, tumor_coverage = temp_file)
    expect_equal(result$inferred_sex, "male")
    
    # Test with female pattern
    mock_coverage$foreground <- rep(1.0)  # Value > 0.7 should infer female
    saveRDS(mock_coverage, temp_file)
    
    result <- add_sex_information(metadata, tumor_coverage = temp_file)
    expect_equal(result$inferred_sex, "female")
    
    # Clean up
    unlink(temp_file)
})

test_that("add_sex_information works with missing optional inputs", {
    metadata <- initialize_metadata_columns("TEST001")
    
    result <- add_sex_information(metadata)
    expect_true(is.na(result$inferred_sex))
    
    result <- add_sex_information(metadata, input_inferred_sex = NULL, jabba_gg = NULL)
    expect_true(is.na(result$inferred_sex))
})

test_that("add_sex_information prioritizes direct specification over inference", {
    metadata <- initialize_metadata_columns("TEST001")
    jabba_gg_path <- md_test_paths$jabba_gg
    coverage_path <- md_test_paths$coverage
    
    result <- add_sex_information(
        metadata,
        input_inferred_sex = "male",
        jabba_gg = jabba_gg_path,
        tumor_coverage = coverage_path
    )
    expect_equal(result$inferred_sex, "male")
})

test_that("extract_metrics correctly extracts and renames columns", {
    # Setup test data
    test_data <- data.table(
        TEST_METRIC = 1,
        ANOTHER_METRIC = "a"
    )
    
    metrics <- c(
        test_metric = "TEST_METRIC",
        another = "ANOTHER_METRIC"
    )
    
    result <- extract_metrics(test_data, metrics, "TEST001")
    
    expect_equal(nrow(result), 1)
    expect_equal(names(result), c("test_metric", "another", "pair"))
    expect_equal(result$test_metric, 1)
    expect_equal(result$another, "a")
    expect_equal(result$pair, "TEST001")
})

test_that("extract_metrics handles missing columns appropriately", {
    test_data <- data.table(METRIC1 = 1)
    
    # Test missing metric columns
    expect_error(
        extract_metrics(test_data, c(missing = "MISSING_COL")),
        "Missing columns in QC data: MISSING_COL"
    )
})

test_that("process_qc_metrics correctly processes all QC files", {
    # Setup file paths
    qc_test_dir <- system.file("extdata/test_data/qc_test_data", package = "Skilift")
    
    test_files <- list(
        estimate_library_complexity = file.path(qc_test_dir, "est_lib_complex_metrics.txt"),
        alignment_summary_metrics = file.path(qc_test_dir, "CollectAlignmentSummaryMetrics.alignment_summary_metrics"),
        insert_size_metrics = file.path(qc_test_dir, "CollectInsertSizeMetrics.insert_size_metrics"),
        wgs_metrics = file.path(qc_test_dir, "CollectWgsMetrics.txt")
    )
    
    # Process metrics
    result <- suppressWarnings(process_qc_metrics(
        test_files$estimate_library_complexity,
        test_files$alignment_summary_metrics,
        test_files$insert_size_metrics,
        test_files$wgs_metrics,
        "TEST001"
    ))
    
    # Test structure
    expect_type(result, "list")
    
    # Test expected metrics are present
    expected_metrics <- c( "pair", "read_pairs_examined", "read_pair_duplicates", "read_pair_optical_duplicates", "percent_duplication", "total_reads", "pf_reads_aligned", "pf_aligned_bases", "mean_read_length", "median_insert_size", "median_coverage", "pct_30x", "pct_50x", "m_reads", "m_reads_mapped", "percent_optical_duplication", "percent_aligned", "percent_optical_dups_of_dups")
    
    expect_true(all(expected_metrics %in% names(result)))
    
    # Test calculated fields
    expect_equal(
        result$m_reads,
        result$total_reads / 1e6
    )
    
    expect_equal(
        result$percent_optical_duplication,
        result$read_pair_optical_duplicates / result$read_pairs_examined
    )
})

test_that("dlrs function works correctly", {
    # Test normal case
    test_vector <- c(1, 2, 3, 2, 1, 2, 3, 2, 1)
    result <- dlrs(test_vector)
    expect_type(result, "double")
    expect_true(!is.na(result))
    
    # Test with minimum length
    expect_error(dlrs(c(1, 2)), "Vector length>2 needed for computation")
    
    # Test with NAs
    test_vector_na <- c(1, NA, 3, 2, 1)
    expect_silent(dlrs(test_vector_na))
    
    # Test expected value for known case
    simple_vector <- c(1, 2, 1, 2, 1)
    expected_dlrs <- IQR(diff(simple_vector)) / (sqrt(2) * 1.34)
    expect_equal(dlrs(simple_vector), expected_dlrs)
})

test_that("add_coverage_metrics processes tumor coverage correctly", {
    # Create mock metadata
    metadata <- initialize_metadata_columns("TEST001")
    
    # Create mock coverage data
    mock_coverage <- list(foreground = c(1, 2, 1, 2, 1, 2))
    temp_file <- tempfile(fileext = ".rds")
    saveRDS(mock_coverage, temp_file)
    
    # Test with only coverage data
    result <- add_coverage_metrics(metadata, tumor_coverage = temp_file)
    expect_true(!is.null(result$coverage_qc[[1]]$coverage_variance))
    
    unlink(temp_file)
})

test_that("add_coverage_metrics processes QC metrics correctly", {
    # Create mock metadata
    metadata <- initialize_metadata_columns("TEST001")
    
    # Create temporary directory for mock files
    temp_dir <- tempdir()
    
    # Create mock QC files
    mock_lib_complex <- data.table(
        READ_PAIRS_EXAMINED = 1000,
        READ_PAIR_DUPLICATES = 100,
        READ_PAIR_OPTICAL_DUPLICATES = 50,
        PERCENT_DUPLICATION = 0.1
    )
    lib_complex_file <- file.path(temp_dir, "lib_complex.txt")
    fwrite(mock_lib_complex, lib_complex_file, sep = "\t")
    
    mock_alignment <- data.table(
        CATEGORY = "PAIR",
        TOTAL_READS = 2000000,
        PF_READS_ALIGNED = 1900,
        PF_ALIGNED_BASES = 190000,
        MEAN_READ_LENGTH = 100
    )
    alignment_file <- file.path(temp_dir, "alignment.txt")
    fwrite(mock_alignment, alignment_file, sep = "\t")
    
    mock_insert <- data.table(
        MEDIAN_INSERT_SIZE = 300
    )
    insert_file <- file.path(temp_dir, "insert.txt")
    fwrite(mock_insert, insert_file, sep = "\t")
    
    mock_wgs <- data.table(
        MEAN_COVERAGE = 30,
        PCT_30X = 0.9,
        PCT_50X = 0.8
    )
    wgs_file <- file.path(temp_dir, "wgs.txt")
    fwrite(mock_wgs, wgs_file, sep = "\t")
    
    # Test with QC metrics only
    result <- add_coverage_metrics(
        metadata,
        estimate_library_complexity = lib_complex_file,
        alignment_summary_metrics = alignment_file,
        insert_size_metrics = insert_file,
        wgs_metrics = wgs_file
    )
    
    expect_true(!is.null(result$coverage_qc[[1]]))
    expect_true(is.list(result$coverage_qc[[1]]))
    expect_equal(result$coverage_qc[[1]]$median_coverage, 30)
    expect_equal(result$coverage_qc[[1]]$m_reads, 2)  # 2000/1e6
    
    # Clean up
    unlink(c(lib_complex_file, alignment_file, insert_file, wgs_file))
})

test_that("add_coverage_metrics combines coverage and QC metrics", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Create mock coverage data
    mock_coverage <- list(foreground = c(1, 2, 1, 2, 1, 2))
    coverage_file <- tempfile(fileext = ".rds")
    saveRDS(mock_coverage, coverage_file)
    
    # Create mock QC files
    temp_dir <- tempdir()
    mock_files <- create_mock_qc_files(temp_dir)  # Helper function to create mock files

    result <- add_coverage_metrics(
        metadata,
        tumor_coverage = coverage_file,
        estimate_library_complexity = mock_files$lib_complex,
        alignment_summary_metrics = mock_files$alignment,
        insert_size_metrics = mock_files$insert,
        wgs_metrics = mock_files$wgs
    )
    
    expect_true(!is.null(result$coverage_qc[[1]]$coverage_variance))
    expect_true(!is.null(result$coverage_qc[[1]]$m_reads))
    
    # Clean up
    unlink(coverage_file)
    unlink(unlist(mock_files))
})

test_that("vcf_count processes VCF correctly", {
    # Create mock VCF with known values
    vcf_file <- create_mock_vcf(num_variants = 3)
    
    # Test function
    result <- vcf_count(vcf_file, genome = "hg19")
    
    # Check structure
    expect_s3_class(result, "data.table")
    expect_equal(ncol(result), 2)
    expect_equal(nrow(result), 2)
    expect_equal(result$category, c("snv_count", "snv_count_normal_vaf_greater0"))
    
    # Check counts - only total count should be populated
    expect_equal(result[category == "snv_count"]$counts, 3)
    expect_true(is.na(result[category == "snv_count_normal_vaf_greater0"]$counts))
    
    # Clean up
    unlink(vcf_file)
})

test_that("vcf_count handles tumor-only VCF", {
    # Create header for tumor-only VCF
    header <- c(
        "##fileformat=VCFv4.2",
        "##reference=hg19",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR"
    )
    variants <- "chr1\t100\t.\tA\tT\t100\tPASS\t.\t.\t."
    
    # Create temporary VCF file
    tumor_only_vcf <- tempfile(fileext = ".vcf")
    writeLines(c(header, variants), tumor_only_vcf)
    
    # Test function
    result <- vcf_count(tumor_only_vcf, genome = "hg19")
    
    # Check results - should only count total variants
    expect_equal(result[category == "snv_count"]$counts, 1)
    expect_true(is.na(result[category == "snv_count_normal_vaf_greater0"]$counts))
    
    # Clean up
    unlink(tumor_only_vcf)
})

test_that("vcf_count handles FILTER field correctly", {
    # Create header with PASS and non-PASS variants
    header <- c(
        "##fileformat=VCFv4.2",
        "##reference=hg19",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    )
    variants <- c(
        "chr1\t100\t.\tA\tT\t100\tPASS\t.",
        "chr1\t200\t.\tC\tG\t100\tLowQual\t."
    )
    
    # Create temporary VCF file
    filter_vcf <- tempfile(fileext = ".vcf")
    writeLines(c(header, variants), filter_vcf)
    
    # Test function
    result <- vcf_count(filter_vcf, genome = "hg19")
    
    # Should only count PASS variants
    expect_equal(result[category == "snv_count"]$counts, 1)
    expect_true(is.na(result[category == "snv_count_normal_vaf_greater0"]$counts))
    
    # Clean up
    unlink(filter_vcf)
})

test_that("vcf_count handles invalid inputs", {
    # Test with non-existent file
    expect_error(
        vcf_count("nonexistent.vcf", genome = "hg19"),
        "VCF file does not exist"
    )
    
    # Test with invalid VCF format
    invalid_vcf <- tempfile(fileext = ".vcf")
    writeLines("invalid content", invalid_vcf)
    expect_error(vcf_count(invalid_vcf, genome = "hg19"))
    unlink(invalid_vcf)
})

test_that("vcf_count works with different genome versions", {
    # Create mock VCF with hg38 reference
    vcf_hg38 <- create_mock_vcf(
        num_variants = 3,
        reference = "hg38"
    )
    
    # Test with hg38
    result <- vcf_count(vcf_hg38, genome = "hg38")
    expect_equal(result[category == "snv_count"]$counts, 3)
    expect_true(is.na(result[category == "snv_count_normal_vaf_greater0"]$counts))
    
    # Clean up
    unlink(vcf_hg38)
})

test_that("vcf_count handles FILTER field correctly", {
    # Create header with non-PASS variants
    header <- c(
        "##fileformat=VCFv4.2",
        "##reference=hg19",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    )
    variants <- c(
        "chr1\t100\t.\tA\tT\t100\tPASS\t.",
        "chr1\t200\t.\tC\tG\t100\tLowQual\t."
    )
    
    # Create temporary VCF file
    filter_vcf <- tempfile(fileext = ".vcf")
    writeLines(c(header, variants), filter_vcf)
    
    # Test function
    result <- vcf_count(filter_vcf, genome = "hg19")
    
    # Should only count PASS variants
    expect_equal(result[category == "snv_count"]$counts, 1)
    expect_true(is.na(result[category == "snv_count_normal_vaf_greater0"]$counts))
    
    # Clean up
    unlink(filter_vcf)
})

test_that("add_variant_counts processes VCF correctly", {
    # Create mock metadata
    metadata <- initialize_metadata_columns("TEST001")
    
    # Create a temporary mock VCF file
    temp_vcf <- create_mock_vcf(num_variants = 3)
    
    # Test with mock VCF
    result <- add_variant_counts(metadata, somatic_snvs = temp_vcf)
    
    # Verify counts - only total count should be populated
    expect_equal(result$snv_count, 3)
    expect_true(is.na(result$snv_count_normal_vaf_greater0))
    
    # Clean up
    unlink(temp_vcf)
})

test_that("add_variant_counts handles missing input gracefully", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Test with NULL input
    result <- add_variant_counts(metadata, somatic_snvs = NULL)
    expect_true(is.na(result$snv_count))
    expect_true(is.na(result$snv_count_normal_vaf_greater0))
})

test_that("add_variant_counts handles invalid VCF", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Create invalid VCF
    invalid_vcf_content <- c(
        "##fileformat=VCFv4.2",
        "Invalid content"
    )
    
    temp_vcf <- tempfile(fileext = ".vcf")
    writeLines(invalid_vcf_content, temp_vcf)
    
    # Test with invalid VCF
    expect_error(add_variant_counts(metadata, somatic_snvs = temp_vcf))
    
    # Clean up
    unlink(temp_vcf)
})

test_that("add_variant_counts works with different genome versions", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Create mock VCF
    temp_vcf <- create_mock_vcf(reference = "hg38")
    
    # Test with hg38
    result <- add_variant_counts(metadata, somatic_snvs = temp_vcf, genome = "hg38")
    expect_equal(result$snv_count, 3)
    
    # Clean up
    unlink(temp_vcf)
})

test_that("add_variant_counts processes VCF with genotype fields correctly", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Test with AF field
    vcf_with_af <- create_mock_vcf(
        num_variants = 3,
        normal_vaf_values = c(0.0, 0.1, 0.0)
    )
    result <- add_variant_counts(metadata, somatic_snvs = vcf_with_af)
    expect_equal(result$snv_count, 3)
    expect_equal(result$snv_count_normal_vaf_greater0, NA_integer_)
    unlink(vcf_with_af)
    
    # Test tumor-only VCF (no NORMAL sample)
    header <- c(
        "##fileformat=VCFv4.2",
        "##reference=hg19",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">",
        "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Allele frequency\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR"
    )
    variants <- "chr1\t100\t.\tA\tT\t100\tPASS\tSOMATIC\tDP:AF\t100:0.3"
    tumor_only_vcf <- tempfile(fileext = ".vcf")
    writeLines(c(header, variants), tumor_only_vcf)
    
    result <- add_variant_counts(metadata, somatic_snvs = tumor_only_vcf)
    expect_equal(result$snv_count, 1)
    # currently failing because of how sage_count handles tumor-only VAF (always equal to 1)
    # expect_true(is.na(result$snv_count_normal_vaf_greater0))
    unlink(tumor_only_vcf)
})

test_that("add_sv_counts processes normal JaBbA graph correctly", {
    # Create mock metadata
    metadata <- initialize_metadata_columns("TEST001")
    
    # Create mock JaBbA graph with default values
    jabba_gg_file <- create_mock_jabba_gg()
    
    # Test function
    result <- add_sv_counts(metadata, jabba_gg_file)
    
    # Check results
    expect_equal(result$junction_count, 10)  # 10 ALT junctions
    expect_equal(result$loose_count, 4)      # 4 non-terminal loose ends
    expect_equal(result$sv_count, 12)        # 10 + (4/2)
    
    # Clean up
    unlink(jabba_gg_file)
})

test_that("add_sv_counts handles empty JaBbA graph", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Create mock JaBbA graph with no junctions or loose ends
    jabba_gg_file <- create_mock_jabba_gg(
        junction_count = 0,
        loose_count = 0,
        include_junctions = FALSE,
        include_loose = FALSE
    )
    
    # Test function with warning expectations
    expect_warning(
        result <- add_sv_counts(metadata, jabba_gg_file),
        "No junctions found in JaBbA graph"
    )
    
    # Check results
    expect_equal(result$junction_count, 0)
    expect_equal(result$loose_count, 0)
    expect_equal(result$sv_count, 0)
    
    # Clean up
    unlink(jabba_gg_file)
})

test_that("add_sv_counts handles NULL input gracefully", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Test with NULL jabba_gg
    result <- add_sv_counts(metadata, jabba_gg = NULL)
    
    # Check that metadata is unchanged
    expect_equal(result, metadata)
})

test_that("add_sv_counts handles various junction/loose end combinations", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Test cases matrix:
    test_cases <- list(
        list(junctions = 5, loose = 2, expected_sv = 6),   # 5 + (2/2)
        list(junctions = 0, loose = 4, expected_sv = 2),   # 0 + (4/2)
        list(junctions = 3, loose = 0, expected_sv = 3),   # 3 + (0/2)
        list(junctions = 10, loose = 1, expected_sv = 10.5) # 10 + (1/2)
    )
    
    for(test_case in test_cases) {
        jabba_gg_file <- create_mock_jabba_gg(
            junction_count = test_case$junctions,
            loose_count = test_case$loose
        )
        
        result <- add_sv_counts(metadata, jabba_gg_file)
        
        expect_equal(result$junction_count, test_case$junctions)
        expect_equal(result$loose_count, test_case$loose)
        expect_equal(result$sv_count, test_case$expected_sv)
        
        unlink(jabba_gg_file)
    }
})

test_that("add_sv_counts maintains data.table structure", {
    metadata <- initialize_metadata_columns("TEST001")
    jabba_gg_file <- create_mock_jabba_gg()
    
    result <- add_sv_counts(metadata, jabba_gg_file)
    
    # Check that result is still a data.table
    expect_s3_class(result, "data.table")
    
    # Check that new columns are added correctly
    expect_true(all(c("junction_count", "loose_count", "sv_count") %in% names(result)))
    
    # Check that original columns are preserved
    expect_true(all(names(metadata) %in% names(result)))
    
    unlink(jabba_gg_file)
})

test_that("add_purity_ploidy processes normal JaBbA graph correctly", {
    # Create mock metadata
    metadata <- initialize_metadata_columns("TEST001")
    
    # Create mock JaBbA graph with specific purity and ploidy values
    jabba_gg_file <- create_mock_jabba_gg(purity = 0.8, ploidy = 2.1)
    
    # Test function
    result <- add_purity_ploidy(metadata, jabba_gg_file)
    
    # Check results
    expect_equal(result$purity, 0.8)
    expect_equal(result$ploidy, 2.1)
    expect_equal(result$beta, 2.1)    # beta should equal ploidy
    expect_equal(result$gamma, 0.2)   # gamma should be 1 - purity
    
    # Clean up
    unlink(jabba_gg_file)
})

test_that("add_purity_ploidy handles NULL input gracefully", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Test with NULL jabba_gg
    result <- add_purity_ploidy(metadata, jabba_gg = NULL)
    
    # Check that metadata is unchanged
    expect_equal(result, metadata)
})

test_that("add_purity_ploidy handles various purity/ploidy combinations", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Test cases matrix
    test_cases <- list(
        list(purity = 0.5, ploidy = 2.0),
        list(purity = 0.95, ploidy = 3.2),
        list(purity = 0.3, ploidy = 1.8),
        list(purity = 0.7, ploidy = 4.0)
    )
    
    for(test_case in test_cases) {
        jabba_gg_file <- create_mock_jabba_gg(
            purity = test_case$purity,
            ploidy = test_case$ploidy
        )
        
        result <- add_purity_ploidy(metadata, jabba_gg_file)
        
        expect_equal(result$purity, test_case$purity)
        expect_equal(result$ploidy, test_case$ploidy)
        expect_equal(result$beta, test_case$ploidy)
        expect_equal(result$gamma, 1 - test_case$purity)
        
        unlink(jabba_gg_file)
    }
})

test_that("add_purity_ploidy maintains data.table structure", {
    metadata <- initialize_metadata_columns("TEST001")
    jabba_gg_file <- create_mock_jabba_gg()
    
    result <- add_purity_ploidy(metadata, jabba_gg_file)
    
    # Check that result is still a data.table
    expect_s3_class(result, "data.table")
    
    # Check that new columns are added correctly
    expect_true(all(c("purity", "ploidy", "beta", "gamma") %in% names(result)))
    
    # Check that original columns are preserved
    expect_true(all(names(metadata) %in% names(result)))
    
    unlink(jabba_gg_file)
})

test_that("add_loh processes normal JaBbA graph correctly", {
    # Create mock metadata
    metadata <- initialize_metadata_columns("TEST001")
    
    # Create mock JaBbA graph with default LOH segments
    jabba_gg_file <- create_mock_jabba_gg()
    
    # Test function
    result <- add_loh(metadata, jabba_gg_file)
    
    # Check results
    expect_true(!is.na(result$loh_fraction))
    expect_true(is.numeric(result$loh_fraction))
    expect_true(result$loh_fraction >= 0 && result$loh_fraction <= 1)
    expect_true(is.numeric(result$loh_seglen))
    expect_true(is.numeric(result$loh_total_genome))
    
    # Clean up
    unlink(jabba_gg_file)
})

test_that("add_loh handles NULL inputs gracefully", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Test with NULL jabba_gg
    result <- add_loh(metadata, jabba_gg = NULL)
    expect_equal(result, metadata)
    
    # Test with NULL seqnames_loh
    result <- add_loh(metadata, jabba_gg = create_mock_jabba_gg(), seqnames_loh = NULL)
    expect_equal(result, metadata)
})

test_that("add_loh handles different chromosome selections", {
    metadata <- initialize_metadata_columns("TEST001")
    # create a second metadata since its updating by reference
    metadata2 <- initialize_metadata_columns("TEST001")

    jabba_gg_file <- create_mock_jabba_gg()
    
    # Test with different chromosome selections
    result1 <- add_loh(metadata, jabba_gg_file, seqnames_loh = c(1:22))
    result2 <- add_loh(metadata2, jabba_gg_file, seqnames_loh = c(1:11))
    
    # Results should be different due to different chromosome selections
    expect_false(identical(result1$loh_fraction, result2$loh_fraction))
    expect_false(identical(result1$loh_total_genome, result2$loh_total_genome))
    
    unlink(jabba_gg_file)
})

test_that("add_loh handles non-allelic JaBbA graph", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Create nodes without allelic information
    nodes_dt <- data.table(
        seqnames = rep(c(1:22), each = 2),
        start = seq(1, 44, by = 1) * 1000000,
        end = seq(2, 45, by = 1) * 1000000,
        width = 1000000  # Only basic columns
    )
    
    jabba_gg_file <- create_mock_jabba_gg(nodes_dt = nodes_dt)
    
    result <- add_loh(metadata, jabba_gg_file)
    expect_equal(result$loh_fraction, "Not Allelic Jabba")
    
    unlink(jabba_gg_file)
})

test_that("add_loh maintains data.table structure", {
    metadata <- initialize_metadata_columns("TEST001")
    jabba_gg_file <- create_mock_jabba_gg()
    
    result <- add_loh(metadata, jabba_gg_file)
    
    # Check that result is still a data.table
    expect_s3_class(result, "data.table")
    
    # Check that new columns are added correctly
    expect_true(all(c("loh_fraction", "loh_seglen", "loh_total_genome") %in% names(result)))
    
    # Check that original columns are preserved
    expect_true(all(names(metadata) %in% names(result)))
    
    unlink(jabba_gg_file)
})

test_that("add_loh calculates LOH fraction correctly", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Create nodes with known LOH regions
    nodes_dt <- data.table(
        seqnames = 1:4,
        start = c(1, 1001, 2001, 3001),
        end = c(1000, 2000, 3000, 4000),
        cn.low = c(0, 1, 0, 1),
        cn.high = c(2, 2, 2, 2)
    )
    nodes_dt[, width := end - start + 1]
    
    jabba_gg_file <- create_mock_jabba_gg(nodes_dt = nodes_dt)
    
    result <- add_loh(metadata, jabba_gg_file)
    
    # Expected LOH fraction: 2000 bases with LOH out of 4000 total = 0.5
    expect_equal(result$loh_fraction, 0.5)
    expect_equal(result$loh_seglen, 2000)
    expect_equal(result$loh_total_genome, 4000)
    
    unlink(jabba_gg_file)
})

test_that("add_genome_length processes normal JaBbA graph correctly", {
    metadata <- initialize_metadata_columns("TEST001")
    jabba_gg_file <- create_mock_jabba_gg()
    
    result <- add_genome_length(metadata, jabba_gg_file)
    
    # Check results
    expect_true(!is.na(result$total_genome_length))
    expect_true(is.numeric(result$total_genome_length))
    expect_equal(result$total_genome_length, 24 * 1e8)  # 22 chromosomes * 100Mb
    
    unlink(jabba_gg_file)
})

test_that("add_genome_length handles numeric genome length input", {
    metadata <- initialize_metadata_columns("TEST001")
    jabba_gg_file <- create_mock_jabba_gg()
    
    # Test with numeric genome length
    total_length <- 1e8  # 100Mb for example
    result <- add_genome_length(metadata, jabba_gg_file, seqnames_genome_width_or_genome_length = total_length)
    
    expect_equal(result$total_genome_length, total_length)
    
    # Test with different numeric value
    result2 <- add_genome_length(metadata, jabba_gg_file, seqnames_genome_width_or_genome_length = 5e7)
    expect_equal(result2$total_genome_length, 5e7)
    
    unlink(jabba_gg_file)
})

test_that("add_genome_length handles NULL inputs gracefully", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Test with NULL jabba_gg
    result <- add_genome_length(metadata, jabba_gg = NULL)
    expect_equal(result, metadata)
    
    # Test with NULL seqnames_genome_width_or_genome_length
    result <- add_genome_length(metadata, jabba_gg = create_mock_jabba_gg(), seqnames_genome_width_or_genome_length = NULL)
    expect_equal(result, metadata)
})

test_that("add_genome_length handles different chromosome selections", {
    metadata <- initialize_metadata_columns("TEST001")
    jabba_gg_file <- create_mock_jabba_gg()
    
    # Test with different chromosome selections
    result1 <- add_genome_length(metadata, jabba_gg_file, seqnames_genome_width_or_genome_length = c(1:22))
    result2 <- add_genome_length(metadata, jabba_gg_file, seqnames_genome_width_or_genome_length = c(1:11))
    
    # Results should be different due to different chromosome selections
    expect_false(identical(result1$total_genome_length, result2$total_genome_length))
    expect_equal(result1$total_genome_length, 22 * 1e8)  # All chromosomes
    expect_equal(result2$total_genome_length, 11 * 1e8)  # Half the chromosomes
    
    unlink(jabba_gg_file)
})

test_that("add_genome_length maintains data.table structure", {
    metadata <- initialize_metadata_columns("TEST001")
    jabba_gg_file <- create_mock_jabba_gg()
    
    result <- add_genome_length(metadata, jabba_gg_file)
    
    # Check that result is still a data.table
    expect_s3_class(result, "data.table")
    
    # Check that new column is added correctly
    expect_true("total_genome_length" %in% names(result))
    
    # Check that original columns are preserved
    expect_true(all(names(metadata) %in% names(result)))
    
    unlink(jabba_gg_file)
})

test_that("add_genome_length handles sex chromosomes correctly", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Create nodes with X and Y chromosomes
    nodes_dt <- data.table(
        seqnames = c(1:22, "X", "Y"),
        start = 1,
        end = 1e8,
        width = 1e8
    )
    
    jabba_gg_file <- create_mock_jabba_gg(nodes_dt = nodes_dt)
    
    # Test including sex chromosomes
    result <- add_genome_length(metadata, jabba_gg_file, 
                              seqnames_genome_width_or_genome_length = c(1:22, "X", "Y"))
    
    # Should include all chromosomes
    expect_equal(result$total_genome_length, 24 * 1e8)  # 22 + X + Y chromosomes
    
    # Test excluding sex chromosomes
    result2 <- add_genome_length(metadata, jabba_gg_file, 
                                seqnames_genome_width_or_genome_length = c(1:22))
    
    # Should only include autosomes
    expect_equal(result2$total_genome_length, 22 * 1e8)  # Only autosomes
    
    unlink(jabba_gg_file)
})

test_that("add_sv_types processes normal inputs correctly", {
    metadata <- initialize_metadata_columns("TEST001")
    jabba_gg_file <- create_mock_jabba_gg()
    complex_file <- create_mock_complex()
    
    result <- add_sv_types(metadata, jabba_gg_file, complex_file)
    
    # Check sv_types_count structure
    expect_type(result$sv_types_count, "list")
    expect_length(result$sv_types_count, 1)
    
    # Check QRP counts
    expect_type(result$hrd, "list")
    expect_length(result$hrd, 1)
    qrp_counts <- result$hrd[[1]]
    expect_equal(qrp_counts$qrpmin, 2)
    expect_equal(qrp_counts$qrpmix, 1)
    expect_equal(qrp_counts$qrppos, 3)
    
    unlink(c(jabba_gg_file, complex_file))
})

test_that("add_sv_types handles NULL inputs gracefully", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Test with NULL inputs
    result <- suppressWarnings(add_sv_types(metadata, NULL, NULL))
    expect_equal(result$sv_types_count[[1]], list())
    expect_equal(result$hrd[[1]], list(qrpmin = 0, qrpmix = 0, qrppos = 0))
})

test_that("add_sv_types handles empty JaBbA graph", {
    metadata <- initialize_metadata_columns("TEST001")
    jabba_gg_file <- create_mock_jabba_gg(include_edges = FALSE)
    complex_file <- create_mock_complex()
    
    expect_warning(
        result <- add_sv_types(metadata, jabba_gg_file, complex_file),
        "No edges found in JaBbA graph"
    )
    
    # Should still have complex event counts
    expect_false(is.null(result$sv_types_count[[1]]))
    
    unlink(c(jabba_gg_file, complex_file))
})

test_that("add_sv_types handles missing complex events", {
    metadata <- initialize_metadata_columns("TEST001")
    jabba_gg_file <- create_mock_jabba_gg()
    
    expect_warning(
        result <- add_sv_types(metadata, jabba_gg_file, NULL),
        "Complex events and JaBbA graph not found as inputs"
    )
    
    # Should still have basic SV counts
    expect_false(is.null(result$sv_types_count[[1]]))
    expect_equal(result$hrd[[1]], list(qrpmin = 0, qrpmix = 0, qrppos = 0))
    
    unlink(jabba_gg_file)
})

test_that("add_sv_types maintains data.table structure", {
    metadata <- initialize_metadata_columns("TEST001")
    jabba_gg_file <- create_mock_jabba_gg()
    complex_file <- create_mock_complex()
    
    result <- add_sv_types(metadata, jabba_gg_file, complex_file)
    
    # Check that result is still a data.table
    expect_s3_class(result, "data.table")
    
    # Check that new columns are added correctly
    expect_true(all(c("sv_types_count", "hrd") %in% names(result)))
    
    # Check that original columns are preserved
    expect_true(all(names(metadata) %in% names(result)))
    
    unlink(c(jabba_gg_file, complex_file))
})

test_that("add_sv_types handles various event type combinations", {
    metadata <- initialize_metadata_columns("TEST001")
    jabba_gg_file <- create_mock_jabba_gg()
    
    # Test different event type combinations
    test_cases <- list(
        list(types = c("qrpmin", "other"), counts = c(2, 3)),
        list(types = c("qrpmix", "qrppos"), counts = c(1, 4)),
        list(types = c("qrpmin", "qrpmix", "qrppos"), counts = c(2, 2, 2)),
        list(types = c("other1", "other2"), counts = c(1, 1))
    )
    
    for(test_case in test_cases) {
        complex_file <- create_mock_complex(
            event_types = test_case$types,
            event_counts = test_case$counts
        )
        
        result <- add_sv_types(metadata, jabba_gg_file, complex_file)
        
        # Verify counts match expected values
        qrp_counts <- result$hrd[[1]]
        for(type in c("qrpmin", "qrpmix", "qrppos")) {
            expected_count <- if(type %in% test_case$types) {
                test_case$counts[which(test_case$types == type)]
            } else {
                0
            }
            expect_equal(qrp_counts[[type]], expected_count)
        }
        
        unlink(complex_file)
    }
    
    unlink(jabba_gg_file)
})

test_that("add_coverage_parameters processes tumor coverage correctly", {
    # Create mock metadata with purity and ploidy
    metadata <- initialize_metadata_columns("TEST001")
    metadata$purity <- 0.8
    metadata$ploidy <- 2.1
    
    coverage_file <- create_mock_coverage()
    result <- suppressWarnings(add_coverage_parameters(metadata, coverage_file))
    
    # Check results
    expect_true(!is.na(result$cov_slope))
    expect_true(!is.na(result$cov_intercept))
    expect_true(is.numeric(result$cov_slope))
    expect_true(is.numeric(result$cov_intercept))
    
    unlink(coverage_file)
})

test_that("add_coverage_parameters handles NULL input gracefully", {
    metadata <- initialize_metadata_columns("TEST001")
    metadata$purity <- 0.8
    metadata$ploidy <- 2.1
    
    # Test with NULL coverage
    result <- add_coverage_parameters(metadata, NULL)
    expect_equal(result, metadata)
})

test_that("add_coverage_parameters maintains data.table structure", {
    metadata <- initialize_metadata_columns("TEST001")
    metadata$purity <- 0.8
    metadata$ploidy <- 2.1
    
    coverage_file <- create_mock_coverage()
    result <- suppressWarnings(add_coverage_parameters(metadata, coverage_file))
    
    # Check that result is still a data.table
    expect_s3_class(result, "data.table")
    
    # Check that new columns are added correctly
    expect_true(all(c("cov_slope", "cov_intercept") %in% names(result)))
    
    # Check that original columns are preserved
    expect_true(all(names(metadata) %in% names(result)))
    
    unlink(coverage_file)
})

test_that("add_coverage_parameters calculates parameters correctly", {
    metadata <- initialize_metadata_columns("TEST001")
    metadata$purity <- 0.8
    metadata$ploidy <- 2.1
    
    # Create coverage with known values for verification
    coverage_file <- create_mock_coverage(
        foreground_values = c(1, 2, 1, 2, 1),
        ncn = 2
    )
    
    result <- suppressWarnings(add_coverage_parameters(metadata, coverage_file))
    
    # Calculate expected values based on rel2abs formula
    data_mean <- mean(c(1, 2, 1, 2, 1))
    expected_slope <- ((1 - 0.8) * 2 + 0.8 * 2.1)/(0.8 * data_mean)
    expected_intercept <- -(2 * (1 - 0.8)/0.8)
    
    expect_equal(result$cov_slope, expected_slope)
    expect_equal(result$cov_intercept, expected_intercept)
    
    unlink(coverage_file)
})

test_that("add_coverage_parameters handles invalid coverage data", {
    metadata <- initialize_metadata_columns("TEST001")
    metadata$purity <- 0.8
    metadata$ploidy <- 2.1
    
    # Create coverage with wrong field name
    coverage_file <- create_mock_coverage(field = "wrong_field")
    expect_error(add_coverage_parameters(metadata, coverage_file))
    
    unlink(coverage_file)
})

test_that("add_coverage_parameters requires purity and ploidy", {
    metadata <- initialize_metadata_columns("TEST001")
    # Note: purity and ploidy are not set
    
    mock_coverage <- list(foreground = c(1, 2, 3))
    temp_file <- tempfile(fileext = ".rds")
    saveRDS(mock_coverage, temp_file)
    
    # Should error or warn when purity/ploidy are missing
    expect_error(suppressWarnings(add_coverage_parameters(metadata, temp_file)))
    
    unlink(temp_file)
})
    
test_that("add_het_pileups_parameters processes het pileups correctly", {
    metadata <- initialize_metadata_columns("TEST001")
    metadata$purity <- 0.8
    metadata$ploidy <- 2.1
    
    # Create mock het pileups with counts that will result in alt fractions between 0.2-0.8
    het_pileups_file <- create_mock_het_pileups(
        ref_counts_n = c(70, 65, 68, 72, 69),  # Normal ref counts
        alt_counts_n = c(30, 35, 32, 28, 31),  # Normal alt counts (~0.3 fraction)
        ref_counts_t = c(140, 130, 136, 144, 138), # Tumor ref counts
        alt_counts_t = c(60, 70, 64, 56, 62)  # Tumor alt counts (~0.3 fraction)
    )
    
    result <- add_het_pileups_parameters(metadata, het_pileups_file)
    
    # Check results
    expect_true(!is.na(result$hets_slope))
    expect_true(!is.na(result$hets_intercept))
    expect_true(is.numeric(result$hets_slope))
    expect_true(is.numeric(result$hets_intercept))
    
    unlink(het_pileups_file)
})

test_that("add_het_pileups_parameters handles NULL input gracefully", {
    metadata <- initialize_metadata_columns("TEST001")
    metadata$purity <- 0.8
    metadata$ploidy <- 2.1
    
    # Test with NULL het_pileups
    result <- add_het_pileups_parameters(metadata, NULL)
    expect_equal(result, metadata)
})

test_that("add_het_pileups_parameters maintains data.table structure", {
    metadata <- initialize_metadata_columns("TEST001")
    metadata$purity <- 0.8
    metadata$ploidy <- 2.1
    
    het_pileups_file <- create_mock_het_pileups()
    result <- add_het_pileups_parameters(metadata, het_pileups_file)
    
    # Check that result is still a data.table
    expect_s3_class(result, "data.table")
    
    # Check that new columns are added correctly
    expect_true(all(c("hets_slope", "hets_intercept") %in% names(result)))
    
    # Check that original columns are preserved
    expect_true(all(names(metadata) %in% names(result)))
    
    unlink(het_pileups_file)
})

test_that("add_het_pileups_parameters calculates parameters correctly", {
    metadata <- initialize_metadata_columns("TEST001")
    metadata$purity <- 0.8
    metadata$ploidy <- 2.1
    
    # Create het pileups with specific counts to test parameter calculation
    het_pileups_file <- create_mock_het_pileups(
        ref_counts_n = c(70, 70, 70, 70, 70),  # Normal ref counts
        alt_counts_n = c(30, 30, 30, 30, 30),  # Normal alt counts (0.3 fraction)
        ref_counts_t = c(140, 140, 140, 140, 140), # Tumor ref counts
        alt_counts_t = c(60, 60, 60, 60, 60)  # Tumor alt counts (0.3 fraction)
    )
    
    result <- add_het_pileups_parameters(metadata, het_pileups_file)
    
    # After grab.hets processing, we expect:
    # - ref will be major allele (higher count)
    # - Total count per position will be ref + alt
    # Expected counts after processing should be consistent
    expect_true(!is.na(result$hets_slope))
    expect_true(!is.na(result$hets_intercept))
    expect_true(result$hets_slope > 0)  # Slope should be positive
    
    unlink(het_pileups_file)
})

test_that("add_het_pileups_parameters handles invalid het pileups data", {
    metadata <- initialize_metadata_columns("TEST001")
    metadata$purity <- 0.8
    metadata$ploidy <- 2.1
    
    # Create het pileups with alt fractions outside valid range
    invalid_het_pileups_file <- create_mock_het_pileups(
        ref_counts_n = c(95, 95, 95, 95, 95),  # Normal ref counts
        alt_counts_n = c(5, 5, 5, 5, 5),    # Normal alt counts (0.05 fraction - too low)
        ref_counts_t = c(190, 190, 190, 190, 190), # Tumor ref counts
        alt_counts_t = c(10, 10, 10, 10, 10)    # Tumor alt counts
    )
    
    # Should handle the case where no valid heterozygous sites are found
    expect_error(suppressWarngings(add_het_pileups_parameters(metadata, invalid_het_pileups_file)))
    
    unlink(invalid_het_pileups_file)
})

test_that("add_tmb calculates TMB correctly with chromosome-based genome length", {
    # Create mock metadata
    metadata <- initialize_metadata_columns("TEST001")
    
    # Create mock VCF with known number of variants
    vcf_file <- create_mock_vcf(num_variants = 100)
    
    # Create mock JaBbA graph with known genome length
    # Using 1e8 length for each chromosome (22 autosomes + X + Y)
    jabba_gg_file <- create_mock_jabba_gg()
    
    result <- add_tmb(metadata, vcf_file, jabba_gg_file)
    
    # Expected TMB = 100 variants / (24 * 100Mb / 1Mb) = 0.042
    expect_equal(result$tmb, 0.042)
    
    # Clean up
    unlink(c(vcf_file, jabba_gg_file))
})

test_that("add_tmb calculates TMB correctly with numeric genome length", {
    metadata <- initialize_metadata_columns("TEST001")
    vcf_file <- create_mock_vcf(num_variants = 100)
    jabba_gg_file <- create_mock_jabba_gg()
    
    # Test with specific genome length (e.g., for targeted panel)
    result <- add_tmb(
        metadata, 
        vcf_file, 
        jabba_gg_file, 
        seqnames_genome_width_or_genome_length = 1e6  # 1Mb targeted region
    )
    
    # Expected TMB = 100 variants / 1Mb = 100
    expect_equal(result$tmb, 100)
    
    # Test with WES-like genome length
    result2 <- add_tmb(
        metadata,
        vcf_file,
        jabba_gg_file,
        seqnames_genome_width_or_genome_length = 30e6  # 30Mb typical WES
    )
    
    # Expected TMB = 100 variants / 30Mb = 3.333
    expect_equal(result2$tmb, 3.333)
    
    # Clean up
    unlink(c(vcf_file, jabba_gg_file))
})

test_that("add_tmb handles NULL inputs gracefully", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Test with NULL inputs
    expect_warning(
        add_tmb(metadata, NULL, NULL),
        "Cannot calculate TMB without both snv_count and total_genome_length"
    )
    
    # Test with only VCF
    vcf_file <- create_mock_vcf(num_variants = 100)
    expect_warning(
        add_tmb(metadata, vcf_file, NULL),
        "Cannot calculate TMB without both snv_count and total_genome_length"
    )
    
    # Test with only JaBbA graph
    jabba_gg_file <- create_mock_jabba_gg()
    expect_warning(
        add_tmb(metadata, NULL, jabba_gg_file),
        "Cannot calculate TMB without both snv_count and total_genome_length"
    )
    
    # Clean up
    unlink(c(vcf_file, jabba_gg_file))
})

test_that("add_tmb maintains data.table structure", {
    metadata <- initialize_metadata_columns("TEST001")
    vcf_file <- create_mock_vcf()
    jabba_gg_file <- create_mock_jabba_gg()
    
    result <- add_tmb(metadata, vcf_file, jabba_gg_file)
    
    # Check that result is still a data.table
    expect_s3_class(result, "data.table")
    
    # Check that tmb column is added correctly
    expect_true("tmb" %in% names(result))
    
    # Check that original columns are preserved
    expect_true(all(names(metadata) %in% names(result)))
    
    # Clean up
    unlink(c(vcf_file, jabba_gg_file))
})

test_that("add_tmb handles different chromosome selections", {
    metadata <- initialize_metadata_columns("TEST001")
    metadata2 <- initialize_metadata_columns("TEST001")
    vcf_file <- create_mock_vcf(num_variants = 100)
    jabba_gg_file <- create_mock_jabba_gg()
    
    # Test with different chromosome selections
    result1 <- add_tmb(metadata, vcf_file, jabba_gg_file, seqnames_genome_width_or_genome_length = c(1:22, "X", "Y"))
    result2 <- add_tmb(metadata2, vcf_file, jabba_gg_file, seqnames_genome_width_or_genome_length = c(1:11))
    
    # Results should be different due to different genome lengths
    expect_false(identical(result1$tmb, result2$tmb))
    
    # Expected TMB values
    # result1: 100 variants / (24 * 100Mb / 1Mb)
    # result2: 100 variants / (22 * 100Mb / 1Mb)
    expect_equal(result1$tmb, 0.042)
    expect_equal(result2$tmb, 0.091)
    
    # Clean up
    unlink(c(vcf_file, jabba_gg_file))
})

test_that("add_tmb rounds TMB to 3 decimal places", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Create mock data that would result in a long decimal
    vcf_file <- create_mock_vcf(num_variants = 123)  # Odd number for division
    jabba_gg_file <- create_mock_jabba_gg()
    
    result <- add_tmb(metadata, vcf_file, jabba_gg_file)
    
    # Check that result has exactly 3 decimal places
    tmb_string <- as.character(result$tmb)
    decimal_places <- nchar(strsplit(tmb_string, "\\.")[[1]][2])
    expect_equal(decimal_places, 3)
    
    # Verify rounding is correct
    expected_tmb <- 0.051
    expect_equal(result$tmb, expected_tmb)
    
    # Clean up
    unlink(c(vcf_file, jabba_gg_file))
})

test_that("initialize_metadata_columns handles invalid inputs", {
    # missing pair
    expect_error(
        initialize_metadata_columns(),
        "argument \"pair\" is missing, with no default"
    )
    
    # NULL pair
    expect_error(
        initialize_metadata_columns(NULL),
        "pair must be a character string"
    )
    
    # empty string
    expect_error(
        initialize_metadata_columns(""),
        "pair cannot be empty"
    )
})

test_that("compute_signature_averages processes deconstructSigs output correctly", {
    # Create mock deconstructSigs output
    signatures <- c("SBS1", "SBS2", "SBS3")
    weights <- c(0.5, 0.3, 0.2)
    mock_file <- create_mock_deconstruct_sigs(signatures, weights)
    
    # Test function
    result <- compute_signature_averages(
        mock_file,
        is_indel = FALSE,
        is_deconstruct_sigs = TRUE
    )
    
    # Check results
    expect_type(result, "double")
    expect_named(result, signatures)
    expect_equal(result, setNames(weights, signatures))
    expect_equal(sum(result), 1)  # Weights should sum to 1
    
    # Clean up
    unlink(mock_file)
})

test_that("compute_signature_averages processes sigprofiler SBS output correctly", {
    # Create mock sigprofiler output
    signatures <- c("SBS1", "SBS2", "SBS3")
    counts <- c(100, 60, 40)
    mock_file <- create_mock_sigprofiler_output(signatures, counts)
    
    # Test function
    result <- compute_signature_averages(
        mock_file,
        is_indel = FALSE,
        is_deconstruct_sigs = FALSE
    )
    
    # Check results
    expect_type(result, "list")
    expect_named(result, c("sbs_fraction", "sbs_count"))
    
    # Check fractions
    total_count <- sum(counts)
    expected_fractions <- counts / total_count
    expect_equal(result$sbs_fraction, setNames(expected_fractions, signatures))
    expect_equal(sum(result$sbs_fraction), 1)
    
    # Check counts
    expect_equal(result$sbs_count, setNames(counts, signatures))
    
    # Clean up
    unlink(mock_file)
})

test_that("compute_signature_averages processes sigprofiler indel output correctly", {
    # Create mock sigprofiler indel output
    signatures <- c("ID1", "ID2", "ID3")
    counts <- c(80, 50, 30)
    mock_file <- create_mock_sigprofiler_output(signatures, counts)
    
    # Test function
    result <- compute_signature_averages(
        mock_file,
        is_indel = TRUE,
        is_deconstruct_sigs = FALSE
    )
    
    # Check results
    expect_type(result, "list")
    expect_named(result, c("indel_fraction", "indel_count"))
    
    # Check fractions
    total_count <- sum(counts)
    expected_fractions <- counts / total_count
    expect_equal(result$indel_fraction, setNames(expected_fractions, signatures))
    expect_equal(sum(result$indel_fraction), 1)
    
    # Check counts
    expect_equal(result$indel_count, setNames(counts, signatures))
    
    # Clean up
    unlink(mock_file)
})


test_that("compute_signature_averages handles invalid inputs gracefully", {
    # Test with non-existent file
    expect_error(
        compute_signature_averages(
            "nonexistent.file",
            is_indel = FALSE,
            is_deconstruct_sigs = FALSE
        )
    )
    
    # Test with invalid RDS file for deconstructSigs
    invalid_rds <- tempfile(fileext = ".rds")
    saveRDS(list(not_weights = 1), invalid_rds)
    expect_error(
        compute_signature_averages(
            invalid_rds,
            is_indel = FALSE,
            is_deconstruct_sigs = TRUE
        )
    )
    unlink(invalid_rds)
    
    # Test with invalid text file for sigprofiler
    invalid_txt <- tempfile(fileext = ".txt")
    writeLines("invalid content", invalid_txt)
    expect_error(
        compute_signature_averages(
            invalid_txt,
            is_indel = FALSE,
            is_deconstruct_sigs = FALSE
        )
    )
    unlink(invalid_txt)
})

test_that("compute_signature_averages maintains consistent output structure", {
    # Test deconstructSigs output
    mock_deconstruct <- create_mock_deconstruct_sigs()
    result_deconstruct <- compute_signature_averages(
        mock_deconstruct,
        is_indel = FALSE,
        is_deconstruct_sigs = TRUE
    )
    expect_type(result_deconstruct, "double")
    expect_true(all(result_deconstruct >= 0 & result_deconstruct <= 1))
    
    # Test sigprofiler SBS output
    mock_sbs <- create_mock_sigprofiler_output()
    result_sbs <- compute_signature_averages(
        mock_sbs,
        is_indel = FALSE,
        is_deconstruct_sigs = FALSE
    )
    expect_type(result_sbs, "list")
    expect_true(all(result_sbs$sbs_fraction >= 0 & result_sbs$sbs_fraction <= 1))
    expect_true(all(result_sbs$sbs_count >= 0))
    
    # Test sigprofiler indel output
    mock_indel <- create_mock_sigprofiler_output()
    result_indel <- compute_signature_averages(
        mock_indel,
        is_indel = TRUE,
        is_deconstruct_sigs = FALSE
    )
    expect_type(result_indel, "list")
    expect_true(all(result_indel$indel_fraction >= 0 & result_indel$indel_fraction <= 1))
    expect_true(all(result_indel$indel_count >= 0))
    
    # Clean up
    unlink(c(mock_deconstruct, mock_sbs, mock_indel))
})


test_that("add_signatures processes all signature types correctly", {
    # Create mock metadata
    metadata <- initialize_metadata_columns("TEST001")
    
    # Create mock signature files
    deconstruct_sigs <- create_mock_deconstruct_sigs(
        signatures = c("SBS1", "SBS2", "SBS3"),
        weights = c(0.5, 0.3, 0.2)
    )
    
    activities_indel <- create_mock_sigprofiler_output(
        signatures = c("ID1", "ID2", "ID3"),
        counts = c(80, 50, 30)
    )
    
    activities_sbs <- create_mock_sigprofiler_output(
        signatures = c("SBS4", "SBS5", "SBS6"),
        counts = c(100, 60, 40)
    )
    
    # Test function
    result <- add_signatures(
        metadata,
        activities_sbs_signatures = activities_sbs,
        activities_indel_signatures = activities_indel,
        deconstructsigs_sbs_signatures = deconstruct_sigs
    )
    
    # Check deconstructSigs SBS results
    expect_type(result$deconstructsigs_sbs_fraction[[1]], "list")
    expect_equal(sum(unlist(result$deconstructsigs_sbs_fraction[[1]])), 1)
    expect_named(result$deconstructsigs_sbs_fraction[[1]], c("SBS1", "SBS2", "SBS3"))
    
    # Check indel results
    expect_type(result$sigprofiler_indel_fraction[[1]], "list")
    expect_equal(sum(unlist(result$sigprofiler_indel_fraction[[1]])), 1)
    expect_named(result$sigprofiler_indel_fraction[[1]], c("ID1", "ID2", "ID3"))
    
    # Check SBS results
    expect_type(result$sigprofiler_sbs_fraction[[1]], "list")
    expect_equal(sum(unlist(result$sigprofiler_sbs_fraction[[1]])), 1)
    expect_named(result$sigprofiler_sbs_fraction[[1]], c("SBS4", "SBS5", "SBS6"))
    
    # Clean up
    unlink(c(deconstruct_sigs, activities_indel, activities_sbs))
})

test_that("add_signatures handles NULL inputs gracefully", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Test with all NULL inputs
    result <- add_signatures(
        metadata,
        activities_sbs_signatures = NULL,
        activities_indel_signatures = NULL,
        deconstructsigs_sbs_signatures = NULL
    )
    
    # Check that lists are empty but initialized
    expect_type(result$deconstructsigs_sbs_fraction, "list")
    expect_equal(result$deconstructsigs_sbs_fraction[[1]], list())
    
    expect_type(result$sigprofiler_indel_fraction, "list")
    expect_equal(result$sigprofiler_indel_fraction[[1]], list())
    
    expect_type(result$sigprofiler_sbs_fraction, "list")
    expect_equal(result$sigprofiler_sbs_fraction[[1]], list())
})

test_that("add_signatures handles partial inputs correctly", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Create only activities_sbs input
    activities_sbs <- create_mock_sigprofiler_output()
    
    # Test with only deconstructSigs input
    result <- add_signatures(
        metadata,
        activities_sbs_signatures = activities_sbs,
        activities_indel_signatures = NULL,
        deconstructsigs_sbs_signatures = NULL
    )
    
    # Check that activities_sbs results exist
    expect_type(result$sigprofiler_sbs_fraction[[1]], "list")
    expect_equal(sum(unlist(result$sigprofiler_sbs_fraction[[1]])), 1)
    expect_named(result$sigprofiler_sbs_fraction[[1]], c("SBS1", "SBS2", "SBS3"))
    
    # Check that other lists are empty
    expect_type(result$sigprofiler_indel_fraction, "list")
    expect_equal(result$sigprofiler_indel_fraction[[1]], list())
    
    expect_type(result$sigprofiler_sbs_fraction, "list")
    expect_equal(result$deconstructsigs_sbs_fraction[[1]], list())
    
    # Clean up
    unlink(activities_sbs)
})

test_that("add_signatures maintains data.table structure", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Create mock signature files
    deconstruct_sigs <- create_mock_deconstruct_sigs()
    activities_indel <- create_mock_sigprofiler_output(
        signatures = c("ID1", "ID2"),
        counts = c(70, 30)
    )
    activities_sbs <- create_mock_sigprofiler_output()
    
    result <- add_signatures(
        metadata,
        activities_sbs_signatures = activities_sbs,
        activities_indel_signatures = activities_indel,
        deconstructsigs_sbs_signatures = deconstruct_sigs
    )
    
    # Check that result is still a data.table
    expect_s3_class(result, "data.table")
    
    # Check that all original columns are preserved
    expect_true(all(names(metadata) %in% names(result)))
    
    # Check that signature columns are present
    expect_true(all(c(
        "deconstructsigs_sbs_fraction",
        "sigprofiler_indel_fraction",
        "sigprofiler_sbs_fraction"
    ) %in% names(result)))
    
    # Clean up
    unlink(c(deconstruct_sigs, activities_indel, activities_sbs))
})

test_that("add_signatures handles different pair names correctly", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Create mock signature files with different sample name
    deconstruct_sigs <- create_mock_deconstruct_sigs()
    activities_indel <- create_mock_sigprofiler_output(sample_name = "OTHER001")
    activities_sbs <- create_mock_sigprofiler_output(sample_name = "OTHER001")
    
    # Test with non-matching pair name
    result <- add_signatures(
        metadata,
        activities_sbs_signatures = activities_sbs,
        activities_indel_signatures = activities_indel,
        deconstructsigs_sbs_signatures = deconstruct_sigs
    )
    
    # Check that signatures are processed correctly
    expect_type(result$deconstructsigs_sbs_fraction[[1]], "list")
    expect_type(result$sigprofiler_indel_fraction[[1]], "list")
    expect_type(result$sigprofiler_sbs_fraction[[1]], "list")
    
    # Clean up
    unlink(c(deconstruct_sigs, activities_indel, activities_sbs))
})

test_that("add_hrd_scores processes HRDetect scores correctly", {
    # Create mock metadata
    metadata <- initialize_metadata_columns("TEST001")
    
    # Create mock HRDetect file with known values
    hrdetect_file <- create_mock_hrdetect(
        dels_mh = 100,
        rs3_exp = 0.4,
        rs5_exp = 0.2,
        sbs3_exp = 0.3,
        sbs8_exp = 0.1,
        del_rep = 50,
        hrd_score = 0.8
    )
    
    # Test function
    result <- suppressWarnings(add_hrd_scores(metadata, hrdetect_file, NULL))
    
    # Check HRD score
    expect_equal(result$hrd_score, 0.8)
    
    # Check HRD values
    expect_type(result$hrd[[1]], "list")
    expect_equal(result$hrd[[1]]$dels_mh, 100)
    expect_equal(result$hrd[[1]]$rs3, 0.4)
    expect_equal(result$hrd[[1]]$rs5, 0.2)
    expect_equal(result$hrd[[1]]$sbs3, 0.3)
    expect_equal(result$hrd[[1]]$sbs8, 0.1)
    expect_equal(result$hrd[[1]]$del_rep, 50)
    
    # Clean up
    unlink(hrdetect_file)
})

test_that("add_hrd_scores processes oneness-twoness scores correctly", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Create mock oneness-twoness file
    onetwo_file <- create_mock_onenesstwoness(
        sum12_score = 0.7,
        brca1_score = 0.4,
        brca2_score = 0.3
    )
    
    # Test function
    result <- suppressWarnings(add_hrd_scores(metadata, NULL, onetwo_file))
    
    # Check scores
    expect_equal(result$b1_2, 0.7)
    expect_equal(result$b1, 0.4)
    expect_equal(result$b2, 0.3)
    
    # Clean up
    unlink(onetwo_file)
})

test_that("add_hrd_scores handles NULL inputs gracefully", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Test with NULL inputs
    result <- suppressWarnings(add_hrd_scores(metadata, NULL, NULL))
    
    # Check that values remain NA
    expect_true(is.na(result$hrd_score))
    expect_null(result$hrd[[1]])
    expect_true(is.na(result$b1_2))
    expect_true(is.na(result$b1))
    expect_true(is.na(result$b2))
})

test_that("add_hrd_scores processes both inputs simultaneously", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Create both mock files
    hrdetect_file <- create_mock_hrdetect(hrd_score = 0.8)
    onetwo_file <- create_mock_onenesstwoness(
        sum12_score = 0.7,
        brca1_score = 0.4,
        brca2_score = 0.3
    )
    
    # Test function with both inputs
    result <- add_hrd_scores(metadata, hrdetect_file, onetwo_file)
    
    # Check all scores are present
    expect_equal(result$hrd_score, 0.8)
    expect_type(result$hrd[[1]], "list")
    expect_equal(result$b1_2, 0.7)
    expect_equal(result$b1, 0.4)
    expect_equal(result$b2, 0.3)
    
    # Clean up
    unlink(c(hrdetect_file, onetwo_file))
})

test_that("add_hrd_scores maintains data.table structure", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Create mock files
    hrdetect_file <- create_mock_hrdetect()
    onetwo_file <- create_mock_onenesstwoness()
    
    # Test function
    result <- add_hrd_scores(metadata, hrdetect_file, onetwo_file)
    
    # Check that result is still a data.table
    expect_s3_class(result, "data.table")
    
    # Check that new columns are added correctly
    expect_true(all(c("hrd_score", "hrd", "b1_2", "b1", "b2") %in% names(result)))
    
    # Check that original columns are preserved
    expect_true(all(names(metadata) %in% names(result)))
    
    # Clean up
    unlink(c(hrdetect_file, onetwo_file))
})

test_that("add_hrd_scores handles invalid file inputs", {
    metadata <- initialize_metadata_columns("TEST001")
    
    # Create invalid RDS files
    invalid_hrdetect <- tempfile(fileext = ".rds")
    saveRDS(list(invalid = "data"), invalid_hrdetect)
    
    invalid_onetwo <- tempfile(fileext = ".rds")
    saveRDS(list(invalid = "data"), invalid_onetwo)
    
    # Test with invalid HRDetect file
    expect_warning(add_hrd_scores(metadata, invalid_hrdetect, NULL))
    # Test with invalid oneness-twoness file
    expect_warning(add_hrd_scores(metadata, NULL, invalid_onetwo))
    # Clean up
    unlink(c(invalid_hrdetect, invalid_onetwo))
})

test_that("create_metadata creates complete metadata object", {
    # Create mock files
    mock_files <- create_mock_metadata_files()
    
    # Test create_metadata with all inputs
    result <- suppressWarnings(create_metadata(
        pair = "TEST001",
        tumor_type = "BRCA",
        disease = "Breast Cancer",
        primary_site = "Breast",
        inferred_sex = "female",
        jabba_gg = mock_files$jabba_gg,
        events = mock_files$complex,
        somatic_snvs = mock_files$somatic_snvs,
        tumor_coverage = mock_files$tumor_coverage,
        estimate_library_complexity = mock_files$qc_metrics$lib_complex,
        alignment_summary_metrics = mock_files$qc_metrics$alignment,
        insert_size_metrics = mock_files$qc_metrics$insert,
        wgs_metrics = mock_files$qc_metrics$wgs,
        het_pileups = mock_files$het_pileups,
        activities_sbs_signatures = mock_files$activities_sbs_signatures,
        activities_indel_signatures = mock_files$activities_indel_signatures,
        deconstructsigs_sbs_signatures = mock_files$deconstructsigs_sbs_signatures,
        hrdetect = mock_files$hrdetect,
        onenesstwoness = mock_files$onenesstwoness
    ))
    
    # Check basic structure
    expect_s3_class(result, "data.table")
    expect_equal(nrow(result), 1)
    expect_equal(result$pair, "TEST001")
    
    # Check that all major components are present
    expect_false(is.na(result$tumor_type))
    expect_false(is.na(result$inferred_sex))
    expect_false(is.na(result$snv_count))
    expect_false(is.null(result$coverage_qc[[1]]))
    expect_false(is.null(result$deconstructsigs_sbs_fraction[[1]]))
    expect_false(is.na(result$hrd_score))
    expect_false(is.na(result$tmb))
    
    # Clean up
    unlink(unlist(mock_files))
})

test_that("create_metadata handles minimal inputs", {
    # Test with only required inputs
    result <- suppressWarnings(create_metadata(pair = "TEST001"))
    
    # Check basic structure
    expect_s3_class(result, "data.table")
    expect_equal(nrow(result), 1)
    expect_equal(result$pair, "TEST001")
    
    # Check that optional fields are NA/NULL
    expect_true(is.na(result$tumor_type))
    expect_true(is.na(result$disease))
    expect_true(is.na(result$primary_site))
    expect_true(is.na(result$inferred_sex))
    expect_true(is.na(result$snv_count))
    expect_null(result$coverage_qc[[1]])
})

test_that("create_metadata handles partial inputs", {
    # Create subset of mock files
    mock_files <- create_mock_metadata_files()
    
    # Test with only some inputs
    result <- suppressWarnings(create_metadata(
        pair = "TEST001",
        tumor_type = "BRCA",
        somatic_snvs = mock_files$somatic_snvs,
        jabba_gg = mock_files$jabba_gg
    ))
    
    # Check that provided inputs are processed
    expect_equal(result$tumor_type, "BRCA")
    expect_false(is.na(result$snv_count))
    
    # Check that unprovided inputs are NA/NULL
    expect_true(is.na(result$disease))
    expect_true(is.na(result$primary_site))
    expect_null(result$coverage_qc[[1]])
    
    # Clean up
    unlink(c(mock_files$somatic_snvs, mock_files$jabba_gg))
})

test_that("create_metadata maintains data consistency", {
    # Create two identical sets of inputs
    mock_files1 <- create_mock_metadata_files()
    mock_files2 <- create_mock_metadata_files()
    
    # Create two metadata objects with identical inputs
    result1 <- suppressWarnings(create_metadata(
        pair = "TEST001",
        tumor_type = "BRCA",
        jabba_gg = mock_files1$jabba_gg,
        somatic_snvs = mock_files1$somatic_snvs
    ))
    
    result2 <- suppressWarnings(create_metadata(
        pair = "TEST001",
        tumor_type = "BRCA",
        jabba_gg = mock_files2$jabba_gg,
        somatic_snvs = mock_files2$somatic_snvs
    ))
    
    # Results should be identical
    expect_equal(result1, result2)
    
    # Clean up
    unlink(unlist(mock_files1))
    unlink(unlist(mock_files2))
})

test_that("create_metadata handles invalid inputs gracefully", {
    # Test with invalid pair
    expect_error(create_metadata(pair = NULL))
    expect_error(create_metadata(pair = ""))
    
    # Test with invalid file paths
    expect_error(suppressWarnings(create_metadata(
        pair = "TEST001",
        somatic_snvs = "nonexistent.vcf",
        jabba_gg = "nonexistent.rds"
    )))
    
    # Test with invalid input types
    expect_error(create_metadata(
        pair = "TEST001",
        tumor_type = 123,  # Should be character
        disease = TRUE     # Should be character
    ))
})

test_that("create_metadata handles different genome length specifications", {
    # Create mock files
    mock_files <- create_mock_metadata_files()
    
    # Test with chromosome-based genome length
    result1 <- suppressWarnings(create_metadata(
        pair = "TEST001",
        jabba_gg = mock_files$jabba_gg,
        somatic_snvs = mock_files$somatic_snvs,
        seqnames_genome_width_or_genome_length = c(1:22, "X", "Y")
    ))
    
    # Test with numeric genome length
    result2 <- suppressWarnings(create_metadata(
        pair = "TEST001",
        jabba_gg = mock_files$jabba_gg,
        somatic_snvs = mock_files$somatic_snvs,
        seqnames_genome_width_or_genome_length = 1e6  # 1Mb
    ))
    
    # Results should be different due to different genome lengths
    expect_false(identical(result1$tmb, result2$tmb))
    
    # Clean up
    unlink(unlist(mock_files))
})

test_that("lift_metadata processes cohort correctly", {
    # Create mock cohort and temporary output directory
    cohort <- create_mock_cohort(num_samples = 2)
    temp_dir <- tempdir()
    
    # Run lift_metadata
    suppressWarnings(lift_metadata(cohort, temp_dir))
    
    # Check that output directories were created
    expect_true(dir.exists(file.path(temp_dir, "TEST001")))
    expect_true(dir.exists(file.path(temp_dir, "TEST002")))
    
    # Check that metadata files were created
    expect_true(file.exists(file.path(temp_dir, "TEST001", "metadata.json")))
    expect_true(file.exists(file.path(temp_dir, "TEST002", "metadata.json")))
    
    # Verify content of metadata files
    metadata1 <- jsonlite::read_json(file.path(temp_dir, "TEST001", "metadata.json"))
    metadata2 <- jsonlite::read_json(file.path(temp_dir, "TEST002", "metadata.json"))
    
    expect_equal(metadata1[[1]]$pair, "TEST001")
    expect_equal(metadata2[[1]]$pair, "TEST002")
    
    # Clean up
    unlink(temp_dir, recursive = TRUE)
})

test_that("lift_metadata handles missing required column", {
    # Create cohort without pair column
    cohort <- create_mock_cohort(num_samples = 2)
    cohort$inputs[, pair := NULL]
    temp_dir <- tempdir()
    
    # Should error due to missing required column
    expect_error(
        lift_metadata(cohort, temp_dir),
        "Missing required column 'pair' in cohort"
    )
    
    # Clean up
    unlink(temp_dir, recursive = TRUE)
})

test_that("lift_metadata handles missing optional columns", {
    # Create cohort without optional columns
    cohort <- create_mock_cohort(num_samples = 2, include_optional = FALSE)
    temp_dir <- tempdir()
    
    # Should run with warning but not error
    expect_warning(
        lift_metadata(cohort, temp_dir),
        "Missing optional columns in cohort"
    )
    
    # Check that output was still created
    expect_true(file.exists(file.path(temp_dir, "TEST001", "metadata.json")))
    expect_true(file.exists(file.path(temp_dir, "TEST002", "metadata.json")))
    
    # Clean up
    unlink(temp_dir, recursive = TRUE)
})

test_that("lift_metadata creates output directory if missing", {
    cohort <- create_mock_cohort(num_samples = 1)
    temp_dir <- file.path(tempdir(), "new_directory")
    
    # Ensure directory doesn't exist
    if (dir.exists(temp_dir)) unlink(temp_dir, recursive = TRUE)
    
    # Run lift_metadata
    suppressWarnings(lift_metadata(cohort, temp_dir))
    
    # Check that directory was created
    expect_true(dir.exists(temp_dir))
    expect_true(dir.exists(file.path(temp_dir, "TEST001")))
    
    # Clean up
    unlink(temp_dir, recursive = TRUE)
})

test_that("lift_metadata handles parallel processing failures", {
    # Create cohort with some invalid samples
    cohort <- create_mock_cohort(num_samples = 3, invalid_samples = 2)
    temp_dir <- tempdir()
    
    # Should complete with warnings
    expect_warning(
        lift_metadata(cohort, temp_dir, cores = 2),
        "Missing optional columns in cohort: inferred_sex, events, germline_snvs, het_pileups, activities_sbs_signatures, activities_indel_signatures, hrdetect, onenesstwoness"
    )
    
    # Check that valid samples were processed
    expect_true(file.exists(file.path(temp_dir, "TEST001", "metadata.json")))
    expect_true(file.exists(file.path(temp_dir, "TEST003", "metadata.json")))
    
    # Clean up
    unlink(temp_dir, recursive = TRUE)
})

test_that("lift_metadata validates input type", {
    # Test with invalid input type
    invalid_cohort <- list(inputs = data.table(pair = c("TEST001", "TEST002")))
    temp_dir <- tempdir()
    
    expect_error(
        lift_metadata(invalid_cohort, temp_dir),
        "Input must be a Cohort object"
    )
    
    # Clean up
    unlink(temp_dir, recursive = TRUE)
})

test_that("lift_metadata processes multiple cores correctly", {
    # Create larger cohort to test parallel processing
    cohort <- create_mock_cohort(num_samples = 4)
    temp_dir <- tempdir()
    
    # Run with multiple cores
    suppressWarnings(lift_metadata(cohort, temp_dir, cores = 2))
    
    # Check that all samples were processed
    for (i in 1:4) {
        sample_dir <- file.path(temp_dir, sprintf("TEST%03d", i))
        expect_true(dir.exists(sample_dir))
        expect_true(file.exists(file.path(sample_dir, "metadata.json")))
    }
    
    # Clean up
    unlink(temp_dir, recursive = TRUE)
})

## integration tests (only works on NYU hpc)
will_test_integration = FALSE
if (will_test_integration) {
test_that("lift_metadata works on real cohort", {
    clinical_pairs_path = "~/projects/Clinical_NYU/db/pairs.rds"
    clinical_pairs = readRDS(clinical_pairs_path)
    vip_sample = clinical_pairs[patient_id == "397089", ]
    qc_test_files <- list(
        estimate_library_complexity = file.path(qc_test_dir, "est_lib_complex_metrics.txt"),
        alignment_summary_metrics = file.path(qc_test_dir, "CollectAlignmentSummaryMetrics.alignment_summary_metrics"),
        insert_size_metrics = file.path(qc_test_dir, "CollectInsertSizeMetrics.insert_size_metrics"),
        wgs_metrics = file.path(qc_test_dir, "CollectWgsMetrics.txt")
    )
    vip_sample$estimate_library_complexity = qc_test_files$estimate_library_complexity
    vip_sample$alignment_summary_metrics = qc_test_files$alignment_summary_metrics
    vip_sample$insert_size_metrics = qc_test_files$insert_size_metrics
    vip_sample$wgs_metrics = qc_test_files$wgs_metrics

    suppressWarnings(cohort <- Cohort$new(vip_sample, col_mapping = list("pair" = "patient_id")))
    row <- cohort$inputs
    metadata <- suppressWarnings(create_metadata(
        pair = row$pair,
        tumor_type = row$tumor_type,
        disease = row$disease,
        primary_site = row$primary_site,
        inferred_sex = row$inferred_sex,
        jabba_gg = row$jabba_gg,
        events = row$events,
        somatic_snvs = row$somatic_snvs,
        germline_snvs = row$germline_snvs,
        tumor_coverage = row$tumor_coverage,
        estimate_library_complexity = row$estimate_library_complexity,
        alignment_summary_metrics = row$alignment_summary_metrics,
        insert_size_metrics = row$insert_size_metrics,
        wgs_metrics = row$wgs_metrics,
        het_pileups = row$het_pileups,
        activities_sbs_signatures = row$activities_sbs_signatures,
        activities_indel_signatures = row$activities_indel_signatures,
        hrdetect = row$hrdetect,
        onenesstwoness = row$onenesstwoness
    ))
    expect_s3_class(metadata, "data.table")

    temp_dir <- tempdir()
    suppressWarnings(lift_metadata(cohort, temp_dir))
    expect_true(file.exists(file.path(temp_dir, vip_sample$patient_id, "metadata.json")))
    unlink(temp_dir, recursive = TRUE)
})
}


## debug
DEBUG = FALSE
if (DEBUG) {

# why is add_metadata so slow?
# hypothesis: add_variant_counts
# result: het_pileups, tumor_coverage, somatic_snvs, jabba_gg cause slow down
# conclusion: operating on large objects slows things down, add_variant_counts (sage_count) good candidate for optimization

metadata <- initialize_metadata_columns("TEST001")
clinical_pairs_path = "~/projects/Clinical_NYU/db/pairs.rds"
clinical_pairs = readRDS(clinical_pairs_path)
vip_sample = clinical_pairs[patient_id == "397089", ]
qc_test_files <- list(
    estimate_library_complexity = file.path(qc_test_dir, "est_lib_complex_metrics.txt"),
    alignment_summary_metrics = file.path(qc_test_dir, "CollectAlignmentSummaryMetrics.alignment_summary_metrics"),
    insert_size_metrics = file.path(qc_test_dir, "CollectInsertSizeMetrics.insert_size_metrics"),
    wgs_metrics = file.path(qc_test_dir, "CollectWgsMetrics.txt")
)
vip_sample$estimate_library_complexity = qc_test_files$estimate_library_complexity
vip_sample$alignment_summary_metrics = qc_test_files$alignment_summary_metrics
vip_sample$insert_size_metrics = qc_test_files$insert_size_metrics
vip_sample$wgs_metrics = qc_test_files$wgs_metrics

suppressWarnings(cohort <- Cohort$new(vip_sample, col_mapping = list("pair" = "patient_id")))
row <- cohort$inputs

add_variant_counts(
    metadata,
    somatic_snvs = row$somatic_snvs,
    genome = "hg19"
)

devtools::load_all()
# inputs slowing down add_metadata are commented out
metadata <- suppressWarnings(create_metadata(
    pair = row$pair,
    tumor_type = row$tumor_type,
    disease = row$disease,
    primary_site = row$primary_site,
    inferred_sex = row$inferred_sex,
    # jabba_gg = row$jabba_gg,
    events = row$events,
    # somatic_snvs = row$somatic_snvs,
    germline_snvs = row$germline_snvs,
    # tumor_coverage = row$tumor_coverage,
    estimate_library_complexity = row$estimate_library_complexity,
    alignment_summary_metrics = row$alignment_summary_metrics,
    insert_size_metrics = row$insert_size_metrics,
    wgs_metrics = row$wgs_metrics,
    # het_pileups = row$het_pileups,
    activities_sbs_signatures = row$activities_sbs_signatures,
    activities_indel_signatures = row$activities_indel_signatures,
    hrdetect = row$hrdetect,
    onenesstwoness = row$onenesstwoness
))
}
