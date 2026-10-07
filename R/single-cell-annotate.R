## Driver, fusion and signature annotation of single-cell patients for gOS.
##
## Runs the lab's nf-gos tools from their containers on hg38:
##   - SnpEff 5.1 with the GRCh38.105 cache (nf-gos VCF_SNPEFF)
##   - the OncoKB annotator in mskilab/unified (nf-gos ONCOKB_ANNOTATOR). SNVs are
##     queried by protein change with -r GRCh38, since the nf-gos wrapper's
##     vcf2maf --inhibit-vep step drops SnpEff's HGVSp and hard-codes GRCh37;
##     fusions and copy-number changes go through nf-gos process_singularity.sh
##   - gGnome fusions via nf-gos bin/Fusions.R in mskilab/unified:0.0.10
##     (nf-gos FUSIONS), with an hg38 "5putr.expanded" GENCODE
##   - SigProfilerAssignment, mskilab/sigprofilerassignment (COSMIC v3.4, GRCh38)
## and skilift's own oncotable() / create_filtered_events() for the drivers tab.
##
## Paths come from sc_annotation_config(); override any of them there.

#' @name sc_annotation_config
#' @title sc_annotation_config
#' @description
#'
#' Paths and settings used by the single-cell annotation functions (nf-gos
#' checkout, container images, SnpEff cache, OncoKB gene list, hg38 GENCODE
#' files, OncoTree code). Pass named values to override the lab defaults.
#'
#' The two hg38 GENCODE files are the hg38 versions of nf-gos's hg19 files:
#' cna_gencode is oncokb/gencode.v19.annotation.gtf.nochr.tsg.onco.rds (the same
#' 399 OncoKB oncogenes / tumor suppressors on GENCODE v32), fusions_gencode is
#' fusions/hg19/gencode.v29lift37.annotation.nochr.5putr.expanded.rds rebuilt
#' on GENCODE v32 (transcripts with a 5' UTR get a 150 kb upstream exon "-1"
#' and UTR, as in the hg19 file).
#'
#' @param ... values to override
#' @return named list
#' @export
#' @author Stanley Clarke
sc_annotation_config <- function(...) {
    images <- "/gpfs/commons/groups/imielinski_lab/data/pipeline/container_images_cache"
    mski_base <- "/gpfs/commons/groups/imielinski_lab/data/nf_gos_files/mskilab_pipeline"
    reference <- "/gpfs/commons/groups/imielinski_lab/projects/single_cell_GBM/gos_datasets/gOS_GBM_sc/reference"
    defaults <- list(
        nf_gos = "/gpfs/commons/home/sclarke/git/nf-gos",
        singularity = "/nfs/sw/easybuild/software/singularity/4.2.2/bin/singularity",
        mski_base = mski_base,
        snpeff_image = file.path(images, "depot.galaxyproject.org-singularity-snpeff-5.1--hdfd78af_2.img"),
        snpeff_cache = file.path(mski_base, "snpeff_cache"),
        snpeff_db = "GRCh38.105",
        unified_image = file.path(images, "mskilab-unified-0.0.11.img"),
        fusions_image = file.path(images, "mskilab-unified-0.0.10.img"),
        sigprofiler_image = file.path(images, "mskilab-sigprofilerassignment-0.0.4.img"),
        oncokb_genes = file.path(mski_base, "oncokb/OncoKB_cancer_genes.rds"),
        oncokb_secrets = "~/.nextflow/secrets/store.json",
        gencode = "/gpfs/commons/groups/imielinski_lab/DB/GENCODE/hg38/v32/gencode.v32.annotation.nochr.rds",
        cna_gencode = file.path(reference, "gencode.v32.annotation.nochr.tsg.onco.rds"),
        fusions_gencode = file.path(reference, "gencode.v32.annotation.nochr.5putr.expanded.rds"),
        oncotree_code = "GBM"
    )
    cfg <- utils::modifyList(defaults, list(...))
    cfg$nf_gos_real <- normalizePath(cfg$nf_gos, mustWork = FALSE)   # containers need the path behind symlinks
    cfg
}

## OncoKB API key from the Nextflow secrets store (as nf-gos uses it)
sc_oncokb_token <- function(cfg) {
    secrets <- jsonlite::fromJSON(cfg$oncokb_secrets)
    token <- secrets$value[secrets$name == "ONCOKB_API_KEY"]
    if (!length(token) || !nzchar(token)) stop("ONCOKB_API_KEY not found in ", cfg$oncokb_secrets)
    token
}

## run a command in a container; secrets (env) go in through SINGULARITYENV_*
sc_singularity_exec <- function(cfg, image, command, binds = character(0), env = character(0)) {
    if (length(env)) {
        keys <- paste0("SINGULARITYENV_", names(env))
        do.call(Sys.setenv, as.list(stats::setNames(env, keys)))
        on.exit(Sys.unsetenv(keys))
    }
    args <- c("exec", "--cleanenv", unlist(lapply(unique(binds), function(b) c("-B", b))), image, "bash", "-c", shQuote(command))
    status <- system2(cfg$singularity, args)
    if (status != 0) stop("container command failed (", status, "): ", substr(command, 1, 200))
    invisible(TRUE)
}

## p.Arg132His -> p.R132H
SC_AA3 <- c(Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C", Gln = "Q", Glu = "E", Gly = "G", His = "H",
            Ile = "I", Leu = "L", Lys = "K", Met = "M", Phe = "F", Pro = "P", Ser = "S", Thr = "T", Trp = "W",
            Tyr = "Y", Val = "V", Ter = "*", Sec = "U", Pyl = "O", Xaa = "X")
sc_hgvsp_short <- function(x) {
    out <- x
    for (k in names(SC_AA3)) out <- gsub(k, SC_AA3[[k]], out, fixed = TRUE)
    out
}

## SnpEff effect -> MAF Variant_Classification
sc_snpeff_class <- function(effect) {
    e <- sub("&.*", "", effect)
    data.table::fcase(
        e == "missense_variant", "Missense_Mutation",
        e %in% c("stop_gained"), "Nonsense_Mutation",
        e %in% c("stop_lost"), "Nonstop_Mutation",
        e %in% c("start_lost", "initiator_codon_variant"), "Translation_Start_Site",
        grepl("splice_acceptor|splice_donor", e), "Splice_Site",
        grepl("splice_region", e), "Splice_Region",
        e %in% c("synonymous_variant", "stop_retained_variant", "start_retained_variant"), "Silent",
        e == "5_prime_UTR_variant", "5'UTR",
        e == "3_prime_UTR_variant", "3'UTR",
        e %in% c("upstream_gene_variant"), "5'Flank",
        e %in% c("downstream_gene_variant"), "3'Flank",
        grepl("intron", e), "Intron",
        e %in% c("non_coding_transcript_exon_variant"), "RNA",
        default = "IGR")
}

#' @name sc_annotate_snv_drivers
#' @title sc_annotate_snv_drivers
#' @description
#'
#' SnpEff + OncoKB on a set of SNV sites (ids like chr1_18258813_G_A). SnpEff's
#' most severe annotation per site gives gene, consequence and protein change;
#' protein-altering and splice sites are queried in OncoKB by protein change
#' (GRCh38). The full OncoKB output is kept as <patient>.oncokb_full.tsv.gz in
#' work_dir (input of sc_cell_snv_oncotable()).
#'
#' @param mutations SNV ids
#' @param patient patient id (file names, OncoKB sample barcode)
#' @param work_dir folder for the VCF, MAF and OncoKB files
#' @param cfg sc_annotation_config()
#' @return data.table, one row per site: gene, consequence, impact, protein,
#'   classification, oncogenic, oncokb_level, mutation_effect, role, driver, ...
#' @export
#' @author Stanley Clarke
sc_annotate_snv_drivers <- function(mutations, patient, work_dir, cfg = sc_annotation_config()) {
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
    work_dir <- normalizePath(work_dir)
    parts <- data.table::tstrsplit(mutations, "_", fixed = TRUE)
    sites <- data.table::data.table(mutation = mutations, chrom = sub("^chr", "", parts[[1]]),
                                    pos = as.integer(parts[[2]]), ref = parts[[3]], alt = parts[[4]])
    sites <- sites[order(match(chrom, c(1:22, "X", "Y", "MT")), pos)]
    vcf <- file.path(work_dir, paste0(patient, ".sites.vcf"))
    writeLines(c("##fileformat=VCFv4.2", paste0("##contig=<ID=", unique(sites$chrom), ">"),
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"), vcf)
    data.table::fwrite(sites[, .(chrom, pos, mutation, ref, alt, ".", "PASS", ".")], vcf,
                       sep = "\t", append = TRUE, col.names = FALSE)

    ## SnpEff (most severe annotation first in ANN)
    ann_vcf <- file.path(work_dir, paste0(patient, ".sites.ann.vcf"))
    sc_singularity_exec(cfg, cfg$snpeff_image,
        sprintf("snpEff -Xmx8g -noStats -dataDir %s %s %s > %s", cfg$snpeff_cache, cfg$snpeff_db, vcf, ann_vcf),
        binds = c(work_dir, cfg$snpeff_cache))
    ann <- data.table::fread(cmd = sprintf("grep -v '^#' %s", ann_vcf), header = FALSE, sep = "\t",
                             select = c(3, 8), col.names = c("mutation", "info"))
    first <- sub(",.*", "", sub(".*ANN=", "", ann$info))
    f <- data.table::tstrsplit(first, "|", fixed = TRUE, keep = c(2, 3, 4, 7, 10, 11))
    ann <- data.table::data.table(mutation = ann$mutation, consequence = f[[1]], impact = f[[2]], gene = f[[3]],
                                  transcript = f[[4]], hgvsc = f[[5]], hgvsp = f[[6]])
    ann[gene == "", gene := NA_character_]
    ann[, protein := ifelse(nzchar(hgvsp), sc_hgvsp_short(hgvsp), NA_character_)]
    ann[, classification := sc_snpeff_class(consequence)]
    sites <- merge(sites, ann, by = "mutation", all.x = TRUE)

    ## OncoKB on protein-altering and splice sites
    query <- sites[classification %in% c("Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
                                          "Translation_Start_Site", "Splice_Site") & !is.na(gene)]
    if (nrow(query)) {
        maf <- query[, .(Hugo_Symbol = gene, Tumor_Sample_Barcode = patient, Chromosome = chrom,
                         Start_Position = pos, End_Position = pos, Reference_Allele = ref,
                         Tumor_Seq_Allele1 = ref, Tumor_Seq_Allele2 = alt, Variant_Type = "SNP",
                         Variant_Classification = classification, HGVSp = hgvsp,
                         HGVSp_Short = protein, ONCOTREE_CODE = cfg$oncotree_code, mutation)]
        maf_in <- file.path(work_dir, paste0(patient, ".oncokb_input.maf"))
        maf_out <- file.path(work_dir, paste0(patient, ".oncokb.maf"))
        data.table::fwrite(maf, maf_in, sep = "\t")
        sc_singularity_exec(cfg, cfg$unified_image,
            sprintf(paste("/opt/conda/envs/pact/bin/python /root/git/oncokb-annotator/MafAnnotator.py",
                          "-i %s -o %s -b \"$ONCOKB_TOKEN\" -t %s -r GRCh38 -q HGVSp_Short -d"),
                    maf_in, maf_out, cfg$oncotree_code),
            binds = work_dir, env = c(ONCOKB_TOKEN = sc_oncokb_token(cfg)))
        oncokb <- data.table::fread(maf_out, sep = "\t", quote = "")
        data.table::fwrite(oncokb, file.path(work_dir, paste0(patient, ".oncokb_full.tsv.gz")), sep = "\t")
        keep <- intersect(c("mutation", "ONCOGENIC", "HIGHEST_LEVEL", "MUTATION_EFFECT", "MUTATION_EFFECT_DESCRIPTION",
                            "GENE_SUMMARY", "VARIANT_SUMMARY", "TUMOR_TYPE_SUMMARY", "DIAGNOSTIC_SUMMARY",
                            "PROGNOSTIC_SUMMARY", "HIGHEST_DX_LEVEL", "HIGHEST_PX_LEVEL"), names(oncokb))
        sites <- merge(sites, oncokb[, ..keep], by = "mutation", all.x = TRUE)
    }
    genes <- data.table::as.data.table(readRDS(cfg$oncokb_genes))[, .(gene = Hugo_Symbol, role = Role)]
    sites <- merge(sites, unique(genes, by = "gene"), by = "gene", all.x = TRUE)
    data.table::setnames(sites, c("ONCOGENIC", "HIGHEST_LEVEL", "MUTATION_EFFECT"),
                         c("oncogenic", "oncokb_level", "mutation_effect"), skip_absent = TRUE)
    for (col in c("oncogenic", "oncokb_level", "mutation_effect")) if (!col %in% names(sites)) sites[[col]] <- NA_character_
    sites[oncogenic %in% c("", "Unknown"), oncogenic := NA_character_]
    sites[oncokb_level == "", oncokb_level := NA_character_]
    sites[, driver := oncogenic %in% c("Oncogenic", "Likely Oncogenic", "Resistance") | !is.na(oncokb_level)]
    sites[]
}

#' @name sc_snv_contexts
#' @title sc_snv_contexts
#' @description
#'
#' Trinucleotide context of each SNV as its SBS96 channel, e.g. "A[C>T]G"
#' (pyrimidine reference, as SigProfiler / COSMIC). NA for non-SNVs or when the
#' reference base does not match.
#'
#' @param mutations SNV ids like chr1_18258813_G_A
#' @param fasta indexed reference fasta
#' @return character vector
#' @export
#' @author Stanley Clarke
sc_snv_contexts <- function(mutations, fasta) {
    if (!requireNamespace("Rsamtools", quietly = TRUE)) stop("sc_snv_contexts needs the Rsamtools package")
    parts <- data.table::tstrsplit(mutations, "_", fixed = TRUE)
    chrom <- parts[[1]]
    pos <- as.integer(parts[[2]])
    ref <- parts[[3]]
    alt <- parts[[4]]
    fa <- Rsamtools::FaFile(fasta)
    seqnames <- as.character(GenomicRanges::seqnames(Rsamtools::scanFaIndex(fa)))
    if (!all(chrom %in% seqnames)) chrom <- ifelse(chrom %in% seqnames, chrom, sub("^chr", "", chrom))
    tri <- as.character(Rsamtools::getSeq(fa, GenomicRanges::GRanges(chrom, IRanges::IRanges(pos - 1L, pos + 1L))))
    ok <- nchar(ref) == 1 & nchar(alt) == 1 & substr(tri, 2, 2) == toupper(ref)
    comp <- c(A = "T", C = "G", G = "C", T = "A", N = "N")
    rc <- function(s) vapply(strsplit(s, ""), function(b) paste(rev(comp[b]), collapse = ""), "")
    flip <- ref %in% c("G", "A")
    tri2 <- ifelse(flip, rc(tri), tri)
    ref2 <- ifelse(flip, comp[ref], ref)
    alt2 <- ifelse(flip, comp[alt], alt)
    out <- paste0(substr(tri2, 1, 1), "[", ref2, ">", alt2, "]", substr(tri2, 3, 3))
    out[!ok | grepl("N", tri)] <- NA_character_
    out
}

SC_SBS96 <- local({
    b <- c("A", "C", "G", "T")
    subs <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    unlist(lapply(subs, function(s) as.vector(outer(b, b, function(x, y) paste0(x, "[", s, "]", y)))))
})

#' @name sc_fit_sbs_signatures
#' @title sc_fit_sbs_signatures
#' @description
#'
#' COSMIC SBS fits of named SNV sets with SigProfilerAssignment (cosmic_fit, in
#' the lab's container). Sets with fewer than min_mutations contexts are skipped.
#'
#' @param sets named list of SBS96 context vectors (see sc_snv_contexts())
#' @param work_dir folder for the SBS96 matrix and SigProfiler output
#' @param cosmic_version COSMIC version
#' @param genome genome build
#' @param min_mutations smallest set fitted
#' @param cfg sc_annotation_config()
#' @return long data.table (set, signature, activity, n), or NULL
#' @export
#' @author Stanley Clarke
sc_fit_sbs_signatures <- function(sets, work_dir, cosmic_version = 3.4, genome = "GRCh38", min_mutations = 20,
                                  cfg = sc_annotation_config()) {
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
    work_dir <- normalizePath(work_dir)
    sets <- sets[vapply(sets, function(x) sum(!is.na(x)), integer(1)) >= min_mutations]
    if (!length(sets)) return(NULL)
    counts <- vapply(sets, function(x) as.numeric(table(factor(x, levels = SC_SBS96))), numeric(96))
    matrix_file <- file.path(work_dir, "sbs96.txt")
    data.table::fwrite(data.table::data.table(MutationType = SC_SBS96, counts), matrix_file, sep = "\t")
    out <- file.path(work_dir, "assignment")
    unlink(out, recursive = TRUE)
    sc_singularity_exec(cfg, cfg$sigprofiler_image,
        sprintf(paste0("python -c \"from SigProfilerAssignment import Analyzer as A; ",
                       "A.cosmic_fit(samples='%s', output='%s', input_type='matrix', genome_build='%s', ",
                       "cosmic_version=%s, make_plots=False, verbose=False)\""),
                matrix_file, out, genome, cosmic_version),
        binds = work_dir)
    act_file <- list.files(out, pattern = "Assignment_Solution_Activities.txt$", recursive = TRUE, full.names = TRUE)[1]
    act <- data.table::fread(act_file)
    data.table::setnames(act, 1, "set")
    long <- data.table::melt(act, id.vars = "set", variable.name = "signature", value.name = "activity")
    long <- long[activity > 0]
    long[, n := colSums(counts)[set]]
    long[, signature := as.character(signature)]
    long[]
}

#' @name sc_cell_snv_oncotable
#' @title sc_cell_snv_oncotable
#' @description
#'
#' A cell's SNV rows in skilift's oncotable layout (as collect_oncokb()
#' returns), for the OncoKB-annotated sites with alt reads in the cell. Tiers
#' come from skilift's parse_oncokb_tier() on the patient's OncoKB MAF;
#' gene_location is the variant's own hg38 position.
#'
#' @param cell cell id
#' @param obs the cell's SNV observations (mutation, ref.count.t, alt.count.t, vaf)
#' @param oncokb_full the patient's <patient>.oncokb_full.tsv.gz (sc_annotate_snv_drivers())
#' @param roles named vector gene -> OncoKB role
#' @return data.table or NULL
#' @export
#' @author Stanley Clarke
sc_cell_snv_oncotable <- function(cell, obs, oncokb_full, roles) {
    hits <- merge(obs[alt.count.t > 0], oncokb_full, by = "mutation")
    if (!nrow(hits)) return(NULL)
    hits <- parse_oncokb_tier(hits, tx_cols = c("LEVEL_1", "LEVEL_2"), rx_cols = c("LEVEL_R1"),
                              dx_cols = c("LEVEL_Dx1"), px_cols = c("LEVEL_Px1"))
    hits[, role := roles[match(Hugo_Symbol, names(roles))]]
    col <- function(x) if (x %in% names(hits)) hits[[x]] else NA_character_
    hits[, .(
        gene = Hugo_Symbol, gene_summary = col("GENE_SUMMARY"), role,
        variant.g = paste0(Chromosome, ":", Start_Position, "-", End_Position, " ", Reference_Allele, ">", Tumor_Seq_Allele2),
        variant.c = NA_character_, variant.p = HGVSp_Short, annotation = Variant_Classification,
        type = data.table::fcase(Variant_Classification == "Missense_Mutation", "missense",
                                 Variant_Classification %in% c("Nonsense_Mutation", "Nonstop_Mutation", "Translation_Start_Site"), "trunc",
                                 Variant_Classification == "Splice_Site", "splice", default = NA_character_),
        tier, tier_description = tier_factor, variant_summary = col("VARIANT_SUMMARY"),
        therapeutics = tx_string, resistances = rx_string, diagnoses = dx_string, prognoses = px_string,
        distance = NA_integer_, effect = col("MUTATION_EFFECT"), effect_description = col("MUTATION_EFFECT_DESCRIPTION"),
        major_count = NA_real_, minor_count = NA_real_, major_snv_copies = NA_real_, minor_snv_copies = NA_real_,
        altered_copies = NA_real_, segment_cn = NA_integer_, ref = as.integer(ref.count.t), alt = as.integer(alt.count.t),
        VAF = vaf, vartype = "SNV", track = "variants", source = "oncokb_maf", is_multi_hit_per_gene = FALSE,
        gene_location = paste0(Chromosome, ":", Start_Position, "-", End_Position),
        id = cell)]
}

#' @name sc_write_cell_filtered_events
#' @title sc_write_cell_filtered_events
#' @description
#'
#' Write a cell's filtered.events.json with create_filtered_events() from its
#' SNV rows and fusion / CNA rows (oncotable layout).
#'
#' @param cell cell id
#' @param snv_rows sc_cell_snv_oncotable() rows (or NULL)
#' @param other_rows fusion / CNA oncotable rows (or NULL)
#' @param jabba_gg the cell's JaBbA graph (path) or NULL
#' @param cell_dir the cell's gOS folder
#' @return create_filtered_events() result
#' @export
#' @author Stanley Clarke
sc_write_cell_filtered_events <- function(cell, snv_rows, other_rows, jabba_gg, cell_dir) {
    ot <- data.table::rbindlist(list(snv_rows, other_rows), fill = TRUE)
    if (!nrow(ot)) ot <- data.table::data.table(type = NA, source = "none")
    rds <- tempfile(fileext = ".rds")
    on.exit(unlink(rds))
    saveRDS(ot, rds)
    create_filtered_events(pair = cell, oncotable = rds, jabba_gg = jabba_gg,
                           out_file = file.path(cell_dir, "filtered.events.json"))
}

#' @name sc_cell_fusions
#' @title sc_cell_fusions
#' @description
#'
#' gGnome fusions of one cell's JaBbA graph as nf-gos FUSIONS does them:
#' bin/Fusions.R in mskilab/unified:0.0.10 with the hg38 fusions GENCODE.
#' Writes fusions.rds and altedge.annotations.tsv to work_dir.
#'
#' @param cell cell id
#' @param jabba_gg path to the cell's JaBbA gGraph rds
#' @param work_dir output folder
#' @param cfg sc_annotation_config()
#' @return path to fusions.rds (invisibly)
#' @export
#' @author Stanley Clarke
sc_cell_fusions <- function(cell, jabba_gg, work_dir, cfg = sc_annotation_config()) {
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
    work_dir <- normalizePath(work_dir)
    jabba_gg <- normalizePath(jabba_gg)
    sc_singularity_exec(cfg, cfg$fusions_image,
        paste("cd", shQuote(work_dir), "&& Rscript", file.path(cfg$nf_gos_real, "bin/Fusions.R"),
              "--id", cell, "--junctions", shQuote(jabba_gg), "--gencode", shQuote(cfg$fusions_gencode),
              "--outdir", shQuote(work_dir), "--cores 1 >", shQuote(file.path(work_dir, "fusions.log")), "2>&1"),
        binds = c(work_dir, dirname(jabba_gg), cfg$nf_gos_real, dirname(cfg$fusions_gencode)))
    invisible(file.path(work_dir, "fusions.rds"))
}

#' @name sc_cell_cna_fusions
#' @title sc_cell_cna_fusions
#' @description
#'
#' Copy-number and fusion drivers of one cell, as nf-gos does them: OncoKB on the
#' cell's fusions (sc_cell_fusions(), run first unless fusions = FALSE) and gene
#' amplifications / deletions via nf-gos process_singularity.sh (SNV part
#' skipped), then skilift oncotable() on the OncoKB output. gene_location is
#' set from the hg38 gene file (skilift's gene_locations.rds is hg19).
#'
#' @param cell cell id
#' @param jabba_gg path to the cell's JaBbA gGraph rds
#' @param work_dir output folder
#' @param out_rds where the oncotable rows are saved
#' @param gencode_gr process_gencode() of cfg$gencode
#' @param cytoband_gr process_cytoband() of the hg38 cytoband
#' @param fusions include fusions
#' @param amp.thresh,del.thresh passed to oncotable()
#' @param cfg sc_annotation_config()
#' @return oncotable rows (data.table)
#' @export
#' @author Stanley Clarke
sc_cell_cna_fusions <- function(cell, jabba_gg, work_dir, out_rds, gencode_gr, cytoband_gr, fusions = TRUE,
                                amp.thresh = 1.5, del.thresh = 0.5, cfg = sc_annotation_config()) {
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
    work_dir <- normalizePath(work_dir)
    jabba_gg <- normalizePath(jabba_gg)
    fusions_rds <- file.path(work_dir, "fusions.rds")
    if (fusions && !file.exists(fusions_rds)) sc_cell_fusions(cell, jabba_gg, work_dir, cfg)
    if (!fusions || !file.exists(fusions_rds)) fusions_rds <- "/dev/null"
    sc_singularity_exec(cfg, cfg$unified_image,
        paste("export HOME=/root; set +u; source /opt/conda/etc/profile.d/conda.sh; conda activate pact;",
              "cd", shQuote(work_dir), "&&",
              "bash", file.path(cfg$nf_gos_real, "bin/oncokb/process_singularity.sh"), "/dev/null", fusions_rds, jabba_gg,
              "/dev/null", "GRCh38", "1", "FALSE", "TRUE", cfg$oncotree_code, "/dev/null", cfg$oncokb_genes, cfg$cna_gencode,
              ">", shQuote(file.path(work_dir, "oncokb.log")), "2>&1"),
        binds = c(work_dir, dirname(jabba_gg), cfg$nf_gos_real, cfg$mski_base, dirname(cfg$cna_gencode)),
        env = c(ONCOKB_TOKEN = sc_oncokb_token(cfg)))
    fusions_tsv <- file.path(work_dir, "merged_oncokb_fusions.tsv")
    cna_tsv <- file.path(work_dir, "merged_oncokb_cna.tsv")
    ot <- oncotable(pair = cell,
                    oncokb_fusions = if (file.exists(fusions_tsv)) fusions_tsv else NA,
                    oncokb_cna = if (file.exists(cna_tsv)) cna_tsv else NA,
                    jabba_gg = jabba_gg, gencode = gencode_gr, cytoband = cytoband_gr,
                    amp.thresh = amp.thresh, del.thresh = del.thresh, verbose = FALSE)
    ot <- ot[!is.na(type) & source %in% c("oncokb_fusions", "oncokb_cna")]
    if (nrow(ot) && "gene" %in% names(ot)) {
        genes <- readRDS(cfg$cna_gencode)
        k <- match(ot$gene, genes$gene_name)
        hg38 <- paste0(as.character(GenomicRanges::seqnames(genes))[k], ":", GenomicRanges::start(genes)[k], "-", GenomicRanges::end(genes)[k])
        ot[, gene_location := ifelse(is.na(k), gene_location, hg38)]
    }
    saveRDS(ot, out_rds)
    ot
}

#' @name sc_drop_chrx_homdels
#' @title sc_drop_chrx_homdels
#' @description
#'
#' Drop chrX homozygous deletions from oncotable rows (chrX copy number is not
#' reliable in single-cell PTA libraries; male cells show CN 0 across chrX).
#' chrX amplifications are kept.
#'
#' @param ot oncotable rows
#' @return filtered rows
#' @export
#' @author Stanley Clarke
sc_drop_chrx_homdels <- function(ot) {
    if (is.null(ot) || !nrow(ot) || !all(c("type", "gene_location") %in% names(ot))) return(ot)
    ot[!(type %in% "homdel" & grepl("^(chr)?X:", gene_location))]
}

#' @name sc_patient_filtered_events
#' @title sc_patient_filtered_events
#' @description
#'
#' Patient-level filtered.events.json from the cells' own files: each event
#' (same gene, fusion genes, vartype, type and variant) once, with how many cells
#' carry it (cells = "n/N", cell_fraction, n_cells, cell_ids), alt / ref reads
#' pooled over those cells and median copy numbers. Counts and fractions are
#' over the tumor cells; carriers among the other cells (normals) are reported
#' separately as n_normal_cells.
#'
#' @param cell_ids the patient's cells
#' @param tumor_cells the tumor cells among them (denominator of cell_fraction)
#' @param data_dir gOS data folder
#' @param patient patient id (output goes to data_dir/patient)
#' @return the events (invisibly), or NULL
#' @export
#' @author Stanley Clarke
sc_patient_filtered_events <- function(cell_ids, data_dir, patient, tumor_cells = cell_ids) {
    events <- data.table::rbindlist(lapply(cell_ids, function(cell) {
        f <- file.path(data_dir, cell, "filtered.events.json")
        if (!file.exists(f) || file.size(f) < 3) return(NULL)
        x <- data.table::as.data.table(jsonlite::fromJSON(f))
        if (!nrow(x)) return(NULL)
        x[, cell := cell]
    }), fill = TRUE)
    if (!nrow(events)) return(invisible(NULL))
    key_cols <- intersect(c("gene", "fusion_genes", "vartype", "type", "Variant"), names(events))
    for (k in key_cols) events[is.na(get(k)), (k) := ""]
    tumor_cells <- intersect(tumor_cells, cell_ids)
    n_total <- length(tumor_cells)
    num <- function(x) suppressWarnings(as.numeric(x))
    out <- events[, c(
        .SD[1, !c("cell", "alt", "ref", "VAF", "estimated_altered_copies", "segment_cn", "id"), with = FALSE],
        list(id = patient,
             alt = if ("alt" %in% names(.SD)) sum(num(alt), na.rm = TRUE) else NA_real_,
             ref = if ("ref" %in% names(.SD)) sum(num(ref), na.rm = TRUE) else NA_real_,
             estimated_altered_copies = if ("estimated_altered_copies" %in% names(.SD)) stats::median(num(estimated_altered_copies), na.rm = TRUE) else NA_real_,
             segment_cn = if ("segment_cn" %in% names(.SD)) stats::median(num(segment_cn), na.rm = TRUE) else NA_real_,
             n_cells = sum(unique(cell) %in% tumor_cells),
             n_normal_cells = sum(!(unique(cell) %in% tumor_cells)),
             cell_ids = paste(sort(unique(cell)), collapse = ","))),
        by = key_cols, .SDcols = setdiff(names(events), key_cols)]
    out[, VAF := ifelse(alt + ref > 0, alt / (alt + ref), NA_real_)]
    out[, cells := paste0(n_cells, "/", n_total)]
    out[, cell_fraction := round(n_cells / n_total, 3)]
    data.table::setorder(out, Tier, -n_cells)
    jsonlite::write_json(out, file.path(data_dir, patient, "filtered.events.json"),
                         dataframe = "rows", na = "null", auto_unbox = TRUE, digits = NA)
    invisible(out)
}
