"_package"

## Global Variables

#' Paired reads factor 
#' 
#' Multiply paired reads by this number for coverage
#' 
#' To get binned coverage into base level coverage,
#' multiply by paired reads factor (2)
PAIRED_READS_FACTOR = 2

#' Read length
#' 
#' Multiply paired reads by this number for read length
#' 
#' To get binned coverage into read length level coverage,
#' multiply by read length (151)
READ_LENGTH = 151

#' Snpeff Protein Coding Annotations
#' 
#' Nuff said.
#' 
#' 
snpeff_protein_coding_annotations = c("frameshift_variant", "stop_lost", "start_lost", "stop_gained", 
"chromosome", "exon_loss_variant", "feature_ablation", "duplication", 
"splice_acceptor_variant", "splice_donor_variant", "splice_region_variant", 
"inframe_insertion", "disruptive_inframe_insertion", "inframe_deletion", 
"disruptive_inframe_deletion", "coding_sequence_variant", "missense_variant", 
"protein_protein_contact", "structural_interaction_variant", 
"rare_amino_acid_variant")


#' QC Flags Tresholds
#'
#' Parse QC flags into strings
#'
#' QC metrics from picard need to be parsed based on coverage, insert size
#' total number of reads, duplicate rate.
#' The strings should be parsable into a form
#' digested in gOS and shown as a single "PASS"/Checkmark", "Warning", or "Fail". 
#' The actual metrics should show up on hover.
qc_flag_thresholds = list(
    list("FAIL", "greater_than_or_equal_to_50x", `<`, 0.99, "Fraction of genome covered at 50X","99%"),
    list("WARN", "purity", `<`, 0.2, "Purity", "20%"),
    list("WARN", "insert_size", `<`, 300, "Insert Size", "300 bp"),
    list("WARN", "percent_duplication", `>`, 0.3, "Duplicate percent", "30%"),
    list("FAIL", "fraction_of_reads_aligned", `<`, 0.9, "Percent of reads aligned", "90%"),
    list("FAIL", "conpair_concordance_metric", `<`, 0.9, "Tumor/Normal SNP concordance", "90%")
)
