#' Skilift: Translating nf-gOS to gOS-readable inputs
#' 
#' @import data.table
#' @import R6
#' @import jsonlite
#' @import gGnome
#' @import ggplot2 
#' @import parallel
#' @import httr
#' @import futile.logger
#' @importFrom pals glasbey
#' @importFrom googlesheets4 read_sheet gs4_auth
#' @importFrom tryCatchLog tryCatchLog
#' @importFrom VariantAnnotation readVcf geno ref alt
#' @importFrom arrow Table write_feather schema float32
#' @importFrom gGnome cov2cov.js
#' @importFrom gUtils gr2dt dt2gr %Q% gr.val %$% %*% gr.chr gr.nochr gr.findoverlaps gr.stripstrand gr.tile
#' @importFrom skitools rel2abs
#' @importFrom magrittr %>%
#' @importFrom parallel mclapply
#' @importFrom VariantAnnotation readVcf geno
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom GenomeInfoDb Seqinfo seqnames seqnames<- seqinfo seqinfo<- seqlengths seqlengths<- seqlevels seqlevels<- seqlevelsStyle seqlevelsStyle<-
#' @importFrom BiocGenerics width
#' @importFrom glue glue
#' @importMethodsFrom BiocGenerics width sort
#' @importMethodsFrom IRanges trim start end
#' @importMethodsFrom GenomeInfoDb seqnames seqnames<- seqinfo seqinfo<- seqlengths seqlengths<- seqlevels seqlevels<- seqlevelsStyle seqlevelsStyle<- sortSeqlevels
#' @importMethodsFrom GenomicRanges start end
#' @importMethodsFrom MatrixGenerics rowRanges
#' @importMethodsFrom S4Vectors split mcols mcols<- values values<-
#' @useDynLib Skilift, .registration=TRUE
registerS3method(genname = "merge", class = "data.table", method = data.table::merge.data.table)
"_PACKAGE"


.onLoad = function(libname, pkgname) {
	skilift_jabba_columns = getOption("skilift_jabba_columns")
	is_skilift_jabba_columns_option = !is.null(skilift_jabba_columns)
	is_skilift_jabba_columns_valid_option = (
		is_skilift_jabba_columns_option && all(skilift_jabba_columns %in% Skilift:::priority_columns_jabba_og)
	)
	if (is_skilift_jabba_columns_valid_option) {
		message("option('skilift_jabba_columns') already set to: ", paste(skilift_jabba_columns, collapse = ", "))
	} else {
		message("option('skilift_jabba_columns') set to: ", paste(Skilift:::priority_columns_jabba_og, collapse = ", "))
		options("skilift_jabba_columns" = Skilift:::priority_columns_jabba_og)
	}
	
}

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