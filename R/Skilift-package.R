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
#' @importFrom tryCatchLog tryCatchLog
#' @importFrom VariantAnnotation readVcf geno ref alt
#' @importFrom arrow Table write_feather schema float32
#' @importFrom gGnome cov2cov.js
#' @importFrom gUtils gr2dt dt2gr %Q% gr.val %$% %*% gr.chr gr.nochr gr.findoverlaps
#' @importFrom skitools rel2abs
#' @importFrom magrittr %>%
#' @importFrom parallel mclapply
#' @importFrom VariantAnnotation readVcf geno
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom GenomeInfoDb Seqinfo seqnames seqnames<- seqinfo seqinfo<- seqlengths seqlengths<- seqlevels seqlevels<- seqlevelsStyle seqlevelsStyle<-
#' @importFrom BiocGenerics width
#' @importMethodsFrom BiocGenerics width sort
#' @importMethodsFrom IRanges trim start end
#' @importMethodsFrom GenomeInfoDb seqnames seqnames<- seqinfo seqinfo<- seqlengths seqlengths<- seqlevels seqlevels<- seqlevelsStyle seqlevelsStyle<- sortSeqlevels
#' @importMethodsFrom GenomicRanges start end
#' @importMethodsFrom MatrixGenerics rowRanges
#' @importMethodsFrom S4Vectors split mcols mcols<- values values<-
#' @useDynLib Skilift, .registration=TRUE
registerS3method(genname = "merge", class = "data.table", method = data.table::merge.data.table)
"_PACKAGE"
