% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variant-qc.R
\name{lift_variant_qc}
\alias{lift_variant_qc}
\title{lift_variant_qc}
\usage{
lift_variant_qc(cohort, output_data_dir, cores = 1)
}
\arguments{
\item{cohort}{Cohort object containing sample information}

\item{output_data_dir}{Base directory for output files}

\item{cores}{Number of cores for parallel processing (default: 1)}
}
\value{
None
}
\description{
Creates variant QC metrics for all samples in a cohort and writes them to JSON files
}
