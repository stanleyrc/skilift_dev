\name{sc_export_cell_signatures}
\alias{sc_export_cell_signatures}
\title{Export per-cell SBS signature fits on root-to-tip mutation sets}
\usage{sc_export_cell_signatures(patient_dir, out_dir = patient_dir,
  work_dir = file.path(out_dir, "cell_signature_work"), min_mutations = 100L,
  counts = NULL, fit_fun = sc_fit_sbs_signatures,
  quality_fun = sc_signature_fit_quality, cfg = sc_annotation_config())}
\arguments{
\item{patient_dir}{Folder containing tree.nwk and snv_matrix.json.}
\item{out_dir}{Output folder, optionally a staging directory.}
\item{work_dir}{Private SigProfiler working directory.}
\item{min_mutations}{Minimum valid SBS96 mutations; at least 100.}
\item{counts}{Optional precomputed SBS96 by cell counts matrix.}
\item{fit_fun}{Fitting function; defaults to sc_fit_sbs_signatures.}
\item{quality_fun}{Fit-quality reader; defaults to sc_signature_fit_quality.}
\item{cfg}{Annotation configuration.}
}
\value{Invisibly, the exported object.}
\description{Fits all eligible cell catalogs with SigProfilerAssignment COSMIC SBS and writes signatures_cells.json with fractions, mutation counts and cosine. Cells below the threshold are marked too_few with no exposures. Existing output is backed up before atomic replacement.}
