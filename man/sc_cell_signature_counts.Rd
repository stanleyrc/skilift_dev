\name{sc_cell_signature_counts}
\alias{sc_cell_signature_counts}
\title{Count SBS96 mutations on each cell's root-to-tip path}
\usage{sc_cell_signature_counts(variants, tree, cell_ids = tree$tip.label)}
\arguments{
\item{variants}{Variant records from snv_matrix.json with id, category, anchor and context.}
\item{tree}{An ape phylo with matching tip identifiers.}
\item{cell_ids}{Cells to return; defaults to the tree tips.}
}
\value{An SBS96 by cell integer matrix.}
\description{Counts mapped truncal, subclonal and private mutations on ancestral paths using MRCA anchors. Ambiguous, unmapped and outside-tumor sites are excluded. Read dropout does not remove ancestral mutations.}
