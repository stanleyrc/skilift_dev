% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{convert_signature_meta_json}
\alias{convert_signature_meta_json}
\title{convert_signature_meta_json}
\usage{
convert_signature_meta_json(jsons.dt, cores = 1)
}
\arguments{
\item{jsons.dt}{json dt after using read_meta_data_json}

\item{cores}{cores for generating each signature}
}
\value{
list object with all of the signature distribution data.tables
}
\description{
converts signatures from the list objects in the meta data to data.tables to write jsons for distributions. Used in create_distributions
}
\author{
Stanley Clarke
}
