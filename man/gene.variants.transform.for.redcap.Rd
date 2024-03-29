% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene.variants.transform.for.redcap.R
\name{gene.variants.transform.for.redcap}
\alias{gene.variants.transform.for.redcap}
\title{Transform a data frame of untidy gene name and variants into a tidy format for import in REDCap}
\usage{
gene.variants.transform.for.redcap(data)
}
\arguments{
\item{data}{A data frame of gene names with their gene variants in an untidy, single cell, e.g. as an output of one of the gene.variants.isolate...() functions. Required.}
}
\value{
a data frame with two columns: variables (the redcap variable names) and results (the data to be imported in to redcap)
}
\description{
`gene.variants.transform.for.redcap()` transforms a data frame from gene_variants.isolate.onco() to be imported into a REDCap Instrument with the "Genomics Instrument"
}
