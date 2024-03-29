% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cnv.fmi.R
\name{cnv.fmi}
\alias{cnv.fmi}
\title{Apply NLP to genomic reports to extract copy number variants (CNVs) from foundation medicine reports}
\usage{
cnv.fmi(data = dplyr::tibble(Results = readr::clipboard()))
}
\arguments{
\item{data}{The data frame of the genomic report of interest. This can be copied to the Clipboard from the EHR report by "selecting all" and then "copy". That is the default. If you don't use the clipboard function, you can use a data frame of the text file. The single column data frame should be labeled "Results". Required.}
}
\value{
a data frame with two columns: variables (the redcap variable names) and results (the data to be imported in to redcap)
}
\description{
`cnv.fmi()` provides natural language processing tools to abstract copy number variants data from a variety of genomic reports
}
