#' Create a single regex of all cancer-related gene symbols
#' @description
#' `genes_regex()` produces a regular expression of over 900 HGNC gene names
#' @param file File of gene symbols. Default imports csv file from data-raw containing over 900 gene symbols. Required
#'
#' @return a single cell string with all gene symbols concatenated by a pipe ("|")
#' @export
#'
#' @examples
#' #Import "genes" data set from package
#' genes %>%
#' genes_regex()
genes_regex <- function(file = readr::read_tsv("data-raw/genes.txt")
){
  gene_symbols <- file %>%
    mutate(genes = base::paste0(
      genes,
      collapse = "|"
      )
      )
  genes_regex <- gene_symbols[1,]


  return(genes_regex)
}

