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
genes_regex <- function(file = genetex::genes){
  ##########################################################################################################################
  # Load the genes look up table
  ##########################################################################################################################
  file <- genetex::genes
  ##########################################################################################################################
  # Concatenate all of the rows separated by "\\b"
  ##########################################################################################################################
  gene_symbols <- file %>%
    mutate(genes = base::paste0(
      genes,
      collapse = "|"
    )
    )
  ##########################################################################################################################
  # Select only the first row
  ##########################################################################################################################
  genes_regex <- gene_symbols[1,]
  ##########################################################################################################################
  # Return the df of interest: genes_regex
  ##########################################################################################################################
  return(genes_regex)
}
