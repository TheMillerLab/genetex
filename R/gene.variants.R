#' Abstract Genes with nucleotide variants from genomic reports
#' @description
#' `gene.variants()` integrates various platform-specific NLP functions to text mine gene names and nucleotide variants from genomic reports and transforms them to structured data for import into REDCap
#' @param data  The data frame of a genomic report. Ideally this is information copied to the Clipboard from the EHR report by selecting the relevant text and then "copy". That is the default. If you don't use the clipboard function, you can use a data frame of the text file. The single column data frame should be labeled "Results". Required.
#'
#' @return a data frame with two columns: variables (the redcap variable names) and results (the data to be imported in to redcap)
#' @export
#'
gene.variants <- function(data = dplyr::tibble(Results = clipboard()),
                                           platform){

  ##########################################################################################################################
  # load data and platform
  ##########################################################################################################################
  dt <- data
  ##########################################################################################################################
  # choose relevant platform
  ##########################################################################################################################
  platform_num <- genetex::platform(platform)
  # BWH == 1
  # MGH == 2
  # MSK == 3
  # Foundation == 4
  # Tempus == 5
  # Guardant == 6
  ##########################################################################################################################
  # final df
  ##########################################################################################################################
  ifelse(test = platform_num ==1,
         yes = gene_variants <- genetex::gene.variants.ready.for.redcap.oncopanel(data = dt),
         no = ifelse(test = platform_num ==2,
                     yes = gene_variants <- genetex::gene.variants.ready.for.redcap.snapshot(data = dt),
                     no = ifelse(test = platform_num == 4,
                                 yes = gene_variants <- genetex::gene.variants.fmi(data = dt),
                                 no = ifelse(test = platform_num == 6,
                                             yes = gene_variants <- genetex::gene.variants.guardant(data = dt),
                                             no = gene_variants <- base::data.frame(variables = c("variant_number"),
                                                                        results = 0)))))
  # change to character vector
  gene_variants$results <- base::as.character(gene_variants$results)

  return(gene_variants)
}

