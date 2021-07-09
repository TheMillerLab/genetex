#' Text Mine a Oncopanel genomic report and prepares the data for import into the Genomics Instrument in REDCap
#' @description
#' `gene.variants.ready.for.redcap.snapshot()` text mines a data frame a Oncopanel report so it can be imported into a REDCap "Genomics Instrument" when combined with the other functions of "genetex_to_redcap"
#' @param data  The data frame of the genomic report of interest, here "Oncopanel". This can be copied to the Clipboard from the EHR report by selecting the relevant text and then "copy". That is the default. If you don't use the clipboard function, you can use a data frame of the text file. The single column data frame should be labeled "Results". Required..
#'
#' @return a data frame with two columns: "variables" (the redcap variable names) and "results" (the data to be imported in to REDCap). This data can then be combined with the outputs of the other functions of "genetex_to_redcap". That end product would then be imported into REDCap
#' @export
#' @examples
#' oncopanel_sample_report %>%
#'  gene.variants.ready.for.redcap.oncopanel()
gene.variants.ready.for.redcap.oncopanel <- function(data = dplyr::tibble(Results = clipboard())){

  ##########################################################################################################################
  # Load Data
  ##########################################################################################################################
  ## Use the gene.variants.isolate.oncopanel() function to text mine the genetic variants from an Oncopanel report and return
  ### a single column, untidy data frame
  dt <- genetex::gene.variants.isolate.oncopanel(data)

  ##########################################################################################################################
  #  Transform Data for Loading Into REDCap
  ##########################################################################################################################
  ## Use the gene.variants.transform.for.redcap() function to return a tidy data frame that has two columns "variables" that
  ### contains the variable names for the Genomics Instrument and "results' which is either the Gene Name, or genetic variant
  gene_variants <- genetex::gene.variants.transform.for.redcap(dt)

  return(gene_variants)
}

