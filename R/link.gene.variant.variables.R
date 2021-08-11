#' Link a column of gene variant variables that are found in the REDCap Genomics Instrument to their appropriate strings
#' @description
#' `link.gene.variant.variables()` creates a column of gene variant variables that are found in the REDCap Genomics Instrument
#' @param data  A Data Frame of Gene Names, Nucleotide Sequences, Amino Acid Sequences and possibly cell-free DNA percentages. This must be a single column data frame, with tidy data (i.e. one name or sequence per cell). These data would like come from a function such as gene.variants.isolate.snapshot()
#'
#' @return a data frame with an additional column with gene variant variables for the REDCap Genomic Instrument. The two columns will be labeled, "genes_variants" and "var". These data could then be used as part of a larger function such as gene.variatns.ready.for.redcap.snapshot()
#' @export
#' @examples
#' # Test with embedded data frame "snapshot_sample_report"
#' snapshot_sample_report %>%
#'   gene.variants.isolate.snapshot() %>%
#'   link.gene.variant.variables()
link.gene.variant.variables <- function(data){

  ##########################################################################################################################
  # Create regexs that we will use to select our data
  ##########################################################################################################################
  genes_boundary_df <- genetex::genes_boundary_regex() # this is a df of unique genes names concatenated with "|"
  genes_boundary_regex <- genes_boundary_df$genes # this is a character string of the genes for our regex
  nuc_regex <- "[ACTG]>[ACTG]|del[ACTG]|([ACTG]{1,}[0-9]{1,}_[0-9]{1,}[ACTG])" # This regex should identify those rows that have nucleotide changes
  aa_regex <- "(\\b([A-Z][0-9]{1,}(([A-Z])|(_[A-Z][1-9]{1,}del)|(fs\\*[1-9]{1,})|(\\*)|(fs)|(del)))|(p\\.[A-Z]))|([0-9]ins[A-Z])|\\b[A-Z][a-z]{1,3}[0-9]{1,}[A-Z][a-z]{1,3}"
  cfdna_regex <- "\\b[0-9]{1,2}\\.[0-9]{1,2}%"

  ##########################################################################################################################
  # Add column with appropriate variable names
  ##########################################################################################################################
  # First rename column to ensure data function will run
  data <- data
  names(data)[1] <- "x"
  # assign the correct variable "variant_gene", "variant_nucleotide", "variant_protein" or "variant_gene_perc_cfdna"
  dt <- data %>%
    dplyr::mutate(var = base::ifelse(test = stringr::str_detect(string = data$x,
                                                                pattern = stringr::regex(genes_boundary_regex)),
                                     yes = "variant_gene",
                                     no = ifelse(test = stringr::str_detect(string = data$x,
                                                                            pattern = stringr::regex(nuc_regex)),
                                                 yes = "variant_nucleotide",
                                                 no = ifelse(test = stringr::str_detect(string = data$x,
                                                                                        pattern = stringr::regex(aa_regex)),
                                                             yes = "variant_protein",
                                                             no = ifelse(test = stringr::str_detect(string = data$x,
                                                                                                    pattern = stringr::regex(cfdna_regex)),
                                                                         yes = "variant_gene_perc_cfdna",
                                                                         no = "")))))
  # rename the first column
  names(dt)[1] <- "genes_variants"
  # change first column to character
  dt$genes_variants <- as.character(dt$genes_variants)

  return(dt)
}

