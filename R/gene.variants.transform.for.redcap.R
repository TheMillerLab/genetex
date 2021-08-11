#' Transform a data frame of untidy gene name and variants into a tidy format for import in REDCap
#' @description
#' `gene.variants.transform.for.redcap()` transforms a data frame from gene_variants.isolate.onco() to be imported into a REDCap Instrument with the "Genomics Instrument"
#' @param data  A data frame of gene names with their gene variants in an untidy, single cell, e.g. as an output of one of the gene.variants.isolate...() functions. Required.
#'
#' @return a data frame with two columns: variables (the redcap variable names) and results (the data to be imported in to redcap)
#' @export
#'
gene.variants.transform.for.redcap <- function(data){

  ##########################################################################################################################
  # Create regexs that we will use to select our data
  ##########################################################################################################################
  genes_boundary_df <- genetex::genes_boundary_regex() # this is a df of unique genes names concatenated with "|"
  genes_boundary_regex <- genes_boundary_df$genes # this is a character string of the genes for our regex

  ##########################################################################################################################
  #  Transform Data for Loading Into REDCap
  ##########################################################################################################################

  # first makes sure first column has the appropriate name
  names(data)[1] <- "genes_variants"

  # tokenize the vector with space as the delimiter
  dt.1 <- splitstackshape::cSplit(
    indt = data,
    splitCols = "genes_variants",
    sep = " ",
    direction = "long"
  )
  # assign the correct variable "variant_gene", "variant_nucleotide", or "variant_protein" using link.gene.variant.variables()
  dt.2 <- genetex::link.gene.variant.variables(data = dt.1)

  # Create a grouping system based on gene names
  dt.3 <- dt.2 %>% dplyr::mutate(keywords = stringr::str_detect(dt.2$genes_variants,
                                                                pattern = regex(genes_boundary_regex)),
                                 group = cumsum(keywords))

  # paste "group" and "var" to create a "variable" column, which are the variable names that are used in the Genomics Instrument
  dt.4 <- dt.3 %>% mutate(variables = paste(var, group, sep = "_"))

  # Select appropriate columns and rename
  dt.5 <- dt.4 %>% select(variables, genes_variants) %>% rename(results = genes_variants)
  # Create a grouping system based on gene names so we can arrange in a moment
  dt.6 <- dt.5 %>%
    dplyr::mutate(keywords = stringr::str_detect(string = dt.5$results,
                                                 pattern = stringr::regex(genes_boundary_regex)),
                  group = base::cumsum(keywords))
  # TERT promoter variants are often reported as a conglomerate of multiple alterations. This can lead to duplicates of
  ## variable names, which cause the algorithm to error. Thus, we should aggregate these in one cell.
  ### Let's do that by transposeing wide via paste0
  dt.7 <- dt.6 %>% dplyr::group_by(variables) %>% dplyr::mutate(merge = base::paste0(results, collapse = "; ")) %>% ungroup()

  # slice to eliminate duplicates
  ## first, we'll make one more merge, so that we don't eliminate two rows with the same gene name that we
  ### actually want to keep (e.g. have two different mutations that we want to capture)
  dt.8a <- dt.7 %>% group_by(variables) %>% dplyr::mutate(merge_2 = base::paste(group, merge, sep = "-"))
  dt.8 <- dt.8a %>% dplyr::group_by(merge_2) %>% dplyr::slice_head() %>% ungroup()

  # Order by group to maintain the order in the report
  dt.9 <- dt.8 %>% dplyr::arrange(group, variables)
  dt.final <- dt.9 %>% select(variables, merge) %>% rename(results = merge)

  # Add variants_number
  variants.number.1 <- stringr::str_detect(string = dt.final$variables, pattern = "variant_gene")
  variants.number.2 <- sum(variants.number.1)
  variants.number <- if(is.na(dt.final[1,2])) 0 else variants.number.2
  variants.number.df <- data.frame(variables = "variant_number",
                              results = variants.number)

  # combine dfs
  gene_variants <- rbind(dt.final, variants.number.df)
  gene_variants$results <- base::as.character(gene_variants$results)

  return(gene_variants)
}

