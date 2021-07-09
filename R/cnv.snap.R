#' Apply NLP to genomic reports to extract copy number variants (CNVs) from MGH SNaPshot reports
#' @description
#' `cnv.snap()` provides natural language processing tools to abstract copy number variants data SNaPshot reports
#' @param data The data frame of the genomic report of interest. This can be copied to the Clipboard from the EHR report by selecting the relevant text and then "copy". That is the default. If you don't use the clipboard function, you can use a data frame of the text file. The single column data frame should be labeled "Results". Required.
#' @return a data frame with two columns: variables (the redcap variable names) and results (the data to be imported in to redcap)
#' @export
#'
cnv.snap <- function(data = dplyr::tibble(Results = readr::clipboard())){
  ##########################################################################################################################
  # load data
  ##########################################################################################################################
  dt <- data

  ##########################################################################################################################
  # Create regexs that we will use to select our data
  ##########################################################################################################################
  # Keep only those rows with the following
  ## call the genes_regex function to create a regex of all the gene names
  genes_boundary_df <- genetex::genes_boundary_regex() # this is a df of unique genes names concatenated with "|"
  genes_boundary_regex <- genes_boundary_df$genes # this is a character string of the genes for our regex
  genes_df <- genetex::genes_regex() # this is a df of genes names concatenated with "|" (no word boundaries; "greedy")
  genes_regex <- genes_df$genes # this is a character string of the genes for our regex
  gain_loss_amplification_regex <- "gain|loss|amplification"
  genes_gain_loss_amplification_regex <- paste(genes_boundary_regex, gain_loss_amplification_regex, sep = "|")


  ##########################################################################################################################
  # SNaPshot
  ##########################################################################################################################
  cnv.snap.1 <- dt %>%
    tidyr::separate(
      col = Results,
      into = c("text","cnv"))
  # use regex to look for those with either gains or losses
  cnv.snap.2 <- cnv.snap.1 %>%dplyr::filter(stringr::str_detect(string = cnv, pattern = ("gain|loss"), negate = FALSE))
  # filter to make sure those rows in "text" are gene names
  cnv.snap.2a <- cnv.snap.2 %>% dplyr::filter(stringr::str_detect(string = cnv.snap.2$text, pattern = stringr::regex(genes_regex)))
  ###### Now let's modify "cnv.snap.3" to set it up with the appropriate variable names for the import_tool
  cnv.snap.3 <- if(is.na(cnv.snap.2a[1,1])) base::data.frame(text = NA, cnv = NA) else cnv.snap.2a # adding blank rows for those reports that have no CNVs
  cnv.snap.3 <- cnv.snap.3 %>% dplyr::mutate(cnv_gene = "cnv_gene")
  cnv.snap.3 <- cnv.snap.3 %>% dplyr::mutate(cnv_gain_or_loss = "cnv_gain_or_loss")
  cnv.snap.3 <- cnv.snap.3 %>% dplyr::mutate(number = 1:nrow(cnv.snap.3))
  cnv.snap.4a <- cnv.snap.3 %>% dplyr::mutate(cnv_gene = paste(cnv_gene, number, sep = "_"))
  cnv.snap.4 <- cnv.snap.4a %>% dplyr::mutate(cnv_gain_or_loss = paste(cnv_gain_or_loss, number, sep = "_"))
  ####### rename variables to set up for combining with import_tool
  ######## cnv.snap.gene is a df that contains the name of the gene side-by-side with the variable name
  cnv.snap.gene <- cnv.snap.4 %>% dplyr::select(text, cnv_gene) %>% dplyr::rename(variables = cnv_gene) %>%
    dplyr::rename(results = text) %>% dplyr::select(variables, results)
  ######## cnv.snap.gene.gain.or.loss is a df that contains the type of cnv side-by-side with the variable name
  cnv.snap.gene.gain.or.loss <- cnv.snap.4 %>% dplyr::select(cnv, cnv_gain_or_loss) %>% dplyr::rename(variables = cnv_gain_or_loss) %>% dplyr::rename(results = cnv) %>% dplyr::select(variables, results)
  # combine the above two dfs to form final df
  cnv.final <- base::rbind(cnv.snap.gene, cnv.snap.gene.gain.or.loss)



  # Add cnv_number to data frame
  cnv.number.1 <- stringr::str_detect(string = cnv.final$variables, pattern = "cnv_gene")
  cnv.number.2 <- sum(cnv.number.1)
  cnv.number <- if(is.na(cnv.final[1,2])) 0 else cnv.number.2
  cnv.number.df <- data.frame(variables = "cnv_number",
                              results = cnv.number)
  # Add amplifications_number to data frame
  amp.number.1 <- stringr::str_detect(string = cnv.final$variables, pattern = "amplifications_gene")
  amp.number.2 <- sum(amp.number.1)
  amp.number <- if(is.na(cnv.final[1,2])) 0 else amp.number.2
  amp.number.df <- data.frame(variables = "amplifications_number",
                              results = amp.number)

  # add cnv.number.df and amp.number.df to cnv.final
  cnv.final <- rbind(cnv.final,
                     cnv.number.df,
                     amp.number.df)
  cnv.final <- cnv.final %>% tidyr::drop_na()
  # change to character vector
  cnv.final$results <- base::as.character(cnv.final$results)

  return(cnv.final)
}
