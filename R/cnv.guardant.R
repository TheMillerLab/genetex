#' Apply NLP to genomic reports to extract copy number variants (CNVs) from Guardant Health reports
#' @description
#' `cnv.guardant()` provides natural language processing tools to abstract copy number variants data Guardant Health reports
#' @param data The data frame of the genomic report of interest. This can be copied to the Clipboard from the EHR report by selecting the relevant text and then "copy". That is the default. If you don't use the clipboard function, you can use a data frame of the text file. The single column data frame should be labeled "Results". Required.
#' @return a data frame with two columns: variables (the redcap variable names) and results (the data to be imported in to redcap)
#' @export
#'
cnv.guardant <- function(data = dplyr::tibble(Results = readr::clipboard())){
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
  genes_amplification_regex <- paste("(",genes_regex,")"," Amplification", sep = "")


  ##########################################################################################################################
  # Text Mine for Amplifications
  ##########################################################################################################################
  dt.1 <- dt %>%
    dplyr::filter(stringr::str_detect(string = dt$Results,
                                      pattern = stringr::regex("gain|loss|amplification", ignore_case = TRUE), negate = FALSE))


  dt.2 <- dt.1 %>%
    dplyr::filter(stringr::str_detect(string = dt.1$Results,
                                      pattern = stringr::regex(genes_amplification_regex)))

  dt.3 <- splitstackshape::cSplit(
    indt = dt.2,
    splitCols = "Results",
    sep = " ",
    direction = "Long"
  )

  dt.4 <- dt.3 %>%
    dplyr::filter(stringr::str_detect(string = dt.3$Results, pattern = stringr::regex(genes_gain_loss_amplification_regex, ignore_case = TRUE)))

  dt.5 <- dt.4 %>%
    dplyr::mutate(divider = stringr::str_detect(string = dt.4$Results, pattern = stringr::regex(genes_regex)),
                  group = cumsum(divider))
  dt.6 <- dt.5 %>% dplyr::count(group) %>% dplyr::filter(n == 2)
  dt.7 <- dt.5 %>% dplyr::filter(dt.5$group %in% dt.6$group)
  dt.8 <- dt.7 %>% dplyr::group_by(group) %>% dplyr::mutate(merge = paste0(Results, collapse = " ")) %>% dplyr::ungroup()
  dt.9 <- dt.8 %>% dplyr::group_by(merge) %>% dplyr::slice_head() %>% dplyr::ungroup()

  # Select only those with amplifications
  dt.10 <- dt.9 %>% dplyr::filter(stringr::str_detect(string = dt.9$merge, pattern = "Amplification"))
  dt.10a <- if(is.na(dt.10[1,1])) data.frame(variables = "amplifications_gene", Results = NA) else dt.10
  dt.11 <- dt.10a %>% dplyr::select(Results) %>% mutate(variables = "amplifications_gene", number = 1:nrow(dt.10a)) %>%
    dplyr::mutate(variables = paste(variables, number, sep = "_")) %>%
    dplyr::rename(results = Results) %>% dplyr::select(variables, results)
  dt.11$results <- as.character(dt.11$results)


  # Add amplifications_number to data frame
  amp.number.1 <- stringr::str_detect(string = dt.11$variables, pattern = "amplifications_gene")
  amp.number.2 <- sum(amp.number.1)
  amp.number <- if(is.na(dt.11[1,2])) 0 else amp.number.2
  amp.number.df <- data.frame(variables = "amplifications_number",
                              results = amp.number)
  amp.number.df$results <- as.character(amp.number.df$results)

  # add df for cnv_gene (which is NA)
  cnv <- data.frame(variables = "cnv_number", results = "0")

  # add cnv.number.df and amp.number.df to cnv.final
  cnv.final <- rbind(dt.11, amp.number.df, cnv)

  cnv.final <- cnv.final %>% tidyr::drop_na()

  return(cnv.final)
}
