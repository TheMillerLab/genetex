#' Apply NLP to genomic reports to extract copy number variants (CNVs) from foundation medicine reports
#' @description
#' `cnv.fmi()` provides natural language processing tools to abstract copy number variants data from a variety of genomic reports
#' @param data The data frame of the genomic report of interest. This can be copied to the Clipboard from the EHR report by electing the relevant text and then "copy". That is the default. If you don't use the clipboard function, you can use a data frame of the text file. The single column data frame should be labeled "Results". Required.
#' @return a data frame with two columns: variables (the redcap variable names) and results (the data to be imported in to redcap)
#' @export
#'
cnv.fmi <- function(data = dplyr::tibble(Results = readr::clipboard())){
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
  gain_loss_amplification_regex <- "gain|loss|amplification"
  genes_gain_loss_amplification_regex <- paste(genes_boundary_regex, gain_loss_amplification_regex, sep = "|")

  ##########################################################################################################################
  # fmi
  ##########################################################################################################################
  fmi.1 <- dt
  fmi.2 <- splitstackshape::cSplit(
    indt = fmi.1,
    splitCols = "Results",
    sep = " ",
    direction = "long"
  )
  fmi.2a <- fmi.2 %>% dplyr::mutate(divider = stringr::str_detect(string = fmi.2$Results, pattern = "appendix|ACTIONABILITY"),
                                    group = cumsum(divider)) %>% dplyr::filter(group == 0)
  fmi.3 <- fmi.2a %>% dplyr::filter(stringr::str_detect(string = fmi.2a$Results, pattern = stringr::regex(genes_gain_loss_amplification_regex, ignore_case = TRUE)))
  fmi.4 <- fmi.3
  fmi.4$Results <- stringr::str_replace(string = fmi.4$Results, pattern = "§", replacement = "")
  fmi.5 <- fmi.4 %>% dplyr::select(Results)
  fmi.5 <- fmi.5 %>% dplyr::mutate(divider = stringr::str_detect(string = fmi.5$Results, pattern = genes_boundary_regex),
                                   group = cumsum(divider))
  ## Now let's collapse these into a another column
  fmi.6 <- fmi.5 %>% dplyr::group_by(group) %>% dplyr::mutate(results_group_paste = paste0(Results, collapse = ", ")) %>%
    dplyr::ungroup()
  # Keep only those with a ","
  fmi.7 <- fmi.6 %>% dplyr::filter(stringr::str_detect(string = fmi.6$results_group_paste, pattern = ","))
  # remove duplicates
  fmi.8 <- fmi.7 %>% dplyr::group_by(group) %>% dplyr::slice_head() %>% dplyr::ungroup()
  # remove the ","
  fmi.9 <- fmi.8
  fmi.9$results_group_paste <- stringr::str_replace(string = fmi.9$results_group_paste, pattern = ",", replacement = "")
  # sort "gain|loss" from "amplification
  fmi.10 <- fmi.9 %>% dplyr::filter(stringr::str_detect(string = fmi.9$results_group_paste, pattern = "gain|loss"))
  fmi.10a <- if(is.na(fmi.10[1,1])) data.frame(Results = NA, divider = NA,group = NA, results_group_paste = NA) else fmi.10
  fmi.11 <- fmi.10a %>% dplyr::mutate(number = 1:nrow(fmi.10a)) %>% dplyr::select(results_group_paste, number) %>%
    dplyr::rename(results = results_group_paste)
  fmi.12 <- splitstackshape::cSplit(
    indt = fmi.11,
    splitCols = "results",
    sep = " ",
    direction = "long"
  )
  fmi.13 <- fmi.12 %>% dplyr::mutate(variables = ifelse(test = stringr::str_detect(string = fmi.12$results,
                                                                                   pattern = stringr::regex(genes_boundary_regex)),
                                                        yes = "cnv_gene",
                                                        no = ifelse(test = stringr::str_detect(string = fmi.12$results,
                                                                                               pattern = "gain|loss"),
                                                                    yes = "cnv_gain_or_loss",
                                                                    no = NA)))
  fmi.14 <- fmi.13 %>% dplyr::mutate(variables = if(is.na(fmi.13$variables)) NA else paste(variables, number, sep = "_")) %>%
    dplyr::select(variables, results)
  # wrangle the amplifications rows
  fmi.15 <- fmi.9 %>% dplyr::filter(stringr::str_detect(string = fmi.9$results_group_paste, pattern = "amplification"))
  fmi.15a <- if(is.na(fmi.15[1,1])) data.frame(Results = NA, divider = NA,group = NA, results_group_paste = NA) else fmi.15
  fmi.16 <- fmi.15a %>% dplyr::mutate(number = 1:nrow(fmi.15a)) %>% dplyr::select(Results, number) %>%
    dplyr::rename(results = Results)
  fmi.17 <- fmi.16 %>% dplyr::mutate(variables = if(is.na(fmi.16$results)) NA else paste("amplifications_gene",number, sep = "_")) %>%
    dplyr::select(variables, results)
  cnv.final <- rbind(fmi.14, fmi.17)

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
