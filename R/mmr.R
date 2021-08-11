#' Abstract Microsatellite data from genomic reports
#' @description
#' `mmr()` provides NLP tools to text mine mismatch repair status from genomic reports and transforms it to structured data for import into REDCap
#' @description
#' `genetex_to_redcap()` provides natural language processing tools to abstract data from a variety of genomic reports and import them to REDCap
#' @param data The data frame of the genomic report of interest. This can be copied to the Clipboard from the EHR report by selecting the relevant text and then "copy". That is the default. If you don't use the clipboard function, you can use a data frame of the text file. The single column data frame should be labeled "Results". Required.
#'
#' @return a data frame with two columns: variables (the redcap variable names) and results (the data to be imported in to redcap)
#' @export
#'
mmr <- function(data = dplyr::tibble(Results = readr::clipboard())){

  ##########################################################################################################################
  # load data
  ##########################################################################################################################
  dt <- data

  ##########################################################################################################################
  # Initial Filter for MMR
  ##########################################################################################################################
  ###1, Proficient/Stable (1)
  ###2, Deficient/Unstable (2)
  ###3, Indeterminate (3)


  ## Filter for only the rows that have regex specific for microsatellite instablity or MMR
  mmr.pre.1 <- dt %>%
    dplyr::filter(stringr::str_detect(string = Results,
                                      pattern = regex("mmr|mss|msi-high|ms-stable", ignore_case = TRUE)))
  ## Continue to clean by removing rows that have words that designate as undesired rows
  mmr.pre.2 <- mmr.pre.1 %>%
    filter(!str_detect(string = mmr.pre.1$Results,
                       pattern = regex("Guardant360|pembrolizumab|samples", ignore_case = TRUE)))
  ## filter for only those rows that have the output of interest
  mmr.pre.3 <- mmr.pre.2 %>%
    dplyr::filter(stringr::str_detect(string = Results,
                                      pattern = regex("mmr-p|detected|not detected|Microsatellite status - MS-Stable|Microsatellite status MSI-High", ignore_case = TRUE)))
  ## Add a row to account for missing data
  mmr.pre.4 <- if(base::is.na(mmr.pre.3[1,1])) data.frame(Results = NA) else mmr.pre.3

  ## create an indicator variable depending on whether or not stable or instable or MSI
  mmr_status <- ifelse(test = str_detect(mmr.pre.4$Results, regex("Stable|MSS|Not Detected", ignore_case = TRUE)),
                       yes = 1,
                       no = ifelse(test = str_detect(mmr.pre.4$Results, regex("Instab|MSI_H|MSI-H|detected|MSI", ingore_case = TRUE)),
                                   yes = 2,
                                    no = ifelse(test = ifelse(test = str_detect(mmr.pre.4$Results, regex("indeterminate", ignore_case = TRUE)),
                                               yes = 3,
                                               no = NA))))
  mmr <- data.frame(variables = "mmr",
                    results = mmr_status
                    )

  mmr <- mmr[1,]

  mmr$results <- base::as.character(mmr$results)

  return(mmr)

}
