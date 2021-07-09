#' Apply NLP to genomic reports to extract out tumor mutation burden
#' @description
#' `tmb()` text mines tumor mutation burden (TMB) data from a variety of genomic reports and transforms it to structured data for import into REDCap
#' @param data The data frame of the genomic report of interest. This can be copied to the Clipboard from the EHR report by selecting the relevant text and then "copy". That is the default. If you don't use the clipboard function, you can use a data frame of the text file. The single column data frame should be labeled "Results". Required.
#' @param platform String of the platform used to generate the genomics report. Acceptable strings include: MGH, SNaPshot, BWH, Oncopanel, Guardant, Guardant 360, Foundation or fmi. Case is not sensitive. Required.
#'
#' @return a data frame with two columns: variables (the redcap variable names) and results (the data to be imported in to redcap)
#' @export
#'
tmb <- function(data = dplyr::tibble(Results = readr::clipboard()),
                platform){

  ##########################################################################################################################
  # load data
  ##########################################################################################################################
  dt <- data

  ##########################################################################################################################
  # Initial Filter for Tumor Mutational Burden
  ##########################################################################################################################
  ## Filter for only the rows that have "tumor mutational burden" in it
  tmb.pre.1 <- dt %>%
    dplyr::filter(stringr::str_detect(string = dt$Results,
                                      pattern = regex("tumor mutational burden", ignore_case = TRUE)))
  ## Filter for only the rows that have numbers in the string
  tmb.pre.2 <- tmb.pre.1 %>%
    dplyr::filter(stringr::str_detect(string = tmb.pre.1$Results,
                                      pattern = "\\d*\\.?\\d+")) # match only the numbers in this string
  ## Continue to clean by removing rows that have words that designate as undesired rows
  tmb.pre.3 <- tmb.pre.2 %>%
    filter(!str_detect(string = tmb.pre.2$Results,
                       pattern = regex("corresponds|microsatellite|et al|version|amplification", ignore_case = TRUE)))
  ## Extract the numbers from this row
  tmb.pre.4 <- stringr::str_extract(string = tmb.pre.3$Results,
                                    pattern = "\\d*\\.?\\d+") # match only the numbers in this string
  ### \d* matches a digit (equivalent to [0-9])
  ### * matches the previous token between zero and unlimited times, as many times as possible, giving back as needed (greedy)
  ### \.? matches the character . literally (case sensitive)
  ### ? matches the previous token between zero and one times, as many times as possible, giving back as needed (greedy)
  ### \d+ matches a digit (equivalent to [0-9]))

  ## Now we have to create a data frame, since now all we have is a character string of the numbers
  tmb.pre.5 <- data.frame(tmb = tmb.pre.4)

  ## Some samples may not have tmb assessed,
  ### Therefore, there is problem of a potential empty
  #### We'll circumvent this problem df by creating a df with NA in tmb if there are no values in tmb.pre.5
  tmb.pre.6 <- if(is.na(tmb.pre.5[1,1])) base::data.frame(tmb = NA) else tmb.pre.5
  #### slice_head
  tmb <- tmb.pre.6 %>% dplyr::slice_head()

  ## Now let's assign the appropriate variable "tmb_abs" or "tmb" based on the platform
  ### The numeric will come from the genetex::platform() which generates a numeric from the string of the platform
  platform_num <- genetex::platform(platform)
  tmb$variables <- if(platform_num == 2) "tmb_abs" else "tmb"
  # add row to fill in field "tmb_yn": is tmb available (0 = No, 1 = yes)
  tmb <- tmb %>% dplyr::add_row()
  tmb[2,1] <- if(is.na(tmb[1,1])) 0 else 1
  tmb[2,2] <- "tmb_yn"
  tmb <- tmb %>% dplyr::rename(results = tmb)
  tmb <- tmb %>% dplyr::select(variables, results)

  return(tmb)

}
