#' Apply NLP to genomic reports to extract copy number variants (CNVs)
#' @description
#' `cnv()` integrates various platform-specific NLP functions to text mine gene names and copy number variants data from a variety of genomic reports and transforms them to structured data for import into REDCap
#' @param data The data frame of the genomic report of interest. This can be copied to the Clipboard from the EHR report by "selecting all" and then "copy". That is the default. If you don't use the clipboard function, you can use a data frame of the text file. The single column data frame should be labeled "Results". Required.
#' @param platform String of the platform used to generate the genomics report. Acceptable strings include: MGH, SNaPshot, BWH, Oncopanel, Guardant, Guardant 360, Foundation or fmi. Case is not sensitive. Required.
#' @return a data frame with two columns: variables (the redcap variable names) and results (the data to be imported in to redcap)
#' @export
#'
cnv <- function(data = dplyr::tibble(Results = readr::clipboard()),
                platform){
  ##########################################################################################################################
  # load data
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
                      yes = cnv.final <- cnv.onco(data = dt),
                      no = ifelse(test = platform_num ==2,
                                  yes = cnv.final <- cnv.snap(data = dt),
                                  no = ifelse(test = platform_num == 4,
                                              yes = cnv.final <- cnv.fmi(data = dt),
                                              no = ifelse(test = platform_num ==6,
                                                          yes = cnv.final <- cnv.guardant(data = dt),
                                                          no = cnv.final <- base::data.frame(variables = c("cnv_number","amplifications_number"),
                                                                                             results = c("0","0"))))))

  # change to character vector
  cnv.final$results <- base::as.character(cnv.final$results)

  return(cnv.final)
}

