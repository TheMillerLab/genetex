#' Use text mining to assing numeric to genomics platforms
#' @description
#' `platform()` applies regular expressions to assign a numerical value for the various platforms used for genomic reports that aligns with the  "genomics_platform" field in the REDCAp Genomics Instrument
#' @param platform String of the platform used to generate the genomics report. Acceptable strings include: MGH, SNaPshot, BWH, Oncopanel, Guardant, Guardant Health, Foundation or fmi. Case is not sensitive. Required.
#'
#' @return numeric value corresponding to the various platforms used to generate genomic reports
#' @export
platform <- function(platform){


  ##########################################################################################################################
  # Key for REDCap Genomics Instrument Field "genomics_platform"
  ##########################################################################################################################
  ## 1, BWH Oncopanel (1)
  ## 2, MGH SNaPshot (2)
  ## 3, MSK-IMPACT (3)
  ## 4, FoundationOne CDx (4)
  ## 5, TEMPUS xT Gene Panel (5)
  ## 6, Guardant 360 (6)
  ##########################################################################################################################
  # Platform
  ##########################################################################################################################

  platform_number <- ifelse(
    test = stringr::str_detect(platform, stringr::regex("BWH OncoPanel|BWH|Oncopanel", ignore_case = TRUE)),
    yes = 1,
    no = ifelse(test = stringr::str_detect(platform, stringr::regex("MGH SNaPShot|mgh|snapshot", ignore_case = TRUE)),
                yes = 2,
                no = ifelse(test = stringr::str_detect(platform, stringr::regex("msk|impact|MSK-IMPACT", ignore_case = TRUE)),
                            yes = 3,
                            no = ifelse(test = stringr::str_detect(platform, stringr::regex("foundation|fmi|FoundationOne CDx", ignore_case =   TRUE)),
                                        yes = 4,
                                        no = ifelse(test = stringr::str_detect(platform, stringr::regex("tempus|TEMPUS xT Gene Panel",   ignore_case = TRUE)),
                                                    yes = 5,
                                                    no = ifelse(test = stringr::str_detect(platform, stringr::regex("guardant|Guardant 360",   ignore_case = TRUE)),
                                                                yes = 6,
                                                                no = "")
                                        )
                            )
                )
    )
  )

  return(platform_number)
}
