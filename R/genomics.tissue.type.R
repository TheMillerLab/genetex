#' Assing the appropriate numeric to a string designating the tissue type sequence according to the MCC Patient Registry
#' @description
#' `genomics.tissue.type()` applies regular expressions to assign a numerical value for the various tissues that may be used for genomic analysis that aligns with the  "genomics_tissue_type" field in the REDCAp Genomics Instrument
#'
#' @param genomics_tissue_type A string of the genomics tissue type used in the report (e.g. "primary cutaneous lesion" or "MCCUP" or "met"). Default is blank. Optional.
#'
#' @return A numeric
#' @export
#'
genomics.tissue.type <- function(genomics_tissue_type = ""){
  ##########################################################################################################################
  # Key for REDCap Genomics Instrument Field "genomics_tissue_type"
  ##########################################################################################################################
  ## 1, Primary Cutaneous Tumor (1)
  ## 2, Metastases (2)
  ## 3, MCC of Unknown Primary (non-cutaneous lesion at initial presentation) (3)
  ## 4, Local Recurrence (4)
  ## 5, Blood/Liquid Biopsy (5)
  ## 98, Unknown/Not Reported (98)
  ##########################################################################################################################
  # Genomics Tissue Type
  ##########################################################################################################################
  gtt <- genomics_tissue_type
  gtt_number <- ifelse(test = str_detect(gtt, regex("unknown primary", ignore_case = TRUE)),
                       yes = 3,
                       no = ifelse(test = str_detect(gtt, regex("pct", ignore_case = TRUE)),
                                   yes = 1,
                                   no = ifelse(test = str_detect(gtt, regex("met", ignore_case = TRUE)),
                                               yes = 2,
                                               no = ifelse(test = str_detect(gtt, regex("MCCUP", ignore_case = TRUE)),
                                                           yes = 3,
                                                           no = ifelse(test = str_detect(gtt, regex("prim", ignore_case = TRUE)),
                                                                       yes = 1,
                                                                       no = ifelse(test = str_detect(gtt, regex("local", ignore_case = TRUE)),
                                                                                   yes = 4,
                                                                                   no = ifelse(test = str_detect(gtt, regex("lr", ignore_case = TRUE)),
                                                                                               yes = 4,
                                                                                               no = ifelse(test = str_detect(gtt, regex("blood", ignore_case = TRUE)),
                                                                                                           yes = 5,
                                                                                                           no = ifelse(test = str_detect(gtt, regex("liquid", ignore_case = TRUE)),
                                                                                                                       yes = 5,
                                                                                                                       no = ""))) ))))))

  return(gtt_number)

}
