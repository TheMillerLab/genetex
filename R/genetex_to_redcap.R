#' Apply NLP to genomic reports to import into the Genomics Instrument in REDCap (https://www.themillerlab.io/post/optimizing_rwd_collection-genomics_instrument/)
#' @description
#' `genetex_to_redcap()` provides natural language processing tools to abstract data from a variety of genomic reports and import them to REDCap
#' @param data The data frame of the genomic report of interest. This can be copied to the Clipboard from the EHR report by "selecting all" and then "copy". That is the default. If you don't use the clipboard function, you can use a data frame of the text file. The single column data frame should be labeled "Results". Required.
#' @param record_id A string of the the record_id of the subject you want to import (e.g. "John Doe" or "1-1"). Required.
#' @param instrument_instance A numeric representing the instrument instance for the genomics instrument (e.g. 1). Required.
#' @param platform String of the platform used to generate the genomics report. Acceptable strings include: MGH, SNaPshot, BWH, Oncopanel, Guardant, Guardant Health, Foundation or fmi. Case is not sensitive. Required.
#' @param lesion_tag A string of the lesion tag for the genomics instrument (e.g. "Right Arm Primary Lesion"). Optional.
#' @param genomics_tissue_type A string of the genomics tissue type used in the report (e.g. "primary cutaneous lesion" or "MCCUP" or "met"). Default is blank. Optional.
#' @param date_collected A string of the date the genomic sample was collected in YYYY-MM-DD format (e.g. "2020-12-31). Optional.
#' @param redcap_uri A string of the url for the REDCap system you are working with.
#' @param redcap_api_token A user-specific alphanumeric string that serves as the password for the REDCap project. The default is blank, which means the data will not be sent to REDCap. Optional.
#' @return A data frame that can be imported into the Genomics Instrument in REDCap (https://www.themillerlab.io/post/optimizing_rwd_collection-genomics_instrument/)
#' @export
#' @examples
#' snapshot_sample_report %>%
#'  genetex_to_redcap(
#'   record_id = "Jane Doe",
#'   instrument_instance = 1,
#'   genomics_tissue_type = "Skin",
#'   platform = "mgh",
#'   date_collected = "2021-01-01",
#'   lesion_tag = "Right Arm Primary",
#'   )
genetex_to_redcap <- function(
  data = dplyr::tibble(Results = readr::clipboard()), # allow clipboard to bring in data as default
  record_id,
  instrument_instance,
  platform,
  lesion_tag = "",
  genomics_tissue_type = "",
  date_collected = "",
  redcap_uri,
  redcap_api_token = ""
  )
  {

  ##########################################################################################################################
  # Tumor Mutational Burden
  ##########################################################################################################################
  ## Use the tmb() function to isolate tmb data from the report
  tmb <- genetex::tmb(data = data, platform = platform)

  ##########################################################################################################################
  # MMR/MSI
  ##########################################################################################################################
  ## Use the mmr() function to isolate mismatch repair status data from the report
  mmr <- genetex::mmr(data = data)

  ##########################################################################################################################
  # Gene Variants
  ##########################################################################################################################
  ## Use the gene.variants() function to isolate those genes with genetic variants from the report
  gene_var <- genetex::gene.variants(data = data, platform = platform)

  ##########################################################################################################################
  # Mutational Signatures
  ##########################################################################################################################
  ## Use the mutational_signatures() function to isolate mutational signature(s) data from the report
  mut_signatures <- genetex::mutational.signatures(data = data)

  ##########################################################################################################################
  # Structural Variants (e.g. CNVs)
  ##########################################################################################################################
  ## Use the cnv() function to isolate any structural changes data from the report
  cnv <- genetex::cnv(data = data, platform = platform)


  ##########################################################################################################################
  # Platform
  ##########################################################################################################################
  ## Use the platform() function to assign the appropriate numeric for the platform used to generate the report
  platform_num <- genetex::platform(platform = platform)
  platform_num <- as.character(platform_num)

  ##########################################################################################################################
  # Genomics Tissue Type
  ##########################################################################################################################
  ## Use the genomics_tissue_type() function to assign the appropriate numeric of the type of tissue used to generate
  ### the report
  gtt_num <- genetex::genomics.tissue.type(genomics_tissue_type = genomics_tissue_type)
  gtt_num <- as.character(gtt_num)

  ##########################################################################################################################
  # Additional Info for Import
  ##########################################################################################################################
  ## Add direct data to the "genomics" instrument
  add_info <- dplyr::tribble(
    ~variables, ~results,
    "redcap_repeat_instrument", "genomics",
    "redcap_repeat_instance", as.character(instrument_instance),
    "lesion_tag_genomics", lesion_tag,
    "genomics_date",date_collected,
    "genomics_platform", platform_num,
    "genomics_tissue_type", gtt_num,
    "genomics_yn", "1",
    "record_id", record_id,
    "genomics_add_notes", "Data imported via genomics package",
    "gen_qcdash___dac", "1"
  )
  add_info$results <- as.character(add_info$results)

  ##########################################################################################################################
  # Load Import Tool (of note, if this changes, one needs to reload the new file into the "data" file)
  ##########################################################################################################################
  import_tool <- genetex::genomics_import_tool

  ##########################################################################################################################
  # Build genetex-to-redcap ("gtr") Data Frames
  ##########################################################################################################################
  gtr.1 <- left_join(import_tool, tmb, by = "variables")
  gtr.2 <- left_join(gtr.1, mmr, by = "variables")
  gtr.3 <- left_join(gtr.2, gene_var, by = "variables")
  gtr.4 <- left_join(gtr.3, mut_signatures, by = "variables")
  gtr.5 <- left_join(gtr.4, cnv, by = "variables")
  gtr.6 <- left_join(gtr.5, add_info, by = "variables")
  ##########################################################################################################################
  # Now let's unite all of the columns that are not "variable" to form just two columns
  ##########################################################################################################################
  dt_unite <- gtr.6 %>%
    tidyr::unite(col = "results",
                 2:length(gtr.6),
                 na.rm = TRUE,
                 remove = TRUE)

  ##########################################################################################################################
  # Reconcile Gene Names
  ##########################################################################################################################
  ## Because there are several different names for genes, and various panels uses some different ones, we need to reconcile
  ### other gene names to match what we have in REDCAp
  dt_unite <- genetex::reconcile_gene_names(data = dt_unite)

  ##########################################################################################################################
  # Now let's pivot wider to form the table that we will import
  dt_wide <- pivot_wider(data = dt_unite,
                         names_from = "variables",
                         values_from = "results")

  ##########################################################################################################################
  # Import Data into REDCap
  ##########################################################################################################################
  ifelse(test = redcap_api_token == "",
         yes = "",
         no = REDCapR::redcap_write_oneshot(ds = dt_wide,
                                            redcap_uri = redcap_uri,
                                            token = redcap_api_token))


  return(dt_unite)
}
