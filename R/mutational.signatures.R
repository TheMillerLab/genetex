#' Apply NLP to genomic reports to extract out mutational signatures
#' @description
#' `mutational.signatures()` text mines mutational signatures data from a variety of genomic reports and transforms it to structured data for import into REDCap
#' @param data The data frame of the genomic report of interest. This can be copied to the Clipboard from the EHR report by "selecting all" and then "copy". That is the default. If you don't use the clipboard function, you can use a data frame of the text file. The single column data frame should be labeled "Results". Required.
#'
#' @return a data frame with two columns: variables (the redcap variable names) and results (the data to be imported in to redcap)
#' @export
#'
mutational.signatures <- function(data = dplyr::tibble(Results = readr::clipboard())){
  ##########################################################################################################################
  # load data
  ##########################################################################################################################
  dt <- data

  ##########################################################################################################################
  # Mutational Signatures
  ##########################################################################################################################
  ## Let's create a df for the different types of mutational signatures that are found Oncopanel
  ### Of note, we will have to continually update these as specimens with these mutations are found, as the key strings
  #### that we use "str_detect()" to find are contained within these reports. To date, we just have UVA and "Too few
  ##### mutations"
  uv.tf <- str_detect(string = dt$Results,   # dt.2 is the starting df that we've pre-processed
                      pattern = "UVA signature") # "UVA signature is a string found in the Oncopanel reports

  ## Eventually we want a df that we can tally so we know how many "mutation_signature_number" we need to select
  ### Let's create an indicator df that tells us if there are any signatures present
  uv_sum <- sum(uv.tf)
  ### this df will serve as that; we will have each signature in a cell and either a 0 or 1 in the "sum" column
  #### we can then sum that "sum" category to get the number for "mutation_signature_number"
  uv_df <- data.frame(
    signature = "uvr",
    sum = if(uv_sum >=1) 1 else 0
    )

  ## Perform an analogous set of steps for the "Too Few Mutations" group, which is an entry in REDCAp
  too_few.tf <- str_detect(string = dt$Results,
                           pattern = "Too few mutations")

  too_few_sum <- sum(too_few.tf)
  too_few_df <- data.frame(
    signature = "too_few",
    sum = if(too_few_sum >=1) 1 else 0
  )

  ## Let's bind the various signature dfs into "signatures_df", which will be the df that we sum "sum from for the
  ### "mutation_signature_number" variable
  ### *******This will require re-visiting as we add new signatures to the package*******
  signatures_df <- rbind(
    uv_df,
    too_few_df   ######## We are going to need to add new dfs in the future, e.g. "alkylating_df", they'll be added here
  )
  ## create a new df so we can add mutation_signature_number or mutation_signature_*
  mut_signature_df <- base::data.frame(results = NA,
                                       variables = c("mutation_signature_number",
                                                     "mutation_signature_1",
                                                     "mutation_signature_2"))
  ## use the sum() function of "signatures_df" for "mutation_signature_number"
  mut_signature_df[1,1] <- sum(signatures_df$sum)

  ## Now we will add values to the "mutation_signature" variables based on the above data cleaning
  ### Given that UV is the most common in MCC, let's have that as the first signature if it is present
  mut_signature_df[2,1] <- ifelse(test = uv_sum >= 1, # so if the uv_sum integer is >=1, the signaure will be UVR
                                         yes = 7, # 7 is the UVR signature
                                         no = ifelse(test = too_few_sum >= 1, # if UVR is not present, check for too-few
                                                     yes = 99, # 99 is the redcap code for "too-few mutations"
                                                     no = NA) # we will need to add additional options as we add more
                                         # signatures
  )
  ### if there is more than one signature present, we can add another. Since UVR will have been filled in if present
  #### we can start with the "Too-Few Mutations" entry
  mut_signature_df[3,1] <- if(too_few_sum >= 1) 99 else NA
  mut_signature_df$results <- base::as.character(mut_signature_df$results)

  return(mut_signature_df)
}
