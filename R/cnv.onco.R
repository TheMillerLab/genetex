#' Apply NLP to genomic reports to extract copy number variants (CNVs) from oncopanel reports
#' @description
#' `cnv()` provides natural language processing tools to abstract copy number variants data from oncopanel genomic reports
#' @param data The data frame of the genomic report of interest. This can be copied to the Clipboard from the EHR report by selecting the relevant text and then "copy". That is the default. If you don't use the clipboard function, you can use a data frame of the text file. The single column data frame should be labeled "Results". Required.
#' @return a data frame with two columns: variables (the redcap variable names) and results (the data to be imported in to redcap)
#' @export
#'
cnv.onco <- function(data = dplyr::tibble(Results = readr::clipboard())){
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
  nuc_regex <- "[ACTG]>[ACTG]|del[ACTG]" # This regex should identify those rows that have nucleotide changes
  aa_regex <- "(\\b([A-Z][0-9]{1,}(([A-Z])|(_[A-Z][1-9]{1,}del)|(fs\\*[1-9]{1,})|(\\*)|(fs)|(del)))|(p\\.[A-Z]))|([0-9]ins[A-Z])"
  # This regex will detect all the possible AA alterations with single initials
  genes_boundary_nuc_aa_regex <- paste(genes_boundary_regex, nuc_regex, aa_regex, sep = "|") # concatenate several regexs
  # Assemble a regex_discard to perform an initial screen
  regex_discard <- "H3K27|Sf3b1|H3B-6545|\\bfusion\\b|\\bthe\\b|\\bincluding\\b|\\bcells\\b|\\bcancer\\b"
  gain_loss_amplification_regex <- "gain|loss|amplification"
  genes_gain_loss_amplification_regex <- paste(genes_boundary_regex, gain_loss_amplification_regex, sep = "|")

  ##########################################################################################################################
  # Initial cleaning and tokenization
  ##########################################################################################################################
  dt.1a <- dt
  dt.1a$Results <- stringr::str_trim(string = dt.1a$Results) # Trim white space
  # We need to eliminate much of this report in order to mine only those gene alterations relevant to each subject
  ## So we will use the word "references" as a divider and take only those rows that are above it
  references_interpretation_definitions <- "references|interpretation|definitions"
  dt.1b <- dt.1a %>%
    mutate(divider = stringr::str_detect(string = dt.1a$Results,
                                         pattern = stringr::regex(references_interpretation_definitions, ignore_case = TRUE)),
           group = base::cumsum(divider)) %>%
    dplyr::filter(group == 0)
  # Perform a clean to remove rows that have problematic strings
  dt.1c <- dt.1b %>% dplyr::filter(!stringr::str_detect(string = dt.1b$Results,
                                                        pattern = stringr::regex(regex_discard, ignore_case = FALSE)))
  # remove ?
  dt.1d <- dt.1c
  dt.1d$Results <- stringr::str_replace_all(string = dt.1c$Results, pattern = "\\?", replacement = "")
  # tokenize the vector with space
  dt.1 <- splitstackshape::cSplit(
    indt = dt.1d,
    splitCols = "Results",
    sep = " ",
    direction = "long"
  )

  # now let's create an object with a divider and a group based on that divider
  ## To do that, we will create two additional vectors, first is a logical vector "divider" that will result TRUE if the
  ### token in the "word" vector is contained in the "keyword" vector. The "group" vector, then will perform a cumulative
  #### sum of the "divider" vector, such that all those tokens that follow a specific keyword will have the same "group",
  ##### we can then use this group to tag those gene names with the specific type of alteration they had
  dt.2 <- dt.1 %>% dplyr::mutate(divider = stringr::str_detect(string = dt.1$Results,
                                                               pattern = stringr::regex(gain_loss_amplification_regex, ignore_case = TRUE)),
                                 group = base::cumsum(divider))
  # Now let's create an important column that will serve as a look up table of sorts for the specific group that is
  ## associated with the gene alteration (e.g. "Gain" or "Loss")
  dt.3 <- dt.2 %>% dplyr::mutate(results_group = paste(Results,group, sep = "_"))
  # To refine our search, we will narrow our range. The oncopanel form has "cytoband" before our relevant section and
  ## "interpretation" after, so we can define our range that way as well
  cytoband_interpretation_regex <- "Cytoband|INTERPRETATION"
  ### We will use the same strategy as above with the keywords look up table and the cumulative sum grouping
  dt.4a <- dt.3 %>% dplyr::mutate(divider_1 = stringr::str_detect(string = dt.3$Results, pattern = stringr::regex(cytoband_interpretation_regex)),
                                 group_1 = cumsum(divider_1)) %>%
    dplyr::filter(group_1 == 1)
  ## drop the "," at the end of the gene names
  dt.4 <- dt.4a
  dt.4$Results <- stringr::str_replace(string = dt.4$Results, pattern = ",$", replacement = "")


  # Now we need a way to identify each of these groups (i.e. "group_1") with the appropriate gene modification
  ## e.g. "Gain" or "Loss"
  ### To do this, we will create an object that is a vector of only those "results_group"s that contain "Gain"; thus
  #### all the numbers in this group that follow the "_" will be the "group_1" numbers associated with "Gain
  gain_idx <- dt.4 %>% dplyr::select(results_group) %>% dplyr::filter(stringr::str_detect(string = dt.4$results_group, pattern = "Gain"))
  # Let's remove the "Gain_" part so we can identify the numbers that will be the look up table numbers for "Gain"
  gain_idx$results_group <- str_remove_all(string = gain_idx$results_group, pattern = "Gain_") %>% as.numeric() # and let's turn them into numerics
  # rename "results_group" "group" so that it can link up with the appropriate column in "divided.1"
  gain_idx <- gain_idx %>% rename(group = results_group)
  # Add a "cnv" column and populate the cells with "gain"
  ## "cnv" will be used for all three types of cnvs
  gain_idx <- gain_idx %>% dplyr::mutate(cnv = "gain")

  # We are now going to execute an analogous set of steps with "Amplification"
  amplification_idx <- dt.4 %>% dplyr::select(results_group)
  amplification_idx <- amplification_idx %>% dplyr::filter(str_detect(string = results_group, pattern = "Amplification"))
  amplification_idx$results_group <- stringr::str_remove_all(string = amplification_idx$results_group, pattern = "Amplification_") %>% as.numeric()

  amplification_idx <- amplification_idx %>% dplyr::rename(group = results_group)

  amplification_idx <- amplification_idx %>%  dplyr::mutate(cnv = "amplification")

  # We are now going to execute an analogous set of steps with "Loss"
  loss_idx <- dt.4 %>% dplyr::select(results_group) %>% dplyr::filter(str_detect(string = results_group, pattern = "Loss"))

  loss_idx$results_group <- stringr::str_remove_all(string = loss_idx$results_group, pattern = "Loss_") %>% as.numeric()

  loss_idx <- loss_idx %>% dplyr::rename(group = results_group)

  loss_idx <- loss_idx %>% dplyr::mutate(cnv = "loss")


  # combine these three tables with "dt.4" to create a "cnv" table that we can use to ID the genes associated
  ## with "Gain" "Loss" "Amplification"
  cnv <- left_join(dt.4, gain_idx)

  cnv <- left_join(cnv,amplification_idx, by = c("group"))

  cnv <- left_join(cnv, loss_idx, by = "group")

  cnv <- cnv %>% unite("cnv_label",
                       cnv,
                       cnv.x,
                       cnv.y,
                       na.rm = TRUE,
                       remove = TRUE
  )


  ##########################################################################################################################
  # Let's create the keyword list of genes that we are going to use to pull out only those rows in "cnv" that have Genes
  ## as tokens

  # Now we will pull out the rows that have gene names in the "Results" column
  cnv.1 <- cnv %>% dplyr::filter(stringr::str_detect(string = cnv$Results, pattern = stringr::regex(genes_regex)))

  # Let's select only "Results" and "cnv_label" and rename "Results" into "gene"
  cnv.2 <- cnv.1 %>% dplyr::select(Results, cnv_label) %>% dplyr::rename(gene = Results)
  ## Thus, cnv.2 is the dataframe that contains those genes with gains, losses or amplifications
  ### But because in our REDCap tool, amplifications are separate from gains or losses, we will remove them from our data
  ####frame
  ##### Let's create a separate table for amplifications since these have their own variables
  amplifications <- cnv.2 %>% dplyr::filter(cnv_label == "amplification")
  ###### Use anti-join to remove "amplifications" from cnv.2
  cnv.3a <- dplyr::anti_join(cnv.2, amplifications)
  cnv.3 <- if(is.na(cnv.3a[1,1])) base::data.frame(gene = NA, cnv_label = NA) else cnv.3a
  ###### Now let's modify "cnv.3" to set it up with the appropriate variable names for the import_tool
  cnv.3 <- cnv.3 %>% dplyr::mutate(cnv_gene = "cnv_gene")
  cnv.3 <- cnv.3 %>% dplyr::mutate(cnv_gain_or_loss = "cnv_gain_or_loss")
  cnv.3 <- cnv.3 %>% dplyr::mutate(number = 1:nrow(cnv.3))
  cnv.4a <- cnv.3 %>% dplyr::mutate(cnv_gene = paste(cnv_gene, number, sep = "_"))
  cnv.4 <- cnv.4a %>% dplyr::mutate(cnv_gain_or_loss = paste(cnv_gain_or_loss, number, sep = "_"))
  ####### rename variables to set up for combining with import_tool
  ######## cnv.gene is a df that contains the name of the gene side-by-side with the variable name
  cnv.gene <- cnv.4 %>% dplyr::select(gene, cnv_gene) %>% dplyr::rename(variables = cnv_gene) %>% dplyr::rename(results = gene) %>% dplyr::select(variables, results)
  ######## cnv.gene.gain.or.loss is a df that contains the type of cnv side-by-side with the variable name
  cnv.gene.gain.or.loss <- cnv.4 %>% dplyr::select(cnv_label, cnv_gain_or_loss) %>% dplyr::rename(variables = cnv_gain_or_loss) %>% dplyr::rename(results = cnv_label) %>% dplyr::select(variables, results)



  ##########################################################################################################################
  # Let's continue with our processing of amplifications
  ## This is a bit more complicated because for some patient there will be no amplifications, so we need to come up with a
  ### way to prevent the code from producing an error if the amplifications are absent from the sample
  #### we will solve this by adding a NA row as a default and then removing that row later
  amplifications.1 <- amplifications %>% dplyr::select(gene)
  ##*******************************************************KEY STEP*******************************************************##
  #### we will solve this problem of a potential empty df by adding a NA row as a default, adding a random column here
  ##### "blank" to then add an actual value so we can run the nrow(). If we don't have a value in at least one of the rows
  ###### then the code will produce an error. We will remove this later, so it doesn't affect the final product
  amplifications.2 <- if(is.na(amplifications.1[1,1])) base::data.frame(gene = NA, amplification_gene = "amplifications_gene") else amplifications.1
  amplifications.3 <- amplifications.2 %>% dplyr::mutate(number = 1:nrow(amplifications.2))
  ##### Similar to what we did above, we will build a variables column to associate these values with the appropriate
  ##### variable that we can import later
  amplifications.3 <- amplifications.3 %>% dplyr::mutate(amplification_gene = "amplifications_gene")

  amplifications.4 <- amplifications.3 %>% dplyr::mutate(variables = paste(amplification_gene, number, sep = "_")) %>%
    dplyr::rename(results = gene) %>% dplyr::select(variables, results)

  cnv.final <- base::rbind(cnv.gene, cnv.gene.gain.or.loss, amplifications.4)

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
  # change to character vector
  cnv.final$results <- base::as.character(cnv.final$results)

  return(cnv.final)
}

