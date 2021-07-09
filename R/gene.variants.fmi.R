#' Apply NLP to genomic reports to extract gene variants (CNVs) from foundation medicine reports
#' @description
#' `gene.variants.fmi()` provides natural language processing tools to abstract gene variants data from foundation medicine genomic reports
#' @param data The data frame of the genomic report of interest. This can be copied to the Clipboard from the EHR report by selecting the relevant text and then "copy". That is the default. If you don't use the clipboard function, you can use a data frame of the text file. The single column data frame should be labeled "Results". Required.
#' @return a data frame with two columns: variables (the redcap variable names) and results (the data to be imported in to redcap)
#' @export
#'
gene.variants.fmi <- function(data = dplyr::tibble(Results = readr::clipboard())){
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
  ENST <- "ENST"
  nuc_c. <- "c\\."
  nuc_regex <- "[ACTG]>[ACTG]|del[ACTG]" # This regex should identify those rows that have nucleotide changes
  aa_p._regex <-"p\\.[A-Z]" # some platforms designate amino acid changes using a "p"".""AA" structure
  #aa_regex1 <- "\\b[A-Z][0-9]{1,}([A-Z]|fs|\\*)"
  aa_regex <- "(\\b([A-Z][0-9]{1,}(([A-Z])|(_[A-Z][1-9]{1,}del)|(fs\\*[1-9]{1,})|(\\*)|(fs)|(del)))|(p\\.[A-Z]))|([0-9]ins[A-Z])"
  # This regex will detect all the possible AA alterations with single initials
  cfdna_regex <- "\\b[0-9]{1,2}\\.[0-9]{1,2}%"
  genes_boundary_nuc_aa_regex <- paste(genes_boundary_regex, nuc_regex, aa_regex, sep = "|") # concatenate several regexs
  genes_boundary_nuc_aa_aa_p._regex <- paste(genes_boundary_regex, nuc_regex, aa_regex, sep = "|") #concat more regexs
  genes_boundary_nuc_aa_regex_cfdna <- paste(genes_boundary_regex, aa_regex, nuc_regex, cfdna_regex, sep = "|")
  # Assemble a regex_discard to perform an initial screen
  regex_discard <- "H3K27|Sf3b1|H3B-6545|\\bfusion\\b|rearrangement|\\bloss\\b|\\bgain\\b|\\bthe\\b|\\bincluding\\b|\\bcells\\b|\\bcancer\\b|amplification"

  ##########################################################################################################################
  # Trim Data
  ##########################################################################################################################
  dt.1a <- dt
  dt.1a$Results <- stringr::str_trim(string = dt.1a$Results) # Trim white space
  ##########################################################################################################################
  # Isolate those findings on the first page under "Genomic Findings"
  ##########################################################################################################################

  genomic_findings_regex <- "Genomic Findings"
  page_1_of_regex <- "PAGE 1 of"

  gf.1b <- dt.1a %>%
    mutate(divider = stringr::str_detect(string = dt.1a$Results,
                                         pattern = stringr::regex(page_1_of_regex, ignore_case = TRUE)),
           group = cumsum(divider)) %>%
    filter(group == 0)

  # Perform a clean to remove rows that have problematic strings
  gf.1c <- gf.1b %>% dplyr::filter(!stringr::str_detect(string = gf.1b$Results,
                                                        pattern = stringr::regex(regex_discard, ignore_case = FALSE)))

  # tokenize the vector with space
  gf.1 <- splitstackshape::cSplit(
    indt = gf.1c,
    splitCols = "Results",
    sep = " ",
    direction = "long"
  )

  # Filter only those rows with the above genes_nuc_aa_regex
  gf.2 <- gf.1 %>% dplyr::filter(stringr::str_detect(string = gf.1$Results,
                                                     pattern = stringr::regex(genes_boundary_nuc_aa_regex)))


  gf.3 <- gf.2 %>% dplyr::filter(!stringr::str_detect(string = gf.2$Results,
                                                      pattern = regex(regex_discard, ignore_case = FALSE)))

  ##########################################################################################################################
  # Isolate those findings on the page with variants of unknown significance
  ##########################################################################################################################

  vus_regex <- "One or more variants of unknown significance"
  page_of_regex <- "PAGE [0-9]{1,} of"
  vus_page_regex <- paste(vus_regex, page_of_regex, sep = "|")


  vus.1 <- dt.1a %>%
    mutate(divider = stringr::str_detect(string = dt.1a$Results,
                                         pattern = stringr::regex(vus_regex, ignore_case = FALSE)),
           group = cumsum(divider)) %>%
    filter(group == 1)

  vus.2 <- vus.1 %>%
    mutate(divider = stringr::str_detect(string = vus.1$Results,
                                         pattern = stringr::regex(page_of_regex, ignore_case = FALSE)),
           group = cumsum(divider)) %>%
    filter(group == 0)

  # Filter only those rows with the above genes_nuc_aa_regex
  vus.3 <- vus.2 %>% dplyr::filter(stringr::str_detect(string = vus.2$Results,
                                                       pattern = stringr::regex(genes_boundary_nuc_aa_regex)))

  # tokenize the vector with space
  vus.4a <- splitstackshape::cSplit(
    indt = vus.3,
    splitCols = "Results",
    sep = " ",
    direction = "long"
  )

  # Filter only those rows with the above genes_nuc_aa_regex
  vus.4 <- vus.4a %>% dplyr::filter(stringr::str_detect(string = vus.4a$Results,
                                                     pattern = stringr::regex(genes_boundary_nuc_aa_regex)))

  ##########################################################################################################################
  # Combine vus and genomic finding data
  ##########################################################################################################################

  ifelse(test = (is.na(gf.3[1,1]) & is.na(vus.4[1,1])),
         yes = dt.3 <- base::data.frame(Results = NA,
                                         variables = NA),
         no = ifelse(test = is.na(gf.3[1,1]),
                     yes = dt.3 <- vus.4,
                     no = ifelse(test = is.na(vus.4[1,1]),
                                 yes = dt.3 <- gf.3,
                                 no = dt.3 <- base::rbind(gf.3, vus.4))))


  dt.3 <- dt.3

  # Create a grouping system based on gene names
  dt.4 <- dt.3 %>% dplyr::mutate(keywords = stringr::str_detect(string = dt.3$Results,
                                                                pattern = stringr::regex(genes_boundary_regex)),
                                 group = base::cumsum(keywords))
  # find those groups that are duplicated, but do not have more than 4
  dt.group <- dt.4 %>% dplyr::count(group) %>% dplyr::filter(n >= 2 & n<= 4)
  # filter only those rows from dt.4 that are in the new dt.group
  dt.5 <- dt.4 %>% dplyr::filter(dt.4$group %in% dt.group$group)
  # drop ; or . or , at the end of the string
  dt.5$Results <- stringr::str_replace(string = dt.5$Results, pattern = ";$", replacement = "")
  dt.5$Results <- stringr::str_replace(string = dt.5$Results, pattern = "\\.$", replacement = "")
  dt.5$Results <- stringr::str_replace(string = dt.5$Results, pattern = ",$", replacement = "")
  dt.5$Results <- stringr::str_replace(string = dt.5$Results, pattern = "\\(", replacement = "")
  dt.5$Results <- stringr::str_replace(string = dt.5$Results, pattern = "\\)", replacement = "")
  dt.5$Results <- stringr::str_replace_all(string = dt.5$Results, pattern = "†", replacement = "")

  # remove the prefix before AA or NT
  dt.5a <- dt.5
  dt.5a$Results <- stringr::str_replace(string = dt.5a$Results, pattern = "p\\.=$", replacement = "")
  # transpose wide via paste0 so we can eliminate duplicates
  dt.6 <- dt.5a %>% dplyr::group_by(group) %>% dplyr::mutate(merge = base::paste0(Results, collapse = " ")) %>% ungroup()
  # slice to eliminate duplicates
  dt.7 <- dt.6 %>% dplyr::group_by(merge) %>% dplyr::slice_head()
  # tokenize the vector with space
  dt.8 <- splitstackshape::cSplit(
    indt = dt.7,
    splitCols = "merge",
    sep = " ",
    direction = "long"
  )
  # assign the correct variable "variant_gene", "variant_nucleotide", or "variant_protein"
  dt.9 <- dt.8 %>%
    dplyr::mutate(var = base::ifelse(test = stringr::str_detect(string = dt.8$merge,
                                                                pattern = stringr::regex(genes_boundary_regex)),
                                     yes = "variant_gene",
                                     no = ifelse(test = stringr::str_detect(string = dt.8$merge,
                                                                            pattern = stringr::regex(nuc_regex)),
                                                 yes = "variant_nucleotide",
                                                 no = ifelse(test = stringr::str_detect(string = dt.8$merge,
                                                                                        pattern = stringr::regex(aa_regex)),
                                                             yes = "variant_protein",
                                                             no = ifelse(test = stringr::str_detect(string = dt.8$merge,
                                                                                                    pattern = stringr::regex(cfdna_regex)),
                                                                         yes = "variant_gene_perc_cfdna",
                                                                         no = "")))))

  # Create a grouping system based on gene names
  dt.10 <- dt.9 %>%
    dplyr::mutate(keywords = stringr::str_detect(string = dt.9$merge,
                                                 pattern = stringr::regex(genes_boundary_regex)),
                  group = base::cumsum(keywords))

  # paste group and var
  dt.11 <- dt.10 %>% mutate(variables = paste(var, group, sep = "_"))
  dt.12a <- dt.11 %>% select(merge, variables, group) %>% rename(results = merge)
  # find those groups that are duplicated, but do not have more than 4
  dt.group.1 <- dt.12a %>% dplyr::count(group) %>% dplyr::filter(n >= 2 & n<= 4)
  # filter only those rows from dt.4 that are in the new dt.group
  dt.12 <- dt.12a %>% dplyr::filter(dt.12a$group %in% dt.group.1$group)

  ##########################################################################################################################
  # Isolate genes with more than one subclone with different genetic alterations
  ##########################################################################################################################
  ## Tumors have various subclones that can be identified with deep sequencing
  ### Thus, tumors may have a multiple mutations in the same gene
  #### We need to find genes that have subclones and devise a system to appropriate pull them out, duplicate the gene name
  ##### and then order them so that they can fit in with our general structured system
  ###### This is necessary, b/c certain reports will have the name of the gene e.g. "TP53" followed by two AAs
  ####### e.g. R121I Q192N
  dt.13 <- dt.12 %>% filter(duplicated(x = dt.12$variables)) # This df has all the rows that are duplicated, which will let
  ## us create a look up table of groups that need to be pulled out from the main df so we can appropriately amplifiy the
  ### gene names to organize them so they can be re-emerged with the df that will be used as the final product of gene
  #### variants
  dt.14 <- dplyr::anti_join(dt.12, dt.13, by = "group") # This df is all the genes with single nuc/aa sequences that can
  ## be used to build our final df; if there are no "subclones" then this df will be the final working df
  dt.15 <- dplyr::semi_join(dt.12, dt.13, by = "group") # This df has all the genes with subclones
  ## Now we must come up with a way to duplicate the gene names just the appropriate amount to correlate to the number of
  ### subclones that exist
  #### first we will count how many times each variable is used in dt.15
  dt.16 <- dt.15 %>% dplyr::count(variables)
  # Now let's filter for only those variables without "variant_gene" (i.e. just amino acids and nucleotides)
  dt.17 <- dt.16 %>% dplyr::filter(!stringr::str_detect(dt.16$variables, "variant_gene"))
  # Now let's create a look up table by pulling out the number associated with the variable and call this "group"
  ## i.e. all the variant_protein_5 will be group "5". That way we can use that number to link it back to variant_gene_5
  dt.18 <- dt.17 %>% dplyr::mutate(group = stringr::str_extract(variables, regex("[0-9]{1,}")))
  # Create a group called "ntimes" that is the sum of the number of times that group was in dt.18
  ## This is important, because a genomic report may have reported two amino acids and only one nucleotide and we need to
  ### account for this
  dt.19 <- dt.18 %>% dplyr::group_by(group) %>% dplyr::mutate(ntimes = sum(n))
  dt.20 <- dt.19 %>% dplyr::select(group, ntimes)
  dt.20$group <- base::as.numeric(dt.20$group)
  dt.20a <- dt.20 %>% dplyr::group_by(group) %>% dplyr::slice_head() # This df tells us how many times each group
  ## needs to be repeated
  # isolate just the gene names from our df with the genes with duplicates (ie. dt.15)
  dt.21 <-  dt.15 %>% dplyr::filter(stringr::str_detect(string = dt.15$variables, pattern = "variant_gene"))
  # left join dt.21 with dt.20a, so we can see how many times each gene name has to be repeated
  dt.22 <- dplyr::left_join(dt.21, dt.20a, by = "group")
  dt.22$results <- base::as.character(dt.22$results)
  # Now we need a script to repeat the gene names based on the variable "ntimes"
  ## However, first we need dt.22 may not exist for certain genomic reports,
  ### thus, "dt.22a[rep(1:nrow(dt.22a), times = dt.22a$ntimes),]" will throw an error; therefore, we need to generate a hold df,
  #### that we will call dt.22a, and ntimes will equal 0
  dt.22a <- if(base::is.na(dt.22[1,1])) data.frame(results = NA, variables = NA, group = NA, ntimes = 0) else dt.22
  # now the below script wll repeat the variant_gene based on the value of "ntimes"
  dt.23 <- dt.22a[rep(1:nrow(dt.22a), times = dt.22a$ntimes),]
  dt.24 <- dt.23 %>% dplyr::select(results, variables, group)
  dt.25 <- dt.15 %>% dplyr::filter(!stringr::str_detect(string = dt.15$variables,
                                                        pattern = "variant_gene"))
  dt.26 <- rbind(dt.24,dt.25)
  # because dt.26 will be blank for certain forms, we need to create a fake df with NAs for values
  dt.26a <- base::data.frame(results = NA, variables = NA, group = NA)
  dt.26b <- if(base::is.na(dt.26[1,1])) dt.26a else dt.26
  # arrange the rows based on groups and then create a new variable "num" which will allow us to renumber the genes
  dt.27 <- dt.26b %>% dplyr::arrange(group) %>% dplyr::mutate(num = 1:base::nrow(dt.26b))
  # Now we need a method that allows us to pull the first and last rows within a group so that we can have a
  ## gene name, followed by either a nucleotide or a amino acid
  dt.28 <- dt.27 %>% dplyr::group_by(group) %>% dplyr::filter(dplyr::row_number() %in% base::c(1L, dplyr::n()))
  # next remove those rows from the larger df with anti-join
  dt.29 <- dplyr::anti_join(dt.27,dt.28)
  # essentially repeat the above process with filtered df (df.29)
  dt.30 <- dt.29 %>% dplyr::group_by(group) %>% dplyr::filter(dplyr::row_number() %in% base::c(1L, dplyr::n()))
  dt.31 <- dplyr::anti_join(dt.29,dt.30)
  # repeat the process one more time in case there are multiple subclones (i.e 3)
  dt.32 <- dt.31 %>% dplyr::group_by(group) %>% dplyr::filter(dplyr::row_number() %in% base::c(1L, dplyr::n()))
  # rbind dfs 28, 30, 32
  dt.33 <- base::rbind(dt.28, dt.30, dt.32) %>% select(results, variables, group)
  # add conditional statement to choose dt.14 for those data frames that don't have duplicates that have been
  ## dealt with. Because dt.33 will just be NA if there are no subclones, then then if statement will be TRUE and
  ### we will use dt.14; but if dt.33 has data in row 1, column 1, we want to keep that, and we will then rbind it
  #### to dt.14
  dt.34 <- if(base::is.na(dt.33[1,1])) dt.14  else base::rbind(dt.14, dt.33)
  dt.35 <- dt.34 %>% select(results)
  # We now need to renumber this df, so let's first assign the correct variable "variant_gene", "variant_nucleotide", or "variant_protein"
  dt.36 <- dt.35 %>%
    dplyr::mutate(var = base::ifelse(test = stringr::str_detect(string = dt.35$results,
                                                                pattern = stringr::regex(genes_boundary_regex)),
                                     yes = "variant_gene",
                                     no = ifelse(test = stringr::str_detect(string = dt.35$results,
                                                                            pattern = stringr::regex(nuc_regex)),
                                                 yes = "variant_nucleotide",
                                                 no = ifelse(test = stringr::str_detect(string = dt.35$results,
                                                                                        pattern = stringr::regex(aa_regex)),
                                                             yes = "variant_protein",
                                                             no = ifelse(test = stringr::str_detect(string = dt.8$merge,
                                                                                                    pattern = stringr::regex(cfdna_regex)),
                                                                         yes = "variant_gene_perc_cfdna",
                                                                         no = "")))))
  # Create a grouping system based on gene names (similar to what we've done before)
  dt.37 <- dt.36 %>% dplyr::mutate(keywords = stringr::str_detect(string = dt.36$results,
                                                                  pattern = stringr::regex(genes_boundary_regex)),
                                   group = base::cumsum(keywords))
  # paste group and var (similar to what we've done before)
  dt.38 <- dt.37 %>% mutate(variables = paste(var, group, sep = "_"))
  dt.39 <- dt.38 %>% select(results, variables)
  # Let's remove ENST/ENSP if they exist (they may not, so this will have no effect on those genomic reports)
  dt.39$results <- stringr::str_replace(string = dt.39$results,
                                        pattern = "(ENST|ENSP)[0-9]{1,}\\.[0-9]:[pc]\\.",
                                        replacement = "")
  # Some reports may still have duplicates after this process
  ## Thus, we will essentially repeat the steps we've used before
  dt.40 <-  dt.39 %>% dplyr::mutate(group = stringr::str_extract(string = dt.39$variables,
                                                                 pattern = "[0-9]{1,}"))
  # transpose wide via paste0 so we can eliminate duplicates
  dt.41 <- dt.40 %>% dplyr::group_by(group) %>% dplyr::mutate(merge = base::paste0(results, collapse = " ")) %>% ungroup()
  # slice to eliminate duplicates
  dt.42 <- dt.41 %>% dplyr::group_by(merge) %>% dplyr::slice_head()
  # There are some rows that have the same gene and AA, but one may have the nucleotides and the other does not, so let's
  ## try and eliminate the duplicates and preserve those rows that have all three
  dt.43 <- dt.42 %>% tidyr::separate(col = merge, into = c("a","b","c"), sep = " ")
  # now paste "a" and "b" together so we can slice
  dt.44 <- dt.43
  dt.44$merge <-  base::paste(dt.43$a, dt.43$b, sep = " ")
  # now group by merge and then sort by "c" and slice_head
  dt.45 <- dt.44 %>% dplyr::group_by(merge) %>% dplyr::arrange(dt.44$c) %>% dplyr::slice_head() %>% dplyr::ungroup()
  # now we should have non-duplicates, so let's merge "a", "b", "c" and tokenize
  dt.46 <- dt.45 %>% select(a, b, c)
  dt.46$merge <- paste(dt.46$a, dt.46$b, dt.46$c, sep = " ")
  dt.47 <- dt.46 %>% select(merge)
  dt.47$merge <- stringr::str_replace(string = dt.47$merge, pattern = "\\bNA\\b", replacement = "")
  # tokenize by space to transpose to long format
  dt.48 <- splitstackshape::cSplit(
    indt = dt.47,
    splitCols = "merge",
    sep = " ",
    direction = "long"
  )
  dt.49 <- dt.48 %>% rename(results = merge)
  dt.50 <- dt.49 %>%
    dplyr::mutate(var = base::ifelse(test = stringr::str_detect(string = dt.49$results,
                                                                pattern = stringr::regex(genes_boundary_regex)),
                                     yes = "variant_gene",
                                     no = ifelse(test = stringr::str_detect(string = dt.49$results,
                                                                            pattern = stringr::regex(nuc_regex)),
                                                 yes = "variant_nucleotide",
                                                 no = ifelse(test = stringr::str_detect(string = dt.49$results,
                                                                                        pattern = stringr::regex(aa_regex)),
                                                             yes = "variant_protein",
                                                             no = ifelse(test = stringr::str_detect(string = dt.8$merge,
                                                                                                    pattern = stringr::regex(cfdna_regex)),
                                                                         yes = "variant_gene_perc_cfdna",
                                                                         no = "")))))

  # Create a grouping system based on gene names
  dt.51 <- dt.50 %>%
    dplyr::mutate(keywords = stringr::str_detect(string = dt.50$results,
                                                 pattern = stringr::regex(genes_boundary_regex)),
                  group = base::cumsum(keywords))
  # paste group and var
  dt.52 <- dt.51 %>% mutate(variables = paste(var, group, sep = "_"))
  dt.53 <- dt.52 %>% select(variables, results)

  # Add variants_number
  variants.number.1 <- stringr::str_detect(string = dt.53$variables, pattern = "variant_gene")
  variants.number.2 <- sum(variants.number.1)
  variants.number <- if(is.na(dt.53[1,2])) 0 else variants.number.2
  variants.number.df <- data.frame(variables = "variant_number",
                                   results = variants.number)

  # combine dfs
  gene_variants <- rbind(dt.53, variants.number.df)
  gene_variants$results <- base::as.character(gene_variants$results)

  return(gene_variants)

}
