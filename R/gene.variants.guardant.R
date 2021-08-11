#' Abstract Genes with nucleotide variants from Guardant Health genomic reports
#' @description
#' `gene.variants.guardant()` provides natural language processing tools to abstract genes with nucleotide variants
#' @param data  The data frame of a genomic report. Ideally this is information copied to the Clipboard from the EHR report by "selecting all" and then "copy". That is the default. If you don't use the clipboard function, you can use a data frame of the text file. The single column data frame should be labeled "Results". Required.
#'
#' @return a data frame with two columns: variables (the redcap variable names) and results (the data to be imported in to redcap)
#' @export
#'
gene.variants.guardant <- function(data = dplyr::tibble(Results = clipboard())){

  ##########################################################################################################################
  # load data and platform
  ##########################################################################################################################
  dt <- data

  ##########################################################################################################################
  # Create regexs that we will use to select our data
  ##########################################################################################################################
  # Keep only those rows with the following
  ## call the genes_regex function to create a regex of all the gene names
  genes_boundary_df <- genetex::genes_boundary_regex() # this is a df of unique genes names concatenated with "|"
  genes_boundary_regex <- genes_boundary_df$genes # this is a character string of the genes for our regex
  nuc_regex <- "[ACTG]>[ACTG]|del[ACTG]|promoter" # This regex should identify those rows that have nucleotide changes
  aa_regex <- "(\\b([A-Z][0-9]{1,}(([A-Z])|(_[A-Z][1-9]{1,}del)|(fs\\*[1-9]{1,})|(\\*)|(fs)|(del)))|(p\\.[A-Z]))|([0-9]ins[A-Z])"
  # This regex will detect all the possible AA alterations with single initials
  cfdna_regex <- "\\b[0-9]{1,2}\\.[0-9]{1,2}%"
  genes_nuc_aa_cfdna_regex <- paste(genes_boundary_regex, aa_regex, nuc_regex, cfdna_regex, sep = "|")
  not_percent <- "(\\b^.+%)"

  ##########################################################################################################################
  # Cleaning The Data
  ##########################################################################################################################
  dt.1a <- dt
  dt.1a$Results <- stringr::str_trim(string = dt.1a$Results) # Trim white space
  # We need to eliminate much of this report in order to mine only those gene alterations relevant to each subject
  ## So we will use the word "Clinical Trial Page" as a divider and take only those rows that are above it
  keyword <- "Clinical Trial Page"
  dt.1b <- dt.1a %>%
    mutate(divider = stringr::str_detect(string = dt.1a$Results,
                                         pattern = stringr::regex(keyword, ignore_case = TRUE)),
           group = cumsum(divider)) %>%
    filter(group == 0)

  # tokenize the vector with space
  dt.1 <- splitstackshape::cSplit(
    indt = dt.1b,
    splitCols = "Results",
    sep = " ",
    direction = "long"
  )


  # Filter only those rows with the above genes_nuc_aa_cfdna_regex
  dt.2 <- dt.1 %>%
    dplyr::filter(stringr::str_detect(string = dt.1$Results,
                                      pattern = stringr::regex(genes_nuc_aa_cfdna_regex, ignore_case = TRUE)))

  # Create a grouping system based on gene names
  dt.3 <- dt.2 %>%
    dplyr::mutate(keywords = stringr::str_detect(string = dt.2$Results,
                                                 pattern = stringr::regex(genes_boundary_regex)),
                  group = base::cumsum(keywords))
  # find those groups that are duplicated, but do not have more than 3
  dt.4 <- dt.3 %>% dplyr::count(group) %>% dplyr::filter(n >= 2 & n<= 3)
  # filter only those rows from dt.4 that are in the new dt.group
  dt.5 <- dt.3 %>% dplyr::filter(dt.3$group %in% dt.4$group)
  # drop ; or . or , at the end of the string
  dt.5a <- dt.5
  dt.5a$Results <- stringr::str_replace(string = dt.5a$Results, pattern = ";$", replacement = "")
  dt.5a$Results <- stringr::str_replace(string = dt.5a$Results, pattern = "\\.$", replacement = "")
  dt.5a$Results <- stringr::str_replace(string = dt.5a$Results, pattern = ",$", replacement = "")
  dt.5a$Results <- stringr::str_replace(string = dt.5a$Results, pattern = "\\(", replacement = "")
  dt.5a$Results <- stringr::str_replace(string = dt.5a$Results, pattern = "\\)", replacement = "")
  dt.5a$Results <- stringr::str_replace_all(string = dt.5a$Results, pattern = "†", replacement = "")
  dt.5a$Results <- stringr::str_replace_all(string = dt.5a$Results, pattern = "p.", replacement = "")

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
                                                                            pattern = stringr::regex(nuc_regex, ignore_case = TRUE)),
                                                 yes = "variant_nucleotide",
                                                 no = ifelse(test = stringr::str_detect(string = dt.8$merge,
                                                                                        pattern = stringr::regex(aa_regex)),
                                                             yes = "variant_protein",
                                                             no = ifelse(test = stringr::str_detect(string = dt.8$merge,
                                                                                                    pattern = stringr::regex(cfdna_regex)),
                                                                         yes = "variant_gene_perc_cfdna",
                                                                         no = "")))))
  dt.9a <- dt.9 %>% dplyr::filter(stringr::str_detect(string = dt.9$merge, pattern = "%")) %>% select(group, merge)
  dt.9b <- dt.9a
  dt.9b$merge <- stringr::str_replace(string = dt.9b$merge, pattern = "%", replacement = "")
  dt.9b$merge <- as.numeric(dt.9b$merge)
  dt.9b <- dt.9b %>% dplyr::arrange(desc(merge)) %>% dplyr::mutate(order = 1:nrow(dt.9b)) %>% select(group, order)


  # Left join dt.9 with dt.9b
  dt.10 <- dplyr::left_join(dt.9, dt.9b, by = "group") %>% dplyr::arrange(order)

  # paste group and var
  dt.11 <- dt.10 %>% mutate(variables = paste(var, order, sep = "_"))
  dt.12 <- dt.11 %>% rename(results = merge) %>%  select(variables, results)

  # Account for "Promoter SNV"
  ## The guardant report reports Promoter SNVs. Let's account for this and replace "promoter" with "Promoter SNV".

  dt.12.promoter_snv <- dt.12
  dt.12.promoter_snv$results <- stringr::str_replace(string = dt.12.promoter_snv$results,
                                                     pattern = "Promoter",
                                                     replacement = "Promoter SNV")
  ### Only replace if the exact string "Promoter SNV" is found in the initial report (since there may be other
  #### promoter variations and we don't want to replace those, only Promoter SNV)
  dt.13 <- if(any(stringr::str_detect(string = dt$Results, pattern = "Promoter SNV"))) dt.12.promoter_snv else dt.12


  # Add variants_number
  variants.number.1 <- stringr::str_detect(string = dt.13$variables, pattern = "variant_gene_[0-9]")
  variants.number.2 <- sum(variants.number.1)
  variants.number <- if(is.na(dt.13[1,2])) 0 else variants.number.2
  variants.number.df <- data.frame(variables = "variant_number",
                              results = variants.number)

  # combine dfs
  gene_variants <- rbind(dt.13, variants.number.df)
  gene_variants$results <- base::as.character(gene_variants$results)

  return(gene_variants)
}

