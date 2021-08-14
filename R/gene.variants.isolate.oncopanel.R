#' Abstract Genes with Nucleotide variants from Oncopanel genomic reports
#' @description
#' `gene.variants.isolate.oncopanel()` provides natural language processing tools to abstract genes with nucleotide variants from Oncopanel reports
#' @param data  The data frame of a genomic report, here an Oncopanel report. Ideally this is information copied to the Clipboard from the EHR report by "selecting all" and then "copy". That is the default. If you don't use the clipboard function, you can use a data frame of the text file. The single column data frame should be labeled "Results". Required.
#'
#' @return a data frame with a column that has the gene names and variants in one cell (i.e. untidy data)
#' @examples
#' # Test with embedded data frame "oncopanel_sample_report"
#' oncopanel_sample_report %>%
#'   gene.variants.isolate.oncopanel()
#' @export
#'
gene.variants.isolate.oncopanel <- function(data = dplyr::tibble(Results = clipboard())){

  ##########################################################################################################################
  # Load Data and Trim
  ##########################################################################################################################
  dt <- data
  dt$Results <- stringr::str_trim(string = dt$Results) # Trim white space

  ##########################################################################################################################
  # Create regexs that we will use to select our data
  ##########################################################################################################################
  ## call the genes_regex function to create a regex of all the gene names
  genes_boundary_df <- genetex::genes_boundary_regex() # this is a df of unique genes names concatenated with "|"
  genes_boundary_regex <- genes_boundary_df$genes # this is a character string of the genes for our regex
  nuc_regex <- "[ACTG]>[ACTG]|del[ACTG]" # This regex should identify those rows that have nucleotide changes
  aa_regex <- "(\\b([A-Z][0-9]{1,}(([A-Z])|(_[A-Z][1-9]{1,}del)|(fs\\*[1-9]{1,})|(\\*)|(fs)|(del)))|(p\\.[A-Z]))|([0-9]ins[A-Z])"
  genes_boundary_nuc_aa_regex <- paste(genes_boundary_regex, nuc_regex, aa_regex, sep = "|") # concatenate several regexs

  ##########################################################################################################################
  # Isolate Gene Variants
  ##########################################################################################################################
  # tokenize the vector with space as the delimiter
  dt.1a <- splitstackshape::cSplit(
    indt = dt,
    splitCols = "Results",
    sep = " ",
    direction = "long"
  )

  # We need to eliminate much of this report in order to mine only those gene alterations relevant to each subject
  ## So we will use the word "cytoband/size" as a divider and take only those rows that are above it
  keyword <- "cytoband/size"
  dt.1b <- dt.1a %>%
    mutate(divider = stringr::str_detect(string = dt.1a$Results,
                                         pattern = stringr::regex(keyword, ignore_case = TRUE)),
           group = cumsum(divider)) %>%
    filter(group == 0)


  # Eliminate erroneous strings, such as reads###
  dt.1c <- dt.1b
  dt.1c$Results <- stringr::str_replace_all(
    string = dt.1c$Results,
    pattern = "reads###",
    replacement = " "
  )

  dt.1c$Results <- stringr::str_replace_all(
    string = dt.1c$Results,
    pattern = "reads##",
    replacement = " "
  )

  dt.1c$Results <- stringr::str_trim(dt.1c$Results)


  # Filter only those rows with the above genes_nuc_aa_regex
  dt.2 <- dt.1c %>%
    dplyr::filter(stringr::str_detect(string = dt.1c$Results,
                                      pattern = stringr::regex(genes_boundary_nuc_aa_regex)))

  # Create a grouping system based on gene names
  dt.3 <- dt.2 %>%
    dplyr::mutate(keywords = stringr::str_detect(string = dt.2$Results,
                                                 pattern = stringr::regex(genes_boundary_regex)),
                  group = base::cumsum(keywords))
  # find those groups that are duplicated, but do not have more than 3
  dt.4 <- dt.3 %>% dplyr::count(group) %>% dplyr::filter(n >= 2 & n<= 3)
  # filter only those rows from dt.3 that are in the new dt.4
  dt.5 <- dt.3 %>% dplyr::filter(dt.3$group %in% dt.4$group)
  # drop ; or . or , at the end of the string
  dt.5$Results <- stringr::str_replace(string = dt.5$Results, pattern = ";$", replacement = "")
  dt.5$Results <- stringr::str_replace(string = dt.5$Results, pattern = "\\.$", replacement = "")
  dt.5$Results <- stringr::str_replace(string = dt.5$Results, pattern = ",$", replacement = "")
  dt.5$Results <- stringr::str_replace(string = dt.5$Results, pattern = "\\(", replacement = "")
  dt.5$Results <- stringr::str_replace(string = dt.5$Results, pattern = "\\)", replacement = "")
  dt.5$Results <- stringr::str_replace_all(string = dt.5$Results, pattern = "†", replacement = "")

  # remove the prefix before AA or NT
  dt.5a <- dt.5
  dt.5a$Results <- stringr::str_replace(string = dt.5a$Results, pattern = "^[pc]\\.", replacement = "")

  # Replace the string "variants:" that is in a few rows
  dt.5b <- dt.5a
  dt.5b$Results <- stringr::str_replace(
    string = dt.5b$Results,
    pattern = "variants:",
    replacement = ""
  )

  dt.5b$Results <- stringr::str_trim(dt.5b$Results)

  # transpose wide via paste0 so we can eliminate duplicates
  dt.6 <- dt.5b %>% dplyr::group_by(group) %>% dplyr::mutate(merge = base::paste0(Results, collapse = " ")) %>% ungroup()

  # slice to eliminate duplicates
  dt.7 <- dt.6 %>% dplyr::group_by(merge) %>% dplyr::slice_head()

  # Order by group to maintain the order in the report
  dt.8 <- dt.7 %>% dplyr::arrange(group)
  dt.9 <- dt.8 %>% select(merge) %>% rename(genes_variants = merge)


  return(dt.9)
}

