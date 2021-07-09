#' Abstract Genes with genetic variants from SNaPshot genomic reports
#' @description
#' `gene.variants.isolate.snapshot()` provides natural language processing tools to abstract genetic variants and gene names from SNaPshot reports
#' @param data  The data frame of a genomic report. Ideally this is information copied to the Clipboard from the EHR report by "selecting all" and then "copy". That is the default. If you don't use the clipboard function, you can use a data frame of the text file. The single column data frame should be labeled "Results". Required.
#'
#' @return a data frame with a column that has the gene names and variants in one cell
#' @examples
#' # Test with embedded data frame "snapshot_sample_report"
#' snapshot_sample_report %>%
#'   gene.variants.isolate.snapshot()
#' @export
#'
gene.variants.isolate.snapshot <- function(data = dplyr::tibble(Results = clipboard())){

  ##########################################################################################################################
  # Load Data and Trim
  ##########################################################################################################################
  dt <- data
  dt$Results <- stringr::str_trim(string = dt$Results) # Trim white space

  ##########################################################################################################################
  # Create regexs that we will use to select our data
  ##########################################################################################################################
  genes_boundary_df <- genetex::genes_boundary_regex() # this is a df of unique genes names concatenated with "|"
  genes_boundary_regex <- genes_boundary_df$genes # this is a character string of the genes for our regex
  nuc_regex <- "[ACTG]>[ACTG]|del[ACTG]" # This regex should identify those rows that have nucleotide changes
  aa_regex <- "(\\b([A-Z][0-9]{1,}(([A-Z])|(_[A-Z][1-9]{1,}del)|(fs\\*[1-9]{1,})|(\\*)|(fs)|(del)))|(p\\.[A-Z]))|([0-9]ins[A-Z])"
  genes_boundary_nuc_aa_regex <- paste(genes_boundary_regex, nuc_regex, aa_regex, sep = "|") # concatenate several regexs


  ##########################################################################################################################
  # Isolate Gene Variants
  ##########################################################################################################################
  # We need to eliminate much of this report in order to mine only those gene alterations relevant to each subject
  ## So we will use the word "test information" as a divider and take only those rows that are above it
  test_information <- "test information"
  dt.1a <- dt %>%
    mutate(divider = stringr::str_detect(string = dt$Results,
                                         pattern = stringr::regex(test_information, ignore_case = TRUE)),
           group = cumsum(divider)) %>%
    filter(group == 0)

  # tokenize the vector with space
  dt.1 <- splitstackshape::cSplit(
    indt = dt.1a,
    splitCols = "Results",
    sep = " ",
    direction = "long"
  )

  # Filter only those rows with the above genes_nuc_aa_regex
  dt.2 <- dt.1 %>%
    dplyr::filter(stringr::str_detect(string = dt.1$Results,
                                      pattern = stringr::regex(genes_boundary_nuc_aa_regex)))

  # Create a grouping system based on gene names
  dt.3 <- dt.2 %>%
    dplyr::mutate(keywords = stringr::str_detect(string = dt.2$Results,
                                                 pattern = stringr::regex(genes_boundary_regex)),
                  group = base::cumsum(keywords))
  # find those groups that are duplicated, but do not have more than 3 as SNaPshot reports should only have a
  ## Gene Name, Amino Acid variant and a Nucleotide variant. Groups with only one component or more than 3
  ### would be erroneous rows from the report
  dt.4 <- dt.3 %>% dplyr::count(group) %>% dplyr::filter(n >= 2 & n<= 3)
  # filter only those rows from dt.3 that are in the new dt.4
  dt.5 <- dt.3 %>% dplyr::filter(dt.3$group %in% dt.4$group)
  # drop ; or . or , at the end of the string
  dt.5$Results <- stringr::str_replace(string = dt.5$Results, pattern = ";$", replacement = "")
  dt.5$Results <- stringr::str_replace(string = dt.5$Results, pattern = "\\(p\\.=\\)", replacement = "")
  dt.5$Results <- stringr::str_replace(string = dt.5$Results, pattern = "\\.$", replacement = "")
  dt.5$Results <- stringr::str_replace(string = dt.5$Results, pattern = ",$", replacement = "")
  dt.5$Results <- stringr::str_replace(string = dt.5$Results, pattern = "\\(", replacement = "")
  dt.5$Results <- stringr::str_replace(string = dt.5$Results, pattern = "\\)", replacement = "")
  dt.5$Results <- stringr::str_replace_all(string = dt.5$Results, pattern = "†", replacement = "")
  dt.5$Results <- stringr::str_replace(string = dt.5$Results, pattern = "\\(p\\.=\\)", replacement = "")


  # remove all the information prior to either the Amino Acid or nucleotide sequence
  dt.5a <- dt.5
  dt.5a$Results <- stringr::str_replace(string = dt.5a$Results, pattern = "^(.+):[pc]\\.", replacement = "")

  # Remove any potential duplicates within groups (which can happen if there is an error in the report)
  dt.5b <- dt.5a %>% dplyr::group_by(group) %>% dplyr::distinct(Results)

  # transpose wide via paste0 so we can eliminate duplicate entries (different from above, this is a duplicate gene name
  ## plus variant vs. the above, which fixed situations in which there is two exact copies of a nucleotide variant under
  ### the same gene name
  dt.6 <- dt.5b %>% dplyr::group_by(group) %>% dplyr::mutate(merge = base::paste0(Results, collapse = " ")) %>% ungroup()
  # slice to eliminate duplicates
  dt.7 <- dt.6 %>% dplyr::group_by(merge) %>% dplyr::slice_head()
  # Order by group to maintain the order in the report
  dt.8 <- dt.7 %>% dplyr::arrange(group)
  dt.9 <- dt.8 %>% select(merge) %>% rename(genes_variants = merge)

  return(dt.9)
}

