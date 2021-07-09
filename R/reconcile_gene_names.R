#' Replace synonyms for gene names with those that correspond to the data dictionary in the MCC Registry
#' @description
#' `reconcile_gene_names()` reconciles synonyms import them to REDCap
#'
#' @param data A data frame from that is the downstream product of the genomics-to-redcap algorithm
#'
#' @return a data frame with gene names replaced
#' @export
#'
reconcile_gene_names <- function(data){
  dt_unite <- data

  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "FAM123B", "AMER1")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "BIRC3API2", "BIRC3")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "PD-L1", "CD274")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "C11orf30", "EMSY")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "FOX01", "FOXO1")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "C17orf39", "GID4")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "MLL", "KMT2A")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "MLL3", "KMT2C")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "MLL2", "KMT2D")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "MEK1", "MAP2K1")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "MEK2", "MAP2K2")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "MYCL1", "MYCL")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "CAN", "NUP214")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "NUT", "NUTM1")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "PD-1", "PDCD1")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "ARIDA", "ARID1A")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "WHSC1L1", "NSD3")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "WHSC1", "NSD2")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "TMEM173", "STING1")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "TCEB1", "ELOC")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "ROS", "ROS1")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "RFWD2", "COP1")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "PVRL4", "NECTIN4")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "PD-L2", "PDCD1LG2")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "PARK2", "PRKN")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "PAK7", "PAK5")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "MRE11A", "MRE11")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "HIST3H3", "H3-4")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "HIST2H3D", "H3C13")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "HIST2H3C", "H3C14")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "BRE", "BABAM2")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "C10orf54", "VSIR")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "C10orf86", "NSMCE4A")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "C17orf70", "FAAP100")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "C19orf40", "FAAP24")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "C1orf86", "FAAP20")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "DIRC2", "SLC49A4")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "FAM175A", "ABRAXAS1")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "FAM46C", "TENT5C")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "H3F3A", "H3-3A")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "H3F3B", "H3-3B")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "H3F3C", "H3-5")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "HIST1H1C", "H1-2")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "HIST1H1E", "H1-4")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "HIST1H2BD", "H2BC5")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "HIST1H3E", "H3C6")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "HIST1H3G", "H3C8")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "HIST1H3F", "H3C7")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "HIST1H3D", "H3C4")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "HIST1H3A", "H3C1")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "HIST1H3B", "H3C2")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "HIST1H3C", "H3C3")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "HIST1H3H", "H3C10")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "HIST1H3I", "H3C11")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "HIST1H3J", "H3C12")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "HIST1H4E", "H4C5")

  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "API2", "BIRC3")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "C15orf55", "NUTM1")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "C1orf4 or SMARCF1", "ARID1A")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "MMSET", "NSD2")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "MCF3", "ROS1")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "PDL2", "PDCD1LG2")

  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "RNF49", "BIRC3")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "SMARCF1", "ARID1A")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "c-ros-1", "ROS1")

  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "MALT2", "BIRC3")
  dt_unite$results <- replace(x = dt_unite$results, list = dt_unite$results == "GPR124", "ADGRA2")

  return(dt_unite)
}
