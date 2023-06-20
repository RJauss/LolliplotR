#' get ClinVar variants from a specified gene via API
#'
#' @param gene The gene from which the ClinVar variants should be retreived, e.g. "KDM2A"
#' @param min_clinsig The minimum clinical significance, with `1` = "(Likely) Benign",
#' `2` = "Uncertain Significance" and `3` = "(Likely) Pathogenic" (default)
#'
#' @return A tibble containing the ClinVar variants with the transcript, gene name,
#' c_code, p_code and clinical significance
#' @export
#'
#' @import httr
#' @import jsonlite
#' @import dplyr
#' @import tidyr
#' @import stringr
#'
#' @examples
#' MBOAT7_clinvar_variants = get_clinvar("MBOAT7", min_clinsig = 2)

get_clinvar = function(gene, min_clinsig = 3){

  # use NCBI's esearch to get a list of IDs per gene
  res = GET(paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=",
                   gene, "[gene]+AND+single_gene[prop]&retmode=json&retmax=10000"))
  cont = content(res)

  idlist = unlist(cont$esearchresult$idlist)

  # get the total number of entries (for a progress bar)
  n_entries = as.numeric(cont$esearchresult$count)

  # empty tibble which will be filled with the variants in the next step
  variants = tibble(variant = character(), clinsig = character())

  # use the ids from the previous step to loop
  message(paste0("Reading ", n_entries, " ClinVar entries"))

  # generate a progress bar
  pb = txtProgressBar(min = 0, max = n_entries, initial = 0, style = 3)
  i = 1

  for(id in idlist){
    # print a convenient progress bar
    setTxtProgressBar(pb, i)

    # here we use NCBI's esummary to get the ClinVar entry per ID from the previous step
    res2 = GET(paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=",
                      id, "&retmode=json"))
    cont2 = content(res2)

    # title contains transcript, gene, c.code and p.code
    title = cont2$result[[id]]$title

    # here we find the clinical significance
    clinsig = cont2$result[[id]]$clinical_significance$description

    # bind each variant with the (empty) variant tibble
    variant_tibble = tibble(variant = title, clinsig = clinsig)
    variants = bind_rows(variants, variant_tibble)

    i = i+1
  }

  # close progress bar
  close(pb)

  # split variant title into transcript, gene, c.code and p.code
  variants_edited = variants %>%
    # remove variants based on gDNA position
    filter(!startsWith(variant, "NC_")) %>%
    # first remove brackets and blanks, replace by semicolon
    mutate(variant = str_remove_all(variant, " ")) %>%
    mutate(variant = str_replace_all(variant, "\\(", ";")) %>%
    mutate(variant = str_replace_all(variant, "\\):", ";")) %>%
    mutate(variant = str_remove_all(variant, "\\)")) %>%
    # then split the columns by the semicolon
    separate(variant, into = c("transcript", "gene", "c_code", "p_code"),
             sep = ";", fill = "right") %>%
    # extract the protein position from the p_code
    mutate(protein_position = as.numeric(str_extract(p_code, "[0-9]+"))) %>%
    # splice variants will be NA, replace that
    mutate(p_code = if_else(is.na(p_code), "p.?", p_code)) %>%
    # splice variants do not have a protein position, so calculate the expected position
    mutate(protein_position = if_else(is.na(protein_position),
                                      # divide c_code by 3 (codon triplet) and round up
                                      ceiling(as.numeric(str_extract(c_code, "[0-9]+"))/3),
                                      protein_position)) %>%
    # filter variants with conflicting and unprovided interpretations
    filter(!startsWith(clinsig, "Conflict")) %>%
    filter(clinsig != "not provided") %>%
    # translate clinsig into numbers for a later filtering step
    mutate(clinsig_numeric = if_else(str_detect(clinsig, "[Pp]atho"), 3, 0)) %>%
    mutate(clinsig_numeric = if_else(str_detect(clinsig, "Uncertain"), 2, clinsig_numeric)) %>%
    mutate(clinsig_numeric = if_else(str_detect(clinsig, "[Bb]enign"), 1, clinsig_numeric)) %>%
    # rename and collapse "likely" variants
    mutate(clinsig_simple = if_else(clinsig_numeric == 1, "(Likely) Benign", NA_character_)) %>%
    mutate(clinsig_simple = if_else(clinsig_numeric == 2, "Uncertain Significance", clinsig_simple)) %>%
    mutate(clinsig_simple = if_else(clinsig_numeric == 3, "(Likely) Pathogenic", clinsig_simple)) %>%
    # filter for minimum clinical significance (variable from function above)
    filter(clinsig_numeric >= min_clinsig)

  if(nrow(variants_edited) == 0){
    warning(paste0("\nBeware, the results table is empty. \nProbably because there are no variants with a minimum clinical significance of ", min_clinsig, "\n
                   Continute with the empty tibble or set the min_clinsig to a lower value"))
  }

  return(variants_edited)
}
