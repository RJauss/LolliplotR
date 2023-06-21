#' Import the cohort data
#'
#' @description Function to import cohort data (.csv or .tsv) and adding variant frequency and protein position
#'
#' @param path_to_cohort filepath to the cohort data, can be .tsv or .csv
#' @param p_code_column the name of the column that contains the variant's p.code (required)
#'
#' @return A tibble containing the cohort data with variant frequency and protein positon added
#' @export
#'
#' @import dplyr
#' @import stringr
#' @importFrom data.table fread
#'
#' @examples
#' MBOAT7_cohort = import_cohort(path_to_cohort = "inst/extdata/MBOAT7_patient_variants.csv",
#' p_code_column = "Variant_p")

import_cohort = function(path_to_cohort, p_code_column){
  cohort = as_tibble(fread(path_to_cohort)) %>%
    # count the variants
    # "across" needed because we use a variable from the function
    add_count(across(p_code_column)) %>%
    # extract amino acid position
    # "get" needed because we use a variable from the function
    mutate(protein_position = as.numeric(str_extract(get(p_code_column), "[0-9]+")))

  return(cohort)
}

