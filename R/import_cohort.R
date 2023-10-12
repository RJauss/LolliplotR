#' Import the cohort data
#'
#' @description Function to import cohort data (.csv or .tsv) and adding variant frequency and protein position
#'
#' @param path_to_cohort filepath to the cohort data, can be .tsv or .csv
#' @param variant_column the name of the column that contains either
#' the variant's p.code or c.code (the code is specified in the next variable)
#' @param code either "c" if you provide the variant's c.code or "p" if you provide the variant's p.code
#'
#' @return A tibble containing the cohort data with variant frequency and protein position added
#' @export
#'
#' @import dplyr
#' @import stringr
#' @importFrom data.table fread
#'
#' @examples
#' examplepath = system.file("extdata", "MBOAT7_patient_variants.csv", package="LolliplotR")
#' MBOAT7_cohort = import_cohort(path_to_cohort = examplepath, variant_column = "Variant_p", code = "p")

import_cohort = function(path_to_cohort, variant_column, code){
  cohort = as_tibble(fread(path_to_cohort)) %>%
    # count the variants
    # "across" needed because we use a variable from the function
    add_count(across(all_of(variant_column)), name = "Cohort variant occurence")
    # extract amino acid position
    # "get" needed because we use a variable from the function
    # if c.code is provided, we calculate the protein position manually
  if(code == "p"){
    cohort = cohort %>%
      mutate(protein_position = as.numeric(str_extract(get(variant_column), "[0-9]+")))
  } else if (code == "c"){
    cohort = cohort %>%
      mutate(protein_position = ceiling(as.numeric(str_extract(get(variant_column), "[0-9]+"))/3))
  } else {
    stop("You must either provide 'c' for the c.code or 'p' for the p.code")
  }



  return(cohort)
}

