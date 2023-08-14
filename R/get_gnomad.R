#' Get gnomAD variants per gene
#'
#' @description Takes a gene name as an input and returns the gnomad missense variants, retrieved via API
#'
#' @param gene The gene name, e.g. "KDM2A"
#' @param account_outliers Boolean. Should we account for outliers, i.e. variants with an extraordinary high allele count?
#' @param outlier_threshold The threshold for the outliers, allele counts larger than this value will be collapsed (divided by `outlier_shrinkfactor`)
#' @param outlier_shrinkfactor Allele counts larger than `outlier_threshold` will be divided by this value
#'
#' @import jsonlite
#' @import dplyr
#' @import ghql
#' @import glue
#' @import stringr
#'
#' @return A tibble which contains the gnomad missense variants of the provided gene, with duplicated rows based on the allele count.
#' This will be needed for the histogram plot
#' @export
#'
#' @examples
#' KDM2A_gnomad = get_gnomad(gene = "KDM2A", account_outliers = FALSE)

get_gnomad = function(gene, account_outliers = TRUE,
                      outlier_threshold = 2000, outlier_shrinkfactor = 10){

  apiUrl = "https://gnomad.broadinstitute.org/api"

  # function to make GraphIql query at the gnomad API
  # inspired by Dayne Filer's gnomadR functions
  # takes the query as input (will be generated later)
  .makeAndEvalQuery = function(qfmt, maxTries = 3) {
    # set API URL as connection
    gmCon = GraphqlClient$new(url = apiUrl)
    # glue the query (see below)
    qryBody = glue_collapse(glue(qfmt), sep = "\n")
    # make new query
    qry = Query$new()$query('q', glue('query {{\n{qryBody}\n}}'))
    # establish maximum number of connections
    tries = 1
    repeat {
      if (tries > maxTries) break
      tryres = try(gmCon$exec(qry$q), silent = TRUE)
      if (!is(tryres, 'try-error')) break
      Sys.sleep(2*tries)
      tries = tries + 1
    }
    if (is(tryres, 'try-error')) {
      tryres = list(errorMessage = tryres[1], query = qry$q$query)
      class(tryres) = "failedQuery"
    }
    tryres
  }

  # this is the query for the API
  # note that the brackets need to be double: {{ and }}
  # with single brackets this could be pasted directly into the API (see URL)

  qfmt =
    '
      gene(gene_symbol: "{gene}", reference_genome: GRCh37) {{
            variants(dataset: gnomad_r2_1) {{
              consequence
              pos
              rsid
              variant_id: variantId
        			exome {{
        			  ac
        			  an
        			  homozygote_count
        			  hemizygote_count
        			  af
        			  ac_hom
        			  ac_hemi
        			}}
        			transcript_consequence {{
        			  gene_version
        			  gene_symbol
        			  hgvs
        			  hgvsc
        			  hgvsp
        			}}
            }}
          }}
    '

  # make the API query
  message("Reading gnomAD Data")
  res = .makeAndEvalQuery(qfmt)
  # convert from JSON
  reslist = fromJSON(res, flatten = TRUE)$data$gene$variants

  # select relevant columns and filter absent variants
  gnomad = reslist %>%
    select(variant_id, consequence, "allele_count" = exome.ac, "homozygote_count" = exome.ac_hom,
           "c_code" = transcript_consequence.hgvsc, "p_code" = transcript_consequence.hgvsp) %>%
    filter(!is.na(allele_count) & allele_count != 0) %>%
    mutate(p_code = if_else(is.na(p_code), "p.?", p_code))

  # subset to only include missense variants
  gnomad_missense = gnomad %>%
    filter(consequence == "missense_variant") %>%
    mutate(protein_position = as.numeric(str_extract(p_code, "[0-9]+")))

  # should we account for outliers? If yes, round them
  if(account_outliers == TRUE){
    gnomad_missense = gnomad_missense %>%
      mutate(allele_count = if_else(allele_count >= outlier_threshold,
             round(allele_count/outlier_shrinkfactor), allele_count))
  }

  # duplicates rows based on the allele count, will be used in density plot
  gnomad_missense_uncount = gnomad_missense %>%
    uncount(allele_count)

  message("Done!")

  return(gnomad_missense_uncount)
}
