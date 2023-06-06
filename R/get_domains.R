#' Get the pfam domain coordinates and descriptions
#'
#' @param gene the name of the gene, e.g. BRCA1
#'
#' @return Domains
#' @export
#'
#' @import httr
#' @import jsonlite
#' @import dplyr
#'
#' @examples
#' domains = get_domains("LRP6")
#'
get_domains = function(gene){
  ## query with gene name to get uniprot ID
  a = GET(paste0("https://www.ebi.ac.uk/ebisearch/ws/rest/uniprot?query=", gene, "&fields=gene_primary_name"))
  b = content(a, as = "text", type = "text/csv")
  c = fromJSON(b)

  # get the accession from the output
  # there are several results, we only want the one with the human entry
  uniprot_ID = (c$entries %>% filter(id == paste0(gene, "_HUMAN")))$acc


  # now with the uniprot ID we query interpro to get the pfam annoations with the API
  res= GET(paste0("https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/UniProt/", uniprot_ID))
  res_text = content(res, as = "text", type = "text/csv")
  Test = fromJSON(res_text)
  n_entries = Test$count
  results = Test$results
  protein_length = results$proteins[[1]]$protein_length

  ## loop over entries and get position for each domain ##
  # make empty tibble which will be filled with the info
  # the first line will be the protein
  Domains = tibble(accession = uniprot_ID, name = gene, name_short = gene,
                   type = "protein", start = 1, end = protein_length)

  # loop entries
  for(i in seq(1:n_entries)){
    # get metadata for each entry
    meta = slice(results$metadata, i)

    # get the fragments, they contain the amino acid positions
    fragments = bind_rows(results$proteins[[i]]$entry_protein_locations[[1]]$fragments) %>%
      select(-`dc-status`)

    # assign the metadata
    fragments$name = meta$name
    fragments$accession = meta$accession
    fragments$type = meta$type

    # get the short name per entry via the API
    res2= GET(paste0("https://www.ebi.ac.uk/interpro/api/entry/pfam/", meta$accession))
    res2_text = content(res2, as = "text", type = "text/csv")
    Test2 = fromJSON(res2_text)
    shortname = Test2$metadata$name$short
    fragments$name_short = shortname

    # each row should get a unique name based on the short name
    # paste short name and row number
    fragments = fragments %>%
      mutate(name_short_numbered = paste(name_short, row_number(), sep = "_"))

    # bind to the empty Domains tibble
    Domains = bind_rows(Domains, fragments)
  }

  return(Domains)
}
