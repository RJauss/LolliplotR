#' Get the pfam domain coordinates and descriptions
#'
#' @description Takes a gene name as an input and returns the Pfam domain coordinates, retrieved via API
#'
#' @param gene The name of the gene, e.g. "LRP6"
#'
#' @return A tibble containing the start/end positions and domain names of the specified gene
#' @export
#'
#' @import httr
#' @import jsonlite
#' @import dplyr
#'
#' @examples
#' LRP6_domains = get_domains("LRP6")

get_domains = function(gene){
  ## query with gene name to get uniprot ID
  message("Retreiving UniProt ID")
  a = GET(paste0("https://www.ebi.ac.uk/ebisearch/ws/rest/uniprot?query=", gene, "&fields=gene_primary_name"))
  b = suppressMessages(content(a, as = "text", type = "text/csv"))
  c = fromJSON(b)

  # get the accession from the output
  # there are several results, we only want the one with the human entry
  uniprot_id = (c$entries %>% filter(id == paste0(gene, "_HUMAN")))$acc

  # sometimes the gene name is not explicitly in the pfam ID
  # so try to catch invalid responses and then use the first line with "_HUMAN"

  if(identical(uniprot_id, character(0))){
    uniprot_id = (c$entries %>%
                     filter(grepl("_HUMAN", id)) %>%
                     filter(row_number() == 1))$acc
  } else {
    uniprot_id = uniprot_id
  }

  # now with the uniprot ID we query interpro to get the pfam annoations with the API
  message(paste0("Reading Interpro Data for ", uniprot_id))
  res= GET(paste0("https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/UniProt/", uniprot_id))
  res_text = suppressMessages(content(res, as = "text", type = "text/csv"))
  res_json = fromJSON(res_text)
  n_entries = res_json$count
  message(paste0("Found ", n_entries, " entries in Pfam"))
  results = res_json$results
  protein_length = results$proteins[[1]]$protein_length

  ## loop over entries and get position for each domain ##
  # make empty tibble which will be filled with the info
  # the first line will be the protein
  domains = tibble(accession = uniprot_id, name = gene, name_short = gene,
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
    message(paste0("Reading Pfam Data for ", meta$accession))
    res2= GET(paste0("https://www.ebi.ac.uk/interpro/api/entry/pfam/", meta$accession))
    res2_text = suppressMessages(content(res2, as = "text", type = "text/csv"))
    Test2 = fromJSON(res2_text)
    shortname = Test2$metadata$name$short
    fragments$name_short = shortname

    # each row should get a unique name based on the short name
    # paste short name and row number
    fragments = fragments %>%
      mutate(unique = paste(name_short, row_number(), sep = "_"))

    # bind to the empty domains tibble
    domains = bind_rows(domains, fragments)
  }

  message("Done!")

  return(domains)
}
