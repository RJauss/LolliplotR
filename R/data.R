##### LRP6 Cohort #####
#' LRP6 cohort data
#'
#' Imported cohort data with \code{\link{import_cohort}}
#'
#' @format A dataframe
#' \describe{
#'   \item{Individual}{Individual ID}
#'   \item{variant_c}{c.Code of the variant}
#'   \item{variant_p}{p.Code of the variant}
#'   \item{origin}{origin, e.g. de novo or inherited}
#'   \item{origin_simple}{simplified origin, i.e. maternal=inherited}
#'   \item{classification}{classification of the variant}
#'   \item{Cohort variant occurence}{How often does this variant occur in the cohort}
#'   \item{protein_position}{amino acid position}
#' }
"LRP6_cohort"

##### MBOAT7 Cohort #####
#' MBOAT7 cohort data
#'
#' Imported cohort data with \code{\link{import_cohort}}
#'
#' @format A dataframe
#' \describe{
#'   \item{Individual}{Individual ID}
#'   \item{Family}{Family ID}
#'   \item{Sex}{Sex of the individual}
#'   \item{Age}{Age of the individual}
#'   \item{Zygosity}{Zygosity of the variant}
#'   \item{Variant_c}{c.Code of the variant}
#'   \item{Variant_p}{p.Code of the variant}
#'   \item{Consequence}{molecular consequence, missense or frameshift etc}
#'   \item{Classification}{classification of the variant}
#'   \item{LoF}{Loss of function variant? True or False}
#'   \item{Cohort variant occurence}{How often does this variant occur in the cohort}
#'   \item{protein_position}{amino acid position}
#' }
"MBOAT7_cohort"

##### LRP6 ClinVar #####
#' LRP6 clinvar data
#'
#' Got clinvar data with \code{\link{get_clinvar}}
#'
#' @format A dataframe
#' \describe{
#'   \item{transcript}{Transcript}
#'   \item{gene}{Gene}
#'   \item{c_code}{c.Code of the variant}
#'   \item{p_code}{p.Code of the variant}
#'   \item{clinsig}{classification of the variant}
#'   \item{protein_position}{amino acid position}
#'   \item{clinsig_numeric}{classification of the variant encoded as numeric}
#'   \item{clinsig_simple}{classification of the variant collapsed}
#' }
"LRP6_clinvar"

##### MBOAT7 ClinVar #####
#' MBOAT7 clinvar data
#'
#' Got clinvar data with \code{\link{get_clinvar}}
#'
#' @format A dataframe
#' \describe{
#'   \item{transcript}{Transcript}
#'   \item{gene}{Gene}
#'   \item{c_code}{c.Code of the variant}
#'   \item{p_code}{p.Code of the variant}
#'   \item{clinsig}{classification of the variant}
#'   \item{protein_position}{amino acid position}
#'   \item{clinsig_numeric}{classification of the variant encoded as numeric}
#'   \item{clinsig_simple}{classification of the variant collapsed}
#' }
"MBOAT7_clinvar"

##### LRP6 gnomAD #####
#' LRP6 clinvar data
#'
#' Got gnomad data with \code{\link{get_gnomad}}
#'
#' @format A dataframe
#' \describe{
#'   \item{variant_id}{Variant ID from gnomad}
#'   \item{consequence}{molecular consequence, missense or frameshift etc}
#'   \item{homozygote_count}{number of homozygotes}
#'   \item{c_code}{c.Code of the variant}
#'   \item{p_code}{p.Code of the variant}
#'   \item{protein_position}{amino acid position}
#' }
"LRP6_gnomad"

##### MBOAT7 gnomAD #####
#' MBOAT7 clinvar data
#'
#' Got gnomad data with \code{\link{get_gnomad}}
#'
#' @format A dataframe
#' \describe{
#'   \item{variant_id}{Variant ID from gnomad}
#'   \item{consequence}{molecular consequence, missense or frameshift etc}
#'   \item{homozygote_count}{number of homozygotes}
#'   \item{c_code}{c.Code of the variant}
#'   \item{p_code}{p.Code of the variant}
#'   \item{protein_position}{amino acid position}
#' }
"MBOAT7_gnomad"

##### LRP6 domains #####
#' LRP6 domain data
#'
#' Got domain data with \code{\link{get_domains}}
#'
#' @format A dataframe
#' \describe{
#'   \item{accession}{Pfam accession ID}
#'   \item{name}{domain name}
#'   \item{name_short}{short domain name}
#'   \item{start}{start coordinates}
#'   \item{end}{end coordinates}
#'   \item{name_short_numbered}{unique name for each domain}
#' }
"LRP6_domains"

##### MBOAT7 domains #####
#' MBOAT7 domain data
#'
#' Got domain data with \code{\link{get_domains}}
#'
#' @format A dataframe
#' \describe{
#'   \item{accession}{Pfam accession ID}
#'   \item{name}{domain name}
#'   \item{name_short}{short domain name}
#'   \item{start}{start coordinates}
#'   \item{end}{end coordinates}
#'   \item{name_short_numbered}{unique name for each domain}
#' }
"MBOAT7_domains"
