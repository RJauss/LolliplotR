
#' Plot Lolliplot
#'
#' @description Plots a Lolliplot and takes the cohort data from \code{\link{import_cohort}},
#' clinvar data from \code{\link{get_clinvar}}, gnomad data from \code{\link{get_gnomad}}
#' and the domain data from \code{\link{get_domains}}
#'
#' @param cohort The cohort dataframe (see \code{\link{import_cohort}} for details)
#' @param domains The domain data (see \code{\link{get_domains}} for details)
#' @param clinvar The clinvar data (see \code{\link{get_clinvar}} for details)
#' @param gnomad The gnomad data (see \code{\link{get_gnomad}} for details)
#' @param density_adjustment Enlarge or shrink the gnomad density plot, multiplied by this numeric value
#' @param color_domains_by How the domains should be coloured.
#' One of "accession" (default), which colors same repeated domains in the same colour
#' (e.g. all LDLa domains, but not LDLb domains), "type", which colours same types like repeats or receptors,
#' or "unique", where all domains will get unique colours
#' @param shorten_domainname Boolean. Should the domain name be shortened?
#' This will strip the text after the first underscore, e.g. "LDL_recept_a" will be "LDL"
#' @param .y The cohort column name for the y aesthetic. Defaults to "Cohort variant occurence"
#' @param .label The cohort column name for the label aesthetic. Example: "variant_p"
#' @param .linetype The cohort column name for the linetype aesthetic. Example: "origin_simple"
#'
#' @return A ggplot object (the Lolliplot)
#'
#' @import ggplot2
#' @import stringr
#' @import dplyr
#' @import grDevices
#'
#' @export
#'
#' @examples
#' LRP6_lolliplot = plot_lolliplot(cohort = LRP6_cohort, domains = LRP6_domains,
#' clinvar = LRP6_clinvar, gnomad = LRP6_gnomad, .linetype = "origin_simple")

plot_lolliplot = function(cohort = NULL,
                          domains = NULL,
                          clinvar = NULL,
                          gnomad = NULL,
                          density_adjustment = 500,
                          color_domains_by = "accession",
                          shorten_domainname = T,
                          .y = "Cohort variant occurence",
                          .label = NULL,
                          .linetype = NULL){

  # create ggplot object
  g = ggplot()

  ##### cohort data #####
  if(!is.null(cohort)){

    # make .y discrete value
    cohort[[.y]] = as.factor(cohort[[.y]])

    g = g +
      geom_segment(data = cohort, inherit.aes = F,
                   aes(x = protein_position, xend = protein_position,
                       y = 0, yend = !!checkifnull(.y), linetype = !!checkifnull(.linetype))) +
      geom_point(data = cohort, inherit.aes = F,
                 aes(x = protein_position, y = !!checkifnull(.y), shape = "Cohort Variants"),
                 size = 4)

    if(!is.null(.label)){
      g = g +
        geom_text(data = cohort, inherit.aes = F,
                  aes(x = protein_position, y = !!checkifnull(.y), label = !!checkifnull(.label)),
                  angle = 45, hjust = 0, vjust = 0.25, size = 3)
    }

  }

  ##### gnomad data #####
  if(!is.null(gnomad)){
    g = g +
      geom_density(data = gnomad, inherit.aes = F,
                   aes(x = protein_position, y = -(after_stat(density)*density_adjustment),
                       fill = "gnomAD Variants"),
                   alpha = 0.5) +
      geom_point(data = gnomad %>% filter(homozygote_count != 0), inherit.aes = F,
                 aes(x = protein_position, y = -0.25, fill = "gnomAD Variants",
                     shape = "gnomAD Variants (homozygous)"),
                 size = 4)

  }

  ##### clinvar data #####
  if(!is.null(clinvar)){
    g = g +
      geom_point(data = clinvar %>% filter(!is.na(protein_position)), inherit.aes = F,
                 aes(x = protein_position, y = -0.25, color = clinsig_simple,
                     shape = "ClinVar Variants"),
                 size = 4)
  }

  ##### domain data #####

  if(!is.null(domains)){
    # adjustments #
    protein = domains %>%
      filter(type == "protein")
    domains = domains %>%
      filter(type != "protein")

    # shorten domain name?
    if(shorten_domainname == T){
      domains = domains %>%
        mutate(shortened_label = str_extract(name_short, "^[^_]+"))
      domainlabel = "shortened_label"
    } else{
      domainlabel = "name_short"
    }

    # determine colors
    domaincolors = ukl_colors(n = length(unique(domains[[color_domains_by]])))

    # add colors as column
    domains = domains %>%
      group_by(!!sym(color_domains_by)) %>%
      mutate(ID = cur_group_id()) %>%
      mutate(color = domaincolors[ID])

    # update ggplot object
    g = g +
      geom_segment(data = protein , inherit.aes = F,
                   aes(x = start, xend = end, y = 0, yend = 0),
                   linewidth = 2.5, color = "grey50") +
      geom_segment(data = domains, inherit.aes = F,
                   aes(x = start, xend = end, y = 0, yend = 0),
                   linewidth = 10,
                   color = domains$color) +
      geom_text(data = domains,
                aes(x = (end-start)/2+start, label = get(domainlabel)),
                y = 0, color = "white", fontface = "bold",
                hjust = 0.5, size = 3)
  }


  ##### theme stuff #####
  g = g +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.title = element_text(face = "bold", size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(face = "bold", size = 14),
          plot.background = element_rect(fill = "white"),
          panel.grid.minor.y = element_blank()) +
    labs(linetype = .linetype, color = "Clinical Significance",
         y = .y, x = "Amino acid position", fill = "Missense variant density",
         shape = "Dataset") +
    scale_fill_manual(values = c("darkslategrey")) +
    scale_shape_manual(values = c(17, 18, 19)) +
    scale_y_discrete()

  return(g)
}
