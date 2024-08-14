#' @author Schuyler D. Smith
#' Create a heatmap of the out_table from a
#' phyloseq-object.
#'
#' Uses a \code{\link[phyloseq]{phyloseq-class}} object as input and
#' creates a ggplot-heatmap of the abundances across samples.
#' The default color choice is the viridis palette.
#' @useDynLib phylosmith
#' @usage abundance_heatmap(phyloseq_obj, classification = NULL, 
#'   treatment =  NULL, subset = NULL, transformation = 'none', 
#'   colors = 'default', treatment_labels = NULL, sample_labels = NULL,
#'   classification_labels= NULL)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object.
#' @param classification The level of taxonomy to examine. Must be a column name
#' from the tax_table of the phyloseq_object.
#' \code{\link[phyloseq:tax_table]{tax_table}}.
#' @param treatment Column name as a string, or vector of strings, from the
#' \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param subset A level within the \code{treatment}. Multiple levels can be 
#' given as a vector.
#' @param transformation Transformation to be used on the data. "none",
#' "relative_abundance", "log", "log10", "log1p", "log2", "asn", "atanh",
#' "boxcox", "exp", "identity", "logit", "probability", "probit",
#' "reciprocal", "reverse" and "sqrt"
#' @param colors This can be either a name of a color set from the
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palette of R-accepted
#' colors. The default is an adaption of the palette from 
#' \url{https://www.nature.com/articles/nmeth.1618}
#' @param treatment_labels a vector of names to be used as labels for
#' treatments/facets, in the order they appear in the figure.
#' @param sample_labels a vector of names to be used as labels for Samples, 
#' in the order they appear in the figure.
#' @param classification_labels a vector of names to be used as labels for the
#' taxonomic classifications, in the order they appear in the figure.
#' @importFrom stats reformulate
#' @importFrom viridis scale_fill_viridis
#' @importFrom stringr str_to_title
#' @export
#' @return ggplot-object
#' @examples data(GlobalPatterns, package="phyloseq")
#' abundance_heatmap(GlobalPatterns, classification = "Phylum",
#'   treatment = "SampleType", transformation = "log2")

abundance_heatmap <- function(
  phyloseq_obj,
  classification = NULL,
  treatment = NULL,
  subset = NULL,
  transformation = "none",
  colors = "default",
  treatment_labels = NULL,
  sample_labels = NULL,
  classification_labels = NULL
) {
  check_args(
    phyloseq_obj     = phyloseq_obj,
    classification   = classification,
    treatment        = treatment,
    subset           = subset,
    transformation   = transformation,
    treatment_labels = treatment_labels,
    sample_labels    = sample_labels,
    classification_labels = classification_labels,
    colors           = colors
  )
  phyloseq_obj <- 
    taxa_filter(phyloseq_obj, treatment = treatment, subset = subset)
  if (!(is.null(classification))) {
    phyloseq_obj <- conglomerate_taxa(phyloseq_obj, classification, 
      hierarchical = FALSE)
  }
  if (transformation == "relative_abundance") {
    phyloseq_obj <- relative_abundance(phyloseq_obj)
  } else if (!(transformation %in% c("none", "relative_abundance"))) {
    eval(parse(
      text = paste0(
        "otu_table(phyloseq_obj) <- ",
        transformation,
        " (phyloseq::access(phyloseq_obj, 'otu_table'))"
      )
    ))
  }
  treatment_name <- paste(treatment, collapse = sep)
  if (is.null(classification)) classification <- "OTU"

  graph_data <- melt_phyloseq(phyloseq_obj)
  data.table::set(graph_data, j = classification,
    value = factor(graph_data[[classification]],
      levels = rev(unique(graph_data[[classification]]))))
  data.table::set(graph_data, which(is.na(graph_data[[classification]])),
    classification, "Unclassified")
  data.table::set(graph_data, j = "Sample",
    value = factor(graph_data[["Sample"]], 
    levels = rownames(phyloseq_obj@sam_data)))
  data.table::setkey(graph_data, "Sample", "Abundance")
  data.table::set(graph_data, which(graph_data[["Abundance"]] == 0), 
    "Abundance", NA)
  graph_data <- change_labels(graph_data, treatment_name, treatment_labels,
                      sample_labels, classification, classification_labels)
  g <- ggplot(graph_data, 
    aes_string("Sample", classification, fill = "Abundance")) +
    geom_tile(color = "white", size = 0.025)
  if(!(is.null(treatment))){
    g <- g + 
      facet_grid(reformulate(treatment_name), 
        scales = "free", space = "free")
  }
  g <- g + 
    scale_x_discrete(
        expand = expansion(mult = 0, add = .53)) +
    if (colors == "default") {
      viridis::scale_fill_viridis()
    } else {
      graph_colors <- create_palette(100, colors)
      scale_fill_gradientn(colors = graph_colors)
    }
  if(transformation != "none"){
    if(transformation == "relative_abundance"){
      g <- g + labs(fill = str_to_title("Relative\nAbundance"))
    } else {
      g <- g  + labs(fill = str_to_title(paste(transformation, "\nAbundance")))
    }
  }
  g <- g + theme_schuy("heat", 35)

  return(g)
}










