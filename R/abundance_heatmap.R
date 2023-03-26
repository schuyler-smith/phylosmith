#' Create a heatmap of the abundance table from a
#' phyloseq object. Function from the phylosmith-package.
#'
#' Takes a \code{\link[phyloseq]{phyloseq-class}} object as input and
#' creates a ggplot-heatmap of the abundances across samples.
#' The default color choice is the viridis palette, which is supposed to
#' be both aesthetic for normal and color-blind viewers.
#' @useDynLib phylosmith
#' @usage abundance_heatmap(phyloseq_obj, treatment, subset = NULL,
#' classification = NULL, transformation = 'none', colors = 'default',
#' treatment_labels = NULL, sample_labels = NULL, classification_labels= NULL)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain.
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about
#' each taxa/gene.
#' @param treatment Column name as a string or number in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any
#' samples that to not contain this factor. This can be a vector of multiple
#' factors to subset on.
#' @param classification Column name as a string or number in the
#' \code{\link[phyloseq:tax_table]{tax_table}} for the factor.
#' @param transformation Transformation to be used on the data. "none",
#' "relative_abundance", "log", "log10", "log1p", "log2", "asn", "atanh",
#' "boxcox", "exp", "identity", "logit", "probability", "probit",
#' "reciprocal", "reverse" and "sqrt"
#' @param colors Name of a color set from the
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palette of R-accepted
#' colors.
#' @param treatment_labels a vector of names to be used as labels for
#' treatments/facets.
#' @param sample_labels a vector of names to be used as labels for Samples.
#' @param classification_labels a vector of names to be used as labels for the
#' taxonomic classifications.
#' @importFrom stats reformulate
#' @importFrom ggraph scale_fill_viridis
#' @importFrom stringr str_to_title
#' @export
#' @return ggplot-object
#' @examples abundance_heatmap(soil_column, classification = 'Phylum',
#' treatment = c('Matrix', 'Treatment'), transformation = 'log')

abundance_heatmap <- function(
  phyloseq_obj,
  treatment = NULL,
  subset = NULL,
  classification = NULL,
  transformation = "none",
  colors = "default",
  treatment_labels = NULL,
  sample_labels = NULL,
  classification_labels = NULL
) {
  check_args(
    phyloseq_obj   = phyloseq_obj,
    classification = classification,
    treatment      = treatment,
    transformation = transformation
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










