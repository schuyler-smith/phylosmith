#' Create a ggplot object of the phylogenic barplots from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and
#' creates phylogenic barplots.
#' @useDynLib phylosmith
#' @usage phylogeny_profile(phyloseq_obj, treatment = NULL, subset = NULL,
#' classification = NULL, merge = TRUE, relative_abundance = FALSE,
#' colors = 'default', grid = FALSE, treatment_labels = NULL,
#' sample_labels = NULL, classification_labels= NULL)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param treatment Column name as a string or number in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any
#' samples that to not contain this factor. This can be a vector of multiple
#' factors to subset on.
#' @param classification Column name as a string or number in the
#' \code{\link[phyloseq:tax_table]{tax_table}} for the factor to use for node
#' colors.
#' @param merge if \code{FALSE}, will show separation of individuals within
#' each \code{classification}.
#' @param relative_abundance If \code{TRUE}, transforms the abundance data
#' into relative abundance by sample.
#' @param colors Name of a color set from the
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palette of R-accepted
#' colors.
#' @param grid Wraps the sub-plots into a grid pattern rather than side-by-side.
#' @param treatment_labels a vector of names to be used as labels for
#' treatments/facets.
#' @param sample_labels a vector of names to be used as labels for Samples.
#' @param classification_labels a vector of names to be used as labels for the
#' taxonomic classifications.
#' @export
#' @return ggplot-object
#' @examples phylogeny_profile(soil_column, classification = 'Phylum',
#' treatment = c('Matrix', 'Treatment'), merge = TRUE,
#' relative_abundance = TRUE)

phylogeny_profile <- function(
  phyloseq_obj,
  treatment = NULL,
  subset = NULL,
  classification = NULL,
  merge = TRUE,
  relative_abundance = FALSE,
  colors = "default",
  grid = FALSE,
  treatment_labels = NULL,
  sample_labels = NULL,
  classification_labels = NULL
) {
  check_args(
    phyloseq_obj   = phyloseq_obj,
    sample_data    = phyloseq_obj,
    classification = classification,
    treatment      = treatment,
    subset         = subset,
    merge          = merge,
    relative_abundance = relative_abundance,
    grid           = grid
  )
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, subset = subset)
  if (!(is.null(classification)) & merge) {
    phyloseq_obj <- 
    conglomerate_taxa(phyloseq_obj, classification, hierarchical = FALSE)
  }
  if (relative_abundance) phyloseq_obj <- relative_abundance(phyloseq_obj)
  treatment_name <- paste(treatment, collapse = sep)
  graph_data <- melt_phyloseq(phyloseq_obj)
  if (is.null(classification)) classification <- "OTU"
  data.table::set(graph_data, j = classification,
    value = factor(graph_data[[classification]],
      levels = rev(unique(graph_data[[classification]]))))
  data.table::set(graph_data, which(is.na(graph_data[[classification]])),
      classification, "Unclassified")
  data.table::set(graph_data, j = "Sample",
    value = factor(graph_data[["Sample"]],
      levels = rownames(phyloseq_obj@sam_data)))
  data.table::setkey(graph_data, "Sample", "Abundance")
  graph_data <- change_labels(graph_data, treatment_name, treatment_labels,
                      sample_labels, classification, classification_labels)

  color_count <- length(unique(graph_data[[classification]]))
  graph_colors <- rev(create_palette(color_count, colors))

  g <- ggplot(graph_data, 
    aes_string(x = "Sample", y = "Abundance", fill = classification))
  g <- g +
    guides(fill = guide_legend(ncol = ceiling(length(
      unique(graph_data[[classification]])) / 50))) +
    scale_fill_manual(values = graph_colors, 
      aesthetics = c("color", "fill"))
  if (!(is.null(treatment))) {
    g <- g + facet_grid(reformulate(treatment_name), 
      scales = "free", space = "free")
  }
  if (merge) {
    g <- g + geom_bar(
        aes_string(fill = classification),
        color = "black", stat = "identity", position = "stack", size = 0.2,
        width = 0.95)
  } else {
    g <- g + geom_bar(
      stat = "identity", position = "stack", size = 0.12, width = 0.95,
      color = "black")
  }
  if (relative_abundance) g <- g + ylab("Relative Abundance")
  g <- g +
    scale_y_continuous(expand = expansion(mult = c(0.0037, 0.003), 
      add = c(0, 0))) +
    scale_x_discrete(expand = expansion(mult = 0, add = 0.51))
  if(grid){
    g <- g + facet_wrap(reformulate(treatment_name), scales = "free") + 
      theme(axis.text.x = element_blank())
  }
  g <- g + theme_schuy("bar", 35)
  return(g)
}