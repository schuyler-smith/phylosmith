#' Create a lineplot ggplot object of the abundance table from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and
#' creates line graphs with points across samples.
#' @useDynLib phylosmith
#' @usage abundance_lines(phyloseq_obj, treatment, subset = NULL,
#' classification = NULL, relative_abundance = FALSE, points = TRUE,
#' colors = 'default', treatment_labels = NULL, sample_labels = NULL,
#' classification_labels= NULL)
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
#' @param relative_abundance If \code{TRUE}, transforms the abundance data
#' into relative abundance by sample.
#' @param points if \code{FALSE}, will not display the data-points.
#' @param colors Name of a color set from the
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palette of R-accepted
#' colors.
#' @param treatment_labels a vector of names to be used as labels for
#' treatments/facets.
#' @param sample_labels a vector of names to be used as labels for Samples.
#' @param classification_labels a vector of names to be used as labels for the
#' taxonomic classifications.
#' @export
#' @return ggplot-object
#' @examples abundance_lines(soil_column, classification = 'Phylum',
#' treatment = c('Matrix', 'Treatment'), relative_abundance = TRUE)

abundance_lines <- function(
  phyloseq_obj,
  treatment = NULL,
  subset = NULL,
  classification = NULL,
  relative_abundance = FALSE,
  points = TRUE,
  colors = "default",
  treatment_labels = NULL,
  sample_labels = NULL,
  classification_labels = NULL
) {
  check_args(
    phyloseq_obj   = phyloseq_obj,
    classification = classification,
    treatment      = treatment,
    relative_abundance = relative_abundance,
    points         = points
  )

  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, subset = subset)
  if (relative_abundance) phyloseq_obj <- relative_abundance(phyloseq_obj)
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
  color_count <- length(unique(graph_data[[classification]]))
  graph_colors <- create_palette(color_count, colors)

  g <- ggplot(graph_data,
    aes_string(x = "Sample", y = "Abundance", group = classification))
  if (points) {
    g <-  g + 
      stat_summary(fun.y="sum", geom="point", size = 1.5,
        aes_string(color = classification), show.legend = FALSE)
  }
  g <- g +
    stat_summary(fun.y = "sum", geom = "line", size = 1.2, 
      aes_string(color = classification)) +
    scale_colour_manual(values = graph_colors) +
    guides(colour = guide_legend(
      ncol = ceiling(length(unique(graph_data[[classification]])) / 25),
      override.aes = list(size = 4)))
  g <- g +
    scale_y_continuous(expand = expansion(mult = c(0, 0.003), 
      add = c(0.0015, 0.001)))
  if (!is.null(treatment)) 
  g <- g +
    facet_grid(reformulate(treatment_name), 
      scales = "free", space = "free") 
  if (relative_abundance == TRUE) {
    g <- g + ylab("Relative Abundance")
  }
  g <- g + theme_schuy(angle = 35)

  return(g)
}