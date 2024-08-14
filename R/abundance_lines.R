#' @author Schuyler D. Smith
#' Create a lineplot ggplot object of the abundance table from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and
#' creates line graphs with points across samples.
#' @useDynLib phylosmith
#' @usage abundance_lines(phyloseq_obj, classification = NULL, treatment = NULL, 
#' subset = NULL, relative_abundance = FALSE, points = TRUE, colors = 'default', 
#' treatment_labels = NULL, sample_labels = NULL, classification_labels= NULL)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object.
#' @param classification The level of taxonomy to examine. Must be a column name
#' from the tax_table of the phyloseq_object.
#' \code{\link[phyloseq:tax_table]{tax_table}}.
#' @param treatment Column name as a string, or vector of strings, from the
#' \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param subset A level within the \code{treatment}. Multiple levels can be 
#' given as a vector.
#' @param relative_abundance If \code{TRUE}, transforms the abundance data
#' into relative abundance by sample.
#' @param points if \code{FALSE}, will not display the data-points.
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
#' @export
#' @return ggplot-object
#' @examples abundance_lines(soil_column, classification = 'Phylum',
#' treatment = c('Matrix', 'Treatment'), relative_abundance = TRUE)

abundance_lines <- function(
  phyloseq_obj,
  classification = NULL,
  treatment = NULL,
  subset = NULL,
  relative_abundance = FALSE,
  points = TRUE,
  colors = "default",
  treatment_labels = NULL,
  sample_labels = NULL,
  classification_labels = NULL
) {
  check_args(
    phyloseq_obj       = phyloseq_obj,
    classification     = classification,
    treatment          = treatment,
    subset             = subset,
    relative_abundance = relative_abundance,
    points             = points,
    treatment_labels   = treatment_labels,
    sample_labels      = sample_labels,
    classification_labels = classification_labels,
    colors             = colors
  )

  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, subset)
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