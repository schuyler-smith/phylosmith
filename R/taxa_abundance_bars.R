#' Create a ggplot object of the abundance barplots from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and
#' creates barplots of taxa by treatment.
#' @useDynLib phylosmith
#' @usage taxa_abundance_bars(phyloseq_obj, treatment = NULL, subset = NULL,
#' classification = NULL, transformation = "none", colors = "default",
#' wrap_by = NULL, treatment_labels = NULL, sample_labels = NULL,
#' classification_labels = NULL)
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
#' @param transformation Transformation to be used on the data. "none",
#' "mean", "median", "sd", "log", "log10"
#' @param colors Name of a color set from the
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palette of R-accepted
#' colors.
#' @param wrap_by Column name as a string or number in the
#' phyloseq \code{\link[phyloseq:sample_data]{sample_data}} to facet the ggplot
#' figure by.
#' @param treatment_labels a vector of names to be used as labels for
#' treatments/facets.
#' @param sample_labels a vector of names to be used as labels for Samples.
#' @param classification_labels a vector of names to be used as labels for the
#' taxonomic classifications.
#' @export
#' @return ggplot-object
#' @examples taxa_abundance_bars(
#' taxa_filter(soil_column, frequency = 0.8),
#' classification = 'Phylum', treatment = c('Matrix', 'Treatment'),
#' subset = 'Unamended', transformation = 'mean')

taxa_abundance_bars <- function(
  phyloseq_obj,
  treatment = NULL,
  subset = NULL,
  classification = NULL,
  transformation = "none",
  colors = "default",
  wrap_by = NULL,
  treatment_labels = NULL,
  sample_labels = NULL,
  classification_labels = NULL
) {
  check_args(
    phyloseq_obj     = phyloseq_obj,
    treatment        = treatment,
    subset           = subset,
    classification   = classification,
    transformation   = transformation,
    treatment_labels = treatment_labels
  )
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, subset)
  if (!(is.null(classification))) {
    phyloseq_obj <-
    conglomerate_taxa(phyloseq_obj, classification, hierarchical = FALSE)
  } else {
    classification <- "OTU"
  }
  treatment_name <- paste(treatment, collapse = sep)
  if(is.null(treatment)) treatment_name <- NULL 
  if(!(is.null(wrap_by)) && !(wrap_by %in% treatment)){
    treatment <- c(treatment, wrap_by)
  }

  graph_data <- melt_phyloseq(phyloseq_obj)
  data.table::set(graph_data, j = classification,
    value = factor(graph_data[[classification]],
      levels = unique(graph_data[[classification]])))
  if (transformation == "none") {
    abundance <- "Abundance"
    graph_data <- graph_data[, sum(Abundance), 
      by = c(treatment_name, treatment, classification)]
  }
  if (transformation == "relative_abundance") {
    abundance <- "Relative_Abundance"
    graph_data[, Relative_Abundance := Abundance/sum(Abundance), by = c(treatment_name)]
  }
  if (transformation == "mean") {
    abundance <- "Mean_Abundance"
    graph_data <- graph_data[, mean(Abundance),
      by = c(treatment_name, treatment, classification)]
  }
  if (transformation == "median") {
    abundance <- "Median_Abundance"
    graph_data <- graph_data[, stats::median(Abundance), 
      by = c(treatment_name, treatment, classification)]
  }
  if (transformation == "sd") {
    abundance <- "StdDev_Abundance"
    graph_data <- graph_data[, stats::sd(Abundance),
      by = c(treatment_name, treatment, classification)]
  }
  if (transformation == "log") {
    abundance <- "log_Abundance"
    graph_data <- graph_data[, log(Abundance), 
      by = c(treatment_name, treatment, classification)]
  }
  if (transformation == "log10") {
    abundance <- "log10_Abundance"
    graph_data <- graph_data[, log10(Abundance), 
      by = c(treatment_name, treatment, classification)]
  }
  graph_data <- graph_data[, data.table::setnames(.SD, "V1",
                                                  abundance, skip_absent = TRUE)]
  data.table::set(graph_data, which(is.na(graph_data[[classification]])),
      classification, "Unclassified")

  if (is.null(treatment)) {
    color_count <- 1
  } else {
    color_count <- length(unique(graph_data[[treatment_name]]))
  }
  graph_colors <- create_palette(color_count, colors)
  n_col <- ceiling(length(unique(graph_data[[classification]])) / 25)

  g <- ggplot(graph_data,
    aes_string(x = classification, y = abundance,
    fill = treatment_name))
  g <- g +  geom_bar(stat = "identity",
      position = position_dodge2(padding = 3.5), size = 0.2,
      color = "black", alpha = 0.85, width = 0.62) +
    guides(colour = guide_legend(ncol = n_col)) +
    scale_fill_manual(values = graph_colors, 
      aesthetics = c("color", "fill"))
  g <- g +  scale_y_continuous(expand = expansion(mult = c(0.0025, 0.002)))
  if(!is.null(wrap_by)){
    g <- g + facet_wrap(reformulate(wrap_by))
  }
  g <- g + theme_schuy("bar", 35)

  return(g)
}