
#' Create a boxplot of the alpha-diversity. Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and
#' creates boxplot of the alpha diversity as a ggplot object.
#' @useDynLib phylosmith
#' @usage alpha_diversity_graph(phyloseq_obj, treatment = NULL, subset = NULL,
#' index = "shannon", colors = "default")
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about
#' each taxa/gene.
#' @param treatment Column name as a string or number in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any
#' samples that to not contain this factor. This can be a vector of multiple
#' factors to subset on.
#' @param index The diversity index to calculate ("shannon", "simpson", "invsimpson")
#' @param colors Name of a color set from the
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palette of R-accepted
#' colors.
#' @export
#' @return ggplot-object
#' @examples alpha_diversity_graph(soil_column, index = "shannon",
#' treatment = c("Matrix", "Treatment"), subset = NULL, colors = "default")

alpha_diversity_graph <- function(
  phyloseq_obj,
  treatment = NULL,
  subset = NULL,
  index = "shannon",
  colors = "default"
) {
  check_args(
    phyloseq_obj    = phyloseq_obj,
    treatment       = treatment,
    subset          = subset,
    diversity_index = index,
    colors          = colors
  )
  phyloseq_obj <-
    taxa_filter(phyloseq_obj, treatment = treatment, subset = subset)
  treatment_name <- paste(treatment, collapse = sep)
  alpha <- data.table::data.table(as(phyloseq_obj@otu_table, "matrix"))
  alpha <- alpha[, lapply(.SD, function(sample) sample / sum(sample))]
  if (index == "shannon") {
    alpha <- -alpha * log(alpha)
  } else {
    alpha <- alpha * alpha
  }
  alpha <- alpha[, lapply(.SD, sum, na.rm = TRUE)]
  if (index == "simpson") {
    alpha <- 1 - alpha
  } else if (index == "invsimpson") {
    alpha <- 1 / alpha
  }
  graph_data <- data.table::data.table(
    Sample = sample_names(phyloseq_obj),
    Alpha = unlist(alpha)
  )
  graph_data <- merge(graph_data,
    data.table::as.data.table(as(phyloseq_obj@sam_data, "data.frame"),
    keep.rownames = "Sample"), by = "Sample")
  color_count <- length(unique(graph_data[[treatment_name]]))
  graph_colors <- create_palette(color_count, colors)

  g <- ggplot(graph_data,
    aes_string(treatment_name, "Alpha", fill = treatment_name))
  g <- g + geom_boxplot(show.legend = FALSE) +
    scale_fill_manual(values = graph_colors) +
    scale_y_continuous(limits = c(0, ceiling(max(g$data$Alpha)))) +
    labs(y = paste("Alpha-Diverity (", stringr::str_to_title(index),
      " Index)", sep = "")) +
    theme_schuy("box", angle = 35)
    
  return(g)
}