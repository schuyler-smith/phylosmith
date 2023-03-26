#' Create a ggplot object using t-SNE from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object to plot the t-SNE of a
#' treatment or set of treatments.
#' @useDynLib phylosmith
#' @usage tsne_phyloseq(phyloseq_obj, treatment, method = 'bray', perplexity = 10,
#' circle = TRUE, labels = NULL, colors = 'default')
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param treatment Column name as a string or number in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @param method the distance measurement algorithm to use, match to
#' "euclidean", "manhattan", "canberra", "clark", "bray", "kulczynski",
#' "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup",
#' "binomial", "chao", "cao" or "mahalanobis".
#' @param perplexity similar to selecting the number of neighbors to consider
#' in decision making (should not be bigger than 3 * perplexity < nrow(X) - 1,
#' see \code{\link[=Rtsne]{Rtsne}} for interpretation)
#' @param circle If TRUE, a \code{\link[ggplot2:stat_ellipse]{stat_ellipse}} around
#' each of the \code{treatment} factors (\code{TRUE}). If numeric between 0 and 1,
#' will add ellipse of confidence interval equal to value given (i.e. 0.95 produces
#' ellipses of 95\% confidence intervals)
#' @param labels Column name as a string or number in the
#' \code{\link[phyloseq:sample_data]{sample_data}} to use to place labels of
#' that factor instead of circle points.
#' @param colors Name of a color set from the
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palette of R-accepted
#' colors.
#' @importFrom Rtsne Rtsne
#' @importFrom vegan vegdist
#' @seealso \code{\link[=Rtsne]{Rtsne}}
#' @export
#' @return ggplot-object
#' @examples tsne_phyloseq(soil_column,
#' treatment = c('Matrix', 'Treatment'), perplexity = 8)

tsne_phyloseq <- function (
  phyloseq_obj,
  treatment,
  method = "bray",
  perplexity = 10,
  circle = TRUE,
  labels = NULL,
  colors = "default"
) {
  check_args(
    phyloseq_obj  = phyloseq_obj,
    treatment     = treatment,
    circle        = circle,
    labels        = labels,
    perplexity    = perplexity,
    dist_method   = method
  )
  phyloseq_obj    <- taxa_filter(phyloseq_obj, treatment)
  treatment_name  <- paste(treatment, collapse = sep)
  metadata        <- as(phyloseq_obj@sam_data, "data.frame")
  color_count     <- length(unique(metadata[[treatment_name]]))
  graph_colors    <- create_palette(color_count, colors)

  tsne <- Rtsne::Rtsne(vegan::vegdist(t(phyloseq_obj@otu_table),
    method = method), dims = 2, theta = 0.0, perplexity = perplexity
  )
  tSNE1 <- tsne$Y[, 1]/100
  tSNE2 <- tsne$Y[, 2]/100
  ord   <- data.table::data.table(tSNE1, tSNE2, metadata)
  ord   <- subset(ord, !is.na(treatment_name))
  if (is.character(labels)) {
    eval(parse(text = paste0("ord[, ", labels,
        " := access(phyloseq_obj, 'sam_data')[[labels]]]"
  )))}

  g <- ggplot(data = ord, aes_string("tSNE1", "tSNE2", group = treatment_name))
  if (circle) {
    g <- g + stat_ellipse(
      geom = "polygon",
      type = "norm",
      size = 0.6,
      linetype = 1,
      alpha = 0.3,
      color = "black",
      aes_string(fill = treatment_name),
      show.legend = FALSE
    ) +
      scale_color_manual(values = graph_colors) + guides(color = FALSE)
  } else if (is.numeric(circle)) {
    ellipse_df <-
      CI_ellipse(ggplot_build(g)$data[[1]], groups = "group", level = circle)
    g <- g + geom_polygon(data = ellipse_df, 
      aes(x = x, y = y, group = group),
      color = "black", fill = graph_colors[ellipse_df$group], 
      alpha = 0.3, size = 0.6, linetype = 1)
  }
  g <- g + geom_point(aes_string(fill = treatment_name), shape = 21,
    color = "black", size = 3, alpha = 1.0) +
  scale_fill_manual(values = graph_colors)
  if (is.character(labels)) {
    g <- g + geom_label(
      aes_string(label = labels, fill = treatment_name),
      label.padding = unit(0.35, "lines"), label.r = unit(0.55, "lines"),
      show.legend = FALSE
    )
  }
  g <- g + guides(
      fill = guide_legend(override.aes = list(size = 5))) +
      theme_schuy("tsne")

  return(g)
}
