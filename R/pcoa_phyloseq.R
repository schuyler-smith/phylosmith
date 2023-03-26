#' Create a ggplot object of the PCoA from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and
#' plots the PCoA of a treatment or set of treatments in space.
#' @useDynLib phylosmith
#' @usage pcoa_phyloseq(phyloseq_obj, treatment, x = 1, y = 2, method = 'bray',
#' circle = 0.95, colors = 'default', labels = NULL)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param treatment Column name as a string or number in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @param x the numerical principle compenent to use as the x-axis
#' @param y the numerical principle compenent to use as the y-axis
#' @param method the distance measurement algorithm to use, match to
#' "euclidean", "manhattan", "canberra", "clark", "bray", "kulczynski",
#' "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup",
#' "binomial", "chao", "cao" or "mahalanobis".
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
#' @importFrom vegan vegdist
#' @importFrom stats cmdscale
#' @export
#' @return ggplot-object
#' @examples pcoa_phyloseq(soil_column, c('Matrix', 'Treatment'),
#' circle = TRUE)

pcoa_phyloseq <- function(
  phyloseq_obj,
  treatment = NULL,
  x = 1,
  y = 2,
  method = "bray",
  circle = 0.95,
  colors = "default",
  labels = NULL
  ) {
  if (!(is.logical(circle) | is.numeric(circle) & circle >
    0 & circle <= 1)) {
    stop("`circle` must be either `TRUE`, `FALSE`, or a confidence interval",
      call. = FALSE)
  }
  check_args(
    phyloseq_obj   = phyloseq_obj,
    treatment      = treatment,
    x              = x,
    y              = y,
    dist_method    = method,
    labels         = labels
  )

  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment)
  treatment_name <- paste(treatment, collapse = sep)
  metadata <- as(phyloseq_obj@sam_data, "data.frame")
  color_count <- length(unique(metadata[[treatment_name]]))
  graph_colors <- create_palette(color_count, colors)
  otu_matrix <- as(phyloseq_obj@otu_table, "matrix")
  distance_matrix <- vegan::vegdist(t(otu_matrix), method, na.rm = TRUE)
  distance_matrix[is.na(distance_matrix)] <- 0
  MDS <- cmdscale(distance_matrix, k = max(c(x, y)), eig = TRUE)
  graph_data <- cbind(
    x = MDS$points[, x], 
    y = MDS$points[, y], 
    metadata
  )
  graph_data <-
    data.table::data.table(graph_data, keep.rownames = "Sample")

  g <- ggplot(data = graph_data, 
    aes_string("x", "y", group = treatment_name))
  if (circle == TRUE) {
    g <- g + stat_ellipse(geom = "polygon", type = "norm", size = 0.6,
      linetype = 1, alpha = 0.3, color = "black",
      aes_string(fill = treatment_name), show.legend = FALSE
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
  g <- g + geom_point(aes_string(fill = treatment_name),
      shape = 21, color = "black", size = 3, alpha = 1.0) +
    scale_fill_manual(values = graph_colors)
  if (is.character(labels)) {
    g <- g + geom_label(
      aes_string(label = labels, fill = treatment_name),
      label.padding = unit(0.35, "lines"),
      label.r = unit(0.55, "lines") ,
      show.legend = FALSE
    )
  }
  g <- g + labs(
      x = paste("PCo ", x, "  (", round(sum(MDS$eig[x]) /
        sum(MDS$eig[MDS$eig > 0]), 3) * 100, "%)", sep = ""),
      y = paste("PCo ", y, "  (", round(sum(MDS$eig[y]) /
        sum(MDS$eig[MDS$eig > 0]), 3) * 100, "%)", sep = "")
    ) +
    guides(
      fill = guide_legend(override.aes = list(size = 5))) +
      theme_schuy("pcoa")

  return(g)
}