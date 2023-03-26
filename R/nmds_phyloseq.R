#' Create a ggplot object of the NMDS from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and
#' plots the NMDS of a treatment or set of treatments in space.
#' @useDynLib phylosmith
#' @usage nmds_phyloseq(phyloseq_obj, treatment, method = 'bray', dimensions = 2, trymax = 100, 
#' circle = 0.95, labels = NULL, colors = 'default', verbose = TRUE)
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
#' @param dimensions Number of dimensions. NB., the number of points n should 
#' be n > 2*k + 1, and preferably higher in global non-metric MDS, and still 
#' higher in local NMDS.
#' @param trymax Maximum numbers of random starts in search of stable 
#' solution. The iteration will stop when two convergent solutions were found or trymax was reached.
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
#' @param verbose Whether or not to print the
#' \code{\link[vegan:metaMDS]{metaMDS}} stress convergence to console
#' (\code{TRUE}) or not (\code{FALSE}).
#' @importFrom vegan metaMDS scores
#' @export
#' @return ggplot-object
#' @examples nmds_phyloseq(soil_column, c('Matrix', 'Treatment'),
#' circle = TRUE, verbose = FALSE)

nmds_phyloseq <- function(
  phyloseq_obj,
  treatment,
  method = "bray",
  dimensions = 2,
  trymax = 100,
  circle = 0.95,
  labels = NULL,
  colors = "default",
  verbose = TRUE
) {
  if (!(is.logical(circle) |
    is.numeric(circle) & circle > 0 & circle <= 1)) {
    stop("`circle` must be either `TRUE`, `FALSE`, or a confidence level",
      call. = FALSE)
  }
  check_args(
    phyloseq_obj   = phyloseq_obj,
    treatment      = treatment,
    dist_method    = method,
    labels         = labels,
    verbose        = verbose
  )
  phyloseq_obj <-
    taxa_filter(phyloseq_obj, treatment, subset)
  treatment_name <- paste(treatment, collapse = sep)
  metadata <- as(phyloseq_obj@sam_data, "data.frame")
  color_count <- length(unique(metadata[[treatment_name]]))
  graph_colors <- create_palette(color_count, colors)
  otu_matrix <- as(phyloseq_obj@otu_table, "matrix")
  MDS <- vegan::metaMDS(t(otu_matrix), autotransform = FALSE, 
    distance = method, k = dimensions, trymax = trymax, trace = verbose)
  NMDS1 <- data.table::data.table(vegan::scores(MDS)$sites)$NMDS1
  NMDS2 <- data.table::data.table(vegan::scores(MDS)$sites)$NMDS2
  ord <- data.table::data.table(NMDS1, NMDS2, metadata)
  ord <- subset(ord, !is.na(treatment_name))
  if (!is.null(labels)) {
    eval(parse(
      text = paste0(
        "ord[, ", labels, " := access(phyloseq_obj, 'sam_data')[[labels]]]"
      )
    ))
  }
  g <- ggplot(data = ord,
    aes_string("NMDS1", "NMDS2", group = treatment_name))
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
  g <- g + labs(x = "NMDS Dimension 1", y = "NMDS Dimension 2") +
      guides(
        fill = guide_legend(override.aes = list(size = 5))) + 
      theme_schuy("nmds")

  return(g)
}