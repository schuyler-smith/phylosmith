
#' Create a boxplot of the alpha-diversity. Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and
#' creates boxplot of the alpha diversity as a ggplot object.
#' @useDynLib phylosmith
#' @usage alpha_diversity_graph(phyloseq_obj, treatment = NULL, subset = NULL,
#' index = 'shannon', colors = 'default')
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
#' @param index The diversity index to calculate ('shannon', 'simpson', 'invsimpson')
#' @param colors Name of a color set from the
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palette of R-accepted
#' colors.
#' @export
#' @return ggplot-object
#' @examples alpha_diversity_graph(soil_column, index = 'shannon',
#' treatment = c('Matrix', 'Treatment'), subset = NULL, colors = 'default')

alpha_diversity_graph <- function(phyloseq_obj, treatment = NULL, subset = NULL,
                                  index = 'shannon', colors = 'default'){
  if (!inherits(phyloseq_obj, "phyloseq")) {
    stop("`phyloseq_obj` must be a
         phyloseq-class object", call. = FALSE)
  }
  indices <- c("shannon", "simpson", "invsimpson")
  index <- match.arg(index, indices)
  treatment <- check_index_treatment(phyloseq_obj, treatment)

  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment = treatment, subset = subset)
  treatment_name <- paste(treatment, collapse = sep)

  alpha <- data.table(as(phyloseq_obj@otu_table, 'matrix'))
  alpha <- alpha[,lapply(.SD,function(sample) sample/sum(sample))]
  if (index == "shannon"){
    alpha <- -alpha * log(alpha)
  } else {
    alpha <- alpha * alpha
  }
  alpha <- alpha[,lapply(.SD, sum, na.rm = TRUE)]
  if (index == "simpson") {
    alpha <- 1 - alpha
  } else if (index == "invsimpson"){
    alpha <- 1/alpha
  }

  graph_data <- data.table(Sample = sample_names(phyloseq_obj),
                           Alpha = unlist(alpha))
  graph_data <- merge(graph_data, as.data.table(as(phyloseq_obj@sam_data, 'data.frame'),
                                                keep.rownames = "Sample"), by = "Sample")
  color_count <- length(unique(graph_data[[treatment_name]]))
  graph_colors <- create_palette(color_count, colors)

  g <- ggplot(graph_data, aes_string(treatment_name, "Alpha", fill = treatment_name))
  g <- g + geom_boxplot(show.legend = FALSE) +
    scale_fill_manual(values = graph_colors) +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 30, hjust = 1,
        size = 10, face = 'bold'
      ),
      axis.text.y = element_text(hjust = 0.95, size = 10),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 10, face = 'bold'),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank()
    ) + scale_y_continuous(limits = c(0,5))
    labs(y = paste('Alpha-Diverity (', str_to_title(index), ' Index)', sep = ''))
  return(g)
}

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

nmds_phyloseq <-
  function(phyloseq_obj,
           treatment,
           method = 'bray',
           dimensions = 2,
           trymax = 100,
           circle = 0.95,
           labels = NULL,
           colors = 'default',
           verbose = TRUE) {
    if (!inherits(phyloseq_obj, "phyloseq")) {
      stop("`phyloseq_obj` must be a phyloseq-class
        object", call. = FALSE)
    }
    if (is.null(access(phyloseq_obj, 'sam_data'))) {
      stop("`phyloseq_obj` must contain
        sample_data() information",
           call. = FALSE)
    }
    treatment <- check_index_treatment(phyloseq_obj, treatment)
    if (any(!(treatment %in% colnames(access(
      phyloseq_obj, 'sam_data'
    ))))) {
      stop(
        "`treatment` must be at least one column
        name, or index, from the sample_data()",
        call. = FALSE
      )
    }
    if (!(is.logical(circle) |
          is.numeric(circle) & circle > 0 & circle <= 1)) {
      stop("`circle` must be either `TRUE`, `FALSE`, or a confidence interval",
           call. = FALSE)
    }
    labels <- check_index_treatment(phyloseq_obj, labels)
    if (!(is.null(labels)) &
        any(!(labels %in% colnames(access(
          phyloseq_obj,
          'sam_data'
        ))))) {
      stop("`labels` must be a column name, or
        index, from the sample_data()",
           call. = FALSE)
    }
    if (!(is.logical(verbose))) {
      stop("`verbose` must be either `TRUE`, or
        `FALSE`", call. = FALSE)
    }
    method <- match.arg(method, methods)
    options(warn = -1)
    phyloseq_obj <-
      taxa_filter(phyloseq_obj, treatment, frequency = 0)
    phyloseq_obj <- relative_abundance(phyloseq_obj)
    treatment_name <- paste(treatment, collapse = sep)
    metadata <- as(access(phyloseq_obj, 'sam_data'), 'data.frame')
    color_count <- length(unique(metadata[[treatment_name]]))
    graph_colors <- create_palette(color_count, colors)
    otu_matrix <- as(access(phyloseq_obj, 'otu_table'), 'matrix')
    MDS <-
      metaMDS(
        t(otu_matrix),
        autotransform = FALSE,
        distance = method,
        k = dimensions,
        trymax = trymax,
        trace = verbose
      )
    NMDS1 <- data.table(scores(MDS))$NMDS1
    NMDS2 <- data.table(scores(MDS))$NMDS2
    ord <- data.table(NMDS1, NMDS2, metadata)
    ord <- subset(ord, !is.na(treatment_name))
    if (is.character(labels)) {
      eval(parse(
        text = paste0(
          'ord[, ',
          labels,
          ' := access(phyloseq_obj, "sam_data")[[labels]]]'
        )
      ))
    }

    g <-
      ggplot(data = ord, aes_string('NMDS1', 'NMDS2', group = treatment_name))
    if (circle) {
      g <- g + stat_ellipse(
        geom = "polygon",
        type = "norm",
        size = 0.6,
        linetype = 1,
        alpha = 0.3,
        color = 'black',
        aes_string(fill = treatment_name),
        show.legend = FALSE
      ) +
        scale_color_manual(values = graph_colors) + guides(color = FALSE)
    } else if (is.numeric(circle)) {
      ellipse_df <-
        CI_ellipse(ggplot_build(g)$data[[1]],
                   groups = 'group',
                   level = circle)
      g <-
        g + geom_polygon(
          data = ellipse_df,
          aes(x = x, y = y, group = group),
          color = 'black',
          fill = graph_colors[ellipse_df$group],
          alpha = 0.3,
          size = 0.6,
          linetype = 1
        )
    }
    g <-
      g + geom_point(
        aes_string(fill = treatment_name),
        shape = 21,
        color = 'black',
        size = 5,
        alpha = 1.0
      ) +
      scale_fill_manual(values = graph_colors)
    if (is.character(labels)) {
      g <- g + geom_label(
        aes_string(label = labels,
                   fill = treatment_name),
        label.padding = unit(0.35, "lines"),
        label.r = unit(0.55, "lines") ,
        show.legend = FALSE
      )
    }
    g <- g + theme_classic() +
      theme(
        axis.line.x = element_line(
          colour = 'black',
          size = 1,
          linetype = 'solid'
        ),
        axis.line.y = element_line(
          colour = 'black',
          size = 1,
          linetype = 'solid'
        ),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 10, face = 'bold'),
        legend.text = element_text(size = 8),
        legend.spacing.x = unit(0.005, 'npc'),
        legend.key.size = unit(4, "mm"),
        legend.background = element_rect(fill = (alpha = 0))
      ) + labs(x = 'NMDS Dimension 1', y = 'NMDS Dimension 2') +
      guides(fill = guide_legend(
        override.aes = list(size = 4)
      ))
    return(g)
  }


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

pcoa_phyloseq <- function(phyloseq_obj,
                          treatment = NULL,
                          x = 1,
                          y = 2,
                          method = 'bray',
                          circle = 0.95,
                          colors = 'default',
                          labels = NULL){

  treatment <- check_index_treatment(phyloseq_obj, treatment)
  if (any(!(treatment %in% colnames(access(phyloseq_obj, "sam_data"))))) {
    stop("`treatment` must be at least one column\n        name, or index, from the sample_data()",
         call. = FALSE)
  }
  if (!(is.logical(circle) | is.numeric(circle) & circle >
        0 & circle <= 1)) {
    stop("`circle` must be either `TRUE`, `FALSE`, or a confidence interval",
         call. = FALSE)
  }
  labels <- check_index_treatment(phyloseq_obj, labels)
  if (!(is.null(labels)) & any(!(labels %in% colnames(access(phyloseq_obj,
                                                             "sam_data"))))) {
    stop("`labels` must be a column name, or\n        index, from the sample_data()",
         call. = FALSE)
  }

  method <- match.arg(method, methods)

  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0)
  treatment_name <- paste(treatment, collapse = sep)
  metadata <- as(access(phyloseq_obj, "sam_data"), "data.frame")
  color_count <- length(unique(metadata[[treatment_name]]))
  graph_colors <- create_palette(color_count, colors)
  otu_matrix <- as(access(phyloseq_obj, "otu_table"), "matrix")
  distance_matrix <- vegdist(t(otu_matrix), method, na.rm = TRUE)
  distance_matrix[is.na(distance_matrix)] <- 0
  MDS <- cmdscale(distance_matrix, k = max(c(x,y)), eig = TRUE)
  graph_data <- cbind(x = MDS$points[,x], y = MDS$points[,y], metadata)
  graph_data <- data.table::data.table(graph_data, keep.rownames = "Sample")

  g <- ggplot(data = graph_data, aes_string('x', 'y', group = treatment_name))
  if (circle) {
    g <- g + stat_ellipse(geom = "polygon", type = "norm",
                          size = 0.6, linetype = 1, alpha = 0.3, color = "black",
                          aes_string(fill = treatment_name), show.legend = FALSE) +
      scale_color_manual(values = graph_colors) + guides(color = FALSE)
  }
  else if (is.numeric(circle)) {
    ellipse_df <- CI_ellipse(ggplot_build(g)$data[[1]],
                             groups = 'group', level = circle)
    g <- g + geom_polygon(data = ellipse_df, aes(x = x,
                                                 y = y, group = group), color = "black", fill = graph_colors[ellipse_df$group],
                          alpha = 0.3, size = 0.6, linetype = 1)
  }
  g <- g + geom_point(aes_string(fill = treatment_name), shape = 21,
                      color = "black", size = 5, alpha = 1) + scale_fill_manual(values = graph_colors)
  if (is.character(labels)) {
    g <- g + geom_label(aes_string(label = labels, fill = treatment_name),
                        label.padding = unit(0.35, "lines"), label.r = unit(0.55,
                                                                            "lines"), show.legend = FALSE)
  }
  g <- g + theme_classic() + theme(axis.line.x = element_line(colour = "black",
                                                              size = 1, linetype = "solid"), axis.line.y = element_line(colour = "black",
                                                                                                                        size = 1, linetype = "solid"), axis.text.x = element_text(size = 10),
                                   axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 10,
                                                                                                      face = "bold"), axis.title.y = element_text(size = 10,
                                                                                                                                                  face = "bold"), legend.title = element_text(size = 10,
                                                                                                                                                                                              face = "bold"), legend.text = element_text(size = 8),
                                   legend.spacing.x = unit(0.005, "npc"), legend.key.size = unit(4,
                                                                                                 "mm"), legend.background = element_rect(fill = (alpha = 0))) +
    labs(x = paste("PCo ", x, "  (", round(sum(MDS$eig[x])/sum(MDS$eig[MDS$eig > 0]),3)*100, "%)", sep = ''),
         y = paste("PCo ", y, "  (", round(sum(MDS$eig[y])/sum(MDS$eig[MDS$eig > 0]),3)*100, "%)", sep = '')) +
    guides(fill = guide_legend(override.aes = list(size = 4)))
  return(g)
}

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

tsne_phyloseq <-
  function (phyloseq_obj,
            treatment,
            method = 'bray',
            perplexity = 10,
            circle = TRUE,
            labels = NULL,
            colors = 'default') {
    if (!inherits(phyloseq_obj, "phyloseq")) {
      stop("`phyloseq_obj` must be a phyloseq-class
        object", call. = FALSE)
    }
    if (is.null(access(phyloseq_obj, 'sam_data'))) {
      stop("`phyloseq_obj` must contain
        sample_data() information",
           call. = FALSE)
    }
    treatment <- check_index_treatment(phyloseq_obj, treatment)
    if (any(!(treatment %in% colnames(access(
      phyloseq_obj, 'sam_data'
    ))))) {
      stop(
        "`treatment` must be at least one column
        name, or index, from the sample_data()",
        call. = FALSE
      )
    }
    if (!(is.numeric(perplexity)) | perplexity <= 1) {
      stop("`perplexity` must be a numeric value
        greater than 1", call. = FALSE)
    }
    if (!(is.logical(circle) |
          is.numeric(circle) & circle > 0 & circle <= 1)) {
      stop("`circle` must be either `TRUE`, `FALSE`, or a confidence interval",
           call. = FALSE)
    }
    labels <- check_index_treatment(phyloseq_obj, labels)
    if (!(is.null(labels)) & any(!(labels %in% colnames(access(
      phyloseq_obj, 'sam_data'
    ))))) {
      stop("`labels` must be a column name, or
        index, from the sample_data()",
           call. = FALSE)
    }
    method <- match.arg(method, methods)
    options(warn = -1)
    phyloseq_obj <-
      taxa_filter(phyloseq_obj, treatment, frequency = 0)
    treatment_name <- paste(treatment, collapse = sep)
    metadata <- as(access(phyloseq_obj, 'sam_data'), 'data.frame')
    color_count <- length(unique(metadata[[treatment_name]]))
    graph_colors <- create_palette(color_count, colors)

    tsne <- Rtsne(
      vegdist(t(access(
        phyloseq_obj, 'otu_table'
      )),
      method = method),
      dims = 2,
      theta = 0.0,
      perplexity = perplexity
    )
    tSNE1 <- tsne$Y[, 1]/100
    tSNE2 <- tsne$Y[, 2]/100
    ord <- data.table(tSNE1, tSNE2, metadata)
    ord <- subset(ord, !is.na(treatment_name))
    if (is.character(labels)) {
      eval(parse(
        text = paste0(
          'ord[, ',
          labels,
          ' := access(phyloseq_obj, "sam_data")[[labels]]]'
        )
      ))
    }

    g <-
      ggplot(data = ord, aes_string('tSNE1', 'tSNE2', group = treatment_name))
    if (circle) {
      g <- g + stat_ellipse(
        geom = "polygon",
        type = "norm",
        size = 0.6,
        linetype = 1,
        alpha = 0.3,
        color = 'black',
        aes_string(fill = treatment_name),
        show.legend = FALSE
      ) +
        scale_color_manual(values = graph_colors) + guides(color = FALSE)
    } else if (is.numeric(circle)) {
      ellipse_df <-
        CI_ellipse(ggplot_build(g)$data[[1]],
                   groups = 'group',
                   level = circle)
      g <-
        g + geom_polygon(
          data = ellipse_df,
          aes(x = x, y = y, group = group),
          color = 'black',
          fill = graph_colors[ellipse_df$group],
          alpha = 0.3,
          size = 0.6,
          linetype = 1
        )
    }
    g <-
      g + geom_point(
        aes_string(fill = treatment_name),
        shape = 21,
        color = 'black',
        size = 5,
        alpha = 1.0
      ) +
      scale_fill_manual(values = graph_colors)
    if (is.character(labels)) {
      g <-
        g + geom_label(
          aes_string(label = labels, fill = treatment_name),
          label.padding = unit(0.35, "lines"),
          label.r = unit(0.55, "lines"),
          show.legend = FALSE
        )
    }
    g <- g + theme_classic() +
      theme(
        axis.line.x = element_line(
          colour = 'black',
          size = 1,
          linetype = 'solid'
        ),
        axis.line.y = element_line(
          colour = 'black',
          size = 1,
          linetype = 'solid'
        ),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = 'bold'),
        axis.title.y = element_text(size = 10, face = 'bold'),
        legend.title = element_text(size = 10, face = 'bold'),
        legend.text = element_text(size = 8),
        legend.spacing.x = unit(0.005, 'npc'),
        legend.background = element_rect(fill = (alpha = 0)),
        legend.key.size = unit(4, "mm")
      ) + labs(x = 't-SNE Dimension 1', y = 't-SNE Dimension 2') +
      guides(fill = guide_legend(
        override.aes = list(size = 4)
      ))
    return(g)
  }
