#' Create a heatmap ggplot object of the abundance table from a
#' phyloseq object. Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and
#' creates heatmaps of the abundances across samples as a ggplot object.
#' The default color choice is the viridis palette, which is supposed to
#' be both aesthetic for normal and color-blind viewers.
#' @useDynLib phylosmith
#' @usage abundance_heatmap(phyloseq_obj, classification = NULL,
#' treatment, subset = NULL,
#' transformation = 'none', colors = 'default')
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about
#' each taxa/gene.
#' @param classification Column name as a string or number in the
#' \code{\link[phyloseq:tax_table]{tax_table}} for the factor.
#' @param treatment Column name as a string or number in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any
#' samples that to not contain this factor. This can be a vector of multiple
#' factors to subset on.
#' @param transformation Transformation to be used on the data. "none",
#' "relative_abundance", "log", "log10", "log1p", "log2", "asn", "atanh",
#' "boxcox", "exp", "identity", "logit", "probability", "probit",
#' "reciprocal", "reverse" and "sqrt"
#' @param colors Name of a color set from the
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palette of R-accepted
#' colors.
#' @importFrom stats reformulate
#' @importFrom ggraph scale_fill_viridis
#' @importFrom stringr str_to_title
#' @export
#' @return ggplot-object
#' @examples abundance_heatmap(soil_column, classification = 'Phylum',
#' treatment = c('Matrix', 'Treatment'), transformation = 'log')

abundance_heatmap <-
  function(phyloseq_obj,
           classification = NULL,
           treatment,
           subset = NULL,
           transformation = 'none',
           colors = 'default') {
    if (!inherits(phyloseq_obj, "phyloseq")) {
      stop("`phyloseq_obj` must be a
        phyloseq-class object", call. = FALSE)
    }
    classification <- check_index_classification(phyloseq_obj,
                                                 classification)
    if (!(is.null(classification)) &
        is.null(access(phyloseq_obj, 'tax_table'))) {
      stop(
        "`phyloseq_obj` must contain
        tax_table() information if `classification` argument is used",
        call. = FALSE
      )
    }
    if (!(is.null(classification)) &
        any(!(classification %in% colnames(access(
          phyloseq_obj, 'tax_table'
        ))))) {
      stop("`classification` must be a column
        from the the tax_table()",
           call. = FALSE)
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
        "`treatment` must be at least one
        column name, or index, from the sample_data()",
        call. = FALSE
      )
    }
    if (!(
      transformation %in% c(
        "none",
        "relative_abundance",
        "log",
        "log10",
        "log1p",
        "log2",
        "asn",
        "atanh",
        "boxcox",
        "exp",
        "identity",
        "logit",
        "probability",
        "probit",
        "reciprocal",
        "reverse",
        "sqrt"
      )
    )) {
      stop(
        "argument given to `transformation`
        not able to be applied by this function, please see help files for
        list of acceptable values",
        call. = FALSE
      )
    }
    options(warn = -1)
    phyloseq_obj <-
      taxa_filter(phyloseq_obj,
                  treatment = treatment,
                  frequency = 0,
                  subset = subset)
    if (!(is.null(classification))) {
      phyloseq_obj <- conglomerate_taxa(phyloseq_obj,
                                        classification, hierarchical = FALSE)
    }
    if (transformation == 'relative_abundance') {
      phyloseq_obj <- relative_abundance(phyloseq_obj)
    } else if (!(transformation %in% c('none', 'relative_abundance'))) {
      eval(parse(
        text = paste0(
          'otu_table(phyloseq_obj) <- ',
          transformation,
          ' (access(phyloseq_obj, "otu_table"))'
        )
      ))
    }
    treatment_name <- paste(treatment, collapse = sep)

    if (is.null(classification)) {
      classification <- 'OTU'
      if(is.null(treatment)){
        graph_data <- phyloseq(access(phyloseq_obj, 'otu_table'))
      } else {
        graph_data <- phyloseq(access(phyloseq_obj, 'otu_table'),
                               access(phyloseq_obj, 'sam_data')[, treatment_name])
      }
    } else {
      if(is.null(treatment)){
        graph_data <- phyloseq(
          access(phyloseq_obj, 'otu_table'),
          access(phyloseq_obj, 'tax_table')[, classification],
          access(phyloseq_obj, 'sam_data')
        )
      } else {
        graph_data <- phyloseq(
          access(phyloseq_obj, 'otu_table'),
          access(phyloseq_obj, 'tax_table')[, classification],
          access(phyloseq_obj, 'sam_data')[, treatment_name]
        )
      }

    }
    graph_data <- melt_phyloseq(graph_data)
    graph_data[[classification]] <-
      factor(graph_data[[classification]],
             levels = rev(unique(graph_data[[classification]])))
    set(graph_data, which(is.na(graph_data[[classification]])),
        classification, 'Unclassified')
    graph_data[['Sample']] <- factor(graph_data[['Sample']],
                                     levels = rownames(access(phyloseq_obj, 'sam_data')))
    setkey(graph_data, 'Sample', 'Abundance')
    set(graph_data, which(graph_data[['Abundance']] == 0), 'Abundance', NA)

    g <-
      ggplot(graph_data,
             aes_string('Sample', classification, fill = 'Abundance')) +
      geom_tile(color = "white", size = 0.25)
    if(!(is.null(treatment))){
      g <- g + facet_grid(reformulate(treatment_name),
                          scales = "free",
                          space = "free")
    }
    g <- g + theme_classic() +
      theme(
        axis.text.x = element_text(
          angle = 30,
          hjust = 1,
          size = 10
        ),
        axis.text.y = element_text(hjust = 0.95, size = 10),
        axis.title.x = element_text(size = 10, face = 'bold'),
        axis.title.y = element_text(size = 10, face = 'bold'),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 10, face = 'bold'),
        legend.text = element_text(size = 8),
        legend.spacing.x = unit(0.005, 'npc'),
        legend.key.size = unit(6, "mm"),
        legend.background = element_rect(
          fill = (alpha = 0),
          color = 'black',
          size = 0.25
        ),
        panel.background = element_rect(color = 'black', size = 1.4),
        strip.text.x = element_text(size = 10, face = 'bold'),
        strip.background = element_rect(colour = 'black', size = 1.4)
      ) +
      scale_x_discrete(expand = expand_scale(mult = 0, add = .53)) +
      if (colors == 'default') {
        scale_fill_viridis()
      } else {
        color_count <- 100
        graph_colors <- create_palette(color_count, colors)
        scale_fill_gradientn(colors = graph_colors)
      }
    if(transformation != 'none'){
      if(transformation == 'relative_abundance'){
        g <- g + labs(fill = str_to_title('Relative\nAbundance'))
      } else {
        g <- g  + labs(fill = str_to_title(paste(transformation, '\nAbundance')))
      }
    }
    return(g)
  }

#' Create a lineplot ggplot object of the abundance table from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and
#' creates line graphs with points across samples.
#' @useDynLib phylosmith
#' @usage abundance_lines(phyloseq_obj, classification = NULL,
#' treatment, subset = NULL,
#' relative_abundance = FALSE, points = TRUE, colors = 'default')
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param classification Column name as a string or number in the
#' \code{\link[phyloseq:tax_table]{tax_table}} for the factor to use for node
#' colors.
#' @param treatment Column name as a string or number in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any
#' samples that to not contain this factor. This can be a vector of multiple
#' factors to subset on.
#' @param relative_abundance If \code{TRUE}, transforms the abundance data
#' into relative abundance by sample.
#' @param points if \code{FALSE}, will not display the data-points.
#' @param colors Name of a color set from the
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palette of R-accepted
#' colors.
#' @export
#' @return ggplot-object
#' @examples abundance_lines(soil_column, classification = 'Phylum',
#' treatment = c('Matrix', 'Treatment'), relative_abundance = TRUE)

abundance_lines <-
  function(phyloseq_obj,
           classification = NULL,
           treatment,
           subset = NULL,
           relative_abundance = FALSE,
           points = TRUE,
           colors = 'default') {
    if (!inherits(phyloseq_obj, "phyloseq")) {
      stop("`phyloseq_obj` must be a
        phyloseq-class object", call. = FALSE)
    }
    classification <- check_index_classification(phyloseq_obj,
                                                 classification)
    if (!(is.null(classification)) &
        is.null(access(phyloseq_obj, 'tax_table'))) {
      stop(
        "`phyloseq_obj` must contain
        tax_table() information if `classification` argument is used",
        call. = FALSE
      )
    }
    if (!(is.null(classification)) & any(!(classification %in%
                                           colnames(access(
                                             phyloseq_obj, 'tax_table'
                                           ))))) {
      stop("`classification` must be a column from
        the the tax_table()",
           call. = FALSE)
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
        "`treatment` must be at least one
        column name, or index, from the sample_data()",
        call. = FALSE
      )
    }
    if (!(is.logical(relative_abundance))) {
      stop("`relative_abundance` must be either
        `TRUE`, or `FALSE`",
           call. = FALSE)
    }
    if (!(is.logical(points))) {
      stop("`points` must be either `TRUE`, or
        `FALSE`", call. = FALSE)
    }
    options(warn = -1)
    phyloseq_obj <-
      taxa_filter(phyloseq_obj,
                  treatment,
                  frequency = 0,
                  subset = subset)
    if (!(is.null(classification))) {
      phyloseq_obj <- conglomerate_taxa(phyloseq_obj, classification,
                                        hierarchical = FALSE)
    }
    if (relative_abundance) {
      phyloseq_obj <- relative_abundance(phyloseq_obj)
    }
    treatment_name <- paste(treatment, collapse = sep)

    if (is.null(classification)) {
      classification <- 'OTU'
      graph_data <- phyloseq(access(phyloseq_obj, 'otu_table'),
                             access(phyloseq_obj, 'sam_data')[, treatment_name])
    } else {
      graph_data <- phyloseq(
        access(phyloseq_obj, 'otu_table'),
        access(phyloseq_obj, 'tax_table')[, classification],
        access(phyloseq_obj, 'sam_data')[, treatment_name]
      )
    }
    graph_data <- data.table(melt_phyloseq(graph_data))
    graph_data[[classification]] <-
      factor(graph_data[[classification]],
             levels = rev(unique(graph_data[[classification]])))
    set(graph_data, which(is.na(graph_data[[classification]])),
        classification, 'Unclassified')
    graph_data[['Sample']] <- factor(graph_data[['Sample']],
                                     levels = rownames(access(phyloseq_obj, 'sam_data')))

    color_count <- length(unique(graph_data[[classification]]))
    graph_colors <- create_palette(color_count, colors)

    g <-
      ggplot(graph_data,
             aes_string(x = 'Sample', y = 'Abundance', group = classification))
    if (points == TRUE) {
      g <-
        g + geom_point(size = 1.5,
                       aes_string(color = classification),
                       show.legend = FALSE)
    }
    g <-
      g + geom_line(size = 1.2, aes_string(color = classification)) +
      facet_grid(reformulate(treatment_name), scales = "free", space = "free") +
      scale_colour_manual(values = graph_colors) +
      guides(colour = guide_legend(
        ncol = ceiling(length(unique(graph_data[[classification]])) / 25),
        override.aes = list(size = 4)
      ))
    if (relative_abundance == TRUE) {
      g <- g + ylab('Relative Abundance')
    }
    g <- g + theme_bw() +
      theme(
        axis.text.x = element_text(
          angle = 30,
          hjust = 1,
          size = 10
        ),
        axis.text.y = element_text(hjust = 0.95, size = 10),
        axis.title.x = element_text(size = 10, face = 'bold'),
        axis.title.y = element_text(size = 10, face = 'bold'),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 10, face = 'bold'),
        legend.text = element_text(size = 8),
        legend.spacing.x = unit(0.005, 'npc'),
        legend.key.size = unit(4, "mm"),
        legend.background = element_rect(fill = (alpha = 0)),
        panel.background = element_rect(
          color = 'black',
          size = 1.5,
          fill = 'white'
        ),
        panel.spacing = unit(.015, 'npc'),
        strip.text.x = element_text(
          size = 10,
          face = 'bold',
          color = 'black'
        ),
        strip.background = element_rect(
          colour = 'black',
          size = 1.4,
          fill = 'white'
        )
      ) +
      scale_y_continuous(expand = expand_scale(mult = c(0, 0.003), add = c(0.0015, 0.001))) +
      scale_x_discrete(expand = expand_scale(mult = 0, add = c(0.3, 0.5)))
    return(g)
  }

#' Create a boxplot of the alpha-diversity. Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and
#' creates boxplot of the alpha diversity as a ggplot object.
#' @useDynLib phylosmith
#' @usage alpha_diversity_graph(phyloseq_obj, index = 'shannon',
#' treatment = NULL, subset = NULL, colors = 'default')
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about
#' each taxa/gene.
#' @param index The diversity index to calculate ('shannon', 'simpson', 'invsimpson')
#' @param treatment Column name as a string or number in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any
#' samples that to not contain this factor. This can be a vector of multiple
#' factors to subset on.
#' @param colors Name of a color set from the
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palette of R-accepted
#' colors.
#' @export
#' @return ggplot-object
#' @examples alpha_diversity_graph(soil_column, index = 'shannon',
#' treatment = c('Matrix', 'Treatment'), subset = NULL, colors = 'default')

alpha_diversity_graph <- function(phyloseq_obj, index = 'shannon',
                                  treatment = NULL, subset = NULL, colors = 'default'){
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
                           Treatment = phyloseq_obj@sam_data[[treatment_name]],
                           Alpha = unlist(alpha))
  color_count <- length(unique(graph_data[['Treatment']]))
  graph_colors <- create_palette(color_count, colors)

  g <- ggplot(graph_data, aes(Treatment, Alpha, fill = Treatment))
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
    ) +
    labs(y = paste('Alpha-Diverity (', str_to_title(index), ' Index)', sep = ''))
  return(g)
}


#' Create a ggplot object of the NMDS from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and
#' plots the NMDS of a treatment or set of treatments in space.
#' @useDynLib phylosmith
#' @usage nmds_phyloseq(phyloseq_obj, treatment, circle = 0.95,
#' labels = NULL, colors = 'default', verbose = TRUE)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param treatment Column name as a string or number in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
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
        distance = "bray",
        k = 3,
        trymax = 100,
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

#' Create a ggplot object of the phylogenic barplots from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and
#' creates phylogenic barplots.
#' @useDynLib phylosmith
#' @usage phylogeny_profile(phyloseq_obj, classification = NULL,
#' treatment = NULL, subset = NULL, merge = TRUE, relative_abundance = FALSE,
#' colors = 'default', grid = FALSE)
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
#' @export
#' @return ggplot-object
#' @examples phylogeny_profile(soil_column, classification = 'Phylum',
#' treatment = c('Matrix', 'Treatment'), merge = TRUE,
#' relative_abundance = TRUE)

phylogeny_profile <-
  function(phyloseq_obj,
           classification = NULL,
           treatment = NULL,
           subset = NULL,
           merge = TRUE,
           relative_abundance = FALSE,
           colors = 'default',
           grid = FALSE) {
    if (!inherits(phyloseq_obj, "phyloseq")) {
      stop("`phyloseq_obj` must be a phyloseq-class
        object", call. = FALSE)
    }
    if (is.null(access(phyloseq_obj, 'sam_data'))) {
      stop("`phyloseq_obj` must contain
        sample_data() information",
           call. = FALSE)
    }
    classification <- check_index_classification(phyloseq_obj,
                                                   classification)
    if (!(is.null(classification)) &
        is.null(access(phyloseq_obj, 'tax_table'))) {
      stop(
        "`phyloseq_obj` must contain tax_table()
        information if `classification` argument is used",
        call. = FALSE
      )
    }
    if (!(is.null(classification)) &
        !(classification %in% colnames(access(phyloseq_obj, 'tax_table')))) {
      stop("`classification` must be a column from
        the the tax_table()",
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
    if (!(is.logical(merge))) {
      stop("`merge` must be either `TRUE`, or
        `FALSE`", call. = FALSE)
    }
    if (!(is.logical(relative_abundance))) {
      stop("`relative_abundance` must be either
        `TRUE`, or `FALSE`",
           call. = FALSE)
    }
    options(warn = -1)
    phyloseq_obj <-
      taxa_filter(phyloseq_obj,
                  treatment,
                  frequency = 0,
                  subset = subset)
    if (!(is.null(classification)) & merge) {
      phyloseq_obj <- conglomerate_taxa(phyloseq_obj, classification,
                                        hierarchical = FALSE)
    }
    if (relative_abundance) {
      phyloseq_obj <- relative_abundance(phyloseq_obj)
    }
    treatment_name <- paste(treatment, collapse = sep)

    if (is.null(classification)) {
      classification <- 'OTU'
      graph_data <- phyloseq(access(phyloseq_obj, 'otu_table'),
                             access(phyloseq_obj, 'sam_data'))
    } else {
      graph_data <- phyloseq(
        access(phyloseq_obj, 'otu_table'),
        access(phyloseq_obj, 'tax_table'),
        access(phyloseq_obj, 'sam_data')
      )
    }
    graph_data <- melt_phyloseq(graph_data)
    graph_data[[classification]] <-
      factor(graph_data[[classification]],
             levels = sort(unique(graph_data[[classification]])))
    set(graph_data, which(is.na(graph_data[[classification]])),
        classification, 'Unclassified')
    graph_data[['Sample']] <- factor(graph_data[['Sample']],
                                     levels = rownames(access(phyloseq_obj, 'sam_data')))

    color_count <- length(unique(graph_data[[classification]]))
    graph_colors <- create_palette(color_count, colors)

    g <-
      ggplot(graph_data,
             aes_string(x = "Sample", y = "Abundance", fill = classification))
    g <- g +
      guides(fill = guide_legend(ncol = ceiling(length(
        unique(graph_data[[classification]])
      ) / 50))) +
      scale_fill_manual(values = graph_colors, aesthetics = c('color', 'fill'))
    if (!(is.null(treatment))) {
      g <- g + facet_grid(reformulate(treatment_name), scales = "free", space = "free")
    }
    if (merge) {
      g <-
        g + geom_bar(
          aes_string(fill = classification),
          color = 'black',
          stat = 'identity',
          position = 'stack',
          size = 0.2,
          width = 0.95
        )
    } else {
      g <-
        g + geom_bar(
          stat = "identity",
          position = "stack",
          size = 0.12,
          width = 0.95,
          color = 'black'
        )
    }
    if (relative_abundance) {
      g <- g + ylab('Relative Abundance')
    }

    g <- g + theme_bw() +
      theme(
        axis.text.x = element_text(
          angle = 30,
          hjust = 1,
          size = 10
        ),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = 'bold'),
        axis.title.y = element_text(size = 10, face = 'bold'),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 10, face = 'bold'),
        legend.text = element_text(size = 8),
        legend.spacing.x = unit(0.005, 'npc'),
        legend.key.size = unit(4, "mm"),
        panel.background = element_rect(
          color = 'black',
          size = 1.5
        ),
        panel.spacing = unit(0.01, 'npc'),
        strip.text.x = element_text(size = 10, face = 'bold'),
        strip.background = element_rect(colour = 'black', size = 1.4, fill = 'white'),
        panel.grid.major.x = element_blank()
      ) +
      scale_y_continuous(expand = expand_scale(mult = c(0.0037, 0.003), add = c(0, 0))) +
      scale_x_discrete(expand = expand_scale(mult = 0, add = 0.51))
    if(grid){
      g <- g + facet_wrap(reformulate(treatment_name), scales = "free") + theme(axis.text.x = element_blank())
    }
    return(g)
  }

#' Create a ggplot object of the abundance barplots from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and
#' creates barplots of taxa by treatment.
#' @useDynLib phylosmith
#' @usage taxa_abundance_bars(phyloseq_obj, classification = NULL,
#' treatment, subset = NULL, transformation = 'none', colors = 'default',
#' wrap_by = NULL)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param classification Column name as a string or number in the
#' \code{\link[phyloseq:tax_table]{tax_table}} for the factor to use for node
#' colors.
#' @param treatment Column name as a string or number in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any
#' samples that to not contain this factor. This can be a vector of multiple
#' factors to subset on.
#' @param transformation Transformation to be used on the data. "none",
#' "mean", "median", "sd", "log", "log10"
#' @param colors Name of a color set from the
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palette of R-accepted
#' colors.
#' @param wrap_by Column name as a string or number in the
#' phyloseq \code{\link[phyloseq:sample_data]{sample_data}} to facet the ggplot
#' figure by.
#' @export
#' @return ggplot-object
#' @examples taxa_abundance_bars(
#' taxa_filter(soil_column, frequency = 0.8),
#' classification = 'Phylum', treatment = c('Matrix', 'Treatment'),
#' subset = 'Unamended', transformation = 'mean')

taxa_abundance_bars <-
  function(phyloseq_obj,
           classification = NULL,
           treatment,
           subset = NULL,
           transformation = 'none',
           colors = 'default',
           wrap_by = NULL) {
    if (!inherits(phyloseq_obj, "phyloseq")) {
      stop("`phyloseq_obj` must be a
        phyloseq-class object", call. = FALSE)
    }
    if (is.null(access(phyloseq_obj, 'sam_data'))) {
      stop("`phyloseq_obj` must contain
        sample_data() information",
           call. = FALSE)
    }
    classification <- check_index_classification(phyloseq_obj,
                                                   classification)
    if (!(is.null(classification)) &
        is.null(access(phyloseq_obj, 'tax_table'))) {
      stop(
        "`phyloseq_obj` must contain
        tax_table() information if `classification` argument is used",
        call. = FALSE
      )
    }
    if (!(is.null(classification)) &
        !(classification %in% colnames(access(phyloseq_obj, 'tax_table')))) {
      stop("`classification` must be a column from
        the the tax_table()",
           call. = FALSE)
    }
    treatment <- check_index_treatment(phyloseq_obj, treatment)
    wrap_by <- check_index_treatment(phyloseq_obj, wrap_by)
    if (any(!(treatment %in% colnames(access(
      phyloseq_obj, 'sam_data'
    ))))) {
      stop(
        "`treatment` must be at least one
        column name, or index, from the sample_data()",
        call. = FALSE
      )
    }
    if (!(transformation %in%
          c("none", "mean", "median", "sd", "log", "log10"))) {
      stop(
        "argument given to `transformation`
        not able to be applied by this function, please see help files for
        list of acceptable values",
        call. = FALSE
      )
    }
    options(warn = -1)
    phyloseq_obj <-
      taxa_filter(phyloseq_obj,
                  treatment,
                  frequency = 0,
                  subset = subset)
    if (!(is.null(classification))) {
      phyloseq_obj <- conglomerate_taxa(phyloseq_obj,
                                        classification, hierarchical = FALSE)
    } else {
      classification <- 'OTU'
    }
    treatment_name <- paste(treatment, collapse = sep)
    if(!(is.null(wrap_by)) && !(wrap_by %in% treatment)){
      treatment <- c(treatment, wrap_by)
    }

    graph_data <- melt_phyloseq(phyloseq_obj)
    graph_data[[classification]] <-
      factor(graph_data[[classification]],
             levels = unique(graph_data[[classification]]))
    if (transformation == 'none') {
      abundance <- 'Abundance'
      graph_data <- graph_data[, sum(Abundance),
                               by = c(treatment_name, treatment, classification)][, setnames(.SD, 'V1',
                                                                                             abundance, skip_absent = TRUE)]
    }
    if (transformation == 'mean') {
      abundance <- 'Mean_Abundance'
      graph_data <-
        graph_data[, mean(Abundance), by = c(treatment_name, treatment,
                                             classification)][, setnames(.SD, 'V1', abundance, skip_absent = TRUE)]
    }
    if (transformation == 'median') {
      abundance <- 'Median_Abundance'
      graph_data <- graph_data[, stats::median(Abundance),
                               by = c(treatment_name, treatment, classification)][, setnames(.SD, 'V1',
                                                                                             abundance, skip_absent = TRUE)]
    }
    if (transformation == 'sd') {
      abundance <- 'StdDev_Abundance'
      graph_data <-
        graph_data[, stats::sd(Abundance), by = c(treatment_name, treatment,
                                                  classification)][, setnames(.SD, 'V1', abundance, skip_absent = TRUE)]
    }
    if (transformation == 'log') {
      abundance <- 'log_Abundance'
      graph_data <-
        graph_data[, log(Abundance), by = c(treatment_name, treatment,
                                            classification)][, setnames(.SD, 'V1', abundance, skip_absent = TRUE)]
    }
    if (transformation == 'log10') {
      abundance <- 'log10_Abundance'
      graph_data <-
        graph_data[, log10(Abundance), by = c(treatment_name, treatment,
                                              classification)][, setnames(.SD, 'V1', abundance, skip_absent = TRUE)]
    }
    set(graph_data, which(is.na(graph_data[[classification]])),
        classification, 'Unclassified')

    color_count <- length(unique(graph_data[[treatment_name]]))
    graph_colors <- create_palette(color_count, colors)

    # graph_data[[classification]] <-
    #   factor(graph_data[[classification]], levels = sort(unique(as.character(graph_data[[classification]]))))

    g <-
      ggplot(graph_data,
             aes_string(x = classification, y = abundance, fill = treatment_name))
    g <-
      g +  geom_bar(
        stat = "identity",
        position = position_dodge2(padding = 3.5),
        size = 0.2,
        color = 'black',
        alpha = 0.85,
        width = 0.62
      ) +
      guides(colour = guide_legend(ncol = ceiling(length(
        unique(graph_data[[classification]])
      ) / 25))) +
      scale_fill_manual(values = graph_colors,
                        aesthetics = c('color', 'fill'))
    g <- g + theme_light() +
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
        axis.text.x = element_text(
          size = 10,
          vjust = 1,
          hjust = 1,
          angle = 30
        ),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10, face = "bold"),
        legend.background = element_rect(fill = (alpha = 0)),
        legend.key.size = unit(4, "mm"),
        legend.spacing.x = unit(0.005, 'npc'),
        strip.text.x = element_text(size = 10, face = 'bold', color = 'black'),
        strip.background = element_rect(colour = 'black', size = 1.4, fill = 'white')
      ) +
      scale_y_continuous(expand = expand_scale(mult = c(0.0025, 0.002)))
    if(!is.null(wrap_by)){
      g <- g + facet_wrap(reformulate(wrap_by))
    }
    return(g)
  }

#' Create graph of the core taxa seen in phyloseq-object over a range of
#' abundance and smaple-frequency values.  Function from the phylosmith-package.
#'
#' Inputs a phyloseq object and finds which taxa are seen in a
#' given proportion of samples at a minimum relative abundance, either in the
#' entire dataset, or by treatment, over a range of values. Then draws the
#' distribution.
#' @useDynLib phylosmith
#' @usage taxa_core_graph(phyloseq_obj, treatment = NULL, subset = NULL,
#' frequencies = seq(0.1, 1, 0.1), abundance_thresholds = seq(0.01, 1, 0.01),
#' colors = 'default')
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param treatment Column name as a \code{string} or \code{numeric} in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any
#' samples that to not contain this factor. This can be a vector of multiple
#' factors to subset on.
#' @param frequencies The range of proportions of samples the taxa are found in.
#' @param abundance_thresholds The range of minimum relative abundances the taxa are found
#' in for each sample.
#' @param colors Name of a color set from the
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palette of R-accepted
#' colors.
#' @export
#' @return phyloseq-object
#' @examples taxa_core_graph(soil_column)

taxa_core_graph <-
  function(phyloseq_obj,
           treatment = NULL,
           subset = NULL,
           frequencies = seq(0.1, 1, 0.1),
           abundance_thresholds = seq(0.01, 1, 0.01),
           colors = 'default'
  ) {
    phyloseq_obj <- taxa_filter(phyloseq_obj, treatment = treatment, subset = subset)
    treatment_name <- paste(treatment, collapse = sep)
    phyloseq_obj <- relative_abundance(phyloseq_obj)
    dataset <- melt_phyloseq(phyloseq(phyloseq_obj@otu_table, phyloseq_obj@sam_data))

    graph_data <- data.table()
    if(is.null(treatment)){
      N <- nsamples(phyloseq_obj)
      for(frequency in frequencies){
        for(abundance_threshold in abundance_thresholds){
          sub_table <- dataset[Abundance >= abundance_threshold]
          taxa <- nrow(sub_table[, .(count = .N), by = OTU][count >= floor(N*frequency)])
          if(taxa >= 0){
            graph_data <- rbind(graph_data,
                                data.table(
                                  freq = factor(frequency),
                                  abundance = abundance_threshold,
                                  taxa = taxa)
            )
          } else {
            graph_data <- rbind(graph_data, data.table(freq = factor(frequency), abundance = abundance_threshold, taxa = 0))
          }
        }
      }
    } else {
      treatments <- levels(dataset[[treatment_name]])
      for(treatment in treatments){
        N <- sum(phyloseq_obj@sam_data[[treatment_name]] == treatment)
        sub_table <- dataset[dataset[[treatment_name]] == treatment, ]
        for(frequency in frequencies){
          for(abundance_threshold in abundance_thresholds){
            sub_table <- sub_table[Abundance >= abundance_threshold]
            taxa <- nrow(sub_table[, .(count = .N), by = OTU][count >= floor(N*frequency)])
            if(taxa >= 0){
              graph_data <- rbind(graph_data,
                                  data.table(
                                    freq = factor(frequency),
                                    abundance = abundance_threshold,
                                    taxa = taxa,
                                    treatment = treatment)
              )
            } else {
              graph_data <- rbind(graph_data,
                                  data.table(freq = factor(frequency),
                                             abundance = abundance_threshold,
                                             taxa = 0,
                                             treatment = treatment))
            }
          }
        }
      }
      graph_data[['treatment']] <- factor(graph_data[['treatment']], levels = treatments)
    }

    graph_colors <- create_palette(length(frequencies), colors)

    g <- ggplot(graph_data, aes(x = abundance, y = taxa, color = freq)) +
      geom_line(size = 1.4) +
      scale_colour_manual(values = graph_colors) +
      guides(colour = guide_legend(
        ncol = ceiling(length(unique(graph_data[['freq']])) / 25),
        override.aes = list(size = 4)
      )) +
      theme_bw() +
      theme(
        axis.text.x = element_text(
          angle = 30,
          hjust = 1,
          size = 10
        ),
        axis.text.y = element_text(hjust = 0.95, size = 10),
        axis.title.x = element_text(size = 10, face = 'bold'),
        axis.title.y = element_text(size = 10, face = 'bold'),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 10, face = 'bold'),
        legend.text = element_text(size = 8),
        legend.spacing.x = unit(0.005, 'npc'),
        legend.background = element_rect(fill = (alpha = 0)),
        legend.key.size = unit(4, "mm"),
        panel.background = element_rect(
          color = 'black',
          size = 1.5,
          fill = 'white'
        ),
        panel.spacing = unit(.015, 'npc'),
        strip.text.x = element_text(
          size = 10,
          face = 'bold',
          color = 'black'
        ),
        strip.background = element_rect(
          colour = 'black',
          size = 1.4,
          fill = 'white'
        )
      ) +
      labs(x = 'Relative Abundance', y = 'Number of OTUs', color = 'Proportion\nof Samples') +
      scale_y_continuous(expand = expand_scale(add = c(0.3, 0.5))) +
      scale_x_continuous(expand = expand_scale(add = c(0.0005, 0.001)))
    if(!is.null(treatment)){
      g <- g + facet_wrap(~treatment, ncol = 3, dir = 'v')
    }
    return(g)
  }


#' Create a ggplot object using t-SNE from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object to plot the t-SNE of a
#' treatment or set of treatments.
#' @useDynLib phylosmith
#' @usage tsne_phyloseq(phyloseq_obj, treatment, perplexity = 10,
#' circle = TRUE, labels = NULL, colors = 'default')
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param treatment Column name as a string or number in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
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
      method = 'bray'),
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
#' @param method the distance measurement algorithm to use
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

  match.arg(method, c("manhattan", "euclidean", "canberra", "bray",
       "kulczynski", "gower", "morisita", "horn", "mountford",
       "jaccard", "raup", "binomial", "chao", "altGower", "cao",
       "mahalanobis", "clark"))

  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0)
  phyloseq_obj <- relative_abundance(phyloseq_obj)
  treatment_name <- paste(treatment, collapse = sep)
  metadata <- as(access(phyloseq_obj, "sam_data"), "data.frame")
  color_count <- length(unique(metadata[[treatment_name]]))
  graph_colors <- create_palette(color_count, colors)
  otu_matrix <- as(access(phyloseq_obj, "otu_table"), "matrix")
  distance_matrix <- vegdist(t(otu_matrix), method, na.rm = TRUE)
  distance_matrix[is.na(distance_matrix)] <- 0
  MDS <- cmdscale(distance_matrix, k = max(c(x,y)), eig = TRUE)
  graph_data <- cbind(x = MDS$points[,x], y = MDS$points[,y], metadata)

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


#' Computes the correlation of numerical variables with taxa
#' Function from the phylosmith-package.
#'
#' Computes the correlation of numerical variables with taxa
#' @useDynLib phylosmith
#' @usage variable_correlation_heatmap(phyloseq_obj, variables,
#' treatment = NULL, subset = NULL, classification = NULL,
#' method = 'spearman', limits = c(-0.8, 0.8),
#' colors = 'default', significance_color = 'white', cores = 1)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object.
#' @param variables Numerical factors within the in the
#' \code{\link[phyloseq:sample_data]{sample_data}} to correlate with the
#' abundance data.
#' @param treatment Column name as a \code{string} or \code{numeric} in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any
#' samples that to not contain this factor. This can be a vector of multiple
#' factors to subset on.
#' @param classification Column name as a \code{string} or \code{numeric} in
#' the \code{\link[phyloseq:tax_table]{tax_table}} for the factor to
#' conglomerate by.
#' @param method Which correlation method to calculate, "pearson", "spearman".
#' @param limits The range for the legend, smaller limits will accentuate smaller
#' correlations.
#' @param colors the palette to use for the heatmap, default is viridis.
#' @param significance_color the color to use for the significance stars.
#' @param cores \code{numeric} Number of CPU cores to use for the pair-wise
#' permutations. Default (1), (0) uses max cores available. Parallelization not
#' available for systems running MacOS without openMP configuration.
#' @importFrom parallel detectCores
#' @keywords nonparametric
#' @seealso \code{\link{permute_rho}} \code{\link{phylosmith}}
#' @export
#' @return data.table
#' @examples
#' variable_correlation_heatmap(soil_column, variables = 'Day',
#' treatment = c('Matrix', 'Treatment'), subset = 'Amended',
#' classification = 'Phylum', method = 'spearman', cores = 1,
#' colors = c("#2C7BB6", "white", "#D7191C"),
#' significance_color = 'black')

variable_correlation_heatmap <-
  function(phyloseq_obj,
           variables,
           treatment = NULL,
           subset = NULL,
           classification = NULL,
           method = 'spearman',
           limits = c(-0.8, 0.8),
           colors = 'default',
           significance_color = 'white',
           cores = 1) {
    correlations <- variable_correlation(
      phyloseq_obj,
      variables,
      treatment,
      subset,
      classification,
      method,
      cores
    )
    correlations$Significance<-cut(correlations$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
    if(is.null(treatment)){
      g <- ggplot(data = correlations, aes(x = Y, y = X, fill = rho))
    } else {
      g <- ggplot(data = correlations, aes(x = Treatment, y = X, fill = rho))
    }
    g <- g + geom_tile() +
      geom_text(aes(label = Significance), color = significance_color, size = 3) +
      labs(y = NULL, x = NULL, fill = method) +
      facet_grid(. ~ Y, drop = TRUE, scales = "free", space = "free_x")
    if('default' %in% colors){
      g <- g + ggraph::scale_fill_viridis(limit = limits)
    } else {
      g <- g + scale_fill_gradient2(low = colors[1], mid = colors[2], high = colors[3], limit = limits)
    }
    g <- g + theme_classic() +
      theme(
        axis.text.x = element_text(
          angle = 30,
          hjust = 1,
          size = 10
        ),
        axis.text.y = element_text(hjust = 0.95, size = 10),
        axis.title.x = element_text(size = 10, face = 'bold'),
        axis.title.y = element_text(size = 10, face = 'bold'),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 10, face = 'bold'),
        legend.text = element_text(size = 8),
        legend.spacing.x = unit(0.005, 'npc'),
        legend.key.size = unit(6, "mm"),
        panel.background = element_rect(color = 'black', size = 1.4),
        strip.text.x = element_text(size = 10, face = 'bold'),
        strip.background = element_rect(colour = 'black', size = 1.4)
      ) +
      scale_x_discrete(expand = expand_scale(mult = 0, add = 0))
    return(g)
  }


