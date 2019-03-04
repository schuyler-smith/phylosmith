#' Create a ggplot object of the NMDS from a phyloseq object. Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and plots the NMDS of a treatment or set of treatments.
#' @useDynLib phylosmith
#' @usage nmds_phyloseq_ggplot(phyloseq_obj, treatment, circle = TRUE, colors = 'default',
#' verbose = TRUE)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package (must contain \code{\link[phyloseq:sample_data]{sample_data()}}).
#' @param treatment Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param circle If TRUE, add elipses around each treatment.
#' @param colors Name of a color set from the \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R accepted colors.
#' @param verbose Whether or not to print the \code{\link[vegan:metaMDS]{metaMDS}} convergion to console or not.
#' @import ggplot2
#' @importFrom vegan metaMDS scores
#' @export

nmds_phyloseq_ggplot <- function(phyloseq_obj, treatment, circle = TRUE, colors = 'default', verbose = TRUE){
  options(warn = -1)
  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0)
  treatment_name <- paste(treatment, collapse = '.')

  if(verbose == TRUE){MDS <- metaMDS(t(phyloseq_obj@otu_table), autotransform = FALSE, distance = "bray", k = 3, trymax = 100)}
  else if(verbose == FALSE){(invisible(capture.output(MDS <- metaMDS(t(phyloseq_obj@otu_table), autotransform = FALSE, distance = "bray", k = 3, trymax = 100))))}
  # plot(MDS, display = c("sites", "species"), choices = c(1,2), type = "p");abline(h=0,lty=2);abline(v=0,lty=2)
  # stressplot(MDS)

  Treatment <- phyloseq_obj@sam_data[[treatment_name]]

  color_count <- length(unique(Treatment))
  graph_colors <- create_palette(color_count, colors)

  NMDS1 <- data.table(scores(MDS))$NMDS1
  NMDS2 <- data.table(scores(MDS))$NMDS2
  ord <- data.table(NMDS1,NMDS2,Treatment)
  ord <- subset(ord, !is.na(Treatment))

  g <- ggplot(data = ord, aes(NMDS1, NMDS2))
  if(circle == TRUE){g <- g + stat_ellipse(geom = "polygon", type = "norm", size = 0.6, linetype = 1, alpha = 0.1, color = 'black', aes(fill = Treatment), show.legend = FALSE) +
      scale_color_manual(values = graph_colors) + guides(color = FALSE)}
  g <- g + geom_point(aes(fill = Treatment), shape = 21, color = 'black', size = 5, alpha = 1.0) +
    scale_fill_manual(values = graph_colors) +
    theme_light() +
    theme(aspect.ratio = 1,
          axis.line.x = element_line(colour = 'black', size = 1, linetype = 'solid'),
          axis.line.y = element_line(colour = 'black', size = 1, linetype = 'solid'),
          axis.text.x=element_text(size = 10, face = "bold"),
          axis.text.y=element_text(size = 10, face = "bold"),
          axis.title.x=element_text(size = 12, face= "bold"),
          axis.title.y=element_text(size = 12, face= "bold"),
          legend.title=element_blank(),
          legend.text=element_text(size = 11, face = "bold"),
          legend.background = element_rect(fill = (alpha = 0))
    )
  return(g)
}

#' Create a ggplot object using t-SNE from a phyloseq object. Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and plots the t-SNE of a treatment or set of treatments.
#' @useDynLib phylosmith
#' @usage tsne_phyloseq_ggplot(phyloseq_obj, treatment, perplexity = 10,
#' circle = TRUE, colors = 'default')
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \code{\link[=phyloseq]{phyloseq}} package (must contain \code{\link[phyloseq:sample_data]{sample_data()}}).
#' @param treatment Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param perplexity similar to selecting the number of neighbors to consider in decision making (should not be bigger than 3 * perplexity < nrow(X) - 1, see \code{\link[=Rtsne]{Rtsne}} for interpretation)
#' @param circle If TRUE, add elipses around each treatment.
#' @param colors Name of a color set from the \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R accepted colors.
#' @import ggplot2
#' @importFrom Rtsne Rtsne
#' @importFrom vegan vegdist
#' @seealso \code{\link[=Rtsne]{Rtsne}}
#' @export

tsne_phyloseq_ggplot <- function (phyloseq_obj, treatment, perplexity = 10, circle = TRUE, colors = 'default'){
  options(warn = -1)
  if (is.numeric(treatment)) {treatment <- colnames(phyloseq_obj@sam_data[, treatment])}
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0)
  treatment_name <- paste(treatment, collapse = ".")

  tsne <- Rtsne(vegdist(t(phyloseq_obj@otu_table), method = 'bray'), dims = 2, theta = 0.0, perplexity = perplexity)

  Treatment <- phyloseq_obj@sam_data[[treatment_name]]

  color_count <- length(unique(Treatment))
  graph_colors <- create_palette(color_count, colors)

  tSNE1 <- tsne$Y[,1]
  tSNE2 <- tsne$Y[,2]
  ord <- data.table(tSNE1, tSNE2, Treatment)
  ord <- subset(ord, !is.na(Treatment))

  g <- ggplot(data = ord, aes(tSNE1, tSNE2))
  if(circle == TRUE){g <- g + stat_ellipse(geom = "polygon", type = "norm", size = 0.6, linetype = 1, alpha = 0.1, color = 'black', aes(fill = Treatment), show.legend = FALSE) +
    scale_color_manual(values = graph_colors) + guides(color = FALSE)}
  g <- g + geom_point(aes(fill = Treatment), shape = 21, color = 'black', size = 5, alpha = 1.0) +
    scale_fill_manual(values = graph_colors) +
    theme_light() +
    theme(aspect.ratio = 1,
          axis.line.x = element_line(colour = 'black', size = 1, linetype = 'solid'),
          axis.line.y = element_line(colour = 'black', size = 1, linetype = 'solid'),
          axis.text.x=element_text(size = 10, face = "bold"),
          axis.text.y=element_text(size = 10, face = "bold"),
          axis.title.x=element_text(size = 12, face= "bold"),
          axis.title.y=element_text(size = 12, face= "bold"),
          legend.title=element_blank(),
          legend.text=element_text(size = 11, face = "bold"),
          legend.background = element_rect(fill = (alpha = 0))
    )
  return(g)
}

#' Create a ggplot object of the phylogenic barplots from a phyloseq object. Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and creates phylogenic barplots.
#' @useDynLib phylosmith
#' @usage phylogeny_bars_ggplot(phyloseq_obj, classification, treatment, subset = NULL,
#' merge = TRUE, relative_abundance = TRUE, colors = 'default')
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package (must contain \code{\link[phyloseq:sample_data]{sample_data()}}).
#' @param treatment Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param subset If taxa not needed to be seen in all \code{treatment}, then will subset to treatments containing this string.
#' @param classification Column name or number in the \code{\link[phyloseq:tax_table]{tax_table}}.
#' @param merge if TRUE, does not show separation of individuals within each \code{classification}. FALSE separates with black lines.
#' @param relative_abundance If TRUE, transforms the abundance data into relative abundance by sample.
#' @param colors Name of a color set from the \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R accepted colors.
#' @import ggplot2
#' @export

phylogeny_bars_ggplot <- function(phyloseq_obj, classification, treatment, subset = NULL, merge = TRUE, relative_abundance = TRUE, colors = 'default'){
  options(warn = -1)
  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  if(is.numeric(classification)){classification <- colnames(phyloseq_obj@tax_table[,classification])}
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0, subset = subset)
  if(relative_abundance == TRUE){phyloseq_obj <- relative_abundance(phyloseq_obj)}
  treatment_name <- paste(treatment, collapse = '.')

  # graph_data <- tax_glom(phyloseq_obj, taxrank = classification)
  graph_data <- phyloseq(phyloseq_obj@otu_table, phyloseq_obj@tax_table[,classification], phyloseq_obj@sam_data[,treatment_name])
  graph_data <- data.table(psmelt(graph_data))

  color_count <- length(unique(graph_data[[classification]]))
  graph_colors <- create_palette(color_count, colors)

  g <- ggplot(graph_data, aes_string(x = "Sample", y = "Abundance", fill = classification)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = -35, hjust = 0)) +
    guides(colour = guide_legend(ncol = ceiling(length(levels(graph_data[[classification]]))/30))) +
    facet_grid(treatment_name, scales = "free", space = "free") +
    scale_fill_manual(values = graph_colors, aesthetics = c('color', 'fill'))
  if(merge == TRUE){g <- g + geom_bar(aes_string(color = classification, fill = classification), stat = 'identity', position = 'stack', size = 0.2)
  } else {g <- g + geom_bar(stat = "identity", position = "stack", size = 0.12, color = 'black')}
  if(relative_abundance == TRUE){g <- g + ylab('Relative Abundance')}

  return(g)
}


#' Create a ggplot object of the abundance table from a phyloseq object. Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and creates line graphs across samples.
#' @useDynLib phylosmith
#' @usage abundance_lines_ggplot(phyloseq_obj, classification, treatment, subset = NULL,
#' relative_abundance = FALSE, points = TRUE, colors = 'default')
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package (must contain \code{\link[phyloseq:sample_data]{sample_data()}}).
#' @param treatment Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param subset If taxa not needed to be seen in all \code{treatment}, then will subset to treatments containing this string.
#' @param classification Column name or number in the \code{\link[phyloseq:tax_table]{tax_table}}
#' @param relative_abundance If TRUE, transforms the abundance data into relative abundance by sample.
#' @param points if TRUE, will diplay the data-points.
#' @param colors Name of a color set from the \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R accepted colors.
#' @import ggplot2
#' @export

abundance_lines_ggplot <- function(phyloseq_obj, classification, treatment, subset = NULL, relative_abundance = FALSE, points = TRUE, colors = 'default'){
  options(warn = -1)
  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  if(is.numeric(classification)){classification <- colnames(phyloseq_obj@tax_table[,classification])}
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0, subset = subset)
  if(relative_abundance == TRUE){phyloseq_obj <- relative_abundance(phyloseq_obj)}
  treatment_name <- paste(treatment, collapse = '.')

  # graph_data <- tax_glom(phyloseq_obj, taxrank = classification)
  graph_data <- phyloseq(phyloseq_obj@otu_table, phyloseq_obj@tax_table[,classification], phyloseq_obj@sam_data[,treatment_name])
  graph_data <- data.table(psmelt(graph_data))

  color_count <- length(unique(graph_data[[classification]]))
  graph_colors <- create_palette(color_count, colors)

  g <- ggplot(graph_data, aes_string(x = 'Sample', y = 'Abundance', group = classification)) +
    geom_line(size = 1.0, aes_string(color=classification))+
    theme_bw() +
    theme(axis.text.x = element_text(angle = -35, hjust = 0)) +
    guides(colour = guide_legend(ncol = ceiling(length(levels(graph_data[[classification]]))/30))) +
    facet_grid(treatment_name, scales = "free", space = "free") +
    scale_colour_manual(values = graph_colors)
  if(points == TRUE){g <- g + geom_point(size = 1.8, aes_string(color = classification))}
  if(relative_abundance == TRUE){g <- g + ylab('Relative Abundance')}

  return(g)
}

#' Create a ggplot object of the heatmap of the abundance table from a phyloseq object. Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and creates heatmaps of the abundances across samples.
#' @useDynLib phylosmith
#' @usage abundance_heatmap_ggplot(phyloseq_obj, classification = 'none', treatment, subset = NULL,
#' transformation = 'none', colors = 'default')
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package (must contain \code{\link[phyloseq:sample_data]{sample_data()}}).
#' @param treatment Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param subset If taxa not needed to be seen in all \code{treatment}, then will subset to treatments containing this string.
#' @param classification Column name or number in the \code{\link[phyloseq:tax_table]{tax_table}}.
#' @param transformation Transforms the data. "none", "relative_abundance", "log", "log10", "log1p", "log2", "asn", "atanh", "boxcox", "exp", "identity", "logit", "probability", "probit", "reciprocal", "reverse" and "sqrt"
#' @param colors Name of a color set from the \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R accepted colors.
#' @import ggplot2
#' @export

abundance_heatmap_ggplot <- function(phyloseq_obj, classification = 'none', treatment, subset = NULL, transformation = 'none', colors = 'default'){
  options(warn = -1)
  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  if(is.numeric(classification)){classification <- colnames(phyloseq_obj@tax_table[,classification])}
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0, subset = subset)
  if(transformation == 'relative_abundance'){phyloseq_obj <- relative_abundance(phyloseq_obj)}
  treatment_name <- paste(treatment, collapse = '.')

  color_count <- 10
  if(colors == 'default'){colors <- 'YlOrRd'}
  graph_colors <- create_palette(color_count, colors)

  if(classification == 'none'){classification <- 'OTU'; graph_data <- phyloseq(phyloseq_obj@otu_table, phyloseq_obj@tax_table[,1], phyloseq_obj@sam_data[,treatment_name])
  } else {graph_data <- phyloseq(phyloseq_obj@otu_table, phyloseq_obj@tax_table[,classification], phyloseq_obj@sam_data[,treatment_name])}
  graph_data <- data.table(psmelt(graph_data))
  graph_data[[classification]] <- factor(graph_data[[classification]], levels = rev(unique(graph_data[[classification]])))

  g <- ggplot(graph_data, aes_string('Sample', classification, fill = 'Abundance')) +
    geom_tile(color = "white", size = 0.25) +
    facet_grid(treatment_name, scales = "free", space = "free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = -35, hjust = 0)) +
    if(transformation %in% c('none', 'relative_abundance')){scale_fill_gradientn(colors = graph_colors)
    } else {scale_fill_gradientn(colors = graph_colors, trans = transformation)}

  return(g)
}

#' Create a ggplot object of the abundance barplots from a phyloseq object. Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and creates barplots of taxa by treatment.
#' @useDynLib phylosmith
#' @usage taxa_abundance_bars_ggplot(phyloseq_obj, classification = 'none', treatment,
#' subset = NULL, transformation = 'none', colors = 'default')
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package (must contain \code{\link[phyloseq:sample_data]{sample_data()}}).
#' @param treatment Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param subset If taxa not needed to be seen in all \code{treatment}, then will subset to treatments containing this string.
#' @param classification Column name or number in the \code{\link[phyloseq:tax_table]{tax_table}} for the x-axis.
#' @param transformation Transforms the data. "none", "mean", "median", "sd", "log", "log10"
#' @param colors Name of a color set from the \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R accepted colors.
#' @import ggplot2
#' @export

taxa_abundance_bars_ggplot <- function(phyloseq_obj, classification = 'none', treatment, subset = NULL, transformation = 'none', colors = 'default'){
  options(warn = -1)
  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  if(is.numeric(classification)){classification <- colnames(phyloseq_obj@tax_table[,classification])}
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0, subset = subset)
  treatment_name <- paste(treatment, collapse = '.')
  abundance <- 'Abundance'

  if(classification == 'none'){classification <- 'OTU'; graph_data <- phyloseq(phyloseq_obj@otu_table, phyloseq_obj@tax_table[,1], phyloseq_obj@sam_data[,treatment_name])
  } else {graph_data <- phyloseq(phyloseq_obj@otu_table, phyloseq_obj@tax_table[,classification], phyloseq_obj@sam_data[,treatment_name])}
  graph_data <- data.table(psmelt(graph_data))
  graph_data[[classification]] <- factor(graph_data[[classification]], levels = unique(graph_data[[classification]]))
  graph_data <- graph_data[, -'Sample']
  if(transformation == 'mean'){abundance <- 'Mean.Abundance'
    graph_data <- graph_data[, mean(Abundance), by = c(treatment_name, classification)][, setnames(.SD, 'V1', abundance, skip_absent = TRUE)]}
  if(transformation == 'median'){abundance <- 'Median.Abundance'
    graph_data <- graph_data[, stats::median(Abundance), by = c(treatment_name, classification)][, setnames(.SD, 'V1', abundance, skip_absent = TRUE)]}
  if(transformation == 'sd'){abundance <- 'StdDev.Abundance'
    graph_data <- graph_data[, stats::sd(Abundance), by = c(treatment_name, classification)][, setnames(.SD, 'V1', abundance, skip_absent = TRUE)]}
  if(transformation == 'log'){abundance <- 'log.Abundance'
    graph_data <- graph_data[, log(Abundance), by = c(treatment_name, classification)][, setnames(.SD, 'V1', abundance, skip_absent = TRUE)]}
  if(transformation == 'log10'){abundance <- 'log10.Abundance'
    graph_data <- graph_data[, log10(Abundance), by = c(treatment_name, classification)][, setnames(.SD, 'V1', abundance, skip_absent = TRUE)]}

  color_count <- length(unique(graph_data[[treatment_name]]))
  graph_colors <- create_palette(color_count, colors)

  g <- ggplot(graph_data, aes_string(x = classification, y = abundance, fill = treatment_name)) +
    geom_bar(stat = "identity", position = "dodge", size = 0.12, color = 'black') +
    theme_light() +
    theme(axis.text.x = element_text(angle = -35, hjust = 0)) +
    guides(colour = guide_legend(ncol = ceiling(length(levels(graph_data[[classification]]))/30))) +
    scale_fill_manual(values = graph_colors, aesthetics = c('color', 'fill'))

  return(g)
}
#
#' Internal function for creating color palettes for graphs. Function from the phylosmith-package.
#'
#' This function creates color palettes for graphs.
#' @useDynLib phylosmith
#' @usage create_palette(color_count, colors)
#' @param color_count Number of colors to choose for palette.
#' @param colors Name of a color set from the \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R accepted colors.
#' @import RColorBrewer
#' @import grDevices
#'

create_palette <- function(color_count, colors){
  options(warn = -1)
  cbcolors <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  if(colors == 'default'){colors <- cbcolors}
  if(any(!(colors %in% colors()))){
    if(any(colors %in% rownames(brewer.pal.info))){
      getPalette <- colorRampPalette(brewer.pal(min(c(color_count, brewer.pal.info[rownames(brewer.pal.info) == colors, 1])), colors))
    } else { getPalette <- colorRampPalette(colors)}
  } else { getPalette <- colorRampPalette(colors)}
  return(getPalette(color_count))
}
