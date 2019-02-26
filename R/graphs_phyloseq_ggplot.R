#' Create a ggplot object of the NMDS from a phyloseq object. Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and plots the NMDS of a treatment or set of treatments.
#' @useDynLib phylosmith
#' @usage nmds_phyloseq_ggplot(phyloseq_obj, treatment, circle = TRUE, colors = 'default')
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package (must contain \code{\link[phyloseq:sample_data]{sample_data()}}).
#' @param treatment Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param circle If TRUE, add elipses around each treatment.
#' @param colors Name of a color set from the \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R accepted colors.
#' @import ggplot2
#' @import RColorBrewer
#' @import vegan
#' @export

nmds_phyloseq_ggplot <- function(phyloseq_obj, treatment, circle = TRUE, colors = 'default'){
  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0)
  treatment <- paste(treatment, collapse = '.')

  MDS <- metaMDS(t(phyloseq_obj@otu_table), autotransform = FALSE, distance = "bray", k = 3, trymax = 100)
  # plot(MDS, display = c("sites", "species"), choices = c(1,2), type = "p");abline(h=0,lty=2);abline(v=0,lty=2)
  # stressplot(MDS)

  Treatment <- phyloseq_obj@sam_data[[treatment]]
  if(colors == 'default'){colors <- cbcolors}
  colorCount = length(unique(Treatment))
  if(any(!(colors %in% grDevices::colors()))){
    if(any(colors %in% rownames(RColorBrewer::brewer.pal.info))){
      getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(min(c(colorCount, RColorBrewer::brewer.pal.info[rownames(RColorBrewer::brewer.pal.info) == colors, 1])), colors))
    } else { getPalette <- grDevices::colorRampPalette(colors)}
  } else { getPalette <- grDevices::colorRampPalette(colors)}
  graph_colors = getPalette(colorCount)

  NMDS1 <- data.table(scores(MDS))$NMDS1
  NMDS2 <- data.table(scores(MDS))$NMDS2
  NMDS <- data.table(NMDS1,NMDS2,Treatment)
  NMDS.narm <- subset(NMDS, !is.na(Treatment))

  p <- ggplot(data = NMDS.narm, aes(NMDS1, NMDS2, color = NMDS.narm$Treatment)) +
    # coord_fixed(xlim = c(floor(min(NMDS.narm[,c(1,2)])), ceiling(max(NMDS.narm[,c(1,2)]))),
    #             ylim = c(floor(min(NMDS.narm[,c(1,2)])), ceiling(max(NMDS.narm[,c(1,2)])))) +
    geom_point(aes(color = Treatment), size = 1.5, alpha = 0.75) +
    scale_color_manual(values = graph_colors) +
    theme_classic() +
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
  if(circle == TRUE){
    p <- p + stat_ellipse(geom = "polygon", type = "norm", size = 0.6, linetype = 1, alpha = 0.0, aes(fill = NMDS.narm$Treatment))
  }
  return(p)
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
#' @import vegan
#' @import Rtsne
#' @seealso \code{\link[=Rtsne]{Rtsne}}
#' @export

tsne_phyloseq_ggplot <- function (phyloseq_obj, treatment, perplexity = 10, circle = TRUE, colors = 'default'){
  if (is.numeric(treatment)) {treatment <- colnames(phyloseq_obj@sam_data[, treatment])}
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0)
  treatment <- paste(treatment, collapse = ".")

  tsne <- Rtsne(vegdist(t(phyloseq_obj@otu_table), method = 'bray'), dims = 2, theta = 0.0, perplexity = perplexity)

  Treatment <- phyloseq_obj@sam_data[[treatment]]
  if(colors == 'default'){colors <- cbcolors}
  colorCount = length(unique(Treatment))
  if(any(!(colors %in% grDevices::colors()))){
    if(any(colors %in% rownames(RColorBrewer::brewer.pal.info))){
      getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(min(c(colorCount, RColorBrewer::brewer.pal.info[rownames(RColorBrewer::brewer.pal.info) == colors, 1])), colors))
    } else { getPalette <- grDevices::colorRampPalette(colors)}
  } else { getPalette <- grDevices::colorRampPalette(colors)}
  graph_colors = getPalette(colorCount)

  tsne1 <- tsne$Y[,1]
  tsne2 <- tsne$Y[,2]
  ord <- data.table(tsne1, tsne2, Treatment)
  ord <- subset(ord, !is.na(Treatment))

  p <- ggplot(data = ord, aes(tsne1, tsne2, color = ord$Treatment)) +
    geom_point(aes(color = ord$Treatment), size = 1.5, alpha = 1) +
    scale_color_manual(values = graph_colors) +
    theme_classic() +
    theme(aspect.ratio = 1,
          axis.line.x = element_line(colour = 'black', size = 1, linetype = 'solid'),
          axis.line.y = element_line(colour = 'black', size = 1, linetype = 'solid'),
          axis.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.title.x = element_text(size = 12, face= "bold"),
          axis.title.y = element_text(size = 12, face= "bold"),
          legend.title = element_blank(),
          legend.text = element_text(size = 11, face = "bold"),
          legend.background = element_rect(fill = (alpha = 0))
    )
  if(circle == TRUE){
    p <- p + stat_ellipse(geom = "polygon", type = "norm", size = 0.6, linetype = 1, alpha = 0.0, aes(fill = ord$Treatment))
  }
  return(p)
}

#' Create a ggplot object of the phylogenic barplots from a phyloseq object. Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and creates phylogenic barplots.
#' @useDynLib phylosmith
#' @usage phylogeny_bars_ggplot(phyloseq_obj, classification_level, treatment, subset = NULL,
#' merge = TRUE, relative_abundance = TRUE, colors = 'default')
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package (must contain \code{\link[phyloseq:sample_data]{sample_data()}}).
#' @param classification_level Column name or number in the \code{\link[phyloseq:tax_table]{tax_table}}.
#' @param treatment Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param subset If taxa not needed to be seen in all \code{treatment}, then can check only one particular treatment subset, this works for multiple treatment inputs.
#' @param merge if TRUE, does not show separation of individuals within each \code{classification_level}. FALSE separates with black lines.
#' @param relative_abundance If TRUE, transforms the abundance data into relative abundance by sample.
#' @param colors Name of a color set from the \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R accepted colors.
#' @import ggplot2
#' @import RColorBrewer
#' @export

phylogeny_bars_ggplot <- function(phyloseq_obj, classification_level, treatment, subset = NULL, merge = TRUE, relative_abundance = TRUE, colors = 'default'){

  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  if(is.numeric(classification_level)){classification_level <- colnames(phyloseq_obj@tax_table[,classification_level])}
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0, subset = subset)
  if(relative_abundance == TRUE){phyloseq_obj <- relative_abundance(phyloseq_obj)}
  treatment <- paste(treatment, collapse = '.')

  # graph_data <- tax_glom(phyloseq_obj, taxrank = classification_level)
  graph_data <- phyloseq(phyloseq_obj@otu_table, phyloseq_obj@tax_table[,classification_level], phyloseq_obj@sam_data[,treatment])
  graph_data <- data.table(psmelt(graph_data))

  colorCount = length(unique(phyloseq_obj@tax_table[,classification_level]))
  if(colors == 'default'){colors <- cbcolors}
  if(any(!(colors %in% grDevices::colors()))){
    if(any(colors %in% rownames(RColorBrewer::brewer.pal.info))){
      getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(min(c(colorCount, RColorBrewer::brewer.pal.info[rownames(RColorBrewer::brewer.pal.info) == colors, 1])), colors))
    } else { getPalette <- grDevices::colorRampPalette(colors)}
  } else { getPalette <- grDevices::colorRampPalette(colors)}
  graph_colors = getPalette(colorCount)

  p <- ggplot(graph_data, aes_string(x = "Sample", y = "Abundance", fill = classification_level)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = -35, hjust = 0)) +
    guides(colour = guide_legend(ncol = ceiling(length(levels(graph_data[[classification_level]]))/30))) +
    facet_grid(treatment, scales = "free", space = "free") +
    scale_fill_manual(values = graph_colors, aesthetics = c('color', 'fill'))
  if(merge == TRUE){p <- p + geom_bar(aes_string(color = classification_level, fill = classification_level), stat = 'identity', position = 'stack', size = 0.2)
  } else {p <- p + geom_bar(stat = "identity", position = "stack", size = 0.12, color = 'black')}
  if(relative_abundance == TRUE){p <- p + ylab('Relative Abundance')}

  return(p)
}


#' Create a ggplot object of the abundance table from a phyloseq object. Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and creates line graphs across samples.
#' @useDynLib phylosmith
#' @usage abundance_lines_ggplot(phyloseq_obj, classification_level, treatment, subset = NULL,
#' relative_abundance = FALSE, points = TRUE, colors = 'default')
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package (must contain \code{\link[phyloseq:sample_data]{sample_data()}}).
#' @param classification_level Column name or number in the \code{\link[phyloseq:tax_table]{tax_table}}.
#' @param treatment Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param subset If taxa not needed to be seen in all \code{treatment}, then can check only one particular treatment subset, this works for multiple treatment inputs.
#' @param relative_abundance If TRUE, transforms the abundance data into relative abundance by sample.
#' @param points if TRUE, will diplay the data-points.
#' @param colors Name of a color set from the \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R accepted colors.
#' @import ggplot2
#' @import RColorBrewer
#' @export

abundance_lines_ggplot <- function(phyloseq_obj, classification_level, treatment, subset = NULL, relative_abundance = FALSE, points = TRUE, colors = 'default'){

  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  if(is.numeric(classification_level)){classification_level <- colnames(phyloseq_obj@tax_table[,classification_level])}
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0, subset = subset)
  if(relative_abundance == TRUE){phyloseq_obj <- relative_abundance(phyloseq_obj)}
  treatment <- paste(treatment, collapse = '.')

  # graph_data <- tax_glom(phyloseq_obj, taxrank = classification_level)
  graph_data <- phyloseq(phyloseq_obj@otu_table, phyloseq_obj@tax_table[,classification_level], phyloseq_obj@sam_data[,treatment])
  graph_data <- data.table(psmelt(graph_data))

  colorCount = length(unique(phyloseq_obj@tax_table[,classification_level]))
  if(colors == 'default'){colors <- cbcolors}
  if(any(!(colors %in% grDevices::colors()))){
    if(any(colors %in% rownames(RColorBrewer::brewer.pal.info))){
      getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(min(c(colorCount, RColorBrewer::brewer.pal.info[rownames(RColorBrewer::brewer.pal.info) == colors, 1])), colors))
    } else { getPalette <- grDevices::colorRampPalette(colors)}
  } else { getPalette <- grDevices::colorRampPalette(colors)}
  graph_colors = getPalette(colorCount)

  p <- ggplot(graph_data, aes_string(x = 'Sample', y = 'Abundance', group = classification_level)) +
    geom_line(aes_string(color=classification_level))+
    theme_bw() +
    theme(axis.text.x = element_text(angle = -35, hjust = 0)) +
    guides(colour = guide_legend(ncol = ceiling(length(levels(graph_data[[classification_level]]))/30))) +
    facet_grid(treatment, scales = "free", space = "free") +
    scale_colour_manual(values = graph_colors)
  if(points == TRUE){p <- p + geom_point(aes_string(color=classification_level))}
  if(relative_abundance == TRUE){p <- p + ylab('Relative Abundance')}

  return(p)
}
