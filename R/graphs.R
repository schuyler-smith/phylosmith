#' Create a ggplot object of the heatmap of the abundance table from a phyloseq object. Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and creates heatmaps of the abundances across samples.
#' @useDynLib phylosmith
#' @usage abundance_heatmap_ggplot(phyloseq_obj, classification = 'none', treatment, subset = NULL,
#' transformation = 'none', colors = 'default')
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with information about each sample, and it must contain \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each taxa/gene.
#' @param treatment Column name as a string or number in the \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any samples that to not contain this factor. This can be a vector of multiple factors to subset on.
#' @param classification Column name as a string or number in the \code{\link[phyloseq:tax_table]{tax_table}} for the factor.
#' @param transformation Transformation to be used on the data. "none", "relative_abundance", "log", "log10", "log1p", "log2", "asn", "atanh", "boxcox", "exp", "identity", "logit", "probability", "probit", "reciprocal", "reverse" and "sqrt"
#' @param colors Name of a color set from the \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R-accepted colors.
#' @import ggplot2
#' @export

abundance_heatmap_ggplot <- function(phyloseq_obj, classification = 'none', treatment, subset = NULL, transformation = 'none', colors = 'default'){
  options(warn = -1)
  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  if(is.numeric(classification)){classification <- colnames(phyloseq_obj@tax_table[,classification])}
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0, subset = subset)
  if(classification != 'none'){phyloseq_obj <- conglomerate_taxa(phyloseq_obj, classification, hierarchical = FALSE)}
  if(transformation == 'relative_abundance'){phyloseq_obj <- relative_abundance(phyloseq_obj)}
  treatment_name <- paste(treatment, collapse = sep)

  color_count <- 10
  if(colors == 'default'){colors <- 'YlOrRd'}
  graph_colors <- create_palette(color_count, colors)

  if(classification == 'none'){classification <- 'OTU'; graph_data <- phyloseq(phyloseq_obj@otu_table, phyloseq_obj@tax_table[,1], phyloseq_obj@sam_data[,treatment_name])
  } else {graph_data <- phyloseq(phyloseq_obj@otu_table, phyloseq_obj@tax_table[,classification], phyloseq_obj@sam_data[,treatment_name])}
  graph_data <- data.table(psmelt(graph_data))
  graph_data[[classification]] <- factor(graph_data[[classification]], levels = rev(unique(graph_data[[classification]])))
  graph_data[['Sample']] <- factor(graph_data[['Sample']], levels = sample_names(phyloseq_obj))

  g <- ggplot(graph_data, aes_string('Sample', classification, fill = 'Abundance')) +
    geom_tile(color = "white", size = 0.25) +
    facet_grid(treatment_name, scales = "free", space = "free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = -35, hjust = 0)) +
    if(transformation %in% c('none', 'relative_abundance')){scale_fill_gradientn(colors = graph_colors)
    } else {scale_fill_gradientn(colors = graph_colors, trans = transformation)}

  return(g)
}

#' Create a ggplot object of the abundance table from a phyloseq object. Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and creates line graphs across samples.
#' @useDynLib phylosmith
#' @usage abundance_lines_ggplot(phyloseq_obj, classification, treatment, subset = NULL,
#' relative_abundance = FALSE, points = TRUE, colors = 'default')
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with information about each sample, and it must contain \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each taxa/gene.
#' @param treatment Column name as a string or number in the \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any samples that to not contain this factor. This can be a vector of multiple factors to subset on.
#' @param classification Column name as a string or number in the \code{\link[phyloseq:tax_table]{tax_table}} for the factor to use for node colors.
#' @param relative_abundance If \code{TRUE}, transforms the abundance data into relative abundance by sample.
#' @param points if \code{FALSE}, will not diplay the data-points.
#' @param colors Name of a color set from the \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R-accepted colors.
#' @import ggplot2
#' @export

abundance_lines_ggplot <- function(phyloseq_obj, classification, treatment, subset = NULL, relative_abundance = FALSE, points = TRUE, colors = 'default'){
  options(warn = -1)
  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  if(is.numeric(classification)){classification <- colnames(phyloseq_obj@tax_table[,classification])}
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0, subset = subset)
  phyloseq_obj <- conglomerate_taxa(phyloseq_obj, classification, hierarchical = FALSE)
  if(relative_abundance){phyloseq_obj <- relative_abundance(phyloseq_obj)}
  treatment_name <- paste(treatment, collapse = sep)

  # graph_data <- tax_glom(phyloseq_obj, taxrank = classification)
  graph_data <- phyloseq(phyloseq_obj@otu_table, phyloseq_obj@tax_table[,classification], phyloseq_obj@sam_data[,treatment_name])
  graph_data <- data.table(psmelt(graph_data))
  graph_data[['Sample']] <- factor(graph_data[['Sample']], levels = sample_names(phyloseq_obj))

  color_count <- length(unique(graph_data[[classification]]))
  graph_colors <- create_palette(color_count, colors)

  g <- ggplot(graph_data, aes_string(x = 'Sample', y = 'Abundance', group = classification)) +
    geom_line(size = 1.0, aes_string(color=classification))+
    theme_bw() +
    theme(axis.text.x = element_text(angle = -35, hjust = 0)) +
    guides(colour = guide_legend(ncol = ceiling(length(unique(graph_data[[classification]]))/30))) +
    facet_grid(treatment_name, scales = "free", space = "free") +
    scale_colour_manual(values = graph_colors)
  if(points == TRUE){g <- g + geom_point(size = 1.8, aes_string(color = classification))}
  if(relative_abundance == TRUE){g <- g + ylab('Relative Abundance')}

  return(g)
}

#' Create an object of the abundance barplots from a phyloseq object. Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and creates barplots of taxa by treatment.
#' @useDynLib phylosmith
#' @usage network_phyloseq(phyloseq_obj, treatment = NULL, subset = NULL,
#' co_occurrence_table = NULL, classification = NULL, node_colors = 'default',
#' cluster = FALSE, cluster_colors = 'default', buffer = 0.5)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with information about each sample, and it must contain \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each taxa/gene.
#' @param treatment Column name as a string or number in the \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any samples that to not contain this factor. This can be a vector of multiple factors to subset on.
#' @param co_occurrence_table Table of the co-occurrence of taxa/genes in the \code{phyloseq_obj}, computed using \code{\link{co_occurrence}}. If no table is given, it will be computed with the \code{phyloseq_obj}, using the given \code{treatment} and \code{p} = 0.05.
#' @param classification Column name as a string or number in the \code{\link[phyloseq:tax_table]{tax_table}} for the factor to use for node colors.
#' @param node_colors Name of a color set from the \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R accepted colors.
#' @param cluster if \code{TRUE}, will use igraph's \code{\link[igraph:cluster_fast_greedy]{cluster_fast_greedy}} method. Alternatively, you may pass a vector of cluster assignments with order corresponding to the order of the \code{taxa_names} in the \code{phyloseq_obj}.
#' @param cluster_colors Name of a color set from the \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R accepted colors to use for the clusters.
#' @param buffer Amount of space beyond the points to extend the cluster (aesthetic argument).
#' @import data.table
#' @import igraph
#' @import ggraph
#' @import graphics
#' @importFrom sf st_as_sf st_buffer
#' @export

network_phyloseq <- function(phyloseq_obj, treatment = NULL, subset = NULL, co_occurrence_table = NULL, classification = NULL, node_colors = 'default',
                             cluster = FALSE, cluster_colors = 'default', buffer = 0.5){
  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  if(is.numeric(classification)){classification <- colnames(phyloseq_obj@tax_table[,classification])}
  if(!(is.null(classification))){node_classes = sort(unique(phyloseq_obj@tax_table[,classification]))}
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0, subset = subset)
  treatment_name <- paste(treatment, collapse = sep)

  if(is.null(co_occurrence_table)){co_occurrence_table <- co_occurrence(phyloseq_obj, treatment)}
  co_occurrence_table <- co_occurrence_table[,c('OTU_1','OTU_2','Treatment','rho','p')]
  if(!is.null(subset)){co_occurrence_table <- co_occurrence_table[co_occurrence_table[['Treatment']] %like% subset]}
  colnames(co_occurrence_table)[colnames(co_occurrence_table) == 'rho'] <- 'weight'

  nodes <- data.table(as(phyloseq_obj@tax_table, 'matrix'))
  nodes <- data.table('Node_Name' = rownames(phyloseq_obj@tax_table),
                      'Abundance' = taxa_sums(relative_abundance(phyloseq_obj)), nodes)
  nodes <- nodes[nodes[['Node_Name']] %in% c(as.character(co_occurrence_table$OTU_1), as.character(co_occurrence_table$OTU_2)),]

  net <- graph_from_data_frame(d=co_occurrence_table, vertices=nodes, directed=F)
  net <- simplify(net, remove.multiple = F, remove.loops = T)
  layout <- create_layout(net, layout = 'igraph', algorithm = 'fr')

  if(cluster == TRUE){cluster_table <- co_occurrence_table; cluster_table[['weight']] <- abs(cluster_table[['weight']])
    cluster <- cluster_fast_greedy(simplify(graph_from_data_frame(d=cluster_table, vertices=nodes, directed=F), remove.multiple = F, remove.loops = T))$membership}
  if(length(cluster) > 1){communities <- data.table(layout[,1:2])
    circles <- as(st_buffer(st_as_sf(communities, coords = c('x','y')), dist = buffer, nQuadSegs = 15), 'Spatial')
    circle_coords <- data.frame()
    for(i in 1:nrow(communities)){
      circle_coords <- rbind(circle_coords, circles@polygons[[i]]@Polygons[[1]]@coords)
    }; colnames(circle_coords) <- c('x', 'y')
    communities[, 'Community' := factor(cluster, levels = sort(unique(cluster)))]
    communities <- data.table(circle_coords, 'Community' = unlist(lapply(communities$Community, rep, times = 62)))
    communities <- communities[, .SD[chull(.SD)], by = 'Community', ]
    hulls <- communities[, .SD[chull(.SD)], by = 'Community', ]
    community_count = length(unique(cluster))
    community_colors <- create_palette(community_count, cluster_colors)}

  if(!(is.null(classification))){node_colors <- create_palette(length(node_classes), node_colors)
  node_colors <- node_colors[node_classes %in% sort(unique(nodes[[classification]]))]}
  if(is.null(classification) & node_colors == 'default'){node_colors <- 'steelblue'}

  g <- ggraph(layout) + theme_graph() + coord_fixed()
  if(length(cluster) > 1){g <- g + geom_polygon(data = hulls, aes_string(x = 'x', y = 'y', alpha = 0.4, group = 'Community'), fill = community_colors[hulls$Community])}
  g <- g + geom_edge_link(color = 'grey70') +
    guides(colour = FALSE, alpha = FALSE)
  if(is.null(classification)){g <- g + geom_point(aes_string(x = 'x', y = 'y', fill = classification), pch=21, color = 'black', fill = node_colors, size=5)
  } else {g <- g + geom_point(aes_string(x = 'x', y = 'y', fill = classification), pch=21, color = 'black', size=5) +
    scale_fill_manual(values = node_colors)}

  return(g)
}

#' Create a ggplot object of the NMDS from a phyloseq object. Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and plots the NMDS of a treatment or set of treatments.
#' @useDynLib phylosmith
#' @usage nmds_phyloseq_ggplot(phyloseq_obj, treatment, circle = TRUE, labels = FALSE,
#' colors = 'default', verbose = TRUE)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with information about each sample, and it must contain \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each taxa/gene.
#' @param treatment Column name as a string or number in the \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of multiple columns and they will be combined into a new column.
#' @param circle Add a \code{\link[ggplot2:stat_ellipse]{stat_ellipse}} around each of the \code{treatment} factors (\code{TRUE}).
#' @param labels Column name as a string or number in the \code{\link[phyloseq:sample_data]{sample_data}} to use to place labels of that factor instead of circle points.
#' @param colors Name of a color set from the \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R-accepted colors.
#' @param verbose Whether or not to print the \code{\link[vegan:metaMDS]{metaMDS}} stress convergence to console (\code{TRUE}) or not (\code{FALSE}).
#' @import ggplot2
#' @importFrom vegan metaMDS scores
#' @export

nmds_phyloseq_ggplot <- function(phyloseq_obj, treatment, circle = TRUE, labels = FALSE, colors = 'default', verbose = TRUE){
  options(warn = -1)
  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  if(is.numeric(labels)){labels <- colnames(phyloseq_obj@sam_data[,labels])}
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0)
  treatment_name <- paste(treatment, collapse = sep)

  MDS <- metaMDS(t(phyloseq_obj@otu_table), autotransform = FALSE, distance = "bray", k = 3, trymax = 100, trace = verbose)

  Treatment <- phyloseq_obj@sam_data[[treatment_name]]

  color_count <- length(unique(Treatment))
  graph_colors <- create_palette(color_count, colors)

  NMDS1 <- data.table(scores(MDS))$NMDS1
  NMDS2 <- data.table(scores(MDS))$NMDS2
  ord <- data.table(NMDS1,NMDS2,Treatment)
  ord <- subset(ord, !is.na(Treatment))
  if(is.character(labels)){eval(parse(text=paste0('ord[, ',labels,' := phyloseq_obj@sam_data[[labels]]]')))}

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
  if(is.character(labels)){g <- g + geom_label(aes_string(label = labels, fill = 'Treatment'), label.padding = unit(0.35, "lines"), label.r = unit(0.55, "lines") , show.legend = FALSE)}
  return(g)
}

#' Create a ggplot object of the phylogenic barplots from a phyloseq object. Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and creates phylogenic barplots.
#' @useDynLib phylosmith
#' @usage phylogeny_bars_ggplot(phyloseq_obj, classification, treatment, subset = NULL,
#' merge = TRUE, relative_abundance = TRUE, colors = 'default')
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with information about each sample, and it must contain \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each taxa/gene.
#' @param treatment Column name as a string or number in the \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any samples that to not contain this factor. This can be a vector of multiple factors to subset on.
#' @param classification Column name as a string or number in the \code{\link[phyloseq:tax_table]{tax_table}} for the factor to use for node colors.
#' @param merge if \code{FALSE}, will show separation of individuals within each \code{classification}.
#' @param relative_abundance If \code{TRUE}, transforms the abundance data into relative abundance by sample.
#' @param colors Name of a color set from the \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R-accepted colors.
#' @import ggplot2
#' @export

phylogeny_bars_ggplot <- function(phyloseq_obj, classification, treatment, subset = NULL, merge = TRUE, relative_abundance = TRUE, colors = 'default'){
  options(warn = -1)
  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  if(is.numeric(classification)){classification <- colnames(phyloseq_obj@tax_table[,classification])}
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0, subset = subset)
  if(merge){phyloseq_obj <- conglomerate_taxa(phyloseq_obj, classification, hierarchical = FALSE)}
  if(relative_abundance){phyloseq_obj <- relative_abundance(phyloseq_obj)}
  treatment_name <- paste(treatment, collapse = sep)

  # graph_data <- tax_glom(phyloseq_obj, taxrank = classification)
  graph_data <- phyloseq(phyloseq_obj@otu_table, phyloseq_obj@tax_table[,classification], phyloseq_obj@sam_data[,treatment_name])
  graph_data <- data.table(psmelt(graph_data))
  graph_data[['Sample']] <- factor(graph_data[['Sample']], levels = sample_names(phyloseq_obj))

  color_count <- length(unique(graph_data[[classification]]))
  graph_colors <- create_palette(color_count, colors)

  g <- ggplot(graph_data, aes_string(x = "Sample", y = "Abundance", fill = classification)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = -35, hjust = 0)) +
    guides(colour = guide_legend(ncol = ceiling(length(unique(graph_data[[classification]]))/30))) +
    facet_grid(treatment_name, scales = "free", space = "free") +
    scale_fill_manual(values = graph_colors, aesthetics = c('color', 'fill'))
  if(merge == TRUE){g <- g + geom_bar(aes_string(color = classification, fill = classification), stat = 'identity', position = 'stack', size = 0.2)
  } else {g <- g + geom_bar(stat = "identity", position = "stack", size = 0.12, color = 'black')}
  if(relative_abundance == TRUE){g <- g + ylab('Relative Abundance')}

  return(g)
}

#' Create a ggplot object of the abundance barplots from a phyloseq object. Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and creates barplots of taxa by treatment.
#' @useDynLib phylosmith
#' @usage taxa_abundance_bars_ggplot(phyloseq_obj, classification = 'none', treatment,
#' subset = NULL, transformation = 'none', colors = 'default')
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with information about each sample, and it must contain \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each taxa/gene.
#' @param treatment Column name as a string or number in the \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any samples that to not contain this factor. This can be a vector of multiple factors to subset on.
#' @param classification Column name as a string or number in the \code{\link[phyloseq:tax_table]{tax_table}} for the factor to use for node colors.
#' @param transformation Transformation to be used on the data. "none", "mean", "median", "sd", "log", "log10"
#' @param colors Name of a color set from the \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R-accepted colors.
#' @import ggplot2
#' @export

taxa_abundance_bars_ggplot <- function(phyloseq_obj, classification = 'none', treatment, subset = NULL, transformation = 'none', colors = 'default'){
  options(warn = -1)
  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  if(is.numeric(classification)){classification <- colnames(phyloseq_obj@tax_table[,classification])}
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0, subset = subset)
  if(classification != 'none'){phyloseq_obj <- conglomerate_taxa(phyloseq_obj, classification, hierarchical = FALSE)}
  treatment_name <- paste(treatment, collapse = sep)
  abundance <- 'Abundance'

  if(classification == 'none'){classification <- 'OTU'; graph_data <- phyloseq(phyloseq_obj@otu_table, phyloseq_obj@tax_table[,1], phyloseq_obj@sam_data[,treatment_name])
  } else {graph_data <- phyloseq(phyloseq_obj@otu_table, phyloseq_obj@tax_table[,classification], phyloseq_obj@sam_data[,treatment_name])}
  graph_data <- data.table(psmelt(graph_data))
  graph_data[[classification]] <- factor(graph_data[[classification]], levels = unique(graph_data[[classification]]))
  graph_data <- graph_data[, -'Sample']
  if(transformation == 'mean'){abundance <- 'Mean_Abundance'
    graph_data <- graph_data[, mean(Abundance), by = c(treatment_name, classification)][, setnames(.SD, 'V1', abundance, skip_absent = TRUE)]}
  if(transformation == 'median'){abundance <- 'Median_Abundance'
    graph_data <- graph_data[, stats::median(Abundance), by = c(treatment_name, classification)][, setnames(.SD, 'V1', abundance, skip_absent = TRUE)]}
  if(transformation == 'sd'){abundance <- 'StdDev_Abundance'
    graph_data <- graph_data[, stats::sd(Abundance), by = c(treatment_name, classification)][, setnames(.SD, 'V1', abundance, skip_absent = TRUE)]}
  if(transformation == 'log'){abundance <- 'log_Abundance'
    graph_data <- graph_data[, log(Abundance), by = c(treatment_name, classification)][, setnames(.SD, 'V1', abundance, skip_absent = TRUE)]}
  if(transformation == 'log10'){abundance <- 'log10_Abundance'
    graph_data <- graph_data[, log10(Abundance), by = c(treatment_name, classification)][, setnames(.SD, 'V1', abundance, skip_absent = TRUE)]}

  color_count <- length(unique(graph_data[[treatment_name]]))
  graph_colors <- create_palette(color_count, colors)

  g <- ggplot(graph_data, aes_string(x = classification, y = abundance, fill = treatment_name)) +
    geom_bar(stat = "identity", position = "dodge", size = 0.12, color = 'black') +
    theme_light() +
    theme(axis.text.x = element_text(angle = -35, hjust = 0)) +
    guides(colour = guide_legend(ncol = ceiling(length(unique(graph_data[[classification]]))/30))) +
    scale_fill_manual(values = graph_colors, aesthetics = c('color', 'fill'))

  return(g)
}

#' Create a ggplot object using t-SNE from a phyloseq object. Function from the phylosmith-package.
#'
#' Uses a \code{\link[phyloseq]{phyloseq-class}} object to plot the t-SNE of a treatment or set of treatments.
#' @useDynLib phylosmith
#' @usage tsne_phyloseq_ggplot(phyloseq_obj, treatment, perplexity = 10, circle = TRUE,
#' labels = FALSE, colors = 'default')
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with information about each sample, and it must contain \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each taxa/gene.
#' @param treatment Column name as a string or number in the \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of multiple columns and they will be combined into a new column.
#' @param perplexity similar to selecting the number of neighbors to consider in decision making (should not be bigger than 3 * perplexity < nrow(X) - 1, see \code{\link[=Rtsne]{Rtsne}} for interpretation)
#' @param circle Add a \code{\link[ggplot2:stat_ellipse]{stat_ellipse}} around each of the \code{treatment} factors (\code{TRUE}).
#' @param labels Column name as a string or number in the \code{\link[phyloseq:sample_data]{sample_data}} to use to place labels of that factor instead of circle points.
#' @param colors Name of a color set from the \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R-accepted colors.
#' @import ggplot2
#' @importFrom Rtsne Rtsne
#' @importFrom vegan vegdist
#' @seealso \code{\link[=Rtsne]{Rtsne}}
#' @export

tsne_phyloseq_ggplot <- function (phyloseq_obj, treatment, perplexity = 10, circle = TRUE, labels = FALSE, colors = 'default'){
  options(warn = -1)
  if (is.numeric(treatment)) {treatment <- colnames(phyloseq_obj@sam_data[, treatment])}
  if(is.numeric(labels)){labels <- colnames(phyloseq_obj@sam_data[,labels])}
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0)
  treatment_name <- paste(treatment, collapse = sep)

  tsne <- Rtsne(vegdist(t(phyloseq_obj@otu_table), method = 'bray'), dims = 2, theta = 0.0, perplexity = perplexity)

  Treatment <- phyloseq_obj@sam_data[[treatment_name]]

  color_count <- length(unique(Treatment))
  graph_colors <- create_palette(color_count, colors)

  tSNE1 <- tsne$Y[,1]
  tSNE2 <- tsne$Y[,2]
  ord <- data.table(tSNE1, tSNE2, Treatment)
  ord <- subset(ord, !is.na(Treatment))
  if(is.character(labels)){eval(parse(text=paste0('ord[, ',labels,' := phyloseq_obj@sam_data[[labels]]]')))}

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
  if(is.character(labels)){g <- g + geom_label(aes_string(label = labels, fill = 'Treatment'), label.padding = unit(0.35, "lines"), label.r = unit(0.55, "lines") , show.legend = FALSE)}
  return(g)
}

#' Internal function for creating color palettes for graphs. Function from the phylosmith-package.
#'
#' This function creates color palettes for graphs.
#' @useDynLib phylosmith
#' @usage create_palette(color_count, colors)
#' @param color_count Number of colors to choose for palette.
#' @param colors Name of a color set from the \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R-accepted colors.
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
