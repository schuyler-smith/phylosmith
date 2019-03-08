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
