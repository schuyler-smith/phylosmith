#' Create an object of the abundance barplots from a phyloseq object. Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and creates barplots of taxa by treatment.
#' @useDynLib phylosmith
#' @usage network_phyloseq(phyloseq_obj, treatment, subset = NULL, co_occurrence_table = NULL,
#' classification = 'none', node_colors = 'default', cluster = FALSE,
#' cluster_colors = 'default', buffer = 0.5)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package (must contain \code{\link[phyloseq:sample_data]{sample_data()}}).
#' @param treatment Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param subset If taxa not needed to be seen in all \code{treatment}, then will subset to treatments containing this string.
#' @param co_occurrence_table Co_Occurence table of the \code{phyloseq_obj}, computed using \code{\link{co_occurrence_table}}
#' @param classification Column name or number in the \code{\link[phyloseq:tax_table]{tax_table}} for node colors.
#' @param node_colors Name of a color set from the \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R accepted colors.
#' @param cluster if TRUE, will use igraph's \code{\link[igraph:cluster_fast_greedy]{cluster_fast_greedy}} method. Alternatively, can pass a vctor of cluster assignments with order corresponding to the order of the \code{taxa_names} in the \code{phyloseq_obj}.
#' @param cluster_colors Name of a color set from the \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R accepted colors to use for the cluster.
#' @param buffer Amount of space beyond the points to extend the cluster.
#' @import igraph
#' @import ggraph
#' @importFrom sf st_as_sf st_buffer
#' @import data.table
#' @import graphics
#' @export

network_phyloseq <- function(phyloseq_obj, treatment, subset = NULL, co_occurrence_table = NULL, classification = 'none', node_colors = 'default',
                             cluster = FALSE, cluster_colors = 'default', buffer = 0.5){
  options(warn = -1)
  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  if(is.numeric(classification)){classification <- colnames(phyloseq_obj@tax_table[,classification])}
  node_classes = sort(unique(phyloseq_obj@tax_table[,classification]))
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0, subset = subset)
  treatment_name <- paste(treatment, collapse = '.')

  if(is.null(co_occurrence_table)){co_occurrence_table <- co_occurrence(phyloseq_obj, treatment)}
  links <- co_occurrence_table[co_occurrence_table[['Treatment']] %like% subset]
  links <- links[,c(2,3,1,4,5)]
  colnames(links)[colnames(links)=='rho'] <- 'weight'

  nodes <- data.table(as(phyloseq_obj@tax_table, 'matrix'))
  nodes[, 'Node_Name' := rownames(phyloseq_obj@tax_table)]
  nodes[, 'Abundance' := taxa_sums(relative_abundance(phyloseq_obj))]
  nodes <- nodes[nodes[['Node_Name']] %in% c(as.character(links$OTU_1), as.character(links$OTU_2)),]

  net <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
  net <- simplify(net, remove.multiple = F, remove.loops = T)
  layout <- create_layout(net, layout = 'igraph', algorithm = 'fr')

  if(cluster == TRUE){cluster <- cluster_fast_greedy(net)$membership}
  if(length(cluster) > 1){communities <- data.table(layout[,1:2])
    circles <- as(st_buffer(st_as_sf(communities, coords = c('x','y')), dist = buffer, nQuadSegs = 15), 'Spatial')
    circle_coords <- data.frame()
    for(i in 1:nrow(communities)){
      circle_coords <- rbind(circle_coords, circles@polygons[[i]]@Polygons[[1]]@coords)
    }; colnames(circle_coords) <- c('x', 'y')
    communities[, 'Community' := factor(cluster, levels = sort(unique(cluster)))]
    communities <- data.table(circle_coords, Community = unlist(lapply(communities$Community, rep, times = 62)))
    communities <- communities[, .SD[chull(.SD)], by = Community, ]
    hulls <- communities[, .SD[chull(.SD)], by = Community, ]
    community_count = length(unique(cluster))
    community_colors <- create_palette(community_count, cluster_colors)}

  node_colors <- create_palette(length(node_classes), node_colors)
  node_colors <- node_colors[node_classes %in% sort(unique(nodes[[classification]]))]

  g <- ggraph(layout) + theme_graph() + coord_fixed()
  if(length(cluster) > 1){g <- g + geom_polygon(data = hulls, aes_string(x = 'x', y = 'y', alpha = 0.4, group = 'Community'), fill = community_colors[hulls$Community])}
  g <- g + geom_edge_link(color = 'grey70') +
    geom_point(aes_string(x = 'x', y = 'y', fill = classification), pch=21, color = 'black', size=5)  +
    scale_fill_manual(values = c(node_colors)) +
    labs(x="", y="") +
    guides(colour = FALSE, alpha = FALSE)

  return(g)
}
