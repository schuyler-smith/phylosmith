#' Create an object of the abundance barplots from a phyloseq object. Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and creates barplots of taxa by treatment.
#' @useDynLib phylosmith
#' @usage network_phyloseq(phyloseq_obj, treatment, subset = NULL, co_occurrence,
#' classification = 'none', colors = 'default')
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package (must contain \code{\link[phyloseq:sample_data]{sample_data()}}).
#' @param treatment Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param subset If taxa not needed to be seen in all \code{treatment}, then will subset to treatments containing this string.
#' @param co_occurrence Co_Occurence table of the \code{phyloseq_obj}, computed using \code{\link{co_occurrence}}
#' @param classification Column name or number in the \code{\link[phyloseq:tax_table]{tax_table}} for node colors.
#' @param colors Name of a color set from the \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R accepted colors.
#' @import igraph
#' @import ggraph
#' @importFrom sf st_as_sf st_buffer
#' @import data.table
#' @import graphics
#' @export

network_phyloseq <- function(phyloseq_obj, treatment, subset = NULL, co_occurrence, classification = 'none', colors = 'default'){
  options(warn = -1)
  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0, subset = subset)
  treatment_name <- paste(treatment, collapse = '.')

  links <- co_occurrence[co_occurrence[['Treatment']] %like% subset]
  links <- links[,c(2,3,1,4,5)]
  colnames(links)[colnames(links)=='rho'] <- 'weight'

  nodes <- data.table(as(phyloseq_obj@tax_table, 'matrix'))
  nodes[, 'Node_Name' := rownames(phyloseq_obj@tax_table)]
  nodes[, 'Abundance' := taxa_sums(relative_abundance(phyloseq_obj))]
  nodes <- nodes[nodes[['Node_Name']] %in% c(as.character(links$OTU_1), as.character(links$OTU_2)),]

  net <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
  net <- simplify(net, remove.multiple = F, remove.loops = T)
  clusters <- cluster_fast_greedy(net)
  layout <- create_layout(net, layout = 'igraph', algorithm = 'fr')

  communities <- data.table(layout[,1:2])
  circles <- as(st_buffer(st_as_sf(communities, coords = c('x','y')), dist = 0.5, nQuadSegs = 15), 'Spatial')
  circle_coords <- data.frame()
  for(i in 1:nrow(communities)){
    circle_coords <- rbind(circle_coords, circles@polygons[[i]]@Polygons[[1]]@coords)
  }; colnames(circle_coords) <- c('x', 'y')
  communities[, 'Community' := factor(clusters$membership, levels = sort(unique(clusters$membership)))]
  communities <- data.table(circle_coords, Community = unlist(lapply(communities$Community, rep, times = 62)))
  communities <- communities[, .SD[chull(.SD)], by = Community, ]

  hulls <- communities[, .SD[chull(.SD)], by = Community, ]

  node_count = length(unique(nodes[[classification]]))
  node_colors <- create_palette(node_count, colors)
  community_count = length(unique(clusters$membership))
  community_colors <- create_palette(community_count, 'Pastel2')

  g <- ggraph(layout) +
    geom_polygon(data = hulls, aes_string(x = 'x', y = 'y', fill = 'Community', alpha = 0.4), show.legend = FALSE) +
    geom_edge_link(color = 'grey70') +
    geom_node_point(aes_string(fill = classification), color = 'black', pch=21, size=5) +
    scale_fill_manual(values = c(community_colors, node_colors), aesthetics = c('color', 'fill'), labels = NULL)

  return(g)
}
