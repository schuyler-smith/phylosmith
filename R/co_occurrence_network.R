#' Create a node network ggplot object of the co-occurrence from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and
#' creates a network from the co-occurrence. The co-occurrence can either be
#' input, or it will be calculated with the Spearman-rank correlation. Also,
#' the layout of the graph can be given as an argument as well for reproducibility.
#' @useDynLib phylosmith
#' @usage co_occurrence_network(phyloseq_obj, classification = NULL, 
#' treatment = NULL, subset = NULL, co_occurrence_table = NULL, layout = NULL,
#' nodes_of_interest = NULL, node_colors = 'default',
#' negative_positive_colors = c('tomato3', 'gray22'),
#' cluster = FALSE, cluster_colors = '#979aaa', buffer = 0.5)
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
#' @param co_occurrence_table Table of the co-occurrence of taxa/genes in the
#' \code{phyloseq_obj}, computed using \code{\link{co_occurrence}}. If no
#' table is given, it will be computed with the \code{phyloseq_obj}, using the
#' given \code{treatment} and \code{p} = 0.05.
#' @param layout (optional) an igraph layout of the network, for reproducibility.
#' Can be created with \code{\link{network_layout_ps}}.
#' @param nodes_of_interest A vector of names of classes within the
#' \code{classification} to be labeled.
#' @param node_colors Name of a color set from the
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palette of R accepted
#' colors.
#' @param negative_positive_colors colors to use for the edges to represent negative and
#' positive correlations.
#' @param cluster if \code{TRUE}, will use igraph's
#' \code{\link[igraph:cluster_fast_greedy]{cluster_fast_greedy}} method.
#' Alternatively, you may pass a vector of cluster assignments with order
#' corresponding to the order of the \code{taxa_names} in the
#' \code{phyloseq_obj}.
#' @param cluster_colors Name of a color set from the
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palette of R accepted
#' colors to use for the clusters.
#' @param buffer Amount of space beyond the points to extend the cluster (
#' aesthetic argument).
#' @importFrom igraph graph_from_data_frame simplify cluster_fast_greedy
#' @importFrom ggraph ggraph create_layout theme_graph geom_edge_link
#' @importFrom sf st_as_sf st_buffer
#' @importFrom grDevices chull
#' @export
#' @return ggplot-object
#' @examples
#' #co_occurrence_network(soil_column, treatment = c('Matrix', 'Treatment'),
#' #subset = 'Soil Amended', co_occurrence_table = NULL, layout = NULL,
#' #classification = 'Phylum')

co_occurrence_network <- function(
  phyloseq_obj,
  classification = NULL,
  treatment = NULL,
  subset = NULL,
  co_occurrence_table = NULL,
  layout = NULL,
  nodes_of_interest = NULL,
  node_colors = "default",
  negative_positive_colors = c("tomato3", "gray22"),
  cluster = FALSE,
  cluster_colors = "#979aaa",
  buffer = 0.5
) {
check_args(
    phyloseq_obj = phyloseq_obj,
    tax_table    = phyloseq_obj,
    treatment    = treatment,
    subset       = subset,
    co_occurrence_table = co_occurrence_table
)

  net <- network_ps(phyloseq_obj, treatment, subset, co_occurrence_table)
  if (length(cluster) > 1 & length(cluster) != length(igraph::V(net))) {
    stop(
"`cluster` must be either `TRUE`,`FALSE`, or a vector of 
memborship for each node", call. = FALSE)
  }
  if (is.null(layout)) {
    layout <- ggraph::create_layout(net, layout = "igraph", algorithm = "fr")
  }
  if (cluster == TRUE) {
    cluster <- igraph::cluster_fast_greedy(
      igraph::simplify(net, 
        remove.multiple = TRUE, remove.loops = TRUE))$membership
  }
  if (length(cluster) > 1) {
    communities <- data.table::data.table(layout[, c(1, 2)])
    circles <- as(sf::st_buffer(
        sf::st_as_sf(communities, coords = c("x", "y")),
        dist = buffer,
        nQuadSegs = 15
      ), "Spatial")
    circle_coords <- data.frame()
    for (i in seq(nrow(communities))) {
      circle_coords <- rbind(
        circle_coords, 
        circles@polygons[[i]]@Polygons[[1]]@coords)
    }
    colnames(circle_coords) <- c("x", "y")
    communities[, "Community" :=
      factor(cluster, levels = sort(unique(cluster)))]
    communities <- data.table::data.table(
      circle_coords, 
      Community = unlist(lapply(communities$Community, rep, 
        times = nrow(circles@polygons[[i]]@Polygons[[1]]@coords)))
    )
    communities <-
      communities[, .SD[chull(.SD)], by = "Community", ]
    hulls <- communities[, .SD[chull(.SD)], by = "Community", ]
    community_count = length(unique(cluster))
    community_colors <- create_palette(community_count, cluster_colors)
  }

  if (!is.null(classification)) {
    eval(parse(text = 
    paste0("igraph::V(net)$", classification, "[is.na(igraph::V(net)$", 
      classification, ")] <- 'Unclassified'"
    )))
    node_classes <- c(sort(unique(phyloseq_obj@tax_table[, classification])), 
      "Unclassified")
    node_colors <- create_palette(length(node_classes), node_colors)
    node_colors <- node_colors[node_classes %in%
      eval(parse(text = paste0("igraph::V(net)$", classification)))]
  } else {
    if ("default" %in% node_colors) node_colors <- "black"
  }
  # edge_colors <- c("tomato3", "gray22")[vapply(igraph::E(net)$edge_sign, rep, numeric(100), 100)]

  g <- ggraph::ggraph(layout) + coord_fixed()
  if (length(cluster) > 1) {
    g <- g + 
      geom_polygon(
        data = hulls,
        aes_string(x = "x", y = "y", group = "Community"),
        alpha = 0.2,
        color = "#414a4c",
        fill = community_colors[hulls$Community]
      )
  }
  g <- g + 
    ggraph::geom_edge_link(aes(color = Edge, width = Weight)) +
    ggraph::scale_edge_color_manual(values = negative_positive_colors) +
    ggraph::scale_edge_width_continuous(range = c(0.2, 2))
  if (!is.null(classification)) {
  g <- g + 
      geom_point(
        aes_string(x = "x", y = "y", fill = classification, 
          size = "`Mean Relative Abundance`"), pch = 21, color = "black")
  } else {
     g <- g + 
      geom_point(fill = node_colors,
        aes_string(x = "x", y = "y", 
          size = "`Mean Relative Abundance`"), pch = 21, color = "black")
  }
  g <- g +
    scale_fill_manual(values = node_colors, aesthetics = "fill") +
    scale_size_discrete(range = c(4, 12))
  if (!is.null(nodes_of_interest)) {
    coi <-
      subset(layout, apply(layout, 1, function(class) {
        any(grepl(paste(nodes_of_interest, collapse = "|"), class))
      }))
    g <- g + 
      ggrepel::geom_label_repel(
        data = coi,
        aes_string(x = "x", y = "y", fill = classification),
        box.padding = unit(2.0, "lines"),
        point.padding = unit(0.5, "lines"),
        size = 4,
        arrow = arrow(length = unit(0.2, "inches"),
                      ends="last", type="closed"),
        show.legend = FALSE,
        label = lapply(unname(apply(coi, 1, function(class) {
          class[which(grepl(paste(nodes_of_interest, collapse="|"), class))]
        })), `[[`, 1),
        max.overlaps = 1000
      )
  }
  n_cols <- tryCatch({
      ceiling(length(unique(layout[[classification]])) / 16)
    }, error = function(e) 1)
  g <- g + guides(
    colour = guide_legend(override.aes = list(size=10)),
    alpha = FALSE,
    fill = guide_legend(
      ncol = n_cols, 
      override.aes = list(size = 5)),
    edge_color = guide_legend(override.aes = list(edge_width = 2))) + 
    theme_schuy("graph")
    
  return(g)
}