#' Create an igraph network object of the co-occurrence from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Create an igraph network object of the co-occurrence from a phyloseq object.
#' This can be input into the co_occurrence_network function, or used for other
#' network creating scripts. The purpose is to be able to create reproducible
#' and comparable graphics.
#' @useDynLib phylosmith
#' @usage network_ps(phyloseq_obj,
#' treatment = NULL, subset = NULL, co_occurrence_table = NULL)
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
#' @param co_occurrence_table Table of the co-occurrence of taxa/genes in the
#' \code{phyloseq_obj}, computed using \code{\link{co_occurrence}}. If no
#' table is given, it will be computed with the \code{phyloseq_obj}, using the
#' given \code{treatment} and \code{p} = 0.05.
#' @importFrom igraph graph_from_data_frame simplify cluster_fast_greedy E
#' @importFrom ggraph create_layout
#' @export
#' @return igraph network object
#' @examples
#' network_ps(soil_column, treatment = c('Matrix', 'Treatment'),
#' subset = 'Soil Amended')


network_ps <-
  function (phyloseq_obj,
            treatment = NULL,
            subset = NULL,
            co_occurrence_table = NULL) {
  if (!inherits(phyloseq_obj, "phyloseq")) {
    stop("`phyloseq_obj` must be a phyloseq-class
          object", call. = FALSE)
  }
  if (is.null(access(phyloseq_obj, 'sam_data'))) {
    stop("`phyloseq_obj` must contain sample_data()
          information",
         call. = FALSE)
  }
  if (is.null(access(phyloseq_obj, 'tax_table'))) {
    stop("`phyloseq_obj` must contain tax_table()
          information",
         call. = FALSE)
  }
  treatment <- check_index_treatment(phyloseq_obj, treatment)
  if (!(is.null(treatment)) &
      any(!(treatment %in% colnames(access(
        phyloseq_obj, 'sam_data'
      ))))) {
    stop(
      "`treatment` must be at least one column
          name, or index, from the sample_data()",
      call. = FALSE
    )
  }
  if (!(is.null(treatment)) &
      is.null(subset)) {
    stop(
      "if `treatment` is declared,
        a `subset` must also be declared.",
      call. = FALSE
    )
  }
  if (!(is.null(co_occurrence_table)) &
      !(is.data.frame(co_occurrence_table))) {
    stop("`co_occurrence_table` must be at data.frame
          object", call. = FALSE)
  }
  phyloseq_obj <- taxa_filter(phyloseq_obj,
                              treatment,
                              frequency = 0,
                              subset = subset)
  treatment_name <- paste(treatment, collapse = sep)

  if (is.null(co_occurrence_table)) {
    co_occurrence_table <- co_occurrence(phyloseq_obj,
                                         treatment, method = 'spearman')[rho >= 0.6 |
                                                                           rho <= -0.6]
  } else { if(!is.null(subset)){
    co_occurrence_table <-
      co_occurrence_table[co_occurrence_table[['Treatment']] %like% subset]
    }
  }
  if(is.null(co_occurrence_table[['Treatment']])){
    co_occurrence_table <- cbind(co_occurrence_table, Treatment = 'NA')
  }
  co_occurrence_table <- co_occurrence_table[, c('X',
                                                 'Y', 'Treatment', 'rho', 'p')]
  if (!is.null(access(phyloseq_obj, 'tax_table'))){
    nodes <- data.table(as(access(phyloseq_obj, 'tax_table'), 'matrix'))
    nodes <-
      data.table('Node_Name' = rownames(access(phyloseq_obj,
                                               'tax_table')), nodes)
  } else {
    nodes <- data.table('Node_Name' = rownames(access(phyloseq_obj,
                                                      'otu_table')))
  }
  nodes <- nodes[nodes[['Node_Name']] %in%
                   c(
                     as.character(co_occurrence_table$X),
                     as.character(co_occurrence_table$Y)
                   ),]
  cluster_table <- co_occurrence_table
  cluster_table[['weight']] <- abs(cluster_table[['rho']])
  clusters <- cluster_fast_greedy(simplify(
    graph_from_data_frame(
      d = cluster_table,
      vertices = nodes,
      directed = FALSE
    ),
    remove.multiple = TRUE,
    remove.loops = TRUE
  ))$membership
  cluster_sizes <- table(clusters)
  nodes <- nodes[clusters %in% names(cluster_sizes[cluster_sizes > 3])]
  co_occurrence_table <- co_occurrence_table[X %in% nodes$Node_Name &
                                               Y %in% nodes$Node_Name]
  edge_sign <- vapply(
    co_occurrence_table$rho,
    FUN = function(x) {
      as.numeric(as.logical(sign(x) + 1) + 1)
    }, numeric(1)
  )
  co_occurrence_table$rho <- abs(co_occurrence_table$rho)
  setnames(co_occurrence_table, 'rho', 'weight')
  co_occurrence_table$edge_sign <- edge_sign
  net <- graph_from_data_frame(d = co_occurrence_table,
                               vertices = nodes,
                               directed = FALSE)
  return(net)
}

#' Create an layout_igraph object of the co-occurrence from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Create an layout_igraph object of the co-occurrence from a phyloseq object.
#' This can be input into the co_occurrence_network function, or used for other
#' network creating scripts. The purpose is to be able to create reproducible
#' and comparable graphics.
#' @useDynLib phylosmith
#' @usage network_layout_ps(phyloseq_obj,
#' treatment = NULL, subset = NULL, co_occurrence_table = NULL,
#' algorithm = 'fr')
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
#' @param co_occurrence_table Table of the co-occurrence of taxa/genes in the
#' \code{phyloseq_obj}, computed using \code{\link{co_occurrence}}. If no
#' table is given, it will be computed with the \code{phyloseq_obj}, using the
#' given \code{treatment} and \code{p} = 0.05.
#' @param algorithm Supported \code{\link[igraph:layout_]{igraph::layout_()}} algorithm.
#' @importFrom igraph graph_from_data_frame simplify cluster_fast_greedy
#' @importFrom ggraph create_layout
#' @export
#' @return layout_igraph object
#' @examples
#' network_layout_ps(soil_column, treatment = c('Matrix', 'Treatment'),
#' subset = 'Soil Amended', algorithm = 'kk')

network_layout_ps <-
  function (phyloseq_obj,
            treatment = NULL,
            subset = NULL,
            co_occurrence_table = NULL,
            algorithm = 'fr') {
    net <- network_ps(phyloseq_obj,
               treatment,
               subset,
               co_occurrence_table)
    layout <-
      create_layout(net, layout = 'igraph', algorithm = algorithm)
    return(layout)
  }

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
#' cluster = FALSE, cluster_colors = 'default', buffer = 0.5)
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

co_occurrence_network <- function(phyloseq_obj,
                                  classification = NULL,
                                  treatment = NULL,
                                  subset = NULL,
                                  co_occurrence_table = NULL,
                                  layout = NULL,
                                  nodes_of_interest = NULL,
                                  node_colors = 'default',
                                  cluster = FALSE,
                                  cluster_colors = 'default',
                                  buffer = 0.5) {
  if (!(is.null(nodes_of_interest))) {
    if (!(is.vector(nodes_of_interest))) {
      stop("`nodes_of_interest` must be at vector of
              strings", call. = FALSE)
    }
  }
  if (!(is.numeric(buffer)) | !(buffer >= 0)) {
    stop("`buffer` must be a numeric value >= 0", call.
         = FALSE)
  }
  classification <- check_index_classification(phyloseq_obj,
                                                 classification)
  if (!(is.null(classification)) &
      any(!(classification %in% colnames(access(
        phyloseq_obj,
        'tax_table'
      ))))) {
    stop("`classification` must be a column name, or
          index, from the tax_table()",
         call. = FALSE)
  }
  net <- network_ps(phyloseq_obj,
                    treatment,
                    subset,
                    co_occurrence_table)
  if (length(cluster) > 1 & length(cluster) != length(igraph::V(net))) {
    stop(
      "`cluster` must be either `TRUE`,`FALSE`, or
        a vector of memborship for each node",
      call. = FALSE
    )
  }
  if (is.null(layout)) {
    layout <- create_layout(net, layout = 'igraph', algorithm = 'fr')
  }
  if (cluster == TRUE) {
    cluster <- cluster_fast_greedy(simplify(net,
      remove.multiple = TRUE,
      remove.loops = TRUE
    ))$membership
  }
  if (length(cluster) > 1) {
    communities <- data.table(layout[, c(1, 2)])
    circles <-
      as(st_buffer(
        st_as_sf(communities, coords = c('x', 'y')),
        dist = buffer,
        nQuadSegs = 15
      ),
      'Spatial')
    circle_coords <- data.frame()
    for (i in seq(nrow(communities))) {
      circle_coords <- rbind(circle_coords,
                             circles@polygons[[i]]@Polygons[[1]]@coords)
    }
    colnames(circle_coords) <- c('x', 'y')
    communities[, 'Community' := factor(cluster,
                                        levels = sort(unique(cluster)))]
    communities <- data.table(circle_coords,
                              'Community' = unlist(lapply(communities$Community, rep,
                                                          times = 62)))
    communities <-
      communities[, .SD[chull(.SD)], by = 'Community', ]
    hulls <- communities[, .SD[chull(.SD)], by = 'Community', ]
    community_count = length(unique(cluster))
    community_colors <-
      create_palette(community_count, cluster_colors)
  }

  if (!(is.null(classification))) {
    eval(parse(text = paste0(
      'igraph::V(net)$', classification, '[is.na(igraph::V(net)$', classification, ')] <- "Unclassified"'
    )))
    node_classes <- c(sort(unique(access(
      phyloseq_obj,
      'tax_table'
    )[, classification])), 'Unclassified')
    node_colors <- create_palette(length(node_classes), node_colors)
    node_colors <- node_colors[node_classes %in%
      eval(parse(text = paste0(
        'igraph::V(net)$',
        classification
      )))]
  }
  if (is.null(classification) & node_colors == 'default') {
    node_colors <- 'steelblue'
  }
  edge_colors <- c('pink1', 'gray22')[vapply(igraph::E(net)$edge_sign, rep, numeric(100), 100)]

  g <- ggraph(layout) + theme_graph() + coord_fixed()
  if (length(cluster) > 1) {
    g <-
      g + geom_polygon(
        data = hulls,
        aes_string(
          x = 'x',
          y = 'y',
          alpha = 0.4,
          group = 'Community'
        ),
        fill = community_colors[hulls$Community]
      )
  }
  g <- g + geom_edge_link(width = 0.8, color = edge_colors) +
    guides(colour = FALSE,
           alpha = FALSE,
           fill = guide_legend(ncol = ceiling(length(unique(
             layout[[classification]]
           )) / 25)), override.aes = list(size = 4))
  if (is.null(classification)) {
    g <-
      g + geom_point(
        aes_string(x = 'x', y = 'y', fill = classification),
        pch = 21,
        color = 'black',
        fill = node_colors,
        size = 5
      )
  } else {
    g <-
      g + geom_point(
        aes_string(x = 'x', y = 'y', fill = classification),
        pch = 21,
        color = 'black',
        size = 5
      ) +
      scale_fill_manual(values = node_colors)
  }
  if (!is.null(nodes_of_interest)) {
    coi <-
      subset(layout, apply(layout, 1, function(class) {
        any(class %in% nodes_of_interest)
      }))
    g <-
      g + ggrepel::geom_label_repel(
        data = coi,
        aes_string(x = 'x', y = 'y', fill = classification),
        label = unname(apply(coi, 1, function(class) {
          class[which(class %in% nodes_of_interest)]
        })),
        box.padding = unit(0.8, "lines"),
        point.padding = unit(0.1, "lines"),
        size = 5,
        show.legend = FALSE
      )
  }
  g <- g + theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10, face = 'bold'),
    legend.key.size = unit(4, "mm"),
    legend.spacing.x = unit(0.005, 'npc')
  )
  return(g)
}


#' Create a node network ggplot object of the correlation of taxa and sample variables
#' from a phyloseq object. Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and
#' creates a network from the co-occurrence. The co-occurrence can either be
#' input, or it will be calculated with the Spearman-rank correlation. Also,
#' the layout of the graph can be given as an argument as well for reproducibility.
#' @useDynLib phylosmith
#' @usage variable_correlation_network(phyloseq_obj, variables, classification = NULL,
#' treatment = NULL, subset = NULL, correlation_table = NULL, method = 'spearman',
#' rho_threshold = c(-0.01, 0.01), p_threshold = 0.05, colors = 'default',
#' negative_positive_colors = c('pink1', 'gray22'))
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param variables Numerical factors within the in the
#' \code{\link[phyloseq:sample_data]{sample_data}} to correlate with the
#' abundance data.
#' @param classification Column name as a string or number in the
#' \code{\link[phyloseq:tax_table]{tax_table}} for the factor to use for node
#' colors.
#' @param treatment Column name as a string or number in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any
#' samples that to not contain this factor. This can be a vector of multiple
#' factors to subset on.
#' @param correlation_table Table of the correlation of taxa/variables in the
#' \code{phyloseq_obj}, computed using \code{\link{variable_correlation}}. If no
#' table is given, it will be computed with the \code{phyloseq_obj}, using the
#' given \code{treatment} and \code{p} = 0.05.
#' @param method Which correlation method to calculate, "pearson", "spearman", "kendall.
#' @param rho_threshold Cutoffs to use to subset the `correlation_table` by correlation values.
#' @param p_threshold Cutoffs to use to subset the `correlation_table` by singnificance values.
#' @param colors Name of a color set from the #' \link[=RColorBrewer]{RColorBrewer}
#' package or a vector palette of R accepted colors.
#' @param negative_positive_colors colors to use for the edges to represent negative and
#' positive correlations.
#' @importFrom igraph graph_from_data_frame simplify cluster_fast_greedy
#' @importFrom ggraph ggraph create_layout theme_graph geom_edge_link
#' @export
#' @return ggplot-object
#' @examples variable_correlation_network(soil_column, 'Day', 'Class', c('Matrix','Treatment'),
#' 'Soil Amended')

variable_correlation_network <- function(
  phyloseq_obj,
  variables,
  classification = NULL,
  treatment = NULL,
  subset = NULL,
  correlation_table = NULL,
  method = 'spearman',
  rho_threshold = c(-0.01, 0.01),
  p_threshold = 0.05,
  colors = 'default',
  negative_positive_colors = c('pink1', 'gray22')
){
  if(!(is.null(treatment)) & is.null(subset)) {
    stop("`subset` must be declared if `treatment`
          is decalred", call. = FALSE)
  }
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment = treatment,
                              subset = subset)
  if (!(is.null(classification))) {
    phyloseq_obj <- conglomerate_taxa(phyloseq_obj, classification,
                                      hierarchical = FALSE)
  }
  if(is.null(correlation_table)){
    correlation_table <- variable_correlation(phyloseq_obj, treatment = treatment,
                                              subset = subset, variables = variables,
                                              classification = classification,
                                              method = method, cores = 1)
  }
  correlation_table <- correlation_table[p <= p_threshold]
  correlation_table <- correlation_table[rho <= rho_threshold[1] |
                                           rho >= rho_threshold[2]]
  if(nrow(correlation_table) == 0){
    stop("None of the supplied variables were sgnificantly correlated to the data.",
         call. = FALSE)
  }
  if(is.null(correlation_table[['Treatment']])){
    correlation_table <- cbind(correlation_table, Treatment = 'NA')}
  correlation_table <- correlation_table[, c('X', 'Y', 'Treatment', 'rho', 'p')]
  colnames(correlation_table)[colnames(correlation_table)
                              == 'rho'] <- 'weight'
  if (!is.null(access(phyloseq_obj, 'tax_table'))){
    nodes <- data.table(as(access(phyloseq_obj, 'tax_table'), 'matrix'))
    nodes <-
      data.table('Node_Name' = rownames(access(phyloseq_obj,
                                               'tax_table')), nodes)
  } else {
    nodes <- data.table('Node_Name' = rownames(access(phyloseq_obj,
                                                      'tax_table')))
  }
  nodes <- nodes[nodes[['Node_Name']] %in%
                   c(
                     as.character(correlation_table$X),
                     as.character(correlation_table$Y)
                   ),]
  vars <- matrix(rep(variables, ncol(nodes)), nrow = length(variables))
  colnames(vars) <- colnames(nodes)
  nodes <- rbind(nodes, vars)
  edge_sign <- vapply(
    correlation_table$weight,
    FUN = function(x) {
      as.numeric(as.logical(sign(x) + 1) + 1)
    }, numeric(1)
  )
  correlation_table$weight <- abs(correlation_table$weight)
  net <- graph_from_data_frame(d = correlation_table,
                               vertices = nodes,
                               directed = FALSE)
  igraph::E(net)$edge_sign <- edge_sign
  net <- simplify(net, remove.multiple = FALSE, remove.loops = TRUE)
  layout <- create_layout(net, layout = 'igraph', algorithm = 'fr')

  edge_colors <- negative_positive_colors[vapply(igraph::E(net)$edge_sign, rep, numeric(100), 100)]
  node_colors <- dcast.data.table(correlation_table, X ~ Y, value.var = 'weight')
  node_colors[is.na(node_colors)] <- 0
  node_colors <- apply(node_colors[,-1], 1,FUN=function(x) paste(rev(variables[variables %in% colnames(node_colors)])[x>0], collapse = '_'))
  vars <- factor(c(node_colors, variables),
                 levels = c(variables, unique(node_colors[!(node_colors %in% variables)])))
  node_colors <- create_palette(length(levels(vars)))[vars]
  node_sizes <- rep(5, length(names(igraph::V(net))))
  node_sizes[names(igraph::V(net)) %in% variables] <- 15
  node_sizes[names(igraph::V(net)) %in% correlation_table[['Y']]] <- 20

  variables_layout <- subset(layout, apply(layout, 1, function(class) {
    any(class %in% variables)
  }))

  g <- ggraph(layout) + theme_graph() + coord_fixed() +
    geom_edge_link(width = 0.8, color = edge_colors) +
    guides(colour = FALSE,
           alpha = FALSE,
           fill = guide_legend(ncol = ceiling(length(levels(
             layout[[classification]]
           )) / 25)), override.aes = list(size = 4)) +
    geom_point(
      aes_string(x = 'x', y = 'y'),
      pch = 21,
      color = 'black',
      fill = node_colors,
      size = node_sizes
    ) + theme(
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10, face = 'bold'),
      legend.key.size = unit(4, "mm"),
      legend.spacing.x = unit(0.005, 'npc')
    )
  g <- g + geom_label(data = variables_layout, aes_string(x = 'x', y = 'y'),
                      label = unname(apply(variables_layout, 1, function(x) {
                        x[which(x %in% variables)]
                      }))[1,],
                      size = 3,
                      alpha = 0,
                      label.size = 0,
                      fontface = "bold",
                      show.legend = FALSE
  )
  return(g)
}
