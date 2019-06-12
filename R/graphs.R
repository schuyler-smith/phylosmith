#' Create a heatmap ggplot object of the abundance table from a
#' phyloseq object. Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and
#' creates heatmaps of the abundances across samples as a ggplot object.
#' The default color choice is the viridis palette, which is supposed to
#' be both aesthetic for normal and color-blind viewers.
#' @useDynLib phylosmith
#' @usage abundance_heatmap_ggplot(phyloseq_obj, classification = NULL,
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
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R-accepted
#' colors.
#' @import ggplot2
#' @export
#' @return ggplot-object
#' @examples abundance_heatmap_ggplot(soil_column, classification = 'Phylum',
#' treatment = c('Matrix', 'Treatment'), transformation = 'log')

abundance_heatmap_ggplot <- function(phyloseq_obj, classification = NULL,
    treatment, subset = NULL, transformation = 'none', colors = 'default'){
    if(!inherits(phyloseq_obj, "phyloseq")){
        stop("`phyloseq_obj` must be a
        phyloseq-class object", call. = FALSE)
    }
    classification <- check_numeric_classification(phyloseq_obj,
        classification)
    if(!(is.null(classification)) &
    is.null(access(phyloseq_obj, 'tax_table'))){
        stop("`phyloseq_obj` must contain
        tax_table() information if `classification` argument is used",
            call. = FALSE)
    }
    if(!(is.null(classification)) &
    any(!(classification %in% colnames(access(phyloseq_obj, 'tax_table'))))){
        stop("`classification` must be a column
        from the the tax_table()", call. = FALSE)
    }
    if(is.null(access(phyloseq_obj, 'sam_data'))){
        stop("`phyloseq_obj` must contain
        sample_data() information", call. = FALSE)
    }
    treatment <- check_numeric_treatment(phyloseq_obj, treatment)
    if(any(!(treatment %in% colnames(access(phyloseq_obj, 'sam_data'))))){
        stop("`treatment` must be at least one
        column name, or index, from the sample_data()", call. = FALSE)
    }
    if(!(transformation %in% c("none", "relative_abundance", "log", "log10",
        "log1p", "log2", "asn", "atanh", "boxcox", "exp", "identity", "logit",
        "probability", "probit", "reciprocal", "reverse", "sqrt"))){
        stop("argument given to `transformation`
        not able to be applied by this function, please see help files for
        list of acceptable values", call. = FALSE)
    }
    options(warn = -1)
    phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0,
        subset = subset)
    if(!(is.null(classification))){
        phyloseq_obj <- conglomerate_taxa(phyloseq_obj,
            classification, hierarchical = FALSE)}
    if(transformation == 'relative_abundance'){
        phyloseq_obj <- relative_abundance(phyloseq_obj)}
    treatment_name <- paste(treatment, collapse = sep)

    if(is.null(classification)){
        classification <- 'OTU'
        graph_data <- phyloseq(
            access(phyloseq_obj, 'otu_table'),
            access(phyloseq_obj, 'sam_data')[,treatment_name])
    } else {graph_data <- phyloseq(
            access(phyloseq_obj, 'otu_table'),
            access(phyloseq_obj, 'tax_table')[,classification],
            access(phyloseq_obj, 'sam_data')[,treatment_name])
    }
    graph_data <- melt_phyloseq(graph_data)
    graph_data[[classification]] <- factor(graph_data[[classification]],
        levels = rev(unique(graph_data[[classification]])))
    set(graph_data, which(is.na(graph_data[[classification]])),
        classification, 'Unclassified')
    graph_data[['Sample']] <- factor(graph_data[['Sample']],
        levels = rownames(access(phyloseq_obj, 'sam_data')))
    setkey(graph_data, 'Sample', 'Abundance')
    set(graph_data, which(graph_data[['Abundance']] == 0), 'Abundance', NA)

    g <- ggplot(graph_data, aes_string('Sample', classification, fill = 'Abundance')) +
        geom_tile(color = "white", size = 0.25) +
        facet_grid(treatment_name, scales = "free", space = "free")
    g <- g + theme_classic() +
      theme(
        axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
        axis.text.y = element_text(hjust = 0.95, size = 12),
        axis.title.x = element_text(size = 14, face = 'bold'),
        axis.title.y = element_text(size = 14, face = 'bold'),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 12),
        legend.spacing.x = unit(0.005, 'npc'),
        legend.background = element_rect(fill = (alpha = 0), color = 'black', size = 0.25),
        panel.background = element_rect(color = 'black', size = 1.4),
        strip.text.x = element_text(size = 14, face = 'bold'),
        strip.background = element_rect(colour = 'black', size = 1.4)
      ) +
      scale_x_discrete(expand = expand_scale(mult = 0, add = .53)) +
      if(transformation %in% c('none', 'relative_abundance')){
        if(colors == 'default'){
          scale_fill_viridis()
        } else {
          color_count <- 100
          graph_colors <- create_palette(color_count, colors)
          scale_fill_gradientn(colors = graph_colors)
        }
      } else {
        if(colors == 'default'){
          scale_fill_viridis(trans = transformation, name = transformation)
        } else {
          color_count <- 100
          graph_colors <- create_palette(color_count, colors)
          scale_fill_gradientn(colors = graph_colors, trans = transformation, name = transformation)
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
#' @usage abundance_lines_ggplot(phyloseq_obj, classification = NULL,
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
#' @param points if \code{FALSE}, will not diplay the data-points.
#' @param colors Name of a color set from the
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R-accepted
#' colors.
#' @import ggplot2
#' @export
#' @return ggplot-object
#' @examples abundance_lines_ggplot(soil_column, classification = 'Phylum',
#' treatment = c('Matrix', 'Treatment'), relative_abundance = TRUE)

abundance_lines_ggplot <- function(phyloseq_obj, classification = NULL,
    treatment, subset = NULL, relative_abundance = FALSE, points = TRUE,
    colors = 'default'){
    if(!inherits(phyloseq_obj, "phyloseq")){
        stop("`phyloseq_obj` must be a
        phyloseq-class object", call. = FALSE)
    }
    classification <- check_numeric_classification(phyloseq_obj,
        classification)
    if(!(is.null(classification)) &
        is.null(access(phyloseq_obj, 'tax_table'))){
        stop("`phyloseq_obj` must contain
        tax_table() information if `classification` argument is used",
        call. = FALSE)
    }
    if(!(is.null(classification)) & any(!(classification %in%
        colnames(access(phyloseq_obj, 'tax_table'))))){
        stop("`classification` must be a column from
        the the tax_table()", call. = FALSE)
    }
    if(is.null(access(phyloseq_obj, 'sam_data'))){
        stop("`phyloseq_obj` must contain
        sample_data() information", call. = FALSE)
    }
    treatment <- check_numeric_treatment(phyloseq_obj, treatment)
    if(any(!(treatment %in% colnames(access(phyloseq_obj, 'sam_data'))))){
        stop("`treatment` must be at least one
        column name, or index, from the sample_data()", call. = FALSE)
    }
    if(!(is.logical(relative_abundance))){
        stop("`relative_abundance` must be either
        `TRUE`, or `FALSE`", call. = FALSE)
    }
    if(!(is.logical(points))){
        stop("`points` must be either `TRUE`, or
        `FALSE`", call. = FALSE)
    }
    options(warn = -1)
    phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0,
        subset = subset)
    if(!(is.null(classification))){
        phyloseq_obj <- conglomerate_taxa(phyloseq_obj, classification,
        hierarchical = FALSE)}
    if(relative_abundance){phyloseq_obj <- relative_abundance(phyloseq_obj)}
    treatment_name <- paste(treatment, collapse = sep)

    if(is.null(classification)){
        classification <- 'OTU'
        graph_data <- phyloseq(
            access(phyloseq_obj, 'otu_table'),
            access(phyloseq_obj, 'sam_data')[,treatment_name])
    } else {
        graph_data <- phyloseq(
            access(phyloseq_obj, 'otu_table'),
            access(phyloseq_obj, 'tax_table')[,classification],
            access(phyloseq_obj, 'sam_data')[,treatment_name])
    }
    graph_data <- data.table(melt_phyloseq(graph_data))
    graph_data[[classification]] <- factor(graph_data[[classification]],
        levels = rev(unique(graph_data[[classification]])))
    set(graph_data, which(is.na(graph_data[[classification]])),
        classification, 'Unclassified')
    graph_data[['Sample']] <- factor(graph_data[['Sample']],
        levels = rownames(access(phyloseq_obj, 'sam_data')))

    color_count <- length(unique(graph_data[[classification]]))
    graph_colors <- create_palette(color_count, colors)

    g <- ggplot(graph_data, aes_string(x = 'Sample', y = 'Abundance', group = classification))
    if(points == TRUE){
      g <- g + geom_point(size = 1.5, aes_string(color = classification), show.legend = FALSE)
    }
    g <- g + geom_line(size = 1.2, aes_string(color = classification)) +
      facet_grid(treatment_name, scales = "free", space = "free") +
      scale_colour_manual(values = graph_colors) +
      guides(
        colour = guide_legend(ncol = ceiling(length(unique(graph_data[[classification]]))/25),
                              override.aes = list(size = 4))
      )
    if(relative_abundance == TRUE){g <- g + ylab('Relative Abundance')}
    g <- g + theme_bw() +
      theme(
        axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
        axis.text.y = element_text(hjust = 0.95, size = 12),
        axis.title.x = element_text(size = 14, face = 'bold'),
        axis.title.y = element_text(size = 14, face = 'bold'),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 12),
        legend.spacing.x = unit(0.005, 'npc'),
        legend.background = element_rect(fill = (alpha = 0)),
        panel.background = element_rect(color = 'black', size = 1.5, fill = 'white'),
        panel.spacing = unit(.015, 'npc'),
        strip.text.x = element_text(size = 14, face = 'bold', color = 'black'),
        strip.background = element_rect(colour = 'black', size = 1.4, fill = 'white')
      ) +
      scale_y_continuous(expand = expand_scale(mult = c(0, 0.003), add = c(0.0015, 0.001))) +
      scale_x_discrete(expand = expand_scale(mult = 0, add = c(0.3,0.5)))
    return(g)
}

#' Create a node network ggplot object of the co-occurrence from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and
#' creates a network from the co-occurrence. The co-occurrence can either be
#' input, or it will be calculated with the Spearman-rank correlation. Also,
#' the layout of the graph can be given as an argument as well for reprodusibility.
#' @useDynLib phylosmith
#' @usage network_ps(phyloseq_obj, classification = NULL,
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
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R accepted
#' colors.
#' @param cluster if \code{TRUE}, will use igraph's
#' \code{\link[igraph:cluster_fast_greedy]{cluster_fast_greedy}} method.
#' Alternatively, you may pass a vector of cluster assignments with order
#' corresponding to the order of the \code{taxa_names} in the
#' \code{phyloseq_obj}.
#' @param cluster_colors Name of a color set from the
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R accepted
#' colors to use for the clusters.
#' @param buffer Amount of space beyond the points to extend the cluster (
#' aesthetic argument).
#' @import data.table
#' @import igraph
#' @import ggraph
#' @import graphics
#' @importFrom sf st_as_sf st_buffer
#' @export
#' @return ggplot-object
#' @examples
#' #network_ps(soil_column, treatment = c('Matrix', 'Treatment'),
#' #subset = 'Soil Amended', co_occurrence_table = NULL, layout = NULL,
#' #classification = 'Phylum')

network_ps <- function(phyloseq_obj, classification = NULL,
    treatment = NULL, subset = NULL, co_occurrence_table = NULL,
    layout = NULL, nodes_of_interest = NULL, node_colors = 'default',
    cluster = FALSE, cluster_colors = 'default', buffer = 0.5){
    if(!inherits(phyloseq_obj, "phyloseq")){
        stop("`phyloseq_obj` must be a phyloseq-class
        object", call. = FALSE)
    }
    if(is.null(access(phyloseq_obj, 'sam_data'))){
        stop("`phyloseq_obj` must contain sample_data()
        information", call. = FALSE)
    }
    if(is.null(access(phyloseq_obj, 'tax_table'))){
        stop("`phyloseq_obj` must contain tax_table()
        information", call. = FALSE)
    }
    classification <- check_numeric_classification(phyloseq_obj,
        classification)
    if(!(is.null(classification)) &
        any(!(classification %in% colnames(access(phyloseq_obj,
            'tax_table'))))){
        stop("`classification` must be a column name, or
        index, from the tax_table()", call. = FALSE)
    }
    treatment <- check_numeric_treatment(phyloseq_obj, treatment)
    if(!(is.null(treatment)) &
        any(!(treatment %in% colnames(access(phyloseq_obj, 'sam_data'))))){
        stop("`treatment` must be at least one column
        name, or index, from the sample_data()", call. = FALSE)
    }
    if(!(is.null(co_occurrence_table)) &
        !(is.data.frame(co_occurrence_table))){
        stop("`co_occurrence_table` must be at data.frame
        object", call. = FALSE)
    }
    if(!(is.null(nodes_of_interest))){
        if(!(is.vector(nodes_of_interest))){
            stop("`nodes_of_interest` must be at vector of
            strings", call. = FALSE)
        }
    }
    if(!(is.numeric(buffer)) | !(buffer >= 0)){
        stop("`buffer` must be a numeric value >= 0", call.
        = FALSE)
    }
    node_classes <- c(sort(unique(access(phyloseq_obj,
        'tax_table')[,classification])), 'Unclassified')
    phyloseq_obj <- taxa_filter(phyloseq_obj, treatment,
        frequency = 0, subset = subset)
    treatment_name <- paste(treatment, collapse = sep)

    if(is.null(co_occurrence_table)){
        co_occurrence_table <- co_occurrence(phyloseq_obj,
            treatment)[rho >= 0.8 | rho <= -0.8]
    }
    co_occurrence_table <- co_occurrence_table[,c('OTU_1',
            'OTU_2','Treatment','rho','p')]
    if(!is.null(subset)){
        co_occurrence_table <- co_occurrence_table[co_occurrence_table[[
            'Treatment']] %like% subset]}
    colnames(co_occurrence_table)[colnames(co_occurrence_table)
        == 'rho'] <- 'weight'
    edge_colors <- c('pink1', 'gray22')[sapply(
      co_occurrence_table$weight,
      FUN = function(x){
        rep(as.numeric(as.logical(sign(x)+1)+1), 100)})]
    co_occurrence_table$weight <- abs(co_occurrence_table$weight)
    if(!(is.null(classification))){
        nodes <- data.table(as(access(phyloseq_obj, 'tax_table'), 'matrix'))
        nodes <- data.table('Node_Name' = rownames(access(phyloseq_obj,
            'tax_table')), nodes)
        set(nodes, which(is.na(nodes[[classification]])), classification,
            'Unclassified')
    } else {nodes <- data.table('Node_Name' = rownames(access(phyloseq_obj,
        'tax_table')))}
    nodes <- nodes[nodes[['Node_Name']] %in%
        c(as.character(co_occurrence_table$OTU_1),
            as.character(co_occurrence_table$OTU_2)),]

    net <- graph_from_data_frame(d = co_occurrence_table,
        vertices = nodes, directed = FALSE)
    net <- simplify(net, remove.multiple = FALSE, remove.loops = TRUE)
    if(is.null(layout)){
      layout <- create_layout(net, layout = 'igraph', algorithm = 'fr')
    }

    if(cluster == TRUE){cluster_table <- co_occurrence_table
        cluster_table[['weight']] <- abs(cluster_table[['weight']])
        cluster <- cluster_fast_greedy(simplify(
            graph_from_data_frame(d = cluster_table, vertices = nodes,
                directed = FALSE), remove.multiple = FALSE,
                remove.loops = TRUE))$membership}
    if(length(cluster) > 1 & length(cluster) != nrow(nodes)){
        stop("`cluster` must be either `TRUE`,`FALSE`, or
        a vector of memborship for each node", call. = FALSE)
    }
    if(length(cluster) > 1){communities <- data.table(layout[, c(1,2)])
        circles <- as(st_buffer(st_as_sf(communities, coords = c('x','y')),
            dist = buffer, nQuadSegs = 15), 'Spatial')
        circle_coords <- data.frame()
        for(i in seq(nrow(communities))){
            circle_coords <- rbind(circle_coords,
                circles@polygons[[i]]@Polygons[[1]]@coords)
        }; colnames(circle_coords) <- c('x', 'y')
        communities[, 'Community' := factor(cluster,
            levels = sort(unique(cluster)))]
        communities <- data.table(circle_coords,
            'Community' = unlist(lapply(communities$Community, rep,
                times = 62)))
        communities <- communities[, .SD[chull(.SD)], by = 'Community', ]
        hulls <- communities[, .SD[chull(.SD)], by = 'Community', ]
        community_count = length(unique(cluster))
        community_colors <- create_palette(community_count, cluster_colors)}

    if(!(is.null(classification))){node_colors <- create_palette(length(
        node_classes), node_colors)
    node_colors <- node_colors[node_classes %in%
        sort(unique(nodes[[classification]]))]}
    if(is.null(classification) & node_colors == 'default'){
        node_colors <- 'steelblue'}

    g <- ggraph(layout) + theme_graph() + coord_fixed()
    if(length(cluster) > 1){
      g <- g + geom_polygon(data = hulls, aes_string(x = 'x', y = 'y', alpha = 0.4, group = 'Community'), fill = community_colors[hulls$Community])}
    g <- g + geom_edge_link(color = edge_colors) +
      guides(colour = FALSE, alpha = FALSE, fill = guide_legend(ncol = ceiling(length(levels(layout[[classification]]))/25)))
    if(is.null(classification)){
      g <- g + geom_point(aes_string(x = 'x', y = 'y', fill = classification), pch=21, color = 'black',fill = node_colors, size=5)
    } else {
      g <- g + geom_point(aes_string(x = 'x', y = 'y', fill = classification), pch=21, color = 'black', size=5) +
        scale_fill_manual(values = node_colors)}
    if(!is.null(nodes_of_interest)){
      coi <- subset(layout, apply(layout, 1, function(class){any(class %in% nodes_of_interest)}))
      coi <- subset(layout, apply(layout, 1, function(class){any(class %in% nodes_of_interest)}))
      g <- g + ggrepel::geom_label_repel(data = coi, aes_string(x = 'x', y = 'y', fill = classification),
                                         label = unname(apply(coi, 1, function(class){class[which(class %in% nodes_of_interest)]})),
                                         box.padding = unit(0.8, "lines"), point.padding = unit(0.1, "lines"), size = 5, show.legend = FALSE)
    }
    g <- g + theme(
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12, face = 'bold'),
      legend.spacing.x = unit(0.005, 'npc')
    )
    return(g)
}

#' Create an layout_igraph object of the co-occurrence from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Create an layout_igraph object of the co-occurrence from a phyloseq object.
#' This can be input into the network_ps function, or used for other
#' network creating scripts. The purpose is to be able to create reproducible
#' and comparable graphics.
#' @useDynLib phylosmith
#' @usage network_layout_ps(phyloseq_obj, classification = NULL,
#' treatment = NULL, subset = NULL, co_occurrence_table = NULL,
#' algorithm = 'fr')
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
#' @param algorithm Supported \code{\link[igraph:layout_]{igraph::layout_()}} algorithm.
#' @import data.table
#' @import igraph
#' @import ggraph
#' @import graphics
#' @export
#' @return layout_igraph object
#' @examples
#' network_layout_ps(soil_column, treatment = c('Matrix', 'Treatment'),
#' subset = 'Soil Amended', algorithm = 'kk')

network_layout_ps <- function (phyloseq_obj, classification = NULL, treatment = NULL,
        subset = NULL, co_occurrence_table = NULL, algorithm = 'fr'){

  if(!inherits(phyloseq_obj, "phyloseq")){
    stop("`phyloseq_obj` must be a phyloseq-class
         object", call. = FALSE)
  }
  if(is.null(access(phyloseq_obj, 'sam_data'))){
    stop("`phyloseq_obj` must contain sample_data()
         information", call. = FALSE)
  }
  if(is.null(access(phyloseq_obj, 'tax_table'))){
    stop("`phyloseq_obj` must contain tax_table()
         information", call. = FALSE)
  }
  classification <- check_numeric_classification(phyloseq_obj,
                                                 classification)
  if(!(is.null(classification)) &
     any(!(classification %in% colnames(access(phyloseq_obj,
                                               'tax_table'))))){
    stop("`classification` must be a column name, or
         index, from the tax_table()", call. = FALSE)
  }
  treatment <- check_numeric_treatment(phyloseq_obj, treatment)
  if(!(is.null(treatment)) &
     any(!(treatment %in% colnames(access(phyloseq_obj, 'sam_data'))))){
    stop("`treatment` must be at least one column
         name, or index, from the sample_data()", call. = FALSE)
  }
  if(!(is.null(co_occurrence_table)) &
     !(is.data.frame(co_occurrence_table))){
    stop("`co_occurrence_table` must be at data.frame
         object", call. = FALSE)
  }

  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0, subset = subset)
  treatment_name <- paste(treatment, collapse = sep)

  if(is.null(co_occurrence_table)){
    co_occurrence_table <- co_occurrence(phyloseq_obj, treatment)[rho >= 0.8 | rho <= -0.8]
  }
  co_occurrence_table <- co_occurrence_table[,c('OTU_1', 'OTU_2', 'Treatment', 'rho', 'p')]
  if(!is.null(subset)){
    co_occurrence_table <- co_occurrence_table[co_occurrence_table[['Treatment']] %like% subset]}
  colnames(co_occurrence_table)[colnames(co_occurrence_table) == 'rho'] <- 'weight'
  co_occurrence_table$weight <- abs(co_occurrence_table$weight)
  if(!(is.null(classification))){
    nodes <- data.table(as(access(phyloseq_obj, 'tax_table'), 'matrix'))
    nodes <- data.table('Node_Name' = rownames(access(phyloseq_obj, 'tax_table')), nodes)
    set(nodes, which(is.na(nodes[[classification]])), classification, 'Unclassified')
  } else {
    nodes <- data.table('Node_Name' = rownames(access(phyloseq_obj, 'tax_table')))
  }
  nodes <- nodes[nodes[['Node_Name']] %in%
                   c(as.character(co_occurrence_table$OTU_1),
                     as.character(co_occurrence_table$OTU_2)),]

  net <- graph_from_data_frame(d = co_occurrence_table, vertices = nodes, directed = FALSE)
  net <- simplify(net, remove.multiple = FALSE, remove.loops = TRUE)
  layout <- create_layout(net, layout = 'igraph', algorithm = algorithm)

  return(layout)
}

#' Create a ggplot object of the NMDS from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and
#' plots the NMDS of a treatment or set of treatments in space.
#' @useDynLib phylosmith
#' @usage nmds_phyloseq_ggplot(phyloseq_obj, treatment, circle = 0.95,
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
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R-accepted
#' colors.
#' @param verbose Whether or not to print the
#' \code{\link[vegan:metaMDS]{metaMDS}} stress convergence to console
#' (\code{TRUE}) or not (\code{FALSE}).
#' @import ggplot2
#' @importFrom vegan metaMDS scores
#' @export
#' @return ggplot-object
#' @examples nmds_phyloseq_ggplot(soil_column, c('Matrix', 'Treatment'),
#' circle = TRUE, verbose = FALSE)

nmds_phyloseq_ggplot <- function(phyloseq_obj, treatment, circle = 0.95,
    labels = NULL, colors = 'default', verbose = TRUE){
    if(!inherits(phyloseq_obj, "phyloseq")){
        stop("`phyloseq_obj` must be a phyloseq-class
        object", call. = FALSE)
    }
    if(is.null(access(phyloseq_obj, 'sam_data'))){
        stop("`phyloseq_obj` must contain
        sample_data() information", call. = FALSE)
    }
    treatment <- check_numeric_treatment(phyloseq_obj, treatment)
    if(any(!(treatment %in% colnames(access(phyloseq_obj, 'sam_data'))))){
        stop("`treatment` must be at least one column
        name, or index, from the sample_data()", call. = FALSE)
    }
    if(!(is.logical(circle) | is.numeric(circle) & circle > 0 & circle <= 1)){
        stop("`circle` must be either `TRUE`, `FALSE`, or a confidence interval",
             call. = FALSE)
    }
    labels <- check_numeric_treatment(phyloseq_obj, labels)
    if(!(is.null(labels)) & any(!(labels %in% colnames(access(phyloseq_obj,
        'sam_data'))))){
        stop("`labels` must be a column name, or
        index, from the sample_data()", call. = FALSE)
    }
    if(!(is.logical(verbose))){
        stop("`verbose` must be either `TRUE`, or
        `FALSE`", call. = FALSE)
    }
    options(warn = -1)
    phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0)
    treatment_name <- paste(treatment, collapse = sep)
    metadata <- as(access(phyloseq_obj, 'sam_data'), 'data.frame')
    color_count <- length(unique(metadata[[treatment_name]]))
    graph_colors <- create_palette(color_count, colors)

    MDS <- metaMDS(t(access(phyloseq_obj, 'otu_table')), autotransform = FALSE,
        distance = "bray", k = 3, trymax = 100, trace = verbose)
    NMDS1 <- data.table(scores(MDS))$NMDS1
    NMDS2 <- data.table(scores(MDS))$NMDS2
    ord <- data.table(NMDS1, NMDS2, metadata)
    ord <- subset(ord, !is.na(treatment_name))
    if(is.character(labels)){
        eval(parse(text = paste0('ord[, ', labels,
            ' := access(phyloseq_obj, "sam_data")[[labels]]]')))}

    g <- ggplot(data = ord, aes_string('NMDS1', 'NMDS2', group = treatment_name))
    if(circle){
      g <- g + stat_ellipse(geom = "polygon", type = "norm",
        size = 0.6, linetype = 1, alpha = 0.3, color = 'black',
        aes_string(fill = treatment_name), show.legend = FALSE) +
        scale_color_manual(values = graph_colors) + guides(color = FALSE)
    } else if(is.numeric(circle)){
      ellipse_df <- CI_ellipse(ggplot_build(g)$data[[1]], groups = 'group', level = circle)
      g <- g + geom_polygon(data = ellipse_df, aes(x = x, y = y, group = group), color = 'black',
          fill = graph_colors[ellipse_df$group], alpha = 0.3, size = 0.6, linetype = 1)
    }
    g <- g + geom_point(aes_string(fill = treatment_name), shape = 21, color = 'black', size = 5, alpha = 1.0) +
      scale_fill_manual(values = graph_colors)
    if(is.character(labels)){
      g <- g + geom_label(aes_string(label = labels,
                                     fill = 'Treatment'), label.padding = unit(0.35, "lines"),
                          label.r = unit(0.55, "lines") , show.legend = FALSE)
    }
    g <- g +theme_classic() +
      theme(
        aspect.ratio = 1,
        axis.line.x = element_line(colour = 'black', size = 1,
                                   linetype = 'solid'),
        axis.line.y = element_line(colour = 'black', size = 1,
                                   linetype = 'solid'),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face= "bold"),
        axis.title.y = element_text(size = 14, face= "bold"),
        legend.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 12),
        legend.spacing.x = unit(0.005, 'npc'),
        legend.background = element_rect(fill = (alpha = 0))
      ) + labs(x = 'NMDS Dimension 1', y = 'NMDS Dimension 2')
    return(g)
}

#' Create a ggplot object of the phylogenic barplots from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and
#' creates phylogenic barplots.
#' @useDynLib phylosmith
#' @usage phylogeny_profile_ggplot(phyloseq_obj, classification = NULL,
#' treatment, subset = NULL, merge = TRUE, relative_abundance = FALSE,
#' colors = 'default')
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
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R-accepted
#' colors.
#' @import ggplot2
#' @export
#' @return ggplot-object
#' @examples phylogeny_profile_ggplot(soil_column, classification = 'Phylum',
#' treatment = c('Matrix', 'Treatment'), merge = TRUE,
#' relative_abundance = TRUE)

phylogeny_profile_ggplot <- function(phyloseq_obj, classification = NULL,
    treatment, subset = NULL, merge = TRUE, relative_abundance = FALSE,
    colors = 'default'){
    if(!inherits(phyloseq_obj, "phyloseq")){
        stop("`phyloseq_obj` must be a phyloseq-class
        object", call. = FALSE)
    }
    if(is.null(access(phyloseq_obj, 'sam_data'))){
        stop("`phyloseq_obj` must contain
        sample_data() information", call. = FALSE)
    }
    classification <- check_numeric_classification(phyloseq_obj,
        classification)
    if(!(is.null(classification)) &
        is.null(access(phyloseq_obj, 'tax_table'))){
        stop("`phyloseq_obj` must contain tax_table()
        information if `classification` argument is used", call. = FALSE)
    }
    if(!(is.null(classification)) & !(classification %in% colnames(
        access(phyloseq_obj, 'tax_table')))){
        stop("`classification` must be a column from
        the the tax_table()", call. = FALSE)
    }
    treatment <- check_numeric_treatment(phyloseq_obj, treatment)
    if(any(!(treatment %in% colnames(access(phyloseq_obj, 'sam_data'))))){
        stop("`treatment` must be at least one column
        name, or index, from the sample_data()", call. = FALSE)
    }
    if(!(is.logical(merge))){
        stop("`merge` must be either `TRUE`, or
        `FALSE`", call. = FALSE)
    }
    if(!(is.logical(relative_abundance))){
        stop("`relative_abundance` must be either
        `TRUE`, or `FALSE`", call. = FALSE)
    }
    options(warn = -1)
    phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0,
        subset = subset)
    if(!(is.null(classification)) & merge){
        phyloseq_obj <- conglomerate_taxa(phyloseq_obj, classification,
            hierarchical = FALSE)}
    if(relative_abundance){phyloseq_obj <- relative_abundance(phyloseq_obj)}
    treatment_name <- paste(treatment, collapse = sep)

    if(is.null(classification)){
        classification <- 'OTU'
        graph_data <- phyloseq(
            access(phyloseq_obj, 'otu_table'),
            access(phyloseq_obj, 'sam_data'))
    } else {
        graph_data <- phyloseq(
            access(phyloseq_obj, 'otu_table'),
            access(phyloseq_obj, 'tax_table'),
            access(phyloseq_obj, 'sam_data'))
    }
    graph_data <- melt_phyloseq(graph_data)
    graph_data[[classification]] <- factor(graph_data[[classification]],
        levels = sort(unique(graph_data[[classification]])))
    set(graph_data, which(is.na(graph_data[[classification]])),
        classification, 'Unclassified')
    graph_data[['Sample']] <- factor(graph_data[['Sample']],
        levels = rownames(access(phyloseq_obj, 'sam_data')))

    color_count <- length(unique(graph_data[[classification]]))
    graph_colors <- create_palette(color_count, colors)

    g <- ggplot(graph_data, aes_string(x = "Sample", y = "Abundance", fill = classification))
    g <- g +
      guides(colour = guide_legend(
        ncol = ceiling(length(unique(graph_data[[classification]]))/25))) +
      scale_fill_manual(values = graph_colors, aesthetics = c('color', 'fill'))
    if(!(is.null(treatment))){
      g <- g + facet_grid(treatment_name, scales = "free", space = "free")
    }
    if(merge){
      g <- g + geom_bar(aes_string(fill = classification), color = 'black', stat = 'identity', position = 'stack', size = 0.2, width = 0.95)
    } else {
      g <- g + geom_bar(stat = "identity", position = "stack", size = 0.12, width = 0.95, color = 'black')}
    if(relative_abundance){g <- g + ylab('Relative Abundance')}

    g <- g + theme_classic() +
      theme(
        axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = 'bold'),
        axis.title.y = element_text(size = 14, face = 'bold'),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 12),
        legend.spacing.x = unit(0.005, 'npc'),
        panel.background = element_rect(color = 'black', size = 1.5, fill = 'black'),
        panel.spacing = unit(0.01, 'npc'),
        strip.text.x = element_text(size = 14, face = 'bold'),
        strip.background = element_rect(colour = 'black', size = 1.4)
      ) +
      scale_y_continuous(expand = expand_scale(mult = c(0.0037, 0.003), add = c(0, 0))) +
      scale_x_discrete(expand = expand_scale(mult = 0, add = 0.51))
    return(g)
}

#' Create a ggplot object of the abundance barplots from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and
#' creates barplots of taxa by treatment.
#' @useDynLib phylosmith
#' @usage taxa_abundance_bars_ggplot(phyloseq_obj, classification = NULL,
#' treatment, subset = NULL, transformation = 'none', colors = 'default')
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
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R-accepted
#' colors.
#' @import ggplot2
#' @export
#' @return ggplot-object
#' @examples taxa_abundance_bars_ggplot(
#' taxa_filter(soil_column, frequency = 0.8),
#' classification = 'Phylum', treatment = c('Matrix', 'Treatment'),
#' subset = 'Unamended', transformation = 'mean')

taxa_abundance_bars_ggplot <- function(phyloseq_obj, classification = NULL,
    treatment, subset = NULL, transformation = 'none', colors = 'default'){
    if(!inherits(phyloseq_obj, "phyloseq")){
        stop("`phyloseq_obj` must be a
        phyloseq-class object", call. = FALSE)
    }
    if(is.null(access(phyloseq_obj, 'sam_data'))){
        stop("`phyloseq_obj` must contain
        sample_data() information", call. = FALSE)
    }
    classification <- check_numeric_classification(phyloseq_obj,
        classification)
    if(!(is.null(classification)) &
        is.null(access(phyloseq_obj, 'tax_table'))){
        stop("`phyloseq_obj` must contain
        tax_table() information if `classification` argument is used",
        call. = FALSE)
    }
    if(!(is.null(classification)) & !(classification %in% colnames(
        access(phyloseq_obj, 'tax_table')))){
        stop("`classification` must be a column from
        the the tax_table()", call. = FALSE)
    }
    treatment <- check_numeric_treatment(phyloseq_obj, treatment)
    if(any(!(treatment %in% colnames(access(phyloseq_obj, 'sam_data'))))){
        stop("`treatment` must be at least one
        column name, or index, from the sample_data()", call. = FALSE)
    }
    if(!(transformation %in%
        c("none", "mean", "median", "sd", "log", "log10"))){
        stop("argument given to `transformation`
        not able to be applied by this function, please see help files for
        list of acceptable values", call. = FALSE)
    }
    options(warn = -1)
    phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0,
        subset = subset)
    if(!(is.null(classification))){
        phyloseq_obj <- conglomerate_taxa(phyloseq_obj,
            classification, hierarchical = FALSE)
    } else {classification <- 'OTU'}
    treatment_name <- paste(treatment, collapse = sep)

    graph_data <- melt_phyloseq(phyloseq_obj)
    graph_data[[classification]] <- factor(graph_data[[classification]],
        levels = unique(graph_data[[classification]]))
    if(transformation == 'none'){abundance <- 'Abundance'
    graph_data <- graph_data[, sum(Abundance),
        by = c(treatment_name, classification)][, setnames(.SD, 'V1',
            abundance, skip_absent = TRUE)]}
    if(transformation == 'mean'){abundance <- 'Mean_Abundance'
    graph_data <- graph_data[, mean(Abundance), by = c(treatment_name,
        classification)][, setnames(.SD, 'V1', abundance, skip_absent = TRUE)]}
    if(transformation == 'median'){
        abundance <- 'Median_Abundance'
        graph_data <- graph_data[, stats::median(Abundance),
        by = c(treatment_name, classification)][, setnames(.SD, 'V1',
            abundance, skip_absent = TRUE)]}
    if(transformation == 'sd'){abundance <- 'StdDev_Abundance'
    graph_data <- graph_data[, stats::sd(Abundance), by = c(treatment_name,
        classification)][, setnames(.SD, 'V1', abundance, skip_absent = TRUE)]}
    if(transformation == 'log'){abundance <- 'log_Abundance'
    graph_data <- graph_data[, log(Abundance), by = c(treatment_name,
        classification)][, setnames(.SD, 'V1', abundance, skip_absent = TRUE)]}
    if(transformation == 'log10'){abundance <- 'log10_Abundance'
    graph_data <- graph_data[, log10(Abundance), by = c(treatment_name,
        classification)][, setnames(.SD, 'V1', abundance, skip_absent = TRUE)]}
    set(graph_data, which(is.na(graph_data[[classification]])),
        classification, 'Unclassified')

    color_count <- length(unique(graph_data[[treatment_name]]))
    graph_colors <- create_palette(color_count, colors)

    graph_data[[classification]] <- factor(graph_data[[classification]], levels = sort(unique(as.character(graph_data[[classification]]))))

    g <- ggplot(graph_data, aes_string(x = classification, y = abundance, fill = treatment_name))
    g <- g +  geom_bar(stat = "identity", position = position_dodge2(padding = 3.5),
                       size = 0.2, color = 'black', alpha = 0.85, width=0.62) +
      guides(colour = guide_legend(
        ncol = ceiling(length(unique(graph_data[[classification]]))/25))) +
      scale_fill_manual(values = graph_colors,
                        aesthetics = c('color', 'fill'))
    g <- g + theme_light() +
      theme(
        axis.line.x = element_line(colour = 'black', size = 1,
                                   linetype = 'solid'),
        axis.line.y = element_line(colour = 'black', size = 1,
                                   linetype = 'solid'),
        axis.text.x = element_text(size = 12, vjust = 1, hjust = 1, angle = 30),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        legend.background = element_rect(fill = (alpha = 0)),
        legend.spacing.x = unit(0.005, 'npc')
      ) +
      scale_y_continuous(expand = expand_scale(mult = c(0.0025, 0.002)))
    return(g)
}

#' Create a ggplot object using t-SNE from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object to plot the t-SNE of a
#' treatment or set of treatments.
#' @useDynLib phylosmith
#' @usage tsne_phyloseq_ggplot(phyloseq_obj, treatment, perplexity = 10,
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
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R-accepted
#' colors.
#' @import ggplot2
#' @importFrom Rtsne Rtsne
#' @importFrom vegan vegdist
#' @seealso \code{\link[=Rtsne]{Rtsne}}
#' @export
#' @return ggplot-object
#' @examples tsne_phyloseq_ggplot(soil_column,
#' treatment = c('Matrix', 'Treatment'), perplexity = 8)

tsne_phyloseq_ggplot <- function (phyloseq_obj, treatment, perplexity = 10,
    circle = TRUE, labels = NULL, colors = 'default'){
    if(!inherits(phyloseq_obj, "phyloseq")){
        stop("`phyloseq_obj` must be a phyloseq-class
        object", call. = FALSE)
    }
    if(is.null(access(phyloseq_obj, 'sam_data'))){
        stop("`phyloseq_obj` must contain
        sample_data() information", call. = FALSE)
    }
    treatment <- check_numeric_treatment(phyloseq_obj, treatment)
    if(any(!(treatment %in% colnames(access(phyloseq_obj, 'sam_data'))))){
        stop("`treatment` must be at least one column
        name, or index, from the sample_data()", call. = FALSE)
    }
    if(!(is.numeric(perplexity)) | perplexity <= 1){
        stop("`perplexity` must be a numeric value
        greater than 1", call. = FALSE)
    }
    if(!(is.logical(circle) | is.numeric(circle) & circle > 0 & circle <= 1)){
      stop("`circle` must be either `TRUE`, `FALSE`, or a confidence interval",
           call. = FALSE)
    }
    labels <- check_numeric_treatment(phyloseq_obj, labels)
    if(!(is.null(labels)) & any(!(labels %in% colnames(
        access(phyloseq_obj, 'sam_data'))))){
        stop("`labels` must be a column name, or
        index, from the sample_data()", call. = FALSE)
    }
    options(warn = -1)
    phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0)
    treatment_name <- paste(treatment, collapse = sep)
    metadata <- as(access(phyloseq_obj, 'sam_data'), 'data.frame')
    color_count <- length(unique(metadata[[treatment_name]]))
    graph_colors <- create_palette(color_count, colors)

    tsne <- Rtsne(vegdist(t(access(phyloseq_obj, 'otu_table')),
        method = 'bray'), dims = 2, theta = 0.0, perplexity = perplexity)
    tSNE1 <- tsne$Y[,1]
    tSNE2 <- tsne$Y[,2]
    ord <- data.table(tSNE1, tSNE2, metadata)
    ord <- subset(ord, !is.na(treatment_name))
    if(is.character(labels)){
        eval(parse(text=paste0('ord[, ', labels,
            ' := access(phyloseq_obj, "sam_data")[[labels]]]')))
    }

    g <- ggplot(data = ord, aes_string('tSNE1', 'tSNE2', group = treatment_name))
    if(circle){
      g <- g + stat_ellipse(geom = "polygon", type = "norm",
                            size = 0.6, linetype = 1, alpha = 0.3, color = 'black',
                            aes_string(fill = treatment_name), show.legend = FALSE) +
        scale_color_manual(values = graph_colors) + guides(color = FALSE)
    } else if(is.numeric(circle)){
      ellipse_df <- CI_ellipse(ggplot_build(g)$data[[1]], groups = 'group', level = circle)
      g <- g + geom_polygon(data = ellipse_df, aes(x = x, y = y, group = group), color = 'black',
                            fill = graph_colors[ellipse_df$group], alpha = 0.3, size = 0.6, linetype = 1)
    }
    g <- g + geom_point(aes_string(fill = treatment_name), shape = 21, color = 'black', size = 7, alpha = 1.0) +
      scale_fill_manual(values = graph_colors)
    if(is.character(labels)){
      g <- g + geom_label(aes_string(label = labels, fill = treatment_name),
                          label.padding = unit(0.35, "lines"), label.r = unit(0.55, "lines"),
                          show.legend = FALSE)}
    g <- g + theme_classic() +
      theme(
        aspect.ratio = 1,
        axis.line.x = element_line(colour = 'black', size = 1, linetype = 'solid'),
        axis.line.y = element_line(colour = 'black', size = 1, linetype = 'solid'),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = 'bold'),
        axis.title.y = element_text(size = 14, face = 'bold'),
        legend.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 12),
        legend.spacing.x = unit(0.005, 'npc'),
        legend.background = element_rect(fill = (alpha = 0))
      ) + labs(x = 't-SNE Dimension 1', y = 't-SNE Dimension 2')
    return(g)
}


