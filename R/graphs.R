#' Create a ggplot object of the heatmap of the abundance table from a
#' phyloseq object. Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and
#' creates heatmaps of the abundances across samples.
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
#' @examples abundance_heatmap_ggplot(soil_column, classification = 'phylum',
#' treatment = c('Matrix', 'Treatment'), transformation = 'log')

abundance_heatmap_ggplot <- function(phyloseq_obj, classification = NULL,
    treatment, subset = NULL, transformation = 'none', colors = 'default'){
    if(!inherits(phyloseq_obj, "phyloseq")){
        stop("abundance_heatmap_ggplot(): `phyloseq_obj` must be a
        phyloseq-class object", call. = FALSE)
    }
    classification <- check_numeric_classification(phyloseq_obj,
        classification)
    if(!(is.null(classification)) &
    is.null(access(phyloseq_obj, 'tax_table'))){
        stop("abundance_heatmap_ggplot(): `phyloseq_obj` must contain
        tax_table() information if `classification` argument is used",
            call. = FALSE)
    }
    if(!(is.null(classification)) &
    any(!(classification %in% colnames(access(phyloseq_obj, 'tax_table'))))){
        stop("abundance_heatmap_ggplot(): `classification` must be a column
        from the the tax_table()", call. = FALSE)
    }
    if(is.null(access(phyloseq_obj, 'sam_data'))){
        stop("abundance_heatmap_ggplot(): `phyloseq_obj` must contain
        sample_data() information", call. = FALSE)
    }
    treatment <- check_numeric_treatment(phyloseq_obj, treatment)
    if(any(!(treatment %in% colnames(access(phyloseq_obj, 'sam_data'))))){
        stop("abundance_heatmap_ggplot(): `treatment` must be at least one
        column name, or index, from the sample_data()", call. = FALSE)
    }
    if(!(transformation %in% c("none", "relative_abundance", "log", "log10",
        "log1p", "log2", "asn", "atanh", "boxcox", "exp", "identity", "logit",
        "probability", "probit", "reciprocal", "reverse", "sqrt"))){
        stop("abundance_heatmap_ggplot(): argument given to `transformation`
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

    color_count <- 10
    if(colors == 'default'){colors <- 'YlOrRd'}
    graph_colors <- create_palette(color_count, colors)

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

    g <- ggplot(graph_data, aes_string('Sample', classification,
        fill = 'Abundance')) +
        geom_tile(color = "white", size = 0.25) +
        facet_grid(treatment_name, scales = "free", space = "free") +
        theme_classic() +
        theme(
          axis.text.x = element_text(angle = -35, hjust = 0, size = 12),
          axis.text.y = element_text(hjust = 0.95, size = 12),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.title=element_blank(),
          legend.text=element_text(size = 16),
          legend.background = element_rect(fill = (alpha = 0))
        ) +
        if(transformation %in% c('none', 'relative_abundance')){
            scale_fill_gradientn(colors = graph_colors)
        } else {scale_fill_gradientn(colors = graph_colors,
            trans = transformation)}
    return(g)
}

#' Create a ggplot object of the abundance table from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and
#' creates line graphs across samples.
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
#' @examples abundance_lines_ggplot(soil_column, classification = 'phylum',
#' treatment = c('Matrix', 'Treatment'), relative_abundance = TRUE)

abundance_lines_ggplot <- function(phyloseq_obj, classification = NULL,
    treatment, subset = NULL, relative_abundance = FALSE, points = TRUE,
    colors = 'default'){
    if(!inherits(phyloseq_obj, "phyloseq")){
        stop("abundance_lines_ggplot(): `phyloseq_obj` must be a
        phyloseq-class object", call. = FALSE)
    }
    classification <- check_numeric_classification(phyloseq_obj,
        classification)
    if(!(is.null(classification)) &
        is.null(access(phyloseq_obj, 'tax_table'))){
        stop("abundance_lines_ggplot(): `phyloseq_obj` must contain
        tax_table() information if `classification` argument is used",
        call. = FALSE)
    }
    if(!(is.null(classification)) & any(!(classification %in%
        colnames(access(phyloseq_obj, 'tax_table'))))){
        stop("abundance_lines_ggplot(): `classification` must be a column from
        the the tax_table()", call. = FALSE)
    }
    if(is.null(access(phyloseq_obj, 'sam_data'))){
        stop("abundance_lines_ggplot(): `phyloseq_obj` must contain
        sample_data() information", call. = FALSE)
    }
    treatment <- check_numeric_treatment(phyloseq_obj, treatment)
    if(any(!(treatment %in% colnames(access(phyloseq_obj, 'sam_data'))))){
        stop("abundance_lines_ggplot(): `treatment` must be at least one
        column name, or index, from the sample_data()", call. = FALSE)
    }
    if(!(is.logical(relative_abundance))){
        stop("abundance_lines_ggplot(): `relative_abundance` must be either
        `TRUE`, or `FALSE`", call. = FALSE)
    }
    if(!(is.logical(points))){
        stop("abundance_lines_ggplot(): `points` must be either `TRUE`, or
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

    if(!(is.null(classification))){
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

    g <- ggplot(graph_data, aes_string(x = 'Sample', y = 'Abundance',
        group = classification))
    if(points == TRUE){
      g <- g + geom_point(size = 1.5, aes_string(color = classification))
    }
    g <- g + geom_line(size = 1.2, aes_string(color=classification))+
      theme_classic() +
      theme(
        axis.text.x = element_text(angle = -35, hjust = 0, size = 12),
        axis.text.y = element_text(hjust = 0.95, size = 12),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size = 16),
        legend.title=element_blank(),
        legend.text=element_text(size = 16),
        legend.background = element_rect(fill = (alpha = 0))
       ) +
      scale_y_continuous(expand = expand_scale(mult = c(0, .002))) +
      guides(colour = guide_legend(
          ncol = ceiling(length(unique(graph_data[[classification]]))/30))) +
      facet_grid(treatment_name, scales = "free", space = "free") +
      scale_colour_manual(values = graph_colors)
    if(relative_abundance == TRUE){g <- g + ylab('Relative Abundance')}
    return(g)
}

#' Create an object of the abundance barplots from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and
#' creates barplots of taxa by treatment.
#' @useDynLib phylosmith
#' @usage network_phyloseq(phyloseq_obj, classification = NULL,
#' treatment = NULL, subset = NULL, co_occurrence_table = NULL,
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
#' #network_phyloseq(soil_column, treatment = c('Matrix', 'Treatment'),
#' #subset = 'Soil_Manure', co_occurrence_table = NULL,
#' #classification = 'phylum')

network_phyloseq <- function(phyloseq_obj, classification = NULL,
    treatment = NULL, subset = NULL, co_occurrence_table = NULL,
    nodes_of_interest = NULL, node_colors = 'default',
    cluster = FALSE, cluster_colors = 'default', buffer = 0.5){
    if(!inherits(phyloseq_obj, "phyloseq")){
        stop("network_phyloseq(): `phyloseq_obj` must be a phyloseq-class
        object", call. = FALSE)
    }
    if(is.null(access(phyloseq_obj, 'sam_data'))){
        stop("network_phyloseq(): `phyloseq_obj` must contain sample_data()
        information", call. = FALSE)
    }
    if(is.null(access(phyloseq_obj, 'tax_table'))){
        stop("network_phyloseq(): `phyloseq_obj` must contain tax_table()
        information", call. = FALSE)
    }
    classification <- check_numeric_classification(phyloseq_obj,
        classification)
    if(!(is.null(classification)) &
        any(!(classification %in% colnames(access(phyloseq_obj,
            'tax_table'))))){
        stop("network_phyloseq(): `classification` must be a column name, or
        index, from the tax_table()", call. = FALSE)
    }
    treatment <- check_numeric_treatment(phyloseq_obj, treatment)
    if(!(is.null(treatment)) &
        any(!(treatment %in% colnames(access(phyloseq_obj, 'sam_data'))))){
        stop("network_phyloseq(): `treatment` must be at least one column
        name, or index, from the sample_data()", call. = FALSE)
    }
    if(!(is.null(co_occurrence_table)) &
        !(is.data.frame(co_occurrence_table))){
        stop("network_phyloseq(): `co_occurrence_table` must be at data.frame
        object", call. = FALSE)
    }
    if(!(is.null(nodes_of_interest))){
        if(!(is.vector(nodes_of_interest))){
            stop("network_phyloseq(): `nodes_of_interest` must be at vector of
            strings", call. = FALSE)
        }
    }
    if(!(is.numeric(buffer)) | !(buffer >= 0)){
        stop("network_phyloseq(): `buffer` must be a numeric value >= 0", call.
        = FALSE)
    }
    node_classes <- sort(unique(access(phyloseq_obj,
        'tax_table')[,classification]))
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
    layout <- create_layout(net, layout = 'igraph', algorithm = 'fr')

    if(cluster == TRUE){cluster_table <- co_occurrence_table
        cluster_table[['weight']] <- abs(cluster_table[['weight']])
        cluster <- cluster_fast_greedy(simplify(
            graph_from_data_frame(d = cluster_table, vertices = nodes,
                directed = FALSE), remove.multiple = FALSE,
                remove.loops = TRUE))$membership}
    if(length(cluster) > 1 & length(cluster) != nrow(nodes)){
        stop("network_phyloseq(): `cluster` must be either `TRUE`,`FALSE`, or
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

    g <- ggraph(layout) + theme_graph() + coord_fixed() +
      theme(
        legend.title=element_blank(),
        legend.text=element_text(size = 16)
      )
    if(length(cluster) > 1){
        g <- g + geom_polygon(data = hulls, aes_string(x = 'x', y = 'y',
            alpha = 0.4, group = 'Community'),
            fill = community_colors[hulls$Community])}
    g <- g + geom_edge_link(color = c('pink1', 'gray22')[sapply(
        E(attributes(layout)$graph)$weight, FUN = function(x){
            rep(as.numeric(as.logical(sign(x)+1)+1), 100)})]) +
        guides(colour = FALSE, alpha = FALSE,
            fill = guide_legend(ncol = ceiling(length(node_classes)/30)))
    if(is.null(classification)){
        g <- g + geom_point(aes_string(x = 'x', y = 'y',
            fill = classification), pch=21, color = 'black',
            fill = node_colors, size=5)
    } else {
        g <- g + geom_point(aes_string(x = 'x', y = 'y',
            fill = classification), pch=21, color = 'black', size=5) +
        scale_fill_manual(values = node_colors)}
    if(!is.null(nodes_of_interest)){
        coi <- subset(layout, apply(layout, 1,
            function(class){any(class %in% nodes_of_interest)}))
        coi <- subset(layout, apply(layout, 1,
            function(class){any(class %in% nodes_of_interest)}))
        g <- g + ggrepel::geom_label_repel(data = coi,
            aes_string(x = 'x', y = 'y', fill = classification),
            label = unname(apply(coi, 1, function(class){
                class[which(class %in% nodes_of_interest)]})),
                box.padding = unit(0.8, "lines"),
                point.padding = unit(0.1, "lines"), size = 5,
                show.legend = FALSE)
    }
    return(g)
}

#' Create a ggplot object of the NMDS from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and
#' plots the NMDS of a treatment or set of treatments.
#' @useDynLib phylosmith
#' @usage nmds_phyloseq_ggplot(phyloseq_obj, treatment, circle = TRUE,
#' labels = NULL, colors = 'default', verbose = TRUE)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param treatment Column name as a string or number in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @param circle Add a \code{\link[ggplot2:stat_ellipse]{stat_ellipse}} around
#' each of the \code{treatment} factors (\code{TRUE}).
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

nmds_phyloseq_ggplot <- function(phyloseq_obj, treatment, circle = TRUE,
    labels = NULL, colors = 'default', verbose = TRUE){
    if(!inherits(phyloseq_obj, "phyloseq")){
        stop("nmds_phyloseq_ggplot(): `phyloseq_obj` must be a phyloseq-class
        object", call. = FALSE)
    }
    if(is.null(access(phyloseq_obj, 'sam_data'))){
        stop("nmds_phyloseq_ggplot(): `phyloseq_obj` must contain
        sample_data() information", call. = FALSE)
    }
    treatment <- check_numeric_treatment(phyloseq_obj, treatment)
    if(any(!(treatment %in% colnames(access(phyloseq_obj, 'sam_data'))))){
        stop("nmds_phyloseq_ggplot(): `treatment` must be at least one column
        name, or index, from the sample_data()", call. = FALSE)
    }
    if(!(is.logical(circle))){
        stop("nmds_phyloseq_ggplot(): `circle` must be either `TRUE`, or
        `FALSE`", call. = FALSE)
    }
    labels <- check_numeric_treatment(phyloseq_obj, labels)
    if(!(is.null(labels)) & any(!(labels %in% colnames(access(phyloseq_obj,
        'sam_data'))))){
        stop("nmds_phyloseq_ggplot(): `labels` must be a column name, or
        index, from the sample_data()", call. = FALSE)
    }
    if(!(is.logical(verbose))){
        stop("nmds_phyloseq_ggplot(): `verbose` must be either `TRUE`, or
        `FALSE`", call. = FALSE)
    }
    options(warn = -1)
    phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0)
    treatment_name <- paste(treatment, collapse = sep)
    Treatment <- access(phyloseq_obj, 'sam_data')[[treatment_name]]
    color_count <- length(unique(Treatment))
    graph_colors <- create_palette(color_count, colors)

    MDS <- metaMDS(t(access(phyloseq_obj, 'otu_table')), autotransform = FALSE,
        distance = "bray", k = 3, trymax = 100, trace = verbose)
    NMDS1 <- data.table(scores(MDS))$NMDS1
    NMDS2 <- data.table(scores(MDS))$NMDS2
    ord <- data.table(NMDS1,NMDS2,Treatment)
    ord <- subset(ord, !is.na(Treatment))
    if(is.character(labels)){
        eval(parse(text=paste0('ord[, ', labels,
            ' := access(phyloseq_obj, "sam_data")[[labels]]]')))}

    g <- ggplot(data = ord, aes(NMDS1, NMDS2))
    if(circle == TRUE){
        g <- g + stat_ellipse(geom = "polygon", type = "norm",
        size = 0.6, linetype = 1, alpha = 0.1, color = 'black',
        aes(fill = Treatment), show.legend = FALSE) +
        scale_color_manual(values = graph_colors) + guides(color = FALSE)
    }
    g <- g + geom_point(aes(fill = Treatment), shape = 21, color = 'black',
        size = 5, alpha = 1.0) +
        scale_fill_manual(values = graph_colors) +
        theme_classic() +
        theme(
          aspect.ratio = 1,
          axis.line.x = element_line(colour = 'black', size = 1,
                                     linetype = 'solid'),
          axis.line.y = element_line(colour = 'black', size = 1,
                                     linetype = 'solid'),
          axis.text.x=element_text(size = 12),
          axis.text.y=element_text(size = 12),
          axis.title.x=element_text(size = 16, face= "bold"),
          axis.title.y=element_text(size = 16, face= "bold"),
          legend.title=element_blank(),
          legend.text=element_text(size = 16),
          legend.spacing.x = unit(0.2, 'cm'),
          legend.background = element_rect(fill = (alpha = 0))
        ) + labs(x = 'NMDS 1', y = 'NMDS 2')
    if(is.character(labels)){
        g <- g + geom_label(aes_string(label = labels,
            fill = 'Treatment'), label.padding = unit(0.35, "lines"),
            label.r = unit(0.55, "lines") , show.legend = FALSE)
    }
    return(g)
}

#' Create a ggplot object of the phylogenic barplots from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and
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
#' @examples phylogeny_profile_ggplot(soil_column, classification = 'phylum',
#' treatment = c('Matrix', 'Treatment'), merge = TRUE,
#' relative_abundance = TRUE)

phylogeny_profile_ggplot <- function(phyloseq_obj, classification = NULL,
    treatment, subset = NULL, merge = TRUE, relative_abundance = FALSE,
    colors = 'default'){
    if(!inherits(phyloseq_obj, "phyloseq")){
        stop("phylogeny_profile_ggplot(): `phyloseq_obj` must be a phyloseq-class
        object", call. = FALSE)
    }
    if(is.null(access(phyloseq_obj, 'sam_data'))){
        stop("phylogeny_profile_ggplot(): `phyloseq_obj` must contain
        sample_data() information", call. = FALSE)
    }
    classification <- check_numeric_classification(phyloseq_obj,
        classification)
    if(!(is.null(classification)) &
        is.null(access(phyloseq_obj, 'tax_table'))){
        stop("phylogeny_profile_ggplot(): `phyloseq_obj` must contain tax_table()
        information if `classification` argument is used", call. = FALSE)
    }
    if(!(is.null(classification)) & !(classification %in% colnames(
        access(phyloseq_obj, 'tax_table')))){
        stop("phylogeny_profile_ggplot(): `classification` must be a column from
        the the tax_table()", call. = FALSE)
    }
    treatment <- check_numeric_treatment(phyloseq_obj, treatment)
    if(any(!(treatment %in% colnames(access(phyloseq_obj, 'sam_data'))))){
        stop("phylogeny_profile_ggplot(): `treatment` must be at least one column
        name, or index, from the sample_data()", call. = FALSE)
    }
    if(!(is.logical(merge))){
        stop("phylogeny_profile_ggplot(): `merge` must be either `TRUE`, or
        `FALSE`", call. = FALSE)
    }
    if(!(is.logical(relative_abundance))){
        stop("phylogeny_profile_ggplot(): `relative_abundance` must be either
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

    g <- ggplot(graph_data, aes_string(x = "Sample", y = "Abundance",
        fill = classification)) +
        theme_classic() +
        theme(
          axis.text.x = element_text(angle = -35, hjust = 0, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 16, face = 'bold'),
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.spacing.x = unit(0.2, 'cm')
        ) +
        scale_y_continuous(expand = expand_scale(mult = c(0, 0.002))) +
        guides(colour = guide_legend(
            ncol = ceiling(length(unique(graph_data[[classification]]))/30))) +
        scale_fill_manual(values = graph_colors,
            aesthetics = c('color', 'fill'))
    if(!(is.null(treatment))){
        g <- g + facet_grid(treatment_name, scales = "free", space = "free")
    }
    if(merge){
        g <- g + geom_bar(aes_string(color = classification,
            fill = classification), stat = 'identity', position = 'stack',
        size = 0.2)
    } else {
        g <- g + geom_bar(stat = "identity", position = "stack",
            size = 0.12, color = 'black')}
    if(relative_abundance){g <- g + ylab('Relative Abundance')}
    return(g)
}

#' Create a ggplot object of the abundance barplots from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and
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
#' classification = 'phylum', treatment = c('Matrix', 'Treatment'),
#' subset = 'Control', transformation = 'mean')

taxa_abundance_bars_ggplot <- function(phyloseq_obj, classification = NULL,
    treatment, subset = NULL, transformation = 'none', colors = 'default'){
    if(!inherits(phyloseq_obj, "phyloseq")){
        stop("taxa_abundance_bars_ggplot(): `phyloseq_obj` must be a
        phyloseq-class object", call. = FALSE)
    }
    if(is.null(access(phyloseq_obj, 'sam_data'))){
        stop("taxa_abundance_bars_ggplot(): `phyloseq_obj` must contain
        sample_data() information", call. = FALSE)
    }
    classification <- check_numeric_classification(phyloseq_obj,
        classification)
    if(!(is.null(classification)) &
        is.null(access(phyloseq_obj, 'tax_table'))){
        stop("taxa_abundance_bars_ggplot(): `phyloseq_obj` must contain
        tax_table() information if `classification` argument is used",
        call. = FALSE)
    }
    if(!(is.null(classification)) & !(classification %in% colnames(
        access(phyloseq_obj, 'tax_table')))){
        stop("phylogeny_profile_ggplot(): `classification` must be a column from
        the the tax_table()", call. = FALSE)
    }
    treatment <- check_numeric_treatment(phyloseq_obj, treatment)
    if(any(!(treatment %in% colnames(access(phyloseq_obj, 'sam_data'))))){
        stop("taxa_abundance_bars_ggplot(): `treatment` must be at least one
        column name, or index, from the sample_data()", call. = FALSE)
    }
    if(!(transformation %in%
        c("none", "mean", "median", "sd", "log", "log10"))){
        stop("taxa_abundance_bars_ggplot(): argument given to `transformation`
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

    g <- ggplot(graph_data, aes_string(x = classification, y = abundance,
        fill = treatment_name)) +
        geom_bar(stat = "identity", position = position_dodge2(padding = 3.5),
                 size = 0.2, color = 'black', alpha = 0.85, width=0.62) +
        theme_light() +
        theme(
          axis.line.x = element_line(colour = 'black', size = 1,
                                     linetype = 'solid'),
          axis.line.y = element_line(colour = 'black', size = 1,
                                     linetype = 'solid'),
          axis.text.x=element_text(size = 12, vjust = 0.7, hjust = 0, angle = -35),
          axis.text.y=element_text(size = 12),
          axis.title.x=element_text(size = 16, face= "bold"),
          axis.title.y=element_text(size = 16, face= "bold"),
          legend.text=element_text(size = 16),
          legend.title = element_text(size = 16, face = "bold"),
          legend.background = element_rect(fill = (alpha = 0))
        ) + scale_y_continuous(expand = expand_scale(mult = c(0, 0.002))) +
        guides(colour = guide_legend(
            ncol = ceiling(length(unique(graph_data[[classification]]))/30))) +
        scale_fill_manual(values = graph_colors,
            aesthetics = c('color', 'fill'))
    return(g)
}

#' Create a ggplot object using t-SNE from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Uses a \code{\link[phyloseq]{phyloseq-class}} object to plot the t-SNE of a
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
#' @param circle Add a \code{\link[ggplot2:stat_ellipse]{stat_ellipse}} around
#' each of the \code{treatment} factors (\code{TRUE}).
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
#' treatment = c('Matrix', 'Treatment'), perplexity = 10)


tsne_phyloseq_ggplot <- function (phyloseq_obj, treatment, perplexity = 10,
    circle = TRUE, labels = NULL, colors = 'default'){
    if(!inherits(phyloseq_obj, "phyloseq")){
        stop("tsne_phyloseq_ggplot(): `phyloseq_obj` must be a phyloseq-class
        object", call. = FALSE)
    }
    if(is.null(access(phyloseq_obj, 'sam_data'))){
        stop("tsne_phyloseq_ggplot(): `phyloseq_obj` must contain
        sample_data() information", call. = FALSE)
    }
    treatment <- check_numeric_treatment(phyloseq_obj, treatment)
    if(any(!(treatment %in% colnames(access(phyloseq_obj, 'sam_data'))))){
        stop("tsne_phyloseq_ggplot(): `treatment` must be at least one column
        name, or index, from the sample_data()", call. = FALSE)
    }
    if(!(is.numeric(perplexity)) | perplexity <= 1){
        stop("tsne_phyloseq_ggplot(): `perplexity` must be a numeric value
        greater than 1", call. = FALSE)
    }
    if(!(is.logical(circle))){
        stop("tsne_phyloseq_ggplot(): `circle` must be either `TRUE`, or
        `FALSE`", call. = FALSE)
    }
    labels <- check_numeric_treatment(phyloseq_obj, labels)
    if(!(is.null(labels)) & any(!(labels %in% colnames(
        access(phyloseq_obj, 'sam_data'))))){
        stop("tsne_phyloseq_ggplot(): `labels` must be a column name, or
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

    g <- ggplot(data = ord, aes(tSNE1, tSNE2))
    if(circle == TRUE){g <- g + stat_ellipse(geom = "polygon", type = "norm",
        size = 0.6, linetype = 1, alpha = 0.1, color = 'black',
        aes_string(fill = treatment_name), show.legend = FALSE) +
        scale_color_manual(values = graph_colors) + guides(color = FALSE)}
    g <- g + geom_point(aes_string(fill = treatment_name), shape = 21, color = 'black',
            size = 7, alpha = 1.0) +
        scale_fill_manual(values = graph_colors) +
        theme_classic() +
        theme(
            aspect.ratio = 1,
            axis.line.x = element_line(colour = 'black', size = 1,
                linetype = 'solid'),
            axis.line.y = element_line(colour = 'black', size = 1,
                linetype = 'solid'),
            axis.text.x=element_text(size = 12),
            axis.text.y=element_text(size = 12),
            axis.title.x=element_text(size = 16, face= "bold"),
            axis.title.y=element_text(size = 16, face= "bold"),
            legend.title=element_blank(),
            legend.text=element_text(size = 16),
            legend.spacing.x = unit(0.2, 'cm'),
            legend.background = element_rect(fill = (alpha = 0))
        ) + labs(x = 't-SNE 1', y = 't-SNE 2')
    if(is.character(labels)){
        g <- g + geom_label(aes_string(label = labels, fill = 'Treatment'),
            label.padding = unit(0.35, "lines"), label.r = unit(0.55, "lines"),
            show.legend = FALSE)}
    return(g)
}


