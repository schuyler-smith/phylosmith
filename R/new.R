#' Computes the correlation of numerical variables with taxa
#' Function from the phylosmith-package.
#'
#' Computes the correlation of numerical variables with taxa
#' @useDynLib phylosmith
#' @usage variable_correlation(phyloseq_obj, variables, treatment = NULL,
#'  subset = NULL, classification = NULL, method = 'spearman', cores = 1)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object.
#' @param variables Numericla factors within the in the
#' \code{\link[phyloseq:sample_data]{sample_data}} to correlate with the
#' abundance data.
#' @param treatment Column name as a \code{string} or \code{numeric} in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any
#' samples that to not contain this factor. This can be a vector of multiple
#' factors to subset on.
#' @param classification Column name as a \code{string} or \code{numeric} in
#' the \code{\link[phyloseq:tax_table]{tax_table}} for the factor to
#' conglomerate by.
#' @param method Which correlation method to calculate, "pearson", "spearman".
#' @param cores \code{numeric} Number of CPU cores to use for the pair-wise
#' permutations. Default (0) uses max cores available. Parallelization not
#' available for systems running MacOS without openMP configuration.
#' @importFrom parallel detectCores
#' @keywords nonparametric
#' @seealso \code{\link{permute_rho}} \code{\link{phylosmith}}
#' @export
#' @return data.table
#' @examples
#' variable_correlation(soil_column, variables = 'Day',
#' treatment = c('Matrix', 'Treatment'), subset = 'Amended',
#' classification = 'Phylum', method = 'spearman', cores = 1)
# sourceCpp("src/correlations_Rcpp.cpp")
variable_correlation <-
function(phyloseq_obj,
         variables,
         treatment = NULL,
         subset = NULL,
         classification = NULL,
         method = 'spearman',
         cores = 1) {
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
  treatment <- check_numeric_treatment(phyloseq_obj, treatment)
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
  tryCatch({
    variables <- check_numeric_treatment(phyloseq_obj, variables)
  }, error = function(e){
    stop(
      "`variables` must be at least one column name, or index, from the sample_data() that contains numeric data.",
      call. = FALSE)})
  if (!(is.null(variables)) &
      any(!(variables %in% colnames(access(
        phyloseq_obj, 'sam_data'
      ))))) {
    stop(
      "`variables` must be at least one column
          name, or index, from the sample_data()",
      call. = FALSE
    )
  }
  match.arg(method, c("pearson", "kendall", "spearman"))
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment = treatment, subset = subset)
  treatment_name <- paste(treatment, collapse = sep)
  treatment_classes <- as.character(unique(phyloseq_obj@sam_data[[treatment_name]]))
  if(is.null(treatment)){treatment_classes <- list(NULL)}
  if(!(is.null(classification))){
    phyloseq_obj <- conglomerate_taxa(phyloseq_obj, classification, hierarchical = FALSE)
  }
  correlations <- data.table()
  for(k in treatment_classes){
    phyloseq_obj_subset <- taxa_filter(phyloseq_obj, treatment, k, drop_samples = TRUE)
    treatment_correlations <- Correlation(
      X = phyloseq_obj_subset@otu_table,
      Y = apply(as.matrix(phyloseq_obj_subset@sam_data[,variables]), 2, as.numeric),
      method = method
    )
    if(length(treatment_classes) > 1){
      treatment_correlations <- cbind(Treatment = k, treatment_correlations)
    }
    correlations <- rbind(correlations, treatment_correlations)
  }
  return(correlations)
}
#' Computes the correlation of numerical variables with taxa
#' Function from the phylosmith-package.
#'
#' Computes the correlation of numerical variables with taxa
#' @useDynLib phylosmith
#' @usage variable_correlation_heatmap(phyloseq_obj, variables,
#' treatment = NULL, subset = NULL, classification = NULL,
#' method = 'spearman', limits = c(-0.8, 0.8),
#' colors = 'default', significance_color = 'white', cores = 1)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object.
#' @param variables Numericla factors within the in the
#' \code{\link[phyloseq:sample_data]{sample_data}} to correlate with the
#' abundance data.
#' @param treatment Column name as a \code{string} or \code{numeric} in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any
#' samples that to not contain this factor. This can be a vector of multiple
#' factors to subset on.
#' @param classification Column name as a \code{string} or \code{numeric} in
#' the \code{\link[phyloseq:tax_table]{tax_table}} for the factor to
#' conglomerate by.
#' @param method Which correlation method to calculate, "pearson", "spearman".
#' @param limits The range for the legend, smaller limits will accentuate smaller
#' correlations.
#' @param colors the palette to use for the heatmap, default is viridis.
#' @param significance_color the color to use for the significance stars.
#' @param cores \code{numeric} Number of CPU cores to use for the pair-wise
#' permutations. Default (0) uses max cores available. Parallelization not
#' available for systems running MacOS without openMP configuration.
#' @importFrom parallel detectCores
#' @keywords nonparametric
#' @seealso \code{\link{permute_rho}} \code{\link{phylosmith}}
#' @export
#' @return data.table
#' @examples
#' variable_correlation_heatmap(soil_column, variables = 'Day',
#' treatment = c('Matrix', 'Treatment'), subset = 'Amended',
#' classification = 'Phylum', method = 'spearman', cores = 1,
#' colors = c("#2C7BB6", "white", "#D7191C"),
#' significance_color = 'black')

variable_correlation_heatmap <-
  function(phyloseq_obj,
           variables,
           treatment = NULL,
           subset = NULL,
           classification = NULL,
           method = 'spearman',
           limits = c(-0.8, 0.8),
           colors = 'default',
           significance_color = 'white',
           cores = 1) {
    correlations <- variable_correlation(
            phyloseq_obj,
            variables,
            treatment,
            subset,
            classification,
            method,
            cores
    )
    correlations$Significance<-cut(correlations$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
    if(is.null(treatment)){
      g <- ggplot(data = correlations, aes(x = Y, y = X, fill = rho))
    } else {
      g <- ggplot(data = correlations, aes(x = Treatment, y = X, fill = rho))
    }
    g <- g + geom_tile() +
      geom_text(aes(label = Significance), color = significance_color, size = 3) +
      labs(y = NULL, x = NULL, fill = method) +
      facet_grid(. ~ Y, drop = TRUE, scales = "free", space = "free_x")
    if('default' %in% colors){
      g <- g + ggraph::scale_fill_viridis(limit = limits)
    } else {
      g <- g + scale_fill_gradient2(low = colors[1], mid = colors[2], high = colors[3], limit = limits)
    }
    g <- g + theme_classic() +
      theme(
        axis.text.x = element_text(
          angle = 30,
          hjust = 1,
          size = 10
        ),
        axis.text.y = element_text(hjust = 0.95, size = 10),
        axis.title.x = element_text(size = 10, face = 'bold'),
        axis.title.y = element_text(size = 10, face = 'bold'),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 10, face = 'bold'),
        legend.text = element_text(size = 8),
        legend.spacing.x = unit(0.005, 'npc'),
        legend.key.size = unit(6, "mm"),
        panel.background = element_rect(color = 'black', size = 1.4),
        strip.text.x = element_text(size = 10, face = 'bold'),
        strip.background = element_rect(colour = 'black', size = 1.4)
      ) +
      scale_x_discrete(expand = expand_scale(mult = 0, add = 0))
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
#' @param variables Numericla factors within the in the
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
  variable_colors <- create_palette(sum(names(igraph::V(net)) %in% variables), colors)
  node_colors <- dcast.data.table(correlation_table, X ~ Y, value.var = 'weight')
  node_colors[is.na(node_colors)] <- 0
  node_colors <- apply(node_colors[,-1], 1,FUN=function(x) paste(rev(variables)[x>0], collapse = '_'))
  node_colors <- create_palette(length(unique(node_colors)))[factor(c(node_colors, variables))]
  node_sizes <- rep(5, length(names(igraph::V(net))))
  node_sizes[names(igraph::V(net)) %in% variables] <- 20
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

