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
#' negative_positive_colors = c('tomato3', 'gray22'))
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

variable_correlation_network <- function(
  phyloseq_obj,
  variables,
  classification = NULL,
  treatment = NULL,
  subset = NULL,
  correlation_table = NULL,
  method = "spearman",
  rho_threshold = c(-0.01, 0.01),
  p_threshold = 0.05,
  colors = "default",
  negative_positive_colors = c("tomato3", "gray22")
){
  check_args(
    phyloseq_obj   = phyloseq_obj,
    treatment      = treatment,
    subset         = subset,
    classification = classification,
    corr_method    = method,
    p              = p_threshold
  )
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, subset)
  if (!(is.null(classification))) {
    phyloseq_obj <- conglomerate_taxa(phyloseq_obj, classification,
      hierarchical = FALSE)
  }
  if(is.null(correlation_table)){
    correlation_table <- variable_correlation(phyloseq_obj, treatment = treatment,
      subset = subset, variables = variables, classification = classification,
      method = method, cores = 1)
  }
  correlation_table <- correlation_table[p <= p_threshold]
  correlation_table <- 
    correlation_table[rho <= rho_threshold[1] | rho >= rho_threshold[2]]
  if(nrow(correlation_table) == 0){
    stop("None of the supplied variables were sgnificantly correlated to the data.",
         call. = FALSE)
  }
  if(is.null(correlation_table[["Treatment"]])){
    correlation_table <- cbind(correlation_table, Treatment = "NA")
  }
  correlation_table <- correlation_table[, c("X", "Y", "Treatment", "rho", "p")]
  correlation_table[, Weight := abs(rho)]
  edge_sign <- vapply(correlation_table$rho, FUN = function(x) {
      as.numeric(as.logical(sign(x) + 1) + 1)
    }, numeric(1)
  )
  correlation_table[, Edge := factor(c("Negative", "Positive")[edge_sign],
    levels = c("Positive", "Negative"))]

  if (!is.null(phyloseq_obj@tax_table)){
    nodes <- data.table(as(phyloseq_obj@tax_table, "matrix"))
    nodes <-
      data.table("Node_Name" = rownames(phyloseq_obj@tax_table), nodes)
  } else {
    nodes <- data.table("Node_Name" = rownames(phyloseq_obj@tax_table))
  }
  phyloseq_obj <- relative_abundance(phyloseq_obj)
  nodes[,`Mean Relative Abundance` := bin(taxa_sums(phyloseq_obj) /
    nsamples(phyloseq_obj), nbins = 9)]
  nodes <- nodes[nodes[["Node_Name"]] %in% c(
    as.character(correlation_table$X),
    as.character(correlation_table$Y)
  ),]
  vars <- matrix(rep(variables, ncol(nodes)), nrow = length(variables))
  colnames(vars) <- colnames(nodes)
  nodes <- rbind(nodes, vars)

  net <- igraph::graph_from_data_frame(d = correlation_table, vertices = nodes,
    directed = FALSE)
  net <- igraph::simplify(net, remove.multiple = FALSE, remove.loops = TRUE)
  layout <- create_layout(net, layout = "igraph", algorithm = "fr")

  node_sizes <- rep(16, length(variables))
  node_sizes[names(node_sizes) %in% correlation_table[["Y"]]] <- 20

  variables_layout <- subset(layout, layout[[classification]] %in% variables)

  g <- ggraph::ggraph(layout) + coord_fixed() +
    ggraph::geom_edge_link(aes(color = Edge, width = Weight)) +
    ggraph::scale_edge_color_manual(values = negative_positive_colors) +
    ggraph::scale_edge_width_continuous(range = c(0.2,2))
  layout <- subset(layout, !(layout[[classification]] %in% variables))
  node_colors <- create_palette(length(unique(layout[,classification])), colors)

  g <- g + geom_point(data = layout, aes_string(x = "x", y = "y",
    size = "`Mean Relative Abundance`", fill = classification),
    pch = 21, color = "black", show.legend = TRUE) +
    scale_size_discrete(range = c(4,12)) +
    scale_fill_manual(values = node_colors)

  g <- g + geom_label(data = variables_layout, aes_string(x = "x", y = "y"),
    label = unname(apply(variables_layout, 1, function(x) {
      x[which(x %in% variables)] 
    }))[1, ], size = 5.5, alpha = 1, label.size = 1, fontface = "bold",
    show.legend = FALSE, label.padding = unit(0.65, "lines"), 
    label.r = unit(1, "lines")
  ) 
  n_cols <- ceiling(length(unique(layout[[classification]])) / 12)
  g <- g + guides(colour = guide_legend(override.aes = list(size=10)),
           alpha = FALSE,
           fill = guide_legend(ncol = n_cols, override.aes = list(size = 5)),
           edge_color = guide_legend(override.aes = list(edge_width = 2)))

  g <- g + theme_schuy("net")

  return(g)
}
