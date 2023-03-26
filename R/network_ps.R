#' Create an igraph network object of the co-occurrence from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Create an igraph network object of the co-occurrence from a phyloseq object.
#' This can be input into the co_occurrence_network function, or used for other
#' network creating scripts. The purpose is to be able to create reproducible
#' and comparable graphics.
#' @useDynLib phylosmith
#' @usage network_ps(phyloseq_obj,
#' treatment = NULL, subset = NULL, co_occurrence_table = NULL, rho = 0.6)
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
#' @param rho \code{numeric} The rho-value cutoff. All returned co-occurrences
#' will have a rho-value less than or equal to \code{rho} or less than or
#' equal to -\code{rho}.
#' @importFrom igraph graph_from_data_frame simplify cluster_fast_greedy E
#' @importFrom ggraph create_layout
#' @export
#' @return igraph network object
#' @examples
#' network_ps(soil_column, treatment = c('Matrix', 'Treatment'),
#' subset = 'Soil Amended')


network_ps <- function(
  phyloseq_obj,
  treatment = NULL,
  subset = NULL,
  co_occurrence_table = NULL,
  rho = 0.6
) {
  check_args(
      phyloseq_obj = phyloseq_obj,
      tax_table    = phyloseq_obj,
      sam_data     = phyloseq_obj,
      treatment    = treatment,
      subset       = subset,
      co_occurrence_table = co_occurrence_table,
      rho          = rho
  )
  phyloseq_obj <- 
    taxa_filter(phyloseq_obj, treatment, frequency = 0, subset = subset)
  treatment_name <- paste(treatment, collapse = sep)
  if (is.null(co_occurrence_table)) {
    co_occurrence_table <- 
      co_occurrence(phyloseq_obj, treatment, 
        method = "spearman")[rho >= rho | rho <= -rho]
  } else { 
    if(!is.null(subset)){
      co_occurrence_table <-
        co_occurrence_table[co_occurrence_table[["Treatment"]] %like% subset]
    }
  }
  if(is.null(co_occurrence_table[["Treatment"]])) {
    co_occurrence_table[, Treatment := "NA"]
  }
  co_occurrence_table <- 
    co_occurrence_table[, c("X", "Y", "Treatment", "rho", "p")]
  if (!is.null(phyloseq_obj@tax_table)){
    nodes <- data.table::data.table(as(phyloseq_obj@tax_table, "matrix"))
    nodes <- data.table::data.table(
      "Node_Name" = rownames(phyloseq_obj@tax_table), nodes)
  } else {
    nodes <- data.table::data.table(
      "Node_Name" = rownames(phyloseq_obj@otu_table))
  }
  phyloseq_obj <- relative_abundance(phyloseq_obj)
  nodes[,`Mean Relative Abundance` :=
    bin(taxa_sums(phyloseq_obj) / nsamples(phyloseq_obj), nbins = 9)]
  nodes <- nodes[nodes[["Node_Name"]] %in%
    c(as.character(co_occurrence_table$X),
      as.character(co_occurrence_table$Y)),]
  co_occurrence_table[, Weight := abs(rho)]
  clusters <- igraph::cluster_fast_greedy(igraph::simplify(
    igraph::graph_from_data_frame(
      d = co_occurrence_table, vertices = nodes, directed = FALSE),
    remove.multiple = TRUE, remove.loops = TRUE))$membership
  cluster_sizes <- table(clusters)
  nodes <- nodes[clusters %in% names(cluster_sizes[cluster_sizes >= 3])]
  co_occurrence_table <- 
    co_occurrence_table[X %in% nodes$Node_Name & Y %in% nodes$Node_Name]
  edge_sign <- vapply(
    co_occurrence_table$rho,
    FUN = function(x) as.numeric(as.logical(sign(x) + 1) + 1), numeric(1)
  )
  co_occurrence_table[,
    Edge := factor(c("Negative", "Positive")[edge_sign], 
      levels = c("Positive", "Negative"))]
  net <- igraph::graph_from_data_frame(
    d = co_occurrence_table,
    vertices = nodes,
    directed = FALSE)
    
  return(net)
}