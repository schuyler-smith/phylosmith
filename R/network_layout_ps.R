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

network_layout_ps <- function (
  phyloseq_obj,
  treatment = NULL,
  subset = NULL,
  co_occurrence_table = NULL,
  algorithm = "fr"
) {
  check_args(phyloseq_obj = phyloseq_obj)
  net <- network_ps(phyloseq_obj, treatment, subset, co_occurrence_table)
  layout <-
    ggraph::create_layout(net, layout = "igraph", algorithm = algorithm)
    
  return(layout)
}