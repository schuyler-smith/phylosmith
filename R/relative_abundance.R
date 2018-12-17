#' Relative Abundance
#'
#' This is used to transform abundance data into relative abundance, i.e. proportional data. This is an alternative normalization.
#' @useDynLib phylosmith
#' @usage relative_abundance(phyloseq_obj)
#' @param phyloseq_obj A \code{phyloseq-class} object created with the \code{\link[phyloseq]} package.
#' @keywords manip
#' @export
#' @import phyloseq
#' @import RcppArmadillo
#' @import RcppParallel
#' @examples
#' data(mock_phyloseq)
#' relative_abundance(mock_phyloseq)


relative_abundance <- function(phyloseq_obj){
  options(warnings=-1)
  phyloseq_obj <- transform_sample_counts(phyloseq_obj, function(sample) sample/sum(sample))
  return(phyloseq_obj)
}
