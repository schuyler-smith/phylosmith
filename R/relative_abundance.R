#' Transform abundance data in an \code{otu_table} to relative abundance, sample-by-sample. phylosmith
#'
#' Transform abundance data into relative abundance, i.e. proportional data. This is an alternative method of normalization and may not be appropriate for all datasets, particularly if your sequencing depth varies between samples.
#' @useDynLib phylosmith
#' @usage relative_abundance(phyloseq_obj)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package.
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
