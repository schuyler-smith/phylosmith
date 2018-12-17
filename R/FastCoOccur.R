#' a very efficient c++ implementation of the pair-wise Spearman rank co-occurrence.
#'
#' This function is a rewrite of the pair-wise Spearman rank co-occurence routine written by Jin Choi. I adapted the routine to integrate with the Rcpp API.
#' @aliases cooccurrence co_occurrence
#' @useDynLib phylosmith
#' @usage FastCoOccur(phyloseq_obj, treatment, p = 0.05)
#' @param phyloseq_obj A \code{phyloseq-class} object created with the \link[=phyloseq]{name} package.
#' @param treatment the column name or number of the treatment to be comapred.
#' @param p the p-value cutoff. all returned co-occurrences must have a p-value less than or equal to p.
#' @keywords nonparametric
#' @export
#' @import data.table
#' @import phyloseq
#' @import RcppArmadillo
#' @import RcppParallel
#' @examples
#' data(mock_phyloseq)
#' FastCoOccur(mock_phyloseq, "day", 0.05)


# sourceCpp("src/FastCoOccur_Rcpp.cpp")

FastCoOccur <- function(phyloseq_obj, treatment, p = 0.05){
  options(warnings=-1)
  treatments <- as.character(unique(phyloseq_obj@sam_data[[treatment]]))
  treatment_indices <- lapply(treatments, FUN = function(trt){which(as.character(phyloseq_obj@sam_data[[treatment]]) %in% trt)-1})
  cooccurrence <- FastCoOccur_Rcpp(phyloseq_obj@otu_table, treatment_indices = treatment_indices, treatment_names = treatments, p_cutoff = p)
  return(as.data.table(cooccurrence))
}
