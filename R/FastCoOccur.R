#' FastCoOccur
#'
#' This function is a rewrite of the pair-wise Spearman rank co-occurence routine written by Jin Choi. I adapted the routine to integrate with the Rcpp API.
#' @useDynLib phylosmith
#' @usage FastCoOccur(phyloseq_obj, treatment, p)
#' @param phyloseq_obj A phyloseq-class object created with the phyloseq package.
#' @param treatment the column name or number of the treatment to be comapred.
#' @param p the p-value cutoff. all returned co-occurrences must have a p-value less than or equal to p.
#' @keywords FastCoOccur Spearman co-occurrence coocurrence phylosmith
#' @export
#' @import phyloseq
#' @import RcppArmadillo
#' @import RcppParallel
#' @import data.table
#' @examples
#' data(mock_phyloseq)
#' FastCoOccur(mock_phyloseq, "day", 0.05)

FastCoOccur <- function(phyloseq_obj, treatment, p = 0.05){
  options(warnings=-1)
  treatments <- as.character(unique(sample_data(phyloseq_obj)[[treatment]]))
  treatment_indices <- lapply(treatments, FUN = function(trt){which(as.character(sample_data(phyloseq_obj)[[treatment]]) %in% trt)-1})
  cooccurrence <- FastCoOccur_Rcpp(otu_table(phyloseq_obj), treatment_indices = treatment_indices, treatment_names = treatments, p_cutoff = p)
  return(as.data.table(cooccurrence))
}
