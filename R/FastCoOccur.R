#' pair-wise Spearman rank co-occurrence, written in efficient c++ code.
#'
#' A rewrite of the pair-wise Spearman rank co-occurence routine written by \href{https://github.com/germs-lab/FastCoOccur}{Jin Choi}. The routine has been adapted to integrate with the \link[=Rcpp]{Rcpp} API.
#' @useDynLib phylosmith
#' @usage FastCoOccur(phyloseq_obj, treatment, p = 0.05)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package.
#' @param treatment Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}}.
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
  # phyloseq_obj = mock_phyloseq; treatment = c("treatment", "day"); p = 0.05
  options(warnings=-1)

  phyloseq_obj <- combine_treatments(phyloseq_obj, treatment)
  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  treatment_name <- paste(treatment, collapse = ".")

  treatments <- as.character(unique(phyloseq_obj@sam_data[[treatment_name]]))
  treatment_indices <- lapply(treatments, FUN = function(trt){which(as.character(phyloseq_obj@sam_data[[treatment_name]]) %in% trt)-1})
  cooccurrence <- FastCoOccur_Rcpp(phyloseq_obj@otu_table, treatment_indices = treatment_indices, treatment_names = treatments, p_cutoff = p)
  return(as.data.table(cooccurrence))
}
