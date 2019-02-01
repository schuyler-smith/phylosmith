#' Bootstraps the pair-wise Spearman rank co-occurrence, to determine a significant rho-cutoff. phylosmith
#'
#' Bootstraps the pair-wise Spearman rank co-occurrence, to determine a significant rho-cutoff.
#' @useDynLib phylosmith
#' @usage bootstrap_rho(phyloseq_obj, treatment, replicates = 'independent', permutations = 10000, p = 0, cores = 0, cooccur_p = 0.05)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package.
#' @param treatment Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param replicates Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}} that indicates which samples are non-independent of each other.
#' @param permutations Number of iterations to compute.
#' @param p Significance cut-off for the mean rho values. If set to 0 (default) returns a vector of all permuted rho values.
#' @param cores Number of CPU cores to use for the pair-wise permutations. Default uses all cores available.
#' @param cooccur_p the p-value cutoff. all returned co-occurrences must have a p-value less than or equal to p.
#' @keywords nonparametric
#' @import data.table
#' @import phyloseq
#' @import RcppArmadillo
#' @import RcppParallel
#' @import RcppProgress
#' @seealso \code{\link{FastCoOccur}}
#' @examples
#' data(mock_phyloseq)
#' bootstrap_rho(mock_phyloseq, treatment = "day", permutations = 10)
#' @export

# sourceCpp('src/FastCoOccur_rho_Rcpp.cpp')

bootstrap_rho <- function(phyloseq_obj, treatment, replicates = 'independent', permutations = 10000, p = 0, cores = 0, cooccur_p = 0.05){
  # phyloseq_obj = mock_phyloseq; treatment = c("treatment", "day"); replicates = 'independent'; permutations = 10; p = 0; cores = 0; cooccur_p = 0.05
  options(warnings=-1)

  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  if(is.numeric(replicates)){replicates <- colnames(phyloseq_obj@sam_data[,replicates])}

  if(replicates == 'independent'){
    replicate_indices <- 1:ncol(phyloseq_obj@otu_table)
  } else {
    phyloseq_obj_reps <- combine_treatments(phyloseq_obj, c(treatment, replicates))
    replicate_name <- paste(c(treatment, replicates), collapse = ".")
    replicates <- as.character(unique(phyloseq_obj_reps@sam_data[[replicate_name]]))
    replicate_indices <- lapply(replicates, FUN = function(trt){which(as.character(phyloseq_obj_reps@sam_data[[replicate_name]]) %in% trt)})
  }

  rhos<-list()
  n <- nrow(phyloseq_obj@otu_table)
  permuted_phyloseq_obj <- phyloseq_obj

  tryCatch({
  for(i in 1:permutations){
    for(indices in replicate_indices){
      permuted_phyloseq_obj@otu_table[,indices] <- phyloseq_obj@otu_table[sample(1:n, n),indices]
    }
	  rhos[[i]] <- FastCoOccur_rho(permuted_phyloseq_obj, treatment, p = cooccur_p, cores = cores)
  }},
  interrupt = function(interrupt){rhos <- rhos[-length(rhos)]; message('Interrupted after ', length(rhos), ' permutations.'); return(rhos)})

  rhos <- unlist(rhos)
  if(p == 0){
    return(rhos)
  } else {
  return(stats::quantile(rhos, 1-p, na.rm = TRUE))}
}

# start_time <- Sys.time()
# for(indices in replicate_indices){
#   permuted_phyloseq_obj@otu_table[,indices] <- phyloseq_obj@otu_table[sample(1:n, n),indices]
# }
# rhos[1] <- mean(FastCoOccur(permuted_phyloseq_obj, treatment, p = cooccur_p, cores = cores)$rho)
# end_time <- Sys.time()
# time_taken <- end_time - start_time
# message('Allow up to ', round((time_taken*permutations)/360, 1), ' hours to complete these calculations.')

