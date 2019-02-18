#' Bootstraps the pair-wise Spearman rank co-occurrence, to determine a significant rho-cutoff. phylosmith
#'
#' Bootstraps the pair-wise Spearman rank co-occurrence, to determine a significant rho-cutoff.
#' @useDynLib phylosmith
#' @usage bootstrap_rho(phyloseq_obj, treatment, replicates = 'independent', permutations = 100,
#' cores = 0)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package.
#' @param treatment Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param replicates Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}} that indicates which samples are non-independent of each other.
#' @param permutations Number of iterations to compute.
#' @param cores Number of CPU cores to use for the pair-wise permutations. Default uses all cores available.
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

bootstrap_rho <- function(phyloseq_obj, treatment, replicates = 'independent', permutations = 100, cores = 0){
  # phyloseq_obj = mock_phyloseq; treatment = c("treatment", "day"); replicates = 'independent'; permutations = 10; p = 0; cores = 0;
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
	  rhos[[i]] <- FastCoOccur_rho(permuted_phyloseq_obj, treatment, cores = cores)
  }},
  interrupt = function(interrupt){rhos <- rhos[-length(rhos)]; message('Interrupted after ', length(rhos), ' permutations.'); return(rhos)})

  rhos <- unlist(rhos)
  # if(p == 0){
    return(rhos)
   } #else {
#   return(stats::quantile(rhos, 1-p, na.rm = TRUE))}
# }

# start_time <- Sys.time()
# for(indices in replicate_indices){
#   permuted_phyloseq_obj@otu_table[,indices] <- phyloseq_obj@otu_table[sample(1:n, n),indices]
# }
# rhos[1] <- mean(FastCoOccur(permuted_phyloseq_obj, treatment, cores = cores)$rho)
# end_time <- Sys.time()
# time_taken <- end_time - start_time
# message('Allow up to ', round((time_taken*permutations)/360, 1), ' hours to complete these calculations.')


#' pair-wise Spearman rank co-occurrence, written in efficient c++ code.
#'
#' A rewrite of the pair-wise Spearman rank co-occurrence routine written by \href{https://github.com/germs-lab/FastCoOccur}{Jin Choi}. The routine has been adapted to integrate with the \code{\link[Rcpp]{Rcpp-package}} API.
#' @useDynLib phylosmith
#' @usage FastCoOccur_rho(phyloseq_obj, treatment, cores = 0)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package.
#' @param treatment Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param cores Number of CPU cores to use for the pair-wise permutations. Default uses all cores available.
#' @import phyloseq
#' @import RcppArmadillo
#' @import RcppParallel
#' @import RcppProgress
#' @keywords nonparametric
#' @seealso \code{\link{bootstrap_rho}} \code{\link{phylosmith}}
#' @examples
#' data(mock_phyloseq)
#' phylosmith:::FastCoOccur_rho(mock_phyloseq, "day", 0.05)

# sourceCpp("src/FastCoOccur_rho_Rcpp.cpp")

FastCoOccur_rho <- function(phyloseq_obj, treatment, cores = 0){
  # phyloseq_obj = mock_phyloseq; treatment = c("treatment", "day"); p = 0.05; cores = 0
  options(warnings=-1)

  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment = treatment)
  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  treatment_name <- paste(treatment, collapse = ".")

  treatments <- as.character(unique(phyloseq_obj@sam_data[[treatment_name]]))
  treatment_indices <- lapply(treatments, FUN = function(trt){which(as.character(phyloseq_obj@sam_data[[treatment_name]]) %in% trt)-1})

  if(cores == 0){cores <- parallel::detectCores()}
  rhos <- FastCoOccur_rho_Rcpp(phyloseq_obj@otu_table, treatment_indices, treatments, cores)
  return(rhos)
}
