#' Pair-wise Spearman rank co-occurrence.
#'
#' A rewrite of the pair-wise Spearman rank co-occurrence routine written by
#' \href{https://github.com/germs-lab/FastCoOccur}{Jin Choi}. The routine has
#' been adapted to integrate with the \code{\link[Rcpp]{Rcpp-package}} API.
#' @useDynLib phylosmith
#' @usage co_occurrence(phyloseq_obj, treatment = NULL, subset = NULL,
#' rho = 0, p = 0.05, method = 'spearman', cores = 1)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object.
#' @param treatment Column name as a string, or vector of strings, from the
#' \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param subset A level within the \code{treatment}. Multiple levels can be 
#' given as a vector.
#' @param rho \code{numeric} The rho-value cutoff. All returned co-occurrences
#' will have a rho-value less than or equal to \code{rho} or less than or
#' equal to -\code{rho}.
#' @param p \code{numeric} The p-value cutoff. All returned co-occurrences
#' will have a p-value less than or equal to \code{p}.
#' @param method Which correlation method to calculate, "pearson", "spearman".
#' @param cores \code{numeric} Number of CPU cores to use for the pair-wise
#' permutations. Default (0) uses max cores available. Parallelization not
#' available for systems running MacOS without openMP configuration.
#' @importFrom parallel detectCores
#' @keywords nonparametric
#' @seealso \code{\link{permute_rho}} \code{\link{phylosmith}}
#' @export
#' @return data.table
#' @examples
#' co_occurrence(soil_column, treatment = c('Matrix', 'Treatment'),
#' subset = 'Amended', rho = 0.8, p = 0.05, cores = 1)

co_occurrence <- function(
  phyloseq_obj,
  treatment    = NULL,
  subset       = NULL,
  rho          = 0,
  p            = 0.05,
  method       = "spearman",
  cores        = 1
) {
  check_args(
    phyloseq_obj = phyloseq_obj,
    treatment    = treatment,
    subset       = subset,
    rho          = rho,
    p            = p,
    corr_method  = method,
    cores        = cores
  )
  if (cores == 0) cores <- (detectCores() - 1)
  phyloseq_obj <- 
    taxa_filter(phyloseq_obj, treatment, subset)
  treatment_name <- paste(treatment, collapse = sep)

  treatment_classes <- 
    as.character(unique(phyloseq_obj@sam_data[[treatment_name]]))
  treatment_indices <- lapply(
    treatment_classes,
    FUN = function(trt) {
        which(as.character(phyloseq_obj@sam_data[[treatment_name]]) %in% trt)
    }
  )
  if (is.null(treatment)) {
    treatment_classes <- "Experiment_Wide"
    treatment_indices <- list(seq(nsamples(phyloseq_obj)))
  }
  too_few <- sapply(treatment_indices, length) < 3
  if (any(too_few)) {
    treatment_indices <- treatment_indices[!too_few]
  }
  phyloseq_obj <- phyloseq_obj@otu_table
  co_occurrence <- data.table::data.table()
  if (!is.vector(rho)) rho <- c(-rho, rho)
  for(i in seq_along(treatment_indices)){
    treatment_co_occurrence <- Correlation(
      X               = phyloseq_obj[,treatment_indices[[i]]],
      lowerrho        = rho[1],
      upperrho        = rho[2],
      p_cutoff        = p,
      method          = method,
      ncores          = cores
    )
    treatment_co_occurrence[["X"]] <- rownames(
      phyloseq_obj[,treatment_indices[[i]]])[treatment_co_occurrence[["X"]]]
    treatment_co_occurrence[["Y"]] <- rownames(
      phyloseq_obj[,treatment_indices[[i]]])[treatment_co_occurrence[["Y"]]]
    if(length(treatment_indices) > 0){
      treatment_co_occurrence <- cbind(
        Treatment = treatment_classes[i], 
        treatment_co_occurrence)
    }
    co_occurrence <- rbind(co_occurrence, treatment_co_occurrence)
  }
  
  return(data.table::as.data.table(co_occurrence))
}