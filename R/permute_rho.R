#' Permutes the pair-wise Spearman rank co-occurrence, to determine a
#' significant rho-cutoff. Function from the phylosmith-package.
#'
#' Permutes the pair-wise Spearman rank co-occurrence, to determine a
#' significant rho-cutoff.
#' @useDynLib phylosmith
#' @usage permute_rho(phyloseq_obj, treatment = NULL, subset = NULL,
#' replicate_samples = NULL, permutations = 10, method = 'spearman',
#' cores = 1)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object.
#' @param treatment Column name as a \code{string} or \code{numeric} in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any
#' samples that to not contain this factor. This can be a vector of multiple
#' factors to subset on.
#' @param replicate_samples Column name as a \code{string} or \code{numeric}
#' in the \code{\link[phyloseq:sample_data]{sample_data}} that indicates which
#' samples are non-independent of each other.
#' @param permutations \code{numeric} Number of iterations to compute.
#' @param method Which correlation method to calculate, "pearson", "spearman".
#' @param cores \code{numeric} Number of CPU cores to use for the pair-wise
#' permutations. Default (0) uses max cores available. Parallelization not
#' available for systems running MacOS without openMP configuration.
#' @keywords nonparametric
#' @importFrom parallel detectCores
#' @seealso \code{\link{co_occurrence}}
#' @export
#' @return table
#' @examples
#' permute_rho(soil_column, treatment = c('Matrix', 'Treatment'),
#' subset = 'Amended', replicate_samples = 'Day', permutations = 1,  cores = 0)

permute_rho <- function(
  phyloseq_obj,
  treatment = NULL,
  subset = NULL,
  replicate_samples = NULL,
  permutations = 10,
  method = "spearman",
  cores = 1
) {
check_args(
  phyloseq_obj  = phyloseq_obj,
  treatment     = treatment,
  subset        = subset,
  replicate_samples = replicate_samples,
  permutations  = permutations,
  corr_method   = method,
  cores         = cores
)

  phyloseq_obj <-
    taxa_filter(phyloseq_obj, treatment = treatment, subset = subset)
  phyloseq_obj <- relative_abundance(phyloseq_obj)
  treatment_name <- paste(treatment, collapse = sep)
  treatment_classes <-
    as.character(unique(phyloseq_obj@sam_data[[treatment_name]]))
  treatment_indices <- lapply(treatment_classes, FUN = function(trt) {
    which(as.character(phyloseq_obj@sam_data[[treatment_name]]) %in% trt) - 1
  })
  if (is.null(treatment)) {
    treatment_classes <- "Experiment_Wide"
    treatment_indices <- list(seq(nsamples(phyloseq_obj)) - 1)
  }
  replicate_sample_classes <- vector()
  if (is.null(replicate_samples)){
    replicate_indices <- seq(ncol(phyloseq_obj@otu_table))
  } else if (!(is.null(replicate_samples))){
    phyloseq_obj_reps <-
      merge_treatments(phyloseq_obj, c(treatment, replicate_samples))
    replicate_name <- paste(c(treatment, replicate_samples), collapse = sep)
    replicate_sample_classes <-
      as.character(unique(phyloseq_obj_reps@sam_data[[replicate_name]]))
    replicate_indices <- lapply(replicate_sample_classes, FUN = function(trt) {
      which(as.character(phyloseq_obj_reps@sam_data[[replicate_name]]) %in% trt)
    })
  }
  rhos <- data.table::data.table()
#   phyloseq_obj <- relative_abundance(phyloseq_obj)
  permuted_counts <- as(phyloseq_obj@otu_table, "matrix")
  tryCatch({
    for (j in seq(permutations)){
      co_occurrence_table <- data.table::data.table()
      if (is.null(replicate_samples)){
        n <- nrow(permuted_counts)
        permuted_counts <- apply(permuted_counts, 2, FUN = function(x){
          sample(x, n)
        })
      } else {for (indices in replicate_indices) {
        n <- nrow(permuted_counts)
        permuted_counts[, indices] <-
          permuted_counts[sample(seq(n), n), indices]
      }}
      for (i in seq_along(treatment_classes)){
        rep_co_occurrence_table <- 
          data.table::data.table(permute_rho_Rcpp(
            phyloseq_obj@otu_table[, treatment_indices[[i]]],
            permuted_counts[, treatment_indices[[i]]],
            method,
            cores
          ))
        rep_co_occurrence_table[, rho := round(rho, 2)]
        rep_co_occurrence_table[, Count := .N, by = .(rho)]
        rep_co_occurrence_table <- unique(rep_co_occurrence_table)
        if(length(treatment_classes) > 1){
          rep_co_occurrence_table <-
            cbind(Treatment = treatment_classes[i], rep_co_occurrence_table)
        }
        co_occurrence_table <-
          rbind(co_occurrence_table, unique(rep_co_occurrence_table))
      }

      if (length(treatment_classes) > 1) {
        rhos <- data.table::rbindlist(
          list(rhos, co_occurrence_table))[, lapply(.SD, sum, na.rm = TRUE),
            by = .(Treatment, rho)]
      } else {
        rhos <- data.table::rbindlist(
          list(rhos, co_occurrence_table))[, lapply(.SD, sum, na.rm = TRUE),
            by = .(rho)]
      }
    }
  },
  interrupt = function(interrupt) {
    message("Interrupted after ", i - 1, " permutations completed.")
    data.table::setkey(rhos, Treatment, rho)
    return(rhos)
  })
  if(length(treatment_classes) > 1){
    data.table::setkey(rhos, Treatment, rho)
  } else {
    data.table::setkey(rhos, rho)
  }
  return(rhos)
}