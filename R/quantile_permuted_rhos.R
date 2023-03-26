#' Calculate quantiles for the permuted rho values from the Spearman-rank
#' co-occurrence. Function from the phylosmith-package.
#'
#' Calculate quantiles for the permuted rho values from the Spearman-rank
#' co-occurrence.
#' @useDynLib phylosmith
#' @usage quantile_permuted_rhos(permuted_rhos, p = 0.05, by_treatment = TRUE)
#' @param permuted_rhos A \code{data.table} output from
#' \code{\link[=permute_rho]{permute_rho}}.
#' @param p The significance threshold for setting cutoffs.
#' @param by_treatment Whether to find the rho cutoffs for each treatment
#' individually or for the entire experiment. Suggested to do by treatment
#' first, to see if there is any treatments that are outliers.
#' @seealso \code{\link[=permute_rho]{permute_rho}}
#' @export
#' @return data.frame
#' @examples
#' permuted_rhos <- permute_rho(soil_column,
#' treatment = c('Matrix', 'Treatment'), replicate_samples = 'Day',
#' permutations = 1,  cores = 0)
#' quantile_permuted_rhos(permuted_rhos)
#' quantile_permuted_rhos(permuted_rhos, by_treatment = FALSE)

quantile_permuted_rhos <- function(
  permuted_rhos,
  p = 0.05,
  by_treatment = TRUE
) {
  check_args(
    p            = p,
    by_treatment = by_treatment
  )
  if (by_treatment) {
    permuted_rhos[, Proportion := Count / sum(Count), by = "Treatment"]
    quantiles <-
      permuted_rhos[, list(lower = rho[sum(cumsum(Proportion) <=
                                             (p / 2))], upper = rho[sum(cumsum(Proportion) <= (1 - (p / 2)))]),
                    by = "Treatment"]
  } else {
    permuted_rhos[, Proportion := Count / sum(Count)]
    quantiles <-
      permuted_rhos[, list(lower = rho[sum(cumsum(Proportion) <=
                                             (p / 2))], upper = rho[sum(cumsum(Proportion) <= (1 - (p / 2)))])]
  }
  return(quantiles)
}