#' set_treatment_levels
#'
#' Reorders the levels of a metadata factor column in a
#' \code{\link[phyloseq]{phyloseq-class}} object
#' \code{\link[phyloseq:sample_data]{sample_data}}.
#' @useDynLib phylosmith
#' @usage set_treatment_levels(phyloseq_obj, treatment, order)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param treatment Column name as a \code{string} or \code{numeric} in the
#' \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param order The order of factors in \code{treatment} column as a vector of
#' strings. If assigned "numeric" will set ascending numerical order.
#' @export
#' @return phyloseq-object
#' @examples set_treatment_levels(soil_column, treatment = "Matrix",
#' order = c("Manure", "Soil", "Effluent"))
#' set_treatment_levels(soil_column, "Day", "numeric")

set_treatment_levels <- function(
    phyloseq_obj,
    treatment,
    order
) {
  check_args(
    phyloseq_obj = phyloseq_obj,
    sam_data     = phyloseq_obj,
    treatment    = treatment,
    order        = order
  )
  if (order[1] == "numeric") {
    order <- as.character(sort(as.numeric(unique(
      as.character(phyloseq_obj@sam_data[[treatment]])
    ))))
  }
  sample_data(phyloseq_obj)[[treatment]] <-
    factor(phyloseq_obj@sam_data[[treatment]], levels = order)
  return(phyloseq_obj)
}
