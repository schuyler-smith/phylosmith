#' Compute proportions for taxa.
#'
#' Computes the proportion of a taxa classification. Function from the
#' phylosmith-package.
#' @useDynLib phylosmith
#' @usage taxa_proportions(phyloseq_obj, classification, treatment = NA)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. it
#' must contain \code{\link[phyloseq:tax_table]{tax_table()}}) with
#' information about each taxa/gene.
#' @param classification Column name as a \code{string} or \code{numeric} in
#' the \code{\link[phyloseq:tax_table]{tax_table}} for the proportions to be
#' reported on.
#' @param treatment Column name as a \code{string} or \code{numeric} in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @export
#' @return data.table
#' @examples taxa_proportions(soil_column, 'Phylum', treatment = NULL)
#' taxa_proportions(soil_column, 'Phylum', treatment = 'sample')
#' taxa_proportions(soil_column, 'Phylum', treatment = c('Matrix', 'Treatment'))

taxa_proportions <- function(
  phyloseq_obj,
  classification,
  treatment = NULL
) {
  check_args(
    phyloseq_obj   = phyloseq_obj,
    tax_table      = phyloseq_obj,
    classification = classification
  )
  if (!is.null(treatment)) if (treatment != "Sample"){
    check_args(treatment = treatment)
  }
  phyloseq_obj <- taxa_filter(phyloseq_obj)
  treatment_name <- paste(treatment, collapse = sep)

  phyloseq_obj <- suppressWarnings(
    conglomerate_taxa(phyloseq_obj, classification, FALSE))
  class_table <- melt_phyloseq(phyloseq_obj)
  if (!is.null(treatment)) {
    class_table <- class_table[, .(Abundance = sum(Abundance)),
      by = c(treatment_name, classification)]
    class_table[, Proportion := round(Abundance / sum(Abundance), 4),
                by = treatment_name]
    data.table::setkeyv(class_table, treatment_name)
  } else {
    class_table <- class_table[, .(Abundance = sum(Abundance)),
      by = c(classification)]
    class_table[, Proportion := round(Abundance / sum(Abundance), 4)]
  }
  class_table[, Abundance := NULL][]
  
  return(class_table)
}