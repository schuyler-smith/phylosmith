#' Normalizes abundance data in an \code{otu_table} using library size.
#' Function from the phylosmith-package.
#'
#' Transform abundance data to equal counts using library size and the geometric mean, 
#' i.e. proportional data.
#' @useDynLib phylosmith
#' @usage library_size(phyloseq_obj)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}} with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @export
#' @return phyloseq-object
#' @examples library_size(soil_column)

library_size <- function(
  phyloseq_obj
) {
  check_args(phyloseq_obj = phyloseq_obj)
  phyloseq_obj <- check_TaR(phyloseq_obj)
  phyloseq_obj <- taxa_filter(phyloseq_obj)
  totals       <- phyloseq::sample_sums(phyloseq_obj)
  norm_factor  <- totals / exp(mean(log(totals)))
  phyloseq::otu_table(phyloseq_obj) <-
    round(phyloseq_obj@otu_table / rep(norm_factor,
    each = (phyloseq::ntaxa(phyloseq_obj))))
    
  return(phyloseq_obj)
}