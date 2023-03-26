#' Transform abundance data in an \code{otu_table} to relative abundance,
#' sample-by-sample. Function from the phylosmith-package.
#'
#' Transform abundance data into relative abundance, i.e. proportional data.
#' This is an alternative method of normalization and may not be appropriate
#' for all datasets, particularly if your sequencing depth varies between
#' samples.
#' @useDynLib phylosmith
#' @usage relative_abundance(phyloseq_obj, sig_fig = 4)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param sig_fig Number of significant figures to round to.
#' @export
#' @return phyloseq-object
#' @examples relative_abundance(soil_column, 3)

relative_abundance <- function(
  phyloseq_obj,
  sig_fig = 4
) {
  check_args(
    phyloseq_obj = phyloseq_obj,
    sig_fig      = sig_fig
    )
  phyloseq_obj <- check_TaR(phyloseq_obj)
  abundance_table <- apply(phyloseq_obj@otu_table, 2, FUN = function(c) {
        round(c / sum(c), sig_fig)
  })
  abundance_table[is.na(abundance_table)] <- 0
  phyloseq::otu_table(phyloseq_obj) <-
    phyloseq::otu_table(abundance_table, taxa_are_rows = TRUE)
  return(phyloseq_obj)
}