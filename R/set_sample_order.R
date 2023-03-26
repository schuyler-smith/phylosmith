#' Re-orders the samples of a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and changes the
#' order of sample index either based on the metadata, or a given order. If
#' metadata columns are used, they will take the factor order.
#' @useDynLib phylosmith
#' @usage set_sample_order(phyloseq_obj, treatment)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param treatment Column name as a \code{string} or \code{numeric} in the
#' \code{\link[phyloseq:sample_data]{sample_data}}, or vector of sample names
#' or indices in particular order.
#' @export
#' @return phyloseq object
#' @examples
#' set_sample_order(soil_column, c("Matrix", "Treatment"))

set_sample_order <- function(
  phyloseq_obj,
  treatment = NULL
) {
  check_args(
    phyloseq_obj = phyloseq_obj,
    treatment    = treatment
  )
  metadata <- as(phyloseq_obj@sam_data, "data.frame")
  metadata <- data.table::data.table(samples = rownames(metadata), metadata)
  if (is.null(treatment)) treatment <- "samples"
  data.table::setkeyv(metadata, treatment)
  phyloseq::otu_table(phyloseq_obj) <-
    phyloseq_obj@otu_table[, metadata$samples]
  phyloseq::sample_data(phyloseq_obj) <- data.frame(metadata, row.names = 1)
  
  return(phyloseq_obj)
}