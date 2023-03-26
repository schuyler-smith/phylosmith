#' Remove specific taxa from all samples based on a given classification level.
#' Function from the phylosmith-package.
#'
#' Prunes taxa from phyloseq objects based on taxanomic names.
#' @useDynLib phylosmith
#' @usage taxa_prune(phyloseq_obj, taxa_to_remove,
#' classification=NULL, na_rm=FALSE)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}} with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param taxa_to_remove A vector of the classification names of the taxon to
#' be removed. If you are using ASVs and want to remove specific sequences, this
#' should be set OTU or NULL. By default, it will search for the names in all
#' taxanomic ranks
#' @param classification Column name as a \code{string} or \code{numeric} in
#' the \code{\link[phyloseq:tax_table]{tax_table}} for the factor to
#' conglomerate by.
#' @param na_rm if TRUE, and classification is specified, will remove taxon at
#' the classificaiton level that have NA values.
#' @seealso \code{\link[phyloseq:tax_glom]{tax_glom()}}
#' @export
#' @return phyloseq-object
#' @examples taxa_prune(soil_column, "Firmicutes", "Phylum")

taxa_prune <- function(
    phyloseq_obj,
    taxa_to_remove,
    classification = NULL,
    na_rm = FALSE) {
  check_args(
    phyloseq_obj   = phyloseq_obj,
    tax_table      = phyloseq_obj,
    taxa_to_remove = taxa_to_remove,
    classification = classification,
    na_rm          = na_rm
  )

  taxa <- as(phyloseq_obj@tax_table, "matrix")
  taxa <- data.table::as.data.table(taxa, keep.rownames = "OTU")
  if (is.null(classification)) classification <- colnames(taxa)
  for (rank in classification) {
    data.table::set(taxa, i = which(taxa[[rank]] %in% taxa_to_remove),
        j = rank, value = "REMOVE")
  }
  if (length(classification) == 1 & na_rm) {
    data.table::set(taxa, i = which(is.na(taxa[[classification]])),
        j = classification, value = "REMOVE")
  }
  taxa <- taxa[!(taxa[, Reduce(`|`, lapply(.SD, `%in%`, "REMOVE")),
    .SDcols = classification])]
  taxa <- as.matrix(taxa, rownames = "OTU")
  phyloseq::tax_table(phyloseq_obj) <- phyloseq::tax_table(taxa)

  rm(list = c("taxa"))
  gc()
  return(phyloseq_obj)
}