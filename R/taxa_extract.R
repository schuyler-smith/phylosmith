#' Extracts specific taxa from all samples based on a given classification level.
#' Function from the phylosmith-package.
#' 
#' Extracts taxa from phyloseq objects based on taxanomic names.
#' @useDynLib phylosmith
#' @usage taxa_extract(phyloseq_obj, taxa_to_extract, classification = NULL)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}} with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param taxa_to_extract A vector of the classification names of the taxon to
#' be extracted. If you are using ASVs and want to extract specific sequences, this
#' should be set OTU or NULL. By default, it will search for the names in all
#' taxanomic ranks
#' @param classification Column name as a \code{string} or \code{numeric} in
#' the \code{\link[phyloseq:tax_table]{tax_table}} for the factor to
#' conglomerate by.
#' @seealso \code{\link[phyloseq:tax_glom]{tax_glom()}}
#' @export
#' @return phyloseq-object
#' @examples taxa_extract(soil_column, "Firmicutes", "Phylum")

taxa_extract <- function(
    phyloseq_obj,
    taxa_to_extract,
    classification = NULL
) {
  check_args(
    phyloseq_obj    = phyloseq_obj,
    tax_table       = phyloseq_obj,
    taxa_to_extract = taxa_to_extract,
    classification  = classification
  )
  taxa <- as(phyloseq_obj@tax_table, "matrix")
  taxa <- data.table::as.data.table(taxa, keep.rownames = "OTU")
  if (is.null(classification)) classification <- colnames(taxa)
  taxa <- cbind(taxa, ID = taxa$OTU)
  for (rank in classification) {
    set(taxa, i = unlist(lapply(tolower(taxa_to_extract), grep, 
        tolower(taxa[[rank]]))), j = "ID", value = "EXTRACT")
  }
  taxa <- taxa[ID == "EXTRACT"][, ID := NULL]
  if( nrow(taxa) == 0) {
    stop("no taxa found matching `taxa_to_extract`",
         call. = FALSE)
  }
  taxa <- as.matrix(taxa, rownames = "OTU")
  phyloseq::tax_table(phyloseq_obj) <- phyloseq::tax_table(taxa)

  return(phyloseq_obj)
}
