#' Conglomerate taxa by sample on a given classification level.
#' Function from the phylosmith-package.
#'
#' Conglomerate taxa by sample on a given classification level from the
#' tax_table.
#' @useDynLib phylosmith
#' @usage conglomerate_taxa(phyloseq_obj, classification, hierarchical = TRUE,
#' use_taxonomic_names = TRUE)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}} with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param classification Column name as a \code{string} or \code{numeric} in
#' the \code{\link[phyloseq:tax_table]{tax_table}} for the factor to
#' conglomerate by.
#' @param hierarchical Whether the order of factors in the tax_table represent
#' a decreasing hierarchy (TRUE) or are independent (FALSE). If FALSE, will
#' only return the factor given by \code{classification}.
#' @param use_taxonomic_names If (TRUE), will use the hierarchical taxonomic
#' name from Domain to classification level. If (FALSE) will uses the
#' classification level with a numerical assignment.
#' @seealso \code{\link[phyloseq:tax_glom]{tax_glom()}}
#' @export
#' @return phyloseq-object
#' @examples conglomerate_taxa(soil_column, 'Phylum')

conglomerate_taxa <- function(
    phyloseq_obj,
    classification,
    hierarchical = TRUE,
    use_taxonomic_names = TRUE) {
  check_args(
    phyloseq_obj   = phyloseq_obj,
    taxa           = phyloseq_obj,
    hierarchical   = hierarchical,
    use_taxonomic_names = use_taxonomic_names
  )
  phyloseq_obj <- check_TaR(phyloseq_obj)
  taxa <- as(phyloseq_obj@tax_table, "matrix")
  taxa <- data.table::data.table(taxa, keep.rownames = "OTU")
  for (j in seq_len(ncol(taxa)))
    data.table::set(taxa, which(is.na(taxa[[j]])), j, "Unclassified")
  for (j in seq_len(ncol(taxa)))
    data.table::set(taxa, which(taxa[[j]] %in% "Incertae Sedis"), 
      j, "Unclassified")
  for (j in seq_len(ncol(taxa)))
    data.table::set(taxa, which(taxa[[j]] %in% ""), j, "Unclassified")

  if (hierarchical) {
    if(which(colnames(taxa) %in% classification) != length(colnames(taxa))) {
      taxa[, `:=`(seq(ncol(taxa))[-c(seq(which(colnames(taxa) %in%
        classification)))], NULL)]
    }
    data.table::set(taxa, 
      i = taxa[, .I[get(classification) %in% "Unclassified"]],
      j = classification, value = paste0("Unclassified ",
      unlist(taxa[taxa[, .I[get(classification) %in% "Unclassified"]],
      colnames(taxa)[which(colnames(taxa) %in% classification) - 1],
      with = F])))
    for (i in seq_len(ncol(taxa) - 1)) {
      data.table::set(taxa,
        i = taxa[, .I[get(classification) %in% "Unclassified Unclassified"]],
        j = classification, value = paste0("Unclassified ",
        unlist(taxa[taxa[, .I[get(classification)
        %in% "Unclassified Unclassified"]],
        colnames(taxa)[which(colnames(taxa) %in% classification) - i],
        with = F])))
    }
  } else {
    taxa[, seq(ncol(taxa))[-c(which(colnames(taxa) %in%
      classification))][-1] := NULL]
  }
  otus <- as(phyloseq_obj@otu_table, "matrix")
  otus <- data.table::data.table(otus, keep.rownames = "OTU")
  otus[, OTU := taxa[[classification]]]
  otus <- otus[, lapply(.SD, sum), by = OTU]
  otus <- as.matrix(otus, rownames = "OTU")
  taxa[, OTU := get(classification)]
  taxa <- as.matrix(unique(taxa, by = "OTU"), rownames = "OTU")

  if (!use_taxonomic_names) {
    rownames(taxa) <- paste0(classification, "_", seq(length(rownames(taxa))))
  }
  if (!(is.null(phyloseq_obj@phy_tree))) {
    warning("trees cannot be preserved after taxa conglomeration")
  }
  if (!(is.null(phyloseq_obj@refseq))) {
    warning("reference sequences cannot be preserved after taxa conglomeration")
  }
  phyloseq_obj <- phyloseq::phyloseq(
    phyloseq::otu_table(otus[, phyloseq::sample_names(phyloseq_obj)], 
      taxa_are_rows = TRUE), 
    phyloseq::tax_table(taxa),
    phyloseq_obj@sam_data
  )
  
  return(phyloseq_obj)
}