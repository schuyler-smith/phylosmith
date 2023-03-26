#' Find unique taxa between treatments of a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and
#' finds which taxa are taxa that are unique to a specific subset of the data.
#' @useDynLib phylosmith
#' @usage unique_taxa(phyloseq_obj, treatment, subset = NULL)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param treatment Column name as a \code{string} or \code{numeric} in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any
#' samples that to not contain this factor. This can be a vector of multiple
#' factors to subset on.
#' @seealso \code{\link{common_taxa}}
#' @export
#' @return list
#' @examples unique_taxa(soil_column, c('Matrix', 'Treatment'))

unique_taxa <- function(
  phyloseq_obj, 
  treatment, 
  subset = NULL
) {
  check_args(
    phyloseq_obj  = phyloseq_obj,
    sam_data      = phyloseq_obj,
    treatment     = treatment,
    subset        = subset
  )
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, subset = subset)
  treatment_name <- paste(treatment, collapse = sep)
  treatment_classes <- unique(phyloseq_obj@sam_data[[treatment_name]])
  seen_taxa <- lapply(treatment_classes, FUN = function(trt) {
    taxa_names(taxa_filter(phyloseq_obj, treatment = treatment_name,
      subset = trt, frequency = 0))
  })
  names(seen_taxa) <- treatment_classes
  taxa_counts <- table(unlist(seen_taxa))
  unique_taxa <- names(taxa_counts[taxa_counts == 1])
  unique_taxa <- lapply(seen_taxa, FUN = function(toi) {
    toi[toi %in% unique_taxa]
  })
  
  return(unique_taxa)
}