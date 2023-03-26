#' Merge samples based on common factor within sample_data.
#'
#' Inputs a phyloseq object and merges the samples that meet the
#' specified criteria into a single sample. This is meant for replicates, or
#' samples statistically proven to not be significantly different. This should
#' be used with caution as it may be a misleading representation of the data.
#' @useDynLib phylosmith
#' @usage conglomerate_samples(phyloseq_obj, treatment, subset = NULL)
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
#' @seealso \code{\link[phyloseq:merge_samples]{merge_samples()}}
#' @export
#' @return phyloseq-object
#' @examples conglomerate_samples(soil_column, 
#' treatment = c('Day', 'Matrix', 'Treatment'))

conglomerate_samples <- function(
    phyloseq_obj,
    treatment,
    subset = NULL) {
  check_args(
    phyloseq_obj   = phyloseq_obj,
    metadata       = phyloseq_obj,
    treatment      = treatment,
    subset         = subset
  )
  treatment_name <- paste(treatment, collapse = sep)
  phyloseq_obj <- merge_treatments(phyloseq_obj, treatment)
  sam_data <- data.table::data.table(as(phyloseq_obj@sam_data, "data.frame"),
    keep.rownames = "Sample")
  data.table::set(sam_data, j = "Sample", value = sam_data[[treatment_name]])
  if (length(treatment) > 1) sam_data[, (treatment_name) := NULL]
  otu_table <- data.table::data.table(t(as(phyloseq_obj@otu_table, "matrix")), 
    keep.rownames = "Sample")
  otu_table[, Sample := sam_data[["Sample"]]]
  if (!is.null(subset)) {
    sam_data <- sam_data[sam_data[, Reduce(`|`, lapply(.SD, `%in%`, subset)),
        .SDcols = c(treatment)]]
    otu_table <- otu_table[Sample %in% sam_data$Sample]
  }
  otu_table <- otu_table[, lapply(.SD, mean, na.rm = TRUE), by = Sample]
  otu_table <- t(as.matrix(otu_table, rownames = "Sample"))
  sam_data  <- data.frame(sam_data[, .SD[1], Sample], row.names = 1)

  phyloseq_obj <- phyloseq::phyloseq(
    phyloseq::otu_table(otu_table, taxa_are_rows = TRUE),
    phyloseq_obj@tax_table,
    phyloseq::sample_data(sam_data),
    phyloseq_obj@refseq,
    phyloseq_obj@phy_tree
  )
  
  return(phyloseq_obj)
}
