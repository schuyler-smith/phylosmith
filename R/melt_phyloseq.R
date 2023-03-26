#' Melt a phyloseq object into a data.table.
#' Function from the phylosmith-package.
#'
#' melt_phyloseq inputs a phyloseq object and melts its otu_table, taxa_tale,
#' and sample_Data into a single into a data.table.
#' @useDynLib phylosmith
#' @usage melt_phyloseq(phyloseq_obj)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}} with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @seealso \code{\link[phyloseq:psmelt]{psmelt()}}
#' @export
#' @return data.table
#' @examples melt_phyloseq(soil_column)

melt_phyloseq <- function(
  phyloseq_obj
) {
  check_args(phyloseq_obj = phyloseq_obj)
  phyloseq_obj <- check_TaR(phyloseq_obj)
  melted_phyloseq <- data.table::data.table(as(phyloseq_obj@otu_table, "matrix"
    ), keep.rownames = "OTU")
  if (!(is.null(phyloseq_obj@tax_table))) {
    taxa <- data.table::data.table(as(phyloseq_obj@tax_table, "matrix"),
        keep.rownames = "OTU")
  } else taxa <- NULL

  if (!(is.null(phyloseq_obj@sam_data))) {
    sample_data <- data.table::data.table(as(phyloseq_obj@sam_data, "data.frame"),
        stringsAsFactors = FALSE, keep.rownames = TRUE)
    suppressWarnings(sample_data[, Sample := NULL])
    data.table::setnames(sample_data, "rn", "Sample")
  } else {
    sample_data <- 
        data.table::data.table(Sample = phyloseq::sample_names(phyloseq_obj))
  }
  if (!(is.null(taxa))) {
    melted_phyloseq <- merge(melted_phyloseq, taxa, by = "OTU")
  }
  melted_phyloseq <- data.table::melt(
    melted_phyloseq,
    variable.name = "Sample",
    value.name = "Abundance",
    id.vars = colnames(taxa)
  )
  melted_phyloseq <-
    merge(melted_phyloseq, sample_data, by = "Sample")

  data.table::setorder(melted_phyloseq, -Abundance)

  return(melted_phyloseq)
}
