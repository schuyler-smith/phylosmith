#' Filter taxa based on proportion of samples they are observed in.
#' Function from the phylosmith-package.
#'
#' Inputs a phyloseq object and finds which taxa are seen in a
#' given proportion of samples, either in the entire dataset, by treatment, or
#' a particular treatment of interest.
#' @useDynLib phylosmith
#' @usage taxa_filter(phyloseq_obj, treatment = NULL, subset = NULL,
#' frequency = 0, below = FALSE, drop_samples = FALSE)
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
#' @param frequency The proportion of samples the taxa is found in.
#' @param below Does frequency define the minimum (\code{FALSE}) or maximum
#' (\code{TRUE}) proportion of samples the taxa is found in.
#' @param drop_samples Should the function remove samples that that are empty
#' after removing taxa filtered by frequency (\code{TRUE}).
#' @export
#' @return phyloseq-object
#' @examples taxa_filter(soil_column, frequency = 0.8)
#' taxa_filter(soil_column, treatment = c("Matrix", "Treatment"),
#' subset = "Soil Amended", frequency = 0.8)

taxa_filter <- function(
  phyloseq_obj,
  treatment = NULL,
  subset = NULL,
  frequency = 0,
  below = FALSE,
  drop_samples = FALSE
) {
  check_args(
    phyloseq_obj = phyloseq_obj,
    treatment    = treatment,
    subset       = subset,
    frequency    = frequency,
    below        = below,
    drop_samples = drop_samples
  )
  phyloseq_obj   <- check_TaR(phyloseq_obj)
  phyloseq_table <- data.table::data.table(as(phyloseq_obj@otu_table, "matrix"),
    keep.rownames = "OTU")
  phyloseq_table <- data.table::melt(phyloseq_table, id.vars = "OTU",
    variable.name = "Sample", value.name = "Abundance")
  phyloseq_table <- phyloseq_table[Abundance > 0]
  if (!(is.null(treatment))) {
    phyloseq_obj <- merge_treatments(phyloseq_obj, treatment)
    treatment_name <- paste(treatment, collapse = sep)
    sam_data <- data.table::data.table(as(phyloseq_obj@sam_data, "data.frame"),
      keep.rownames = "Sample")
    sam_data <- sam_data[, c("Sample", treatment, treatment_name), with = FALSE]
    if (!(is.null(subset))) {
      sub_rows <- sam_data[, Reduce(`|`, lapply(.SD, `%in%`, subset)),
        .SDcols = c(treatment, treatment_name)]
      sam_data <- sam_data[sub_rows]
      phyloseq_table <- phyloseq_table[Sample %in% sam_data$Sample]
      phyloseq_obj <-
        phyloseq::prune_samples((sample_names(phyloseq_obj) %in%
        sam_data$Sample), phyloseq_obj)
    }
    sample_counts <-
      phyloseq_table[, .(n_samples = data.table::uniqueN(Sample))]
    taxa_counts <- phyloseq_table[, .(count = .N), by = c("OTU")]
    taxa_counts[, proportion := count / unlist(sample_counts)]
    phyloseq_table <-
      phyloseq_table[OTU %in% taxa_counts[proportion > frequency]$OTU]

    treatment_classes <- unique(sam_data[[treatment_name]])
    taxa <- data.table::data.table()
    for (trt in sam_data) {
      sub_table <- phyloseq_table[Sample %in% 
        sam_data[get(treatment_name) %in% trt]$Sample]
      sample_counts <- sub_table[, .(n_samples = data.table::uniqueN(Sample))]
      taxa_counts <- sub_table[, .(count = .N), by = c("OTU")]
      taxa_counts[, proportion := count/unlist(sample_counts)]
      taxa <- rbind(taxa, taxa_counts)
    }
    rm("sub_table", "sam_data")
  } else {
    sample_counts <-
      phyloseq_table[, .(n_samples = data.table::uniqueN(Sample))]
    taxa_counts <- phyloseq_table[, .(count = .N), by = c("OTU")]
    taxa_counts[, n_samples := sample_counts$n_samples]
    taxa <- taxa_counts[, .(proportion = count / n_samples), by = c("OTU")]
  }
  rm("phyloseq_table", "sample_counts", "taxa_counts")
  if (below) {
    taxa <- taxa[proportion <= frequency]
  } else {
    taxa <- taxa[proportion >= frequency]
  }
  taxa <- unique(taxa[["OTU"]])
  phyloseq_obj <- taxa_prune(phyloseq_obj,
    phyloseq::taxa_names(phyloseq_obj)[!(phyloseq::taxa_names(
      phyloseq_obj) %in% taxa)])
  if (drop_samples == TRUE) {
    phyloseq_obj <-
      phyloseq::prune_samples(phyloseq::sample_sums(phyloseq_obj) > 0,
      phyloseq_obj)
  }
  rm("taxa")
  gc()
  return(phyloseq_obj)
}
