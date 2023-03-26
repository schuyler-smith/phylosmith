#' Filter taxa in phyloseq-object to only include core taxa.
#' Function from the phylosmith-package.
#'
#' Inputs a phyloseq object and finds which taxa are seen in a
#' given proportion of samples at a minimum relative abundance, either in the
#' entire dataset, by treatment, or a particular treatment of interest.
#' @useDynLib phylosmith
#' @usage taxa_core(phyloseq_obj, treatment = NULL, subset = NULL,
#' frequency = 0.5, abundance_threshold = 0.01)
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
#' @param abundance_threshold The minimum relative abundance the taxa is found
#' in for each sample.
#' @export
#' @return phyloseq-object
#' @examples taxa_core(soil_column, frequency = 0.2, abundance_threshold = 0.01)

taxa_core <- function(
  phyloseq_obj,
  treatment = NULL,
  subset = NULL,
  frequency = 0.5,
  abundance_threshold = 0.01
) {
  check_args(
    phyloseq_obj  = phyloseq_obj,
    treatment     = treatment,
    subset        = subset,
    frequency     = frequency,
    abundance_threshold = abundance_threshold
  )
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, subset)
  treatment_name <- paste(treatment, collapse = sep)
  core_taxa <- relative_abundance(phyloseq_obj)
  core_taxa <- melt_phyloseq(core_taxa)

  if(is.null(treatment)){
    N <- nsamples(phyloseq_obj)
    core_taxa <- core_taxa[Abundance >= abundance_threshold]
    core_taxa <- core_taxa[, .(count = .N),
      by = OTU][count >= floor(N * frequency)][["OTU"]]
  } else {
    taxa <- vector()
    treatments <- levels(core_taxa[[treatment_name]])
    for(treatment in treatments){
      N <- sum(phyloseq_obj@sam_data[[treatment_name]] == treatment)
      sub_table <- core_taxa[core_taxa[[treatment_name]] == treatment, ]
      sub_table <- sub_table[Abundance >= abundance_threshold]
      taxa <- c(taxa, sub_table[, .(count = .N),
        by = OTU][count >= floor(N * frequency)][["OTU"]])
    }
    core_taxa <- taxa
  }
  phyloseq_obj <- prune_taxa(core_taxa, phyloseq_obj)
  
  return(phyloseq_obj)
}