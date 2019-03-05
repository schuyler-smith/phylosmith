#' Find unique taxa between treatments of a phyloseq object. Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and finds which taxa are taxa that are unique to a specific subset of the data.
#' @useDynLib phylosmith
#' @usage find_unique_taxa(phyloseq_obj, treatment, subset = NULL)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with information about each sample, and it must contain \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each taxa/gene.
#' @param treatment Column name as a string or number in the \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any samples that to not contain this factor. This can be a vector of multiple factors to subset on.
#' @keywords manip
#' @export
#' @seealso \code{\link{find_common_taxa}}

find_unique_taxa <- function(phyloseq_obj, treatment, subset = NULL){
  #phyloseq_obj=mock_phyloseq; treatment=c(2,3); subset = "control"

  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  phyloseq_obj <- phyloseq(phyloseq_obj@otu_table, phyloseq_obj@tax_table, phyloseq_obj@sam_data)
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, subset = subset)
  treatment_name <- paste(treatment, collapse = sep)
  treatment_classes <- unique(phyloseq_obj@sam_data[[treatment_name]])

  seen_taxa <- lapply(treatment_classes, FUN = function(trt){
    taxa_names(taxa_filter(phyloseq_obj, treatment = treatment_name, subset = trt, frequency = 0))
  }); names(seen_taxa) <- treatment_classes

  taxa_counts <- table(unlist(seen_taxa))
  unique_taxa <- names(taxa_counts[taxa_counts == 1])

  unique_taxa <- lapply(seen_taxa, FUN = function(toi){
    toi[toi %in% unique_taxa]
  })

  return(unique_taxa)
}


#' Find taxa shared between treatments of a phyloseq object. Function from the phylosmith-package.
#'
#' Takes a \code{\link[phyloseq]{phyloseq-class}} object and finds which taxa are shared between all of the specified treatments.
#' @useDynLib phylosmith
#' @usage find_common_taxa(phyloseq_obj, treatment, subset = NULL, n = 'all')
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with information about each sample, and it must contain \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each taxa/gene.
#' @param treatment Column name as a string or number in the \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any samples that to not contain this factor. This can be a vector of multiple factors to subset on.
#' @param n Number of treatment groups that need to share the taxa to be considered a common taxa.
#' @keywords manip
#' @seealso \code{\link{find_unique_taxa}}
#' @export

find_common_taxa <- function(phyloseq_obj, treatment, subset = NULL, n = 'all'){
  #phyloseq_obj=mock_phyloseq; treatment=c(2,3); subset = "control"; n = 2

  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  phyloseq_obj <- phyloseq(phyloseq_obj@otu_table, phyloseq_obj@tax_table, phyloseq_obj@sam_data)
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, subset = subset)
  treatment_name <- paste(treatment, collapse = sep)
  treatment_classes <- unique(phyloseq_obj@sam_data[[treatment_name]])

  seen_taxa <- lapply(treatment_classes, FUN = function(trt){
    taxa_names(taxa_filter(phyloseq_obj, treatment = treatment_name, subset = trt, frequency = 0))
  }); names(seen_taxa) <- treatment_classes

  taxa_counts <- table(unlist(seen_taxa))
  if(n != 'all'){
    if(!(is.numeric(n))){
      message("n must be an integer.")
      return(NA)
    }
    shared_taxa <- names(taxa_counts[taxa_counts == round(n)])
  } else {
    shared_taxa <- names(taxa_counts[taxa_counts == length(treatment_classes)])
  }
  return(shared_taxa)
}
