#' Find unique taxa between treatments of a phyloseq object. Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and finds which taxa are taxa that are unique to a specific subset of the data.
#' @useDynLib phylosmith
#' @usage find_unique_taxa(phyloseq_obj, treatment, subset = NULL)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package (must contain \code{\link[phyloseq:sample_data]{sample_data()}}).
#' @param treatment Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param subset Keyword for a subset of the treatments, can be a substring within the treatment names (e.g. 'control').
#' @keywords manip
#' @export
#' @import phyloseq
#' @examples
#' data(mock_phyloseq)
#' find_unique_taxa(mock_phyloseq, treatment = 2)
#' find_unique_taxa(mock_phyloseq, treatment = c("treatment", "day"), subset = "control")

find_unique_taxa <- function(phyloseq_obj, treatment, subset = NULL){
  #phyloseq_obj=mock_phyloseq; treatment=c(2,3); subset = "control"

  phyloseq_obj <- combine_treatments(phyloseq_obj, treatment)
  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  treatment_name <- paste(treatment, collapse = ".")

  treatments <- eval(parse(text=paste0("unique(phyloseq_obj@sam_data$", paste0(treatment_name), ")")))
  if(!(is.null(subset))){
    treatments <- eval(parse(text=paste0('treatments[grepl("', paste0(subset), '", treatments)]')))
  }

  seen_taxa <- lapply(treatments, FUN = function(trt){
    taxa_names(find_generalists(phyloseq_obj, treatment = treatment_name, subset = trt, frequency = 0))
  }); names(seen_taxa) <- treatments

  taxa_counts <- table(unlist(seen_taxa))
  unique_taxa <- names(taxa_counts[taxa_counts == 1])

  unique_taxa <- lapply(seen_taxa, FUN = function(toi){
    toi[toi %in% unique_taxa]
  })

  return(unique_taxa)
}
