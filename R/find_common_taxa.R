#' Find Unique Taxa
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and finds which taxa are shared between all of the specified treatments.
#' @aliases find_common_genes
#' @useDynLib phylosmith
#' @usage find_common_taxa(phyloseq_obj, column, keyword = NULL)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package (must contain \code{\link[phyloseq:sample_data]{sample_data()}}).
#' @param column Name or column number in the \code{\link[phyloseq:sample_data]{sample_data()}}.
#' @param keyword Looks for a subset of the column, can be a substring within the treatments names (e.g. 'control').
#' @keywords manip
#' @export
#' @import phyloseq
#' @examples
#' data(mock_phyloseq)
#' find_common_taxa(mock_phyloseq, column = 2)
#' find_common_taxa(mock_phyloseq, column = "day")

find_common_taxa <- function(phyloseq_obj, column, keyword = NULL){
  if(is.numeric(column)){column <- colnames(phyloseq_obj@sam_data[,column])
  } else {column <- eval(parse(text=paste0("colnames(phyloseq_obj@sam_data[,c('", paste0(column, collapse = "', '"), "')])")))
  }
  if(is.null(keyword)){
    treatments <- eval(parse(text=paste0("unique(phyloseq_obj@sam_data$", paste0(column), ")")))
  } else {
    treatments <- eval(parse(text=paste0("unique(phyloseq_obj@sam_data$", paste0(column), ")")))
    treatments <- eval(parse(text=paste0('treatments[grepl("', paste0(keyword), '", treatments)]')))
  }

  seen_taxa <- lapply(treatments, FUN = function(treatment){
    taxa_names(find_generalists(phyloseq_obj, treatment = column, subset = treatment, frequency = 0))
  })
  names(seen_taxa) <- treatments

  taxa_counts <- table(unlist(seen_taxa))
  shared_taxa <- names(taxa_counts[taxa_counts == length(treatments)])

  return(shared_taxa)
}
