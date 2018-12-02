#' Find Unique Taxa
#'
#' This function takes a phyloseq object and finds which taxa are taxa that are unique to a specific subset of the data.
#' @useDynLib phylosmith
#' @usage find_unique_taxa(phyloseq_obj, column, keyword = NULL)
#' @param phyloseq_obj A phyloseq-class object created with the phyloseq package (must contain sample_data()).
#' @param column Name or column number in the sample_data(). Function then checks if taxa seen in frequency in each treatment.
#' @param keyword Looks for a subset of the column, can be a substring within the treatments names (e.g. 'control').
#' @keywords unique taxa phyloseq phylosmith
#' @export
#' @import phyloseq
#' @examples
#' data(mock_phyloseq)
#' find_unique_taxa(mock_phyloseq, column = 2)
#' find_unique_taxa(mock_phyloseq, column = "day")

find_unique_taxa <- function(phyloseq_obj, column, keyword = NULL){
  if(is.numeric(column)){column <- colnames(phyloseq::sample_data(phyloseq_obj)[,column])
    } else {column <- eval(parse(text=paste0("colnames(phyloseq::sample_data(phyloseq_obj)[,c('", paste0(column, collapse = "', '"), "')])")))
  }
  if(is.null(keyword)){
    treatments <- eval(parse(text=paste0("unique(phyloseq::sample_data(phyloseq_obj)$", paste0(column), ")")))
  } else {
    treatments <- eval(parse(text=paste0("unique(phyloseq::sample_data(phyloseq_obj)$", paste0(column), ")")))
    treatments <- eval(parse(text=paste0('treatments[grepl("', paste0(keyword), '", treatments)]')))
  }

  seen_taxa <- lapply(treatments, FUN = function(treatment){
    phyloseq::taxa_names(find_generalists(phyloseq_obj, treatments = column, subset = treatment, frequency = 0))
  })
  names(seen_taxa) <- treatments

  taxa_counts <- table(unlist(seen_taxa))
  unique_taxa <- names(taxa_counts[taxa_counts < 2])

  unique_taxa <- lapply(seen_taxa, FUN = function(toi){
    toi[toi %in% unique_taxa]
  })

  return(unique_taxa)
}
