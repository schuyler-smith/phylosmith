#' Curate co-occurrence data. Function from the phylosmith-package.
#'
#' Used to curate a co-occurrence table from the \code{\link{FastCoOccur}} function. Takes a list of taxa and finds all pairs containing any of those taxa.
#' @useDynLib phylosmith
#' @usage curate_cooccurrence(cooccurrence_table, taxa_of_interest, number_of_treatments = 1)
#' @param cooccurrence_table co-occurrence table generated with \code{\link{FastCoOccur}}, or formatted in the same way.
#' @param taxa_of_interest a list or vector of taxa names, like those generated with \code{\link{find_unique_taxa}}.
#' @param number_of_treatments how many treatments should the taxa of interest be seen in? require \code{integer} or 'all' (default = 1).
#' @keywords manip
#' @export
#' @import data.table
#' @import phyloseq
#' @import RcppArmadillo
#' @import RcppParallel
#' @examples
#' data(mock_phyloseq)

curate_cooccurrence <- function(cooccurrence_table, taxa_of_interest, number_of_treatments = 1){
  sub_cooccurrence <- cooccurrence_table[(cooccurrence_table[[2]] %in% taxa_of_interest | cooccurrence_table[[3]] %in% taxa_of_interest),]
  toi_table <- unique(cbind(rbindlist(list(sub_cooccurrence[,1], sub_cooccurrence[,1])), rbindlist(list(sub_cooccurrence[,2], sub_cooccurrence[,3]))))
  toi_table <- toi_table[toi_table$gene_1 %in% taxa_of_interest]
  if(number_of_treatments == 'all'){number_of_treatments <- length(unique(sub_cooccurrence$Treatment))
  } else {number_of_treatments <- number_of_treatments}
  toi <- names(table(toi_table$gene_1)[table(toi_table$gene_1) >= number_of_treatments])
  sub_cooccurrence <- cooccurrence_table[(cooccurrence_table[[2]] %in% toi | cooccurrence_table[[3]] %in% toi),]

  # sourceCpp("src/arrange_cooccurrence_table_tbb.cpp")
  arranged_coocurrence <- as.data.table(arrange_cooccurr_table(sub_cooccurrence, toi))

  setorder(arranged_coocurrence, Treatment)
  return(arranged_coocurrence)
}
