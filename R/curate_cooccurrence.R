#' Curate Co-Occurrence
#'
#' This is used to to curate a co-occurrence table from the FastCoOccur function. It will take a list of taxa, find all pairs with those taxa.
#' @useDynLib phylosmith
#' @usage curate_cooccurrence(cooccurrence_table, taxa_of_interest, number_of_treatments = 'all')
#' @param cooccurrence_table co-occurrence table generated with FastCoOccur().
#' @param taxa_of_interest a list or vector of taxa names.
#' @param number_of_treatments how many treatments should the taxa of interest be seen in? require `int` or 'all' (default).
#' @keywords cooccurrence taxa phyloseq phylosmith
#' @export
#' @import phyloseq
#' @import data.table
#' @examples
#' data(mock_phyloseq)


curate_cooccurrence <- function(cooccurrence_table, taxa_of_interest, number_of_treatments = 'all'){
  sub_cooccurrence <- cooccurrence_table[(cooccurrence_table[[2]] %in% taxa_of_interest | cooccurrence_table[[3]] %in% taxa_of_interest),]
  toi_table <- unique(cbind(rbindlist(list(sub_cooccurrence[,1], sub_cooccurrence[,1])), rbindlist(list(sub_cooccurrence[,2], sub_cooccurrence[,3]))))
  toi_table <- toi_table[toi_table$gene_1 %in% taxa_of_interest]
  if(number_of_treatments == 'all'){number_of_treatments <- length(unique(sub_cooccurrence$Treatment))
  } else {number_of_treatments <- number_of_treatments}
  toi <- names(table(toi_table$gene_1)[table(toi_table$gene_1) >= number_of_treatments])
  sub_cooccurrence <- cooccurrence_table[(cooccurrence_table[[2]] %in% toi | cooccurrence_table[[3]] %in% toi),]
  n_pairs <- nrow(sub_cooccurrence)
  for(row in 1:n_pairs){
    pair <- sub_cooccurrence[row,]
    if(unlist(pair[,3]) %in% toi & !(unlist(pair[,2]) %in% toi)){
      sub_cooccurrence[row,2] <- pair[,3]; sub_cooccurrence[row,3] <- pair[,2]}
    if(unlist(pair[,2]) %in% toi & unlist(pair[,3]) %in% toi){
      new_row <- pair; new_row[,2] <- pair[,3]; new_row[,3] <- pair[,2]
      sub_cooccurrence <- rbind(sub_cooccurrence, new_row)
    }
  }
  return(sub_cooccurrence)
}

