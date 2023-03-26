#' Curate co-occurrence data. Function from the phylosmith-package.
#'
#' Used to curate a co-occurrence table from the \code{\link{co_occurrence}}
#' function. Takes a list of taxa and finds all pairs containing any of those
#' taxa.
#' @useDynLib phylosmith
#' @usage curate_co_occurrence(co_occurrence_table, taxa_of_interest,
#' number_of_treatments = 1)
#' @param co_occurrence_table Table of the co-occurrence of taxa/genes in the
#' \code{phyloseq_obj}, computed using \code{\link{co_occurrence}}.
#' @param taxa_of_interest A list or vector of taxa names, like those
#' generated with \code{\link{unique_taxa}}.
#' @param number_of_treatments How many treatments should the taxa of interest
#' be seen in? Requires \code{integer} or 'all' (default = 1).
#' @keywords manip
#' @seealso \code{\link{co_occurrence}}
#' @return data.table

curate_co_occurrence <- function(
  co_occurrence_table,
  taxa_of_interest,
  number_of_treatments = 1
) {
  sub_co_occurrence <-
    co_occurrence_table[(
      co_occurrence_table[[2]] %in%
        taxa_of_interest |
        co_occurrence_table[[3]] %in% taxa_of_interest
    ), ]
  toi_table <- unique(cbind(rbindlist(
    list(sub_co_occurrence[, 1], sub_co_occurrence[, 1])
  ), rbindlist(
    list(sub_co_occurrence[, 2], sub_co_occurrence[, 3])
  )))
  toi_table <- toi_table[toi_table[[2]] %in% taxa_of_interest]
  if (number_of_treatments == "all") {
    number_of_treatments <- length(unique(sub_co_occurrence[[1]]))
  } else {
    number_of_treatments <- number_of_treatments
  }
  toi <- names(table(toi_table[[2]])[table(toi_table[[2]]) >=
                                        number_of_treatments])
  sub_co_occurrence <-
    co_occurrence_table[(co_occurrence_table[[2]] %in%
                            toi | co_occurrence_table[[3]] %in% toi), ]
  arranged_co_ocurrence <-
    data.table::data.table(arrange_co_occurrence_table(sub_co_occurrence, toi))

  data.table::setorder(arranged_co_ocurrence, "Treatment")
  
  return(arranged_co_ocurrence)
}
