#' Merge amplicon sequence variants (asv).
#' Function from the phylosmith-package.
#'
#' Takes any number of phyloseq objects that have
#' \code{\link[phyloseq:taxa_names]{taxa_names}} declared as true biological
#' sequences or Amplicon Sequence Variants (ASVs). If a dataset has longer
#' reads than the others, it will look for the shorter sequence within the
#' longer and redeclare the ASV.
#' @useDynLib phylosmith
#' @usage merge_asvs(...)
#' @param ... Any number of \code{\link[phyloseq]{phyloseq-class}} objects
#' created with the \link[=phyloseq]{phyloseq} package (must have Amplicon
#' Sequence Variants for calling taxa).
#' @keywords manip
#' @import RcppArmadillo
#' @import RcppParallel
#' @return phyloseq-object
#' @examples
#' phylosmith:::merge_asvs(mock_phyloseq, mock_phyloseq_2)

# sourceCpp("src/rcpp_seq_match.cpp")

merge_asvs <- function(...){
    options(warn=1)
    phyloseq_objects <- list(...)
    p_objects_names <- sapply(substitute(list(...))[-1], deparse)
    names(phyloseq_objects) <- p_objects_names

    asvs <- lapply(phyloseq_objects, otu_table)
    read_size_order <- order(unlist(lapply(lapply(asvs,
        FUN = function(x){rownames(x)[1]}),
        FUN = function(seq){mean(length(strsplit(seq, "")[[1]]))})))
    asvs <- asvs[read_size_order]

    pairs <- utils::combn(names(asvs), m = 2)
    for(run in seq(dim(pairs)[2])){
        asvs[[pairs[,run][2]]] <- match_sequences(asvs[[pairs[,run][1]]],
            asvs[[pairs[,run][2]]])
    }
    all_taxa <- unique(eval(parse(text=paste0("rbind(",
        paste0("tax_table(",p_objects_names,")", collapse = ", "), ")"))))
    for(name in names(asvs)){
        phyloseq_objects[[name]] <- phyloseq(
            otu_table(asvs[[name]], taxa_are_rows = TRUE),
            tax_table(as.matrix(all_taxa[rownames(all_taxa) %in%
                rownames(asvs[[name]]),])),
            sample_data(sample_data(phyloseq_objects[[name]])))
    }
    return(phyloseq_objects)
}












