# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @author Schuyler D. Smith
#' @title Assign rank values to a matrix.
#' @description Arranges the co-occurence table so that taxa of interest are
#' on the left.
#' @param count_matrix A \code{data.frame}
#' @return A rank matrix
assign_rank <- function(count_matrix) {
    .Call('_phylosmith_assign_rank', PACKAGE = 'phylosmith', count_matrix)
}

#' @author Schuyler D. Smith
#' @title Co-occurrence calculation
#' @description Calculate the pair-wise Spearman rank correlation.
#' @usage Correlation(X, Y,
#' cor_coef_cutoff, p_cutoff, method, ncores)
#' @param X An \code{otu_table} in the format from
#' \code{\link[phyloseq:otu_table]{phyloseq}}
#' @param Y An \code{otu_table} in the format from
#' \code{\link[phyloseq:otu_table]{phyloseq}}
#' @param cor_coef_cutoff \code{double} representing the minimum \code{rho-value}
#' accepted for the correlation to be returned.
#' @param p_cutoff \code{double} representing the maximum \code{p-value}
#' accepted for the correlation to be returned.
#' @param method Pearson, Spearman
#' @param ncores An \code{int} for how many cores to use to multithread the
#' calculations.
#' @return A \code{data.frame} with treatment, otu_1, otu_2, rho, p values.
#' @seealso \code{\link{co_occurrence}}
Correlation <- function(X, Y = matrix(1), cor_coef_cutoff = 0, p_cutoff = 1, method = "pearson", ncores = 1L) {
    .Call('_phylosmith_Correlation', PACKAGE = 'phylosmith', X, Y, cor_coef_cutoff, p_cutoff, method, ncores)
}

#' @author Schuyler D. Smith
#' @title Co-occurrence rho calculations
#' @description Calculates the pair-wise Spearman rank correlation without
#' testing for significance.
#' @usage permute_rho_Rcpp(count_matrix, permuted_matrix, method, ncores)
#' @param count_matrix An \code{otu_table} in the format from
#' \code{\link[phyloseq:otu_table]{phyloseq}}
#' @param permuted_matrix An \code{otu_table} in the format from
#' \code{\link[phyloseq:otu_table]{phyloseq}}
#' @param method pearson, spearman
#' @param ncores An \code{int} for how many cores to use to multithread the
#' calculations.
#' @return A \code{vector} with rho values for each pair-wise correlation.
#' @seealso \code{\link{permute_rho}}
permute_rho_Rcpp <- function(count_matrix, permuted_matrix, method = "pearson", ncores = 1L) {
    .Call('_phylosmith_permute_rho_Rcpp', PACKAGE = 'phylosmith', count_matrix, permuted_matrix, method, ncores)
}

#' @author Schuyler D. Smith
#' @title Arrange co-occurence table
#' @description Arranges the co-occurence table so that taxa of interest are
#' on the left.
#' @param co_occurrence_table A \code{data.frame} in the format from
#' \code{\link{co_occurrence}}.
#' @param taxa_of_interest A \code{vector} containing names of taxa to be
#' found in either OTU_1 or OTU_2.
#' @return A \code{data.frame} with treatment, otu_1, otu_2, rho, p values.
#' @seealso \code{\link{co_occurrence}}
arrange_co_occurrence_table <- function(co_occurrence_table, taxa_of_interest) {
    .Call('_phylosmith_arrange_co_occurrence_table', PACKAGE = 'phylosmith', co_occurrence_table, taxa_of_interest)
}

