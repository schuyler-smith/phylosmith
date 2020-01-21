#' Computes the correlation of numerical variables with taxa
#' Function from the phylosmith-package.
#'
#' Computes the correlation of numerical variables with taxa
#' @useDynLib phylosmith
#' @usage variable_correlation(phyloseq_obj, variables, treatment = NULL,
#'  subset = NULL, classification = NULL, method = 'spearman', cores = 1)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object.
#' @param variables Numericla factors within the in the
#' \code{\link[phyloseq:sample_data]{sample_data}} to correlate with the
#' abundance data.
#' @param treatment Column name as a \code{string} or \code{numeric} in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any
#' samples that to not contain this factor. This can be a vector of multiple
#' factors to subset on.
#' @param classification Column name as a \code{string} or \code{numeric} in
#' the \code{\link[phyloseq:tax_table]{tax_table}} for the factor to
#' conglomerate by.
#' @param method Which correlation method to calculate, "pearson", "spearman".
#' @param cores \code{numeric} Number of CPU cores to use for the pair-wise
#' permutations. Default (0) uses max cores available. Parallelization not
#' available for systems running MacOS without openMP configuration.
#' @importFrom parallel detectCores
#' @keywords nonparametric
#' @seealso \code{\link{permute_rho}} \code{\link{phylosmith}}
#' @export
#' @return data.table
#' @examples
#' variable_correlation(soil_column, variables = 'Day',
#' treatment = c('Matrix', 'Treatment'), subset = 'Amended',
#' classification = 'Phylum', method = 'spearman', cores = 1)
# sourceCpp("src/correlations_Rcpp.cpp")
variable_correlation <-
  function(phyloseq_obj,
           variables,
           treatment = NULL,
           subset = NULL,
           classification = NULL,
           method = 'spearman',
           cores = 1) {
    if (!inherits(phyloseq_obj, "phyloseq")) {
      stop("`phyloseq_obj` must be a phyloseq-class
          object", call. = FALSE)
    }
    if (is.null(access(phyloseq_obj, 'sam_data'))) {
      stop("`phyloseq_obj` must contain sample_data()
          information",
           call. = FALSE)
    }
    if (is.null(access(phyloseq_obj, 'tax_table'))) {
      stop("`phyloseq_obj` must contain tax_table()
          information",
           call. = FALSE)
    }
    treatment <- check_index_treatment(phyloseq_obj, treatment)
    if (!(is.null(treatment)) &
        any(!(treatment %in% colnames(access(
          phyloseq_obj, 'sam_data'
        ))))) {
      stop(
        "`treatment` must be at least one column
          name, or index, from the sample_data()",
        call. = FALSE
      )
    }
    tryCatch({
      variables <- check_index_treatment(phyloseq_obj, variables)
    }, error = function(e){
      stop(
        "`variables` must be at least one column name, or index, from the sample_data() that contains numeric data.",
        call. = FALSE)})
    if (!(is.null(variables)) &
        any(!(variables %in% colnames(access(
          phyloseq_obj, 'sam_data'
        ))))) {
      stop(
        "`variables` must be at least one column
          name, or index, from the sample_data()",
        call. = FALSE
      )
    }
    match.arg(method, c("pearson", "kendall", "spearman"))
    phyloseq_obj <- taxa_filter(phyloseq_obj, treatment = treatment, subset = subset)
    treatment_name <- paste(treatment, collapse = sep)
    treatment_classes <- as.character(unique(phyloseq_obj@sam_data[[treatment_name]]))
    if(is.null(treatment)){treatment_classes <- list(NULL)}
    if(!(is.null(classification))){
      phyloseq_obj <- conglomerate_taxa(phyloseq_obj, classification, hierarchical = FALSE)
    }
    correlations <- data.table()
    for(k in treatment_classes){
      phyloseq_obj_subset <- taxa_filter(phyloseq_obj, treatment, k, drop_samples = TRUE)
      treatment_correlations <- Correlation(
        X = phyloseq_obj_subset@otu_table,
        Y = apply(as.matrix(phyloseq_obj_subset@sam_data[,variables]), 2, as.numeric),
        method = method
      )
      treatment_correlations[['X']] <- rownames(phyloseq_obj_subset@otu_table)[treatment_correlations[['X']]]
      treatment_correlations[['Y']] <- colnames(phyloseq_obj_subset@sam_data[,variables])[treatment_correlations[['Y']]]
      if(length(treatment_classes) > 1){
        treatment_correlations <- cbind(Treatment = k, treatment_correlations)
      }
      correlations <- rbind(correlations, treatment_correlations)
    }
    return(correlations)
  }
