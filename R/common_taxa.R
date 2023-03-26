#' Find taxa shared between treatments of a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and finds which taxa
#' are shared between all of the specified treatments from data in the
#' \code{\link[phyloseq:sample_data]{sample_data()}}), or every sample in the
#' dataset.
#' @useDynLib phylosmith
#' @usage common_taxa(phyloseq_obj, treatment = NULL, subset = NULL, n = "all")
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param treatment Column name as a \code{string} or \code{numeric} in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any
#' samples that to not contain this factor. This can be a vector of multiple
#' factors to subset on.
#' @param n Number of treatment groups that need to share the taxa to be
#' considered a common taxa.
#' @seealso \code{\link{unique_taxa}}
#' @export
#' @return vector
#' @examples common_taxa(soil_column, treatment = "Treatment",
#' subset = "Amended", n = "all")

common_taxa <- function(
    phyloseq_obj,
    treatment = NULL,
    subset = NULL,
    n = "all"
) {
check_args(
  phyloseq_obj = phyloseq_obj,
  sam_data     = phyloseq_obj,
  treatment    = treatment,
  subset       = subset,
  n            = n
)
  phyloseq_obj <-
    taxa_filter(phyloseq_obj, treatment, subset = subset)
  treatment_name <- paste(treatment, collapse = sep)
  treatment_classes <- unique(phyloseq_obj@sam_data[[treatment_name]])
  if (n == "all") n <- length(treatment_classes)
  if (n > length(treatment_classes)) {
    stop(
"`n` must be either 'all' or a numeric value less than the number of treatments 
  being compared", call. = FALSE
    )
  }
  seen_taxa <- lapply(
    treatment_classes,
    FUN = function(trt) {
        phyloseq::taxa_names(taxa_filter(
        phyloseq_obj,
        treatment = treatment,
        subset = trt,
        frequency = 0
      ))
    }
  )
  if (n == 0) {
    seen_taxa <-
      tryCatch(
        phyloseq::taxa_names(taxa_filter(phyloseq_obj, frequency = .99)),
      error = function(e) {
        stop(
"No taxa seen in every sample.
  To get a list of taxa seen in a certain proportion of samples use:
  taxa_names(taxa_filter(phyloseq_obj, frequency = x))"
          )
        }
      )
    n <- 1
  }
  names(seen_taxa) <- treatment_classes
  taxa_counts <- table(unlist(seen_taxa))
  shared_taxa <- names(taxa_counts[taxa_counts == round(n)])
  
  return(shared_taxa)
}
