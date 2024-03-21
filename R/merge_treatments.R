#' Combine meta-data columns. Function from the phylosmith-package.
#'
#' Combines multiple columns of a \code{\link[phyloseq]{phyloseq-class}}
#' object \code{\link[phyloseq:sample_data]{sample_data}} into a
#' single-variable column.
#' @useDynLib phylosmith
#' @usage merge_treatments(phyloseq_obj, treatment)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample.
#' @param treatment A vector of any number of column names as \code{string}s or \code{numeric}s
#' in the \code{\link[phyloseq:sample_data]{sample_data}} that are to be
#' combined.
#' @export
#' @return phyloseq-object
#' @examples merge_treatments(soil_column, c("Matrix", "Treatment", "Day"))

merge_treatments <- function(
  phyloseq_obj, 
  treatment
) {
  if(length(treatment) <= 1) return(phyloseq_obj)
  check_args(
    phyloseq_obj = phyloseq_obj,
    sam_data     = phyloseq_obj,
    treatment    = treatment
  )
  treatment_classes <- data.table::setDT(as(phyloseq_obj@sam_data[,
        colnames(phyloseq_obj@sam_data) %in% treatment], "data.frame"))
  treatment_name <- paste(treatment, collapse = sep)
  order <- apply(
    eval(parse(text = paste0(
      "expand.grid(",
      paste0(
        paste0(
          "levels(factor(phyloseq::access(phyloseq_obj, 'sam_data')[['",
          treatment,
          "']]))",
          collapse = ", "
        )
      ), ")"
    ))),
    1,
    FUN = function(combination) {
      paste0(combination, collapse = " ")
    }
  )
  eval(parse(
    text = paste0(
      "treatment_classes[, '",
      treatment_name,
      "' := as.character(paste(",
      paste(treatment, collapse = ", "),
      ", sep = ' '), by = treatment_classes)]"
    )
  ))
  phyloseq::sample_data(phyloseq_obj)[[treatment_name]] <- 
    factor(treatment_classes[[treatment_name]], levels = order)
    
  return(phyloseq_obj)
}