#' Combine meta-data columns.
#'
#' Combines multiple columns of a \code{\link[phyloseq]{phyloseq-class}} object \code{\link[phyloseq:sample_data]{sample_data}} into a single-variable column.
#' @aliases combine_treatments
#' @useDynLib phylosmith
#' @usage concatenate_treatments(phyloseq_obj, treatments, sep = '.')
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package.
#' @param treatments A vector of column names or numbers in the \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param sep Delimiter to separate the treatments.
#' @keywords manip
#' @import data.table
#' @import phyloseq
#' @export
#' @examples
#' data(mock_phyloseq)
#' concatenate_treatments(mock_phyloseq, treatments = c("treatment", "day"))@sam_data
#' concatenate_treatments(mock_phyloseq, treatments = c("treatment", "day"), sep = '_')@sam_data

concatenate_treatments <- function(phyloseq_obj, treatments, sep = '.'){
  # data("mock_phyloseq")
  # phyloseq_obj = mock_phyloseq; sep='.'; treatments = c("treatment", "day");
  if(!(is.null(treatments))){
    if(is.numeric(treatments)){treatments <- colnames(phyloseq_obj@sam_data[,treatments])}

    Treatment_Groups <- setDT(as(phyloseq_obj@sam_data[,colnames(phyloseq_obj@sam_data) %in% treatments], "data.frame"))
    treatment_name <- paste(treatments, collapse = sep)
    eval(parse(text=paste0("Treatment_Groups[, '",treatment_name,"' := as.character(paste(", paste(treatments, collapse = ", "), ", sep = '",sep,"'), by = Treatment_Groups)]")))
    phyloseq_obj@sam_data[[treatment_name]] <- Treatment_Groups[[treatment_name]]
  }
  return(phyloseq_obj)
}
