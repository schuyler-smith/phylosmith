#' Converts numeric values to column names in sample_data.
#'
#' Converts numeric values to column names in sample_data.
#' @useDynLib phylosmith
#' @usage check_numeric_treatment(phyloseq_obj, treatment)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with information about each sample, and it must contain \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each taxa/gene.
#' @param treatment Column name as a \code{string} or \code{numeric} in the \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of multiple columns and they will be combined into a new column.
#'

check_numeric_treatment <- function(phyloseq_obj, treatments){
  sapply(treatments, FUN = function(treatment){
    if(is.numeric(treatment)){
      return(colnames(phyloseq_obj@sam_data[,treatment]))
    } else {return(treatment)}
  })
}

#' Converts numeric values to column names in tax_table
#'
#' Converts numeric values to column names in tax_table.
#' @useDynLib phylosmith
#' @usage check_numeric_treatment(phyloseq_obj, treatment)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with information about each sample, and it must contain \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each taxa/gene.
#' @param classification Column name as a \code{string} or \code{numeric} in the \code{\link[phyloseq:tax_table]{tax_table}} for the factor to conglomerate by.
#'

check_numeric_classification <- function(phyloseq_obj, classification){
  sapply(classification, FUN = function(class){
    if(is.numeric(class)){
      return(colnames(phyloseq_obj@tax_table[,class]))
    } else {return(class)}
  })
}
