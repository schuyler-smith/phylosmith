#' Filter taxa based on proportion of samples they are observed in. Function from the phylosmith-package.
#'
#' This function takes a phyloseq object and finds which taxa are seen in a given proportion of samples, either in the entire dataset, by treatment, or a particular treatment of interest.
#' @useDynLib phylosmith
#' @usage taxa_filter(phyloseq_obj, treatment = NULL, frequency = 0,
#' subset = NULL, below = FALSE, drop_samples = FALSE)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package.
#' @param treatment Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}}. Function then checks if taxa seen frequency number of times in each treatment. If multiple \code{sample_data} columns are given, they will be appended to the \code{sample_data} as one column with '.' separating each.
#' @param frequency The minimum proportion of samples the taxa is seen in.
#' @param subset If taxa not needed to be seen in all \code{treatment}, then can check only one particular treatment subset, this works for multiple treatment inputs.
#' @param below Does frequency define the minimum or maximum, should the presence fall below frequency or not.
#' @param drop_samples Should the function remove samples that that are empty after removing taxa filtered by frequency.
#' @keywords manip
#' @import data.table
#' @import phyloseq
#' @examples
#' data(mock_phyloseq)
#' taxa_filter(mock_phyloseq, frequency = 0.3)
#' taxa_filter(mock_phyloseq, frequency = 0.3, treatment = "day")
#' taxa_filter(mock_phyloseq, frequency = 0.3, treatment = 3)
#' taxa_filter(mock_phyloseq, frequency = 0.3, treatment = c("day", "treatment"))
#' taxa_filter(mock_phyloseq, frequency = 0.3, treatment = c("day", "treatment"),
#' subset = "soil")
#' taxa_filter(mock_phyloseq, frequency = 0.3, treatment = c("day", "treatment"),
#' subset = c("5","soil"))
#' taxa_filter(mock_phyloseq, frequency = 0.3, treatment = c("day", "treatment"),
#' subset = "5.soil")
#' @export

taxa_filter <- function(phyloseq_obj, treatment = NULL, frequency = 0, subset = NULL, below = FALSE, drop_samples = FALSE){
  # data("mock_phyloseq")
  # phyloseq_obj = mock_phyloseq; frequency = 0; treatment = c("treatment", "day"); subset = "5"; below = FALSE; drop_samples = FALSE
  if(!(is.null(treatment))){
    if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}

    phyloseq_obj <- combine_treatments(phyloseq_obj, treatment)
    treatment_name <- paste(treatment, collapse = ".")
    Treatment_Groups <- unique(phyloseq_obj@sam_data[[treatment_name]])

    if(below == TRUE){
      phyloseq_obj <- do.call(merge_phyloseq,
                              apply(array(Treatment_Groups), 1, FUN = function(group){
                                sub_phy <- eval(parse(text=paste0("subset_samples(phyloseq_obj, ",treatment_name," == '",group,"')")))
                                cutoff <- floor(ncol(sub_phy@otu_table) * frequency)
                                sub_phy <- filter_taxa(sub_phy, function(x){sum(x != 0, na.rm = TRUE) <= cutoff}, TRUE)
                                return(sub_phy)}))
    } else {
      phyloseq_obj <- do.call(merge_phyloseq,
                              apply(array(Treatment_Groups), 1, FUN = function(group){
                                sub_phy <- eval(parse(text=paste0("subset_samples(phyloseq_obj, ",treatment_name," == '",group,"')")))
                                cutoff <- floor(ncol(sub_phy@otu_table) * frequency)
                                sub_phy <- filter_taxa(sub_phy, function(x){sum(x != 0, na.rm = TRUE) >= cutoff}, TRUE)
                                if(sum(taxa_sums(sub_phy), na.rm = TRUE) != 0){sub_phy <- filter_taxa(sub_phy, function(x){sum(x, na.rm = TRUE) != 0}, TRUE)}
                                return(sub_phy)}))
    }
  } else {
    cutoff <- floor(ncol(phyloseq_obj@otu_table) * frequency)
    if(below == TRUE){
      phyloseq_obj <- filter_taxa(phyloseq_obj, function(x){sum(x != 0, na.rm = TRUE) <= cutoff}, TRUE)
    } else {
      phyloseq_obj <- filter_taxa(phyloseq_obj, function(x){sum(x != 0, na.rm = TRUE) >= cutoff}, TRUE)
      phyloseq_obj <- filter_taxa(phyloseq_obj, function(x){sum(x, na.rm = TRUE) != 0}, TRUE)
      }
  }
  if(!(is.null(subset))){
    treatments <- eval(parse(text=paste0("unique(phyloseq_obj@sam_data$", paste0(treatment_name), ")")))
    treatments <- eval(parse(text=paste0('treatments[grepl("', paste0(subset), '", treatments)]')))
    phyloseq_obj <- prune_samples(apply(phyloseq_obj@sam_data,1,FUN = function(x){any(x[c(treatment, treatment_name)] %in% treatments)}), phyloseq_obj)
    # apply(sample_data(phyloseq_obj),1,FUN = function(x){(x[c(treatment, treatment_name)] %in% subset)}) == length(subset)
    phyloseq_obj <- filter_taxa(phyloseq_obj, function(x){sum(x, na.rm = TRUE) > 0}, TRUE)
  }
  if(drop_samples == TRUE){
    phyloseq_obj <- prune_samples(sample_sums(phyloseq_obj) > 0, phyloseq_obj)
  }
  return(phyloseq_obj)
}


#' Combine meta-data columns. Function from the phylosmith-package.
#'
#' Combines multiple columns of a \code{\link[phyloseq]{phyloseq-class}} object \code{\link[phyloseq:sample_data]{sample_data}} into a single-variable column.
#' @useDynLib phylosmith
#' @usage combine_treatments(phyloseq_obj, treatments, sep = '.')
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package.
#' @param treatments A vector of column names or numbers in the \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param sep Delimiter to separate the treatments.
#' @keywords manip
#' @import data.table
#' @import phyloseq
#' @examples
#' data(mock_phyloseq)
#' combine_treatments(mock_phyloseq, treatments = c("treatment", "day"))@sam_data
#' combine_treatments(mock_phyloseq, treatments = c("treatment", "day"), sep = '_')@sam_data
#' @export

combine_treatments <- function(phyloseq_obj, treatments, sep = '.'){
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


#' Transform abundance data in an \code{otu_table} to relative abundance, sample-by-sample. Function from the phylosmith-package.
#'
#' Transform abundance data into relative abundance, i.e. proportional data. This is an alternative method of normalization and may not be appropriate for all datasets, particularly if your sequencing depth varies between samples.
#' @useDynLib phylosmith
#' @usage relative_abundance(phyloseq_obj)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package.
#' @keywords manip
#' @import phyloseq
#' @import RcppArmadillo
#' @import RcppParallel
#' @examples
#' data(mock_phyloseq)
#' relative_abundance(mock_phyloseq)
#' @export

relative_abundance <- function(phyloseq_obj){
  options(warnings=-1)
  phyloseq_obj <- transform_sample_counts(phyloseq_obj, function(sample) sample/sum(sample))
  return(phyloseq_obj)
}
