#' Find Generalists
#'
#' This function takes a phyloseq object and finds which taxa are seen in a given proportion of samples, either in the entire dataset, by treatment, or a particular treatment of interest.
#' @useDynLib phylosmith
#' @usage find_generalists(phyloseq_obj, frequency = 0, treatment = NULL,
#' subset = NULL, below = FALSE, drop_samples = FALSE)
#' @param phyloseq_obj A phyloseq-class object created with the phyloseq package.
#' @param frequency The minimum proportion of samples the taxa is seen in.
#' @param treatment Name(s) or column number(s) in the sample_data(). Function then checks if taxa seen in frequency in each treatment. If multiple sample_data() columns are given, they will be appended to the sample_data() as one column with '.' separating each.
#' @param subset If taxa not needed to be seen in all `treatment``, then can check only one particular treatment subset, this works for multiple treatment inputs.
#' @param below Does frequency define the minimum or maximum, should the presence fall below frequency or not.
#' @param drop_samples Should the function remove samples that that are empty after removing taxa filtered by frequency.
#' @keywords generalists frequency phyloseq phylosmith
#' @import phyloseq
#' @import data.table
#' @export
#' @examples
#' data(mock_phyloseq)
#' find_generalists(mock_phyloseq, frequency = 0.3)
#' find_generalists(mock_phyloseq, frequency = 0.3, treatment = "day")
#' find_generalists(mock_phyloseq, frequency = 0.3, treatment = 3)
#' find_generalists(mock_phyloseq, frequency = 0.3, treatment = c("day", "treatment"))
#' find_generalists(mock_phyloseq, frequency = 0.3, treatment = c("day", "treatment"), subset = "soil")
#' find_generalists(mock_phyloseq, frequency = 0.3, treatment = c("day", "treatment"), subset = c("5","soil"))
#' find_generalists(mock_phyloseq, frequency = 0.3, treatment = c("day", "treatment"), subset = "5.soil")

find_generalists <- function(phyloseq_obj, frequency = 0, treatment = NULL, subset = NULL, below = FALSE, drop_samples = FALSE){
  # data("mock_phyloseq")
  # phyloseq_obj = mock_phyloseq; frequency = 0; treatment = c("treatment", "day"); subset = "5"; below = FALSE; drop_samples = FALSE
  if(!(is.null(treatment))){
    if(is.numeric(treatment)){treatment <- colnames(sample_data(phyloseq_obj)[,treatment])}

    Treatment_Groups <- setDT(as(sample_data(phyloseq_obj)[,colnames(sample_data(phyloseq_obj)) %in% treatment], "data.frame"))
    # Treatment_Groups[, Treatment_Group := .GRP, by = Treatment_Groups]
    treatment_name <- paste(treatment, collapse = ".")
    eval(parse(text=paste0("Treatment_Groups[, '",treatment_name,"' := as.character(paste(", paste(treatment, collapse = ", "), ", sep = '.'), by = Treatment_Groups)]")))
    sample_data(phyloseq_obj)[[treatment_name]] <- Treatment_Groups[[treatment_name]]
    Treatment_Groups <- unique(Treatment_Groups[[treatment_name]])
    # phyloseq_obj <- phyloseq(otu_table(phyloseq_obj), tax_table(phyloseq_obj), sample_data(phyloseq_obj))

    if(below == TRUE){
      phyloseq_obj <- do.call(merge_phyloseq,
                              apply(array(Treatment_Groups), 1, FUN = function(group){
                                sub_phy <- eval(parse(text=paste0("subset_samples(phyloseq_obj, ",treatment_name," == '",group,"')")))
                                cutoff <- floor(ncol(otu_table(sub_phy)) * frequency)
                                sub_phy <- filter_taxa(sub_phy, function(x){sum(x != 0) <= cutoff}, TRUE)
                                return(sub_phy)}))
    } else {
      phyloseq_obj <- do.call(merge_phyloseq,
                              apply(array(Treatment_Groups), 1, FUN = function(group){
                                sub_phy <- eval(parse(text=paste0("subset_samples(phyloseq_obj, ",treatment_name," == '",group,"')")))
                                cutoff <- floor(ncol(otu_table(sub_phy)) * frequency)
                                sub_phy <- filter_taxa(sub_phy, function(x){sum(x != 0) >= cutoff}, TRUE)
                                if(sum(taxa_sums(sub_phy)) != 0){sub_phy <- filter_taxa(sub_phy, function(x){sum(x) != 0}, TRUE)}
                                return(sub_phy)}))
    }
  } else {
    cutoff <- floor(ncol(otu_table(phyloseq_obj)) * frequency)
    if(below == TRUE){
      phyloseq_obj <- filter_taxa(phyloseq_obj, function(x){sum(x != 0) <= cutoff}, TRUE)
    } else {
      phyloseq_obj <- filter_taxa(phyloseq_obj, function(x){sum(x != 0) >= cutoff}, TRUE)
      phyloseq_obj <- filter_taxa(phyloseq_obj, function(x){sum(x) != 0}, TRUE)
      }
  }
  if(!(is.null(subset))){
    phyloseq_obj <- prune_samples(apply(sample_data(phyloseq_obj),1,FUN = function(x){any(x[c(treatment, treatment_name)] %in% subset)}), phyloseq_obj)
    # apply(sample_data(phyloseq_obj),1,FUN = function(x){(x[c(treatment, treatment_name)] %in% subset)}) == length(subset)
    phyloseq_obj <- filter_taxa(phyloseq_obj, function(x){sum(x) > 0}, TRUE)
  }
  if(drop_samples == TRUE){
    phyloseq_obj <- prune_samples(sample_sums(phyloseq_obj) > 0, phyloseq_obj)
  }
  return(phyloseq_obj)
}
