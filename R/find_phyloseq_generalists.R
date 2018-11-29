#' Find Generalists
#'
#' This function takes a phyloseq object and finds which taxa are seen in a given proportion of samples, either in the entire dataset, by treatment, or a particular treatment of interest.
#' @useDynLib phyloschuyler
#' @usage find_generalists(phyloseq_obj, frequency = 0, treatments = NULL,
#' subset = NULL, below = FALSE, drop_samples = FALSE)
#' @param phyloseq_obj An object created with the phyloseq package.
#' @param frequency The minimum proportion of samples the taxa is seen in.
#' @param treatments Name(s) or column number(s) in the sample_data(). Function then checks if taxa seen in frequency in each treatment.
#' @param subset If taxa not needed to be seen in all treatments, then can check only one particular treatment subset.
#' @param below Does frequency define the minimum or maximum, should the presence fall below frequency or not.
#' @param drop_samples Should the function remove samples that that are empty after removing taxa filtered by frequency.
#' @keywords generalists frequency phyloseq phyloschuyler
#' @import phyloseq
#' @import data.table
#' @export
#' @examples
#' data(mock_phyloseq)
#' find_generalists(mock_phyloseq, frequency = 0.3)
#' find_generalists(mock_phyloseq, frequency = 0.3, treatments = "day")
#' find_generalists(mock_phyloseq, frequency = 0.3, treatments = 3)
#' find_generalists(mock_phyloseq, frequency = 0.3, treatments = c("day", "treatment"))
#' find_generalists(mock_phyloseq, frequency = 0.3, treatments = c("day", "treatment"), subset = "5-soil")

find_generalists <- function(phyloseq_obj, frequency = 0, treatments = NULL, subset = NULL, below = FALSE, drop_samples = FALSE){
  # phyloseq_obj <- prune_samples(sample_sums(phyloseq_obj) > 0, phyloseq_obj)
  if(!(is.null(treatments))){
    if(is.numeric(treatments)){treatment <- colnames(sample_data(phyloseq_obj)[,treatments])
    } else {treatment <- eval(parse(text=paste0("colnames(sample_data(phyloseq_obj)[,c('", paste0(treatments, collapse = "', '"), "')])")))
    }

    if(!(is.null(subset)) & length(treatment) == 1){
      phyloseq_obj <- eval(parse(text=paste0("subset_samples(phyloseq_obj, ", paste0(treatment), " %in% c('", paste0(subset, collapse = "', '"), "'))")))
      phyloseq_obj <- filter_taxa(phyloseq_obj, function(x){sum(x) > 0}, TRUE)
      }

    Treatment_Groups <- setDT(as(sample_data(phyloseq_obj)[,colnames(sample_data(phyloseq_obj)) %in% treatment], "data.frame"))
    # Treatment_Groups[, Treatment_Group := .GRP, by = Treatment_Groups]
    eval(parse(text=paste0("Treatment_Groups[, Treatment_Group := paste(", paste(treatment, collapse = ", "), ", sep = '-'), by = Treatment_Groups]")))
    sample_data(phyloseq_obj)[,'Treatment_Group'] <- Treatment_Groups[,"Treatment_Group"]
    Treatment_Groups <- unique(Treatment_Groups[,"Treatment_Group"])
    # phyloseq_obj <- phyloseq(otu_table(phyloseq_obj), tax_table(phyloseq_obj), sample_data(phyloseq_obj))
    if(below == TRUE){
      phyloseq_obj <- do.call(merge_phyloseq,
                              apply(array(Treatment_Groups), 1, FUN = function(group){
                                sub_phy <- eval(parse(text=paste0("subset_samples(phyloseq_obj, Treatment_Group == '",group,"')")))
                                cutoff <- floor(ncol(otu_table(sub_phy)) * frequency)
                                sub_phy <- filter_taxa(sub_phy, function(x){sum(x != 0) <= cutoff}, TRUE)
                                return(sub_phy)}))
    } else {
      phyloseq_obj <- do.call(merge_phyloseq,
                              apply(array(Treatment_Groups), 1, FUN = function(group){
                                sub_phy <- eval(parse(text=paste0("subset_samples(phyloseq_obj, Treatment_Group == '",group,"')")))
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
  if(drop_samples == TRUE){
    phyloseq_obj <- prune_samples(sample_sums(phyloseq_obj) > 0, phyloseq_obj)
  }
  return(phyloseq_obj)
}
