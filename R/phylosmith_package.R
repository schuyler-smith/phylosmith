#' phylosmith-package
#'
#' phylosmith is a conglomeration of functions written to process and analyze \code{\link[phyloseq]{phyloseq-class}} objects.
#'
#'
#' \code{\link{FastCoOccur}} calculate Spearman-rank co-occurrence between taxa
#'
#' \code{\link{bootstrap_rho}} runs permutations of the otu_table to calculate a significant rho value
#'
#' \code{\link{combine_treatments}} combines multiple columns in meta-data into a single column
#'
#' \code{\link{curate_cooccurrence}} subsets the co-occurence table to specific taxa
#'
#' \code{\link{taxa_filter}} filter taxa by proportion of samples they are seen in
#'
#' \code{\link{find_common_taxa}} find taxa common to each treatment
#'
#' \code{\link{find_unique_taxa}} find taxa unique to each treatment
#'
#' \code{\link{merge_asvs}} combine ASVs to lowest common biological sequence
#'
#' \code{\link{nmds_phyloseq_ggplot}} create a ggplot object of the NMDS from a phyloseq object
#'
#' \code{\link{relative_abundance}} transform abundance data to relative abundance
#'
#'
#'
#' @name phylosmith
#' @author Schuyler Smith \email{sdsmith@@iastate.edu}
#' @docType package
#' @keywords package
NA
