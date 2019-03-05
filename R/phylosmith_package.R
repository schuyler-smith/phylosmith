#' phylosmith-package
#'
#' phylosmith is a conglomeration of functions written to process and analyze \code{\link[phyloseq]{phyloseq-class}} objects.
#'
#' \code{\link{tsne_phyloseq_ggplot}} create a ggplot object of the t-SNE from a phyloseq object
#'
#' \code{\link{taxa_abundance_bars_ggplot}} create a ggplot object of the abundance of taxa in each sample
#'
#' \code{\link{phylogeny_bars_ggplot}} create a ggplot barplot object of the compositons of each sample at a taxonomic level
#'
#' \code{\link{network_phyloseq}} creates a ggplot object of the co-occurrence network of taxa
#'
#' \code{\link{abundance_lines_ggplot}} create a ggplot object of the abundance data as a line graph
#'
#' \code{\link{abundance_heatmap_ggplot}} create a ggplot object of the heatmaps of the abndance table

#' \code{\link{order_phyloseq_metadata}} sets the orders of the factors in a sample_data column (for ordering graphs)
#'
#' \code{\link{co_occurrence}} calculate Spearman-rank co-occurrence between taxa
#'
#' \code{\link{bootstrap_rho}} runs permutations of the otu_table to calculate a significant rho value
#'
#' \code{\link{combine_treatments}} combines multiple columns in meta-data into a single column
#'
#' \code{\link{curate_co_occurrence}} subsets the co-occurence table to specific taxa
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
#' @author Schuyler Smith \email{schuyler.smith@@iastate.edu}
#' @docType package
#' @keywords package
NA


# OBO_filepath <- '~/Downloads/card-ontology/aro.obo'
# CARD <- processOBO(OBO_filepath)
# load('R/sysdata.rda')
# sep = '_'
# save(CARD, sep, file = 'R/sysdata.rda')



