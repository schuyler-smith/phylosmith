#' Normalizes abundance data in an \code{otu_table} using library size.
#' Function from the phylosmith-package.
#'
#' Transform abundance data to equal counts using library size and the geometric mean, 
#' i.e. proportional data.
#' @useDynLib phylosmith
#' @usage library_size(phyloseq_obj)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}} with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @export
#' @return phyloseq-object
#' @examples library_size(soil_column)

library_size <- function(phyloseq_obj) {
  if (!inherits(phyloseq_obj, "phyloseq")) {
    stop("`phyloseq_obj` must be a phyloseq-class
            object", call. = FALSE)
  }
  phyloseq_obj <- check_TaR(phyloseq_obj)
  phyloseq_obj <- taxa_filter(phyloseq_obj)
  phyloseq_obj <- phylosmith::taxa_filter(phyloseq_obj)
  totCounts <- sample_sums(phyloseq_obj)
  normFactor <- totCounts/exp(mean(log(totCounts)))
  otu_table(phyloseq_obj) <- round(otu_table(phyloseq_obj)/rep(normFactor, each = (ntaxa(phyloseq_obj))))
  return(phyloseq_obj)
}

#' Create a ggplot object of the dendrogram from a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and
#' plots the dendrogram of samples colored by treatment.
#' @useDynLib phylosmith
#' @usage dendrogram_phyloseq(phyloseq_obj, treatment, method = 'bray',
#' colors = 'default')
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param treatment Column name as a string or number in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @param method the distance measurement algorithm to use, match to
#' "euclidean", "manhattan", "canberra", "clark", "bray", "kulczynski",
#' "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup",
#' "binomial", "chao", "cao" or "mahalanobis".
#' @param colors Name of a color set from the
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palette of R-accepted
#' colors.
#' @importFrom vegan vegdist
#' @importFrom stats hclust
#' @export
#' @return ggplot-object
#' @examples dendrogram_phyloseq(soil_column, c('Matrix', 'Treatment'))

dendrogram_phyloseq <- function(phyloseq_obj,
                          treatment = NULL,
                          method = "bray",
                          colors = 'default'
                          ){
  
  treatment <- check_index_treatment(phyloseq_obj, treatment)
  if (any(!(treatment %in% colnames(access(phyloseq_obj, "sam_data"))))) {
    stop("`treatment` must be at least one column\n        name, or index, from the sample_data()",
         call. = FALSE)
  }
  method <- match.arg(method, methods)
  
  graph_data <- taxa_filter(phyloseq_obj, treatment, frequency = 0)
  metadata <- data.table(as(graph_data@sam_data, 'data.frame'), keep.rownames = 'Sample')
  if(!is.null(treatment)){
    metadata[, paste(treatment, collapse = " x ") := do.call(paste,c(.SD)), .SDcols=treatment]
    treatment <- paste(treatment, collapse = " x ") #create column for treatment(s) selected for facet
  } else {treatment <- "Sample"}
  graph_data <- vegan::vegdist(t(graph_data@otu_table@.Data),method, na.rm = TRUE)
  graph_data[is.na(graph_data)] <- 0
  graph_data <- hclust(as.dist(graph_data), method = 'complete')
  
  graph_data <- ggdendro::dendro_data(graph_data)
  sample_colors <- create_palette(length(unique(metadata[[treatment]])))
  labels <- merge(ggdendro::label(graph_data), metadata, by.x='label', by.y='Sample')
  
  ggplot(data = ggdendro::segment(graph_data)) + 
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend, color=y), show.legend=FALSE, size=0.8) +
    geom_label(data = labels, aes(x=x, y=y, label=label, fill=get(treatment)), 
               color = "white", label.padding = unit(0.15, "lines"), fontface = "bold", 
               hjust=1.05, size=2.5) +
    scale_fill_manual(values = sample_colors, name = treatment) +
    guides(fill = guide_legend(override.aes = aes(label = ""))) +
    scale_color_gradientn(colors =  viridis::viridis(10)) + 
    coord_flip() +
    scale_y_continuous(limits = c(-max(sapply(sapply(labels$label, strsplit, ""), length))/150, 
                                  max(ggdendro::segment(graph_data)$yend))) + 
    theme_minimal() + 
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.border = element_blank(),
          legend.position = "top")
}
