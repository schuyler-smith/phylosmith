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

dendrogram_phyloseq <- function(
  phyloseq_obj,
  treatment = NULL,
  method = "bray",
  colors = "default"
){
  check_args(
    phyloseq_obj = phyloseq_obj,
    treatment    = treatment,
    dist_method  = method
  )
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0)
  metadata <- data.table::data.table(as(phyloseq_obj@sam_data, "data.frame"), 
    keep.rownames = "Sample")
  if(!is.null(treatment)){
    metadata[, paste(treatment, collapse = " x ") := do.call(paste, c(.SD)), 
      .SDcols = treatment]
    treatment_name <- paste(treatment, collapse = " x ")
  } else {treatment_name <- "Sample"}
  graph_data <-
    vegan::vegdist(t(phyloseq_obj@otu_table@.Data), method, na.rm = TRUE)
  graph_data[is.na(graph_data)] <- 0
  graph_data <- hclust(as.dist(graph_data), method = "complete")
  
  graph_data <- dendextend::as.ggdend(as.dendrogram(graph_data))
  graph_data$labels <- 
    merge(graph_data$labels, metadata, by.x="label", by.y="Sample")

  g <- ggplot(data = graph_data$segment) + 
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend, color=y), 
      show.legend=FALSE, size=0.8)
  if (is.null(treatment)){
    g <- g + 
      geom_text(data = graph_data$labels, 
        aes_string(x = "x", y = "y", label = "label"), color = "black",
        fontface = "bold", hjust = 1.05, size = 2.8)

  } else {
    sample_colors <- create_palette(length(unique(metadata[[treatment_name]])))
    g <- g + 
      geom_label(data = graph_data$labels, 
        aes_string(x = "x", y = "y", label = "label",
        fill = treatment_name), color = "white",
        label.padding = unit(0.2, "lines"), 
        fontface = "bold", hjust = 1.05, size = 2.5) +
      scale_fill_manual(values = sample_colors, 
        name = treatment_name) + 
      guides(
        fill = guide_legend(override.aes = aes(label = "")))
  }
  g <- g +
    scale_color_gradientn(colors =  viridis::viridis(10)) + 
    scale_y_continuous(
      limits = c(
        -max(sapply(graph_data$labels$label, FUN = function(x) {
          length(strsplit(as.character(x), "")[[1]]) 
        })) / 150, 
        max(graph_data$segment$yend))
    ) +  
    scale_x_discrete(
      limits = c(
        min(graph_data$segment$xend), 
        max(graph_data$segment$xend))
    )
  g <- g +
    coord_flip() + 
    labs(y = paste(tools::toTitleCase(method), "Distance", sep = " ")) + 
    theme_schuy("dend")

  return(g)
}
