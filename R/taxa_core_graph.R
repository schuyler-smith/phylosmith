#' Create graph of the core taxa seen in phyloseq-object over a range of
#' abundance and smaple-frequency values.  Function from the phylosmith-package.
#'
#' Inputs a phyloseq object and finds which taxa are seen in a
#' given proportion of samples at a minimum relative abundance, either in the
#' entire dataset, or by treatment, over a range of values. Then draws the
#' distribution.
#' @useDynLib phylosmith
#' @usage taxa_core_graph(phyloseq_obj, treatment = NULL, subset = NULL,
#' frequencies = seq(0.1, 1, 0.1), abundance_thresholds = seq(0.01, 1, 0.01),
#' colors = 'default',
#' treatment_labels = NULL, sample_labels = NULL, classification_labels= NULL)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param treatment Column name as a \code{string} or \code{numeric} in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any
#' samples that to not contain this factor. This can be a vector of multiple
#' factors to subset on.
#' @param frequencies The range of proportions of samples the taxa are found in.
#' @param abundance_thresholds The range of minimum relative abundances the taxa are found
#' in for each sample.
#' @param colors Name of a color set from the
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palette of R-accepted
#' colors.
#' @param treatment_labels a vector of names to be used as labels for
#' treatments/facets.
#' @param sample_labels a vector of names to be used as labels for Samples.
#' @param classification_labels a vector of names to be used as labels for the
#' taxonomic classifications.
#' @export
#' @return phyloseq-object
#' @examples taxa_core_graph(soil_column)

taxa_core_graph <- function(
  phyloseq_obj,
  treatment = NULL,
  subset = NULL,
  frequencies = seq(0.1, 1, 0.1),
  abundance_thresholds = seq(0.01, 1, 0.01),
  colors = "default",
  treatment_labels = NULL,
  sample_labels = NULL,
  classification_labels = NULL
) {
  check_args(
    phyloseq_obj  = phyloseq_obj,
    treatment     = treatment,
    treatment_labels = treatment_labels
  )
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment = treatment, subset = subset)
  treatment_name <- paste(treatment, collapse = sep)
  phyloseq_obj <- relative_abundance(phyloseq_obj)
  dataset <- melt_phyloseq(phyloseq_obj)

  graph_data <- data.table::data.table()
  if (is.null(treatment)) {
    N <- phyloseq::nsamples(phyloseq_obj)
    for (frequency in frequencies){
      for (abundance_threshold in abundance_thresholds){
        sub_table <- dataset[Abundance >= abundance_threshold]
        taxa <- nrow(
          sub_table[, .(count = .N),
            by = OTU][count >= floor(N * frequency)])
        if (taxa >= 0) {
          graph_data <-
            rbind(graph_data, 
              data.table::data.table(
                freq = factor(frequency),
                abundance = abundance_threshold,
                taxa = taxa
            ))
        } else {
          graph_data <- rbind(
            graph_data, 
            data.table::data.table(
              freq = factor(frequency), 
              abundance = abundance_threshold, 
              taxa = 0
            ))
        }
      }
    }
  } else {
    treatments <- levels(dataset[[treatment_name]])
    for(treatment in treatments){
      N <- sum(phyloseq_obj@sam_data[[treatment_name]] == treatment)
      sub_table <- dataset[dataset[[treatment_name]] == treatment, ]
      for(frequency in frequencies){
        for(abundance_threshold in abundance_thresholds){
          sub_table <- sub_table[Abundance >= abundance_threshold]
          taxa <- 
            nrow(sub_table[, .(count = .N), 
              by = OTU][count >= floor(N * frequency)])
          if(taxa >= 0){
            graph_data <- rbind(graph_data,
                                data.table::data.table(
                                  freq = factor(frequency),
                                  abundance = abundance_threshold,
                                  taxa = taxa,
                                  treatment = treatment)
            )
          } else {
            graph_data <- rbind(graph_data,
                                data.table::data.table(freq = factor(frequency),
                                            abundance = abundance_threshold,
                                            taxa = 0,
                                            treatment = treatment))
          }
        }
      }
    }
  data.table::set(graph_data, j = "treatment",
      value = factor(graph_data[["treatment"]], levels = treatments))
  }

  graph_data <- change_labels(graph_data, treatment_name, treatment_labels,
                      sample_labels)
  graph_colors <- create_palette(length(frequencies), colors)

  g <- ggplot(graph_data, aes(x = abundance, y = taxa, color = freq)) +
    geom_line(size = 1.4) +
    scale_colour_manual(values = graph_colors) +
    guides(colour = guide_legend(
      ncol = ceiling(length(unique(graph_data[["freq"]])) / 25),
      override.aes = list(size = 4)
    )) +
    theme_schuy("bar", 35) +
    labs(x = "Relative Abundance", y = "Number of OTUs", color = "Proportion\nof Samples") +
    scale_y_continuous(expand = expansion(add = c(0.3, 0.5))) +
    scale_x_continuous(expand = expansion(add = c(0.0005, 0.001)))
  if(!is.null(treatment)){
    g <- g + facet_wrap(~treatment, ncol = 3, dir = "v")
  }
  return(g)
}
