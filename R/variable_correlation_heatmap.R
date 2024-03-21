#' Computes the correlation of numerical variables with taxa
#' Function from the phylosmith-package.
#'
#' Computes the correlation of numerical variables with taxa
#' @useDynLib phylosmith
#' @usage variable_correlation_heatmap(phyloseq_obj, variables, treatment = NULL,
#' subset = NULL, classification = NULL, method = "spearman", limits = c(-0.8, 0.8),
#' colors = "default", significance_color = "white", cores = 1, treatment_labels = NULL,
#' sample_labels = NULL, classification_labels = NULL)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object.
#' @param treatment Column name as a \code{string} or \code{numeric} in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @param variables Numerical factors within the in the
#' \code{\link[phyloseq:sample_data]{sample_data}} to correlate with the
#' abundance data.
#' @param subset A factor within the \code{treatment}. This will remove any
#' samples that to not contain this factor. This can be a vector of multiple
#' factors to subset on.
#' @param classification Column name as a \code{string} or \code{numeric} in
#' the \code{\link[phyloseq:tax_table]{tax_table}} for the factor to
#' conglomerate by.
#' @param method Which correlation method to calculate, "pearson", "spearman".
#' @param limits The range for the legend, smaller limits will accentuate smaller
#' correlations.
#' @param colors the palette to use for the heatmap, default is viridis.
#' @param significance_color the color to use for the significance stars.
#' @param cores \code{numeric} Number of CPU cores to use for the pair-wise
#' permutations. Default (1), (0) uses max cores available. Parallelization not
#' available for systems running MacOS without openMP configuration.
#' @param treatment_labels a vector of names to be used as labels for
#' treatments/facets.
#' @param sample_labels a vector of names to be used as labels for Samples.
#' @param classification_labels a vector of names to be used as labels for the
#' taxonomic classifications.
#' @importFrom parallel detectCores
#' @keywords nonparametric
#' @seealso \code{\link{permute_rho}} \code{\link{phylosmith}}
#' @export
#' @return data.table
#' @examples
#' variable_correlation_heatmap(soil_column, variables = 'Day',
#' treatment = c('Matrix', 'Treatment'), subset = 'Amended',
#' classification = 'Phylum', method = 'spearman', cores = 1,
#' colors = c("#2C7BB6", "white", "#D7191C"),
#' significance_color = 'black')

variable_correlation_heatmap <- function(
  phyloseq_obj,
  variables,
  treatment = NULL,
  subset = NULL,
  classification = NULL,
  method = "spearman",
  limits = c(-0.8, 0.8),
  colors = "default",
  significance_color = "white",
  cores = 1,
  treatment_labels = NULL,
  sample_labels = NULL,
  classification_labels = NULL
) {
  check_args(
    phyloseq_obj   = phyloseq_obj,
    treatment      = treatment,
    subset         = subset,
    classification = classification,
    corr_method    = method,
    cores          = cores
  )
  correlations <- variable_correlation(
    phyloseq_obj,
    variables,
    treatment,
    subset,
    classification,
    method,
    cores
  )
  correlations <- change_labels(correlations, treatment_name, treatment_labels,
                      sample_labels, classification, classification_labels)
  correlations$Significance <- cut(correlations$p, 
    breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),
    label=c("***", "**", "*", ""))
  data.table::set(correlations, j = "X", value = factor(correlations[["X"]],
    levels = rev(unique(sort(correlations[["X"]])))
  ))
  if(is.null(treatment)){
    g <- ggplot(data = correlations, aes(x = Y, y = X, fill = rho))
  } else {
    g <- ggplot(data = correlations, aes(x = Treatment, y = X, fill = rho))
  }
  g <- g + geom_tile() + geom_text(aes(label = Significance), 
    color = significance_color, size = 3) +
    labs(y = NULL, x = NULL, fill = method) +
    facet_grid(. ~ Y, drop = TRUE, scales = "free", space = "free_x")
  if("default" %in% colors){
    g <- g + viridis::scale_fill_viridis(limit = limits)
  } else {
    g <- g + scale_fill_gradient2(low = colors[1], mid = colors[2],
      high = colors[3], limit = limits)
  } + scale_x_discrete(expand = expansion(mult = 0, add = 0))
  g <- g + theme_schuy("heatmap", 35)

  return(g)
}
