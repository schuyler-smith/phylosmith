#' Create a ggplot object of the distribution of rho values from
#' permute_rho(). Function from the phylosmith-package.
#'
#' Plots the output of \code{\link[=permute_rho]{permute_rho}} into a
#' histogram with the distributions shown by treatment. This is a
#' visualization tool to help show how the permutation worked, and to see
#' where the cutoffs lie.
#' @useDynLib phylosmith
#' @usage histogram_permuted_rhos(permuted_rhos, p = NULL,
#' x_breaks = 0.25, colors = 'default')
#' @param permuted_rhos A \code{data.table} output from
#' \code{\link[=permute_rho]{permute_rho}}.
#' @param p The significance threshold for setting cutoffs.
#' @param x_breaks What intervals to set the ticks on the x-axis.
#' @param colors Name of a color set from the
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R-accepted
#' colors.
#' @seealso \code{\link[=permute_rho]{permute_rho}}
#' @export
#' @return ggplot
#' @examples
#' permuted_rhos <- permute_rho(soil_column, treatment = c('Matrix', 'Treatment'),
#' replicate_samples = 'Day', permutations = 1,  cores = 0)
#' histogram_permuted_rhos(permuted_rhos)
#' histogram_permuted_rhos(permuted_rhos, p = 0.01)

histogram_permuted_rhos <- function(
  permuted_rhos,
  p = NULL,
  x_breaks = 0.25,
  colors = "default"
) {
  check_args(
    p   = p
  )
  if (!(is.data.frame(permuted_rhos))) {
    stop("`permuted_rhos` must be at data.frame
        object", call. = FALSE)
  }
  if (!(is.numeric(x_breaks)) |
      !(x_breaks >= -1 & x_breaks <= 1)) {
    stop("`x_breaks` must be a numeric value
        between -1 and 1", call. = FALSE)
  }
  color_count <- length(unique(permuted_rhos[["Treatment"]]))
  graph_colors <- create_palette(color_count, colors)

  if(color_count > 0){
    permuted_rhos[, Proportion := Count / sum(Count),by = Treatment]
    quantiles <- permuted_rhos[, 
      list(lower = rho[sum(cumsum(Proportion) <= (p /2))],
      upper = rho[sum(cumsum(Proportion) <= (1 - (p / 2)))]),
      by = Treatment]

    permuted_rhos[,
      bin := findInterval(rho, 
        seq(-1, 1, (max(permuted_rhos$rho, na.rm = TRUE) -
          min(permuted_rhos$rho, na.rm = TRUE)) / 100)
      ), by = Treatment]
    permuted_rhos <- permuted_rhos[, list(
      rho = mean(rho),
      Count = sum(Count),
      Proportion = sum(Proportion)
    ), by = .(Treatment, bin)]

    g <-
      ggplot(permuted_rhos, aes(x = rho, y = Proportion, fill = Treatment))
  } else {
    permuted_rhos[, Proportion := Count / sum(Count)]
    quantiles <- permuted_rhos[, list(lower = rho[sum(cumsum(Proportion) <= (p /
                                                                               2))],
                                      upper = rho[sum(cumsum(Proportion) <= (1 - (p / 2)))])]

    permuted_rhos[, bin := findInterval(rho, seq(-1, 1, (max(permuted_rhos$rho) - min(permuted_rhos$rho))/100))]
    permuted_rhos <- permuted_rhos[, list(
      rho = mean(rho),
      Count = sum(Count),
      Proportion = sum(Proportion)
    ), by = .(bin)]

    g <-
      ggplot(permuted_rhos, aes(x = rho, y = Proportion, fill = graph_colors))
  }
  if (!(is.null(p))) {
    g <-
      g + geom_vline(
        data = quantiles,
        xintercept = c(quantiles$lower, quantiles$upper),
        color = c(graph_colors, graph_colors),
        size = 1.5,
        alpha = 0.8
      )
  }
  g <- g + scale_x_continuous(breaks = seq(-1, 1, x_breaks)) +
    scale_fill_manual(values = graph_colors) +
    geom_bar(stat = "identity",
             position = position_dodge(width = .02),
             width = .05)
  if (!(is.null(p))) {
    for (i in seq(nrow(quantiles))) {
      g <- g + geom_text(
        data = quantiles,
        x = quantiles$lower[i] - .15,
        label = paste0(round(quantiles$lower[i], 2)),
        y = (.015 + (i * .005)),
        color = graph_colors[i],
        size = 5
      ) +
        geom_text(
          data = quantiles,
          x = quantiles$upper[i] + .15,
          label = paste0(round(quantiles$upper[i], 2)),
          y = (.015 + (i * .005)),
          color = graph_colors[i],
          size = 5
        )
    }
  }
  g <- g + theme_schuy("bar", 35) + 
    scale_y_continuous(expand = expansion(mult = c(0.0025, 0.002)))
  return(g)
}