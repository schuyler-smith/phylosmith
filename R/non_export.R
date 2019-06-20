#' Converts numeric values to column names in sample_data.
#'
#' Converts numeric values to column names in sample_data.
#' @useDynLib phylosmith
#' @usage check_numeric_treatment(phyloseq_obj, ...)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param ... Column name as a \code{string} or \code{numeric} in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be any number of
#' multiple columns and they will be combined into a new column.
#' @return string

check_numeric_treatment <- function(phyloseq_obj, ...) {
  treatments <- list(...)
  if (any(unlist(lapply(treatments, is.null))) |
      any(unlist(lapply(treatments, is.na)))) {
    return(NULL)
  } else {
    return(tryCatch({
      unlist(lapply(
        treatments,
        FUN = function(treatment) {
          colnames(access(phyloseq_obj, 'sam_data')[, treatment])
        }
      ))
    }, error = function(e) {
      stop(
        "`treatment` must be at least one column name, or index, from the sample_data()",
        call. = FALSE
      )
    }))
  }
}

#' Converts numeric values to column names in tax_table
#'
#' Converts numeric values to column names in tax_table.
#' @useDynLib phylosmith
#' @usage check_numeric_classification(phyloseq_obj, ...)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param ... Column name as a \code{string} or \code{numeric} in the
#' \code{\link[phyloseq:tax_table]{tax_table}} for the factor to conglomerate
#' by.
#' @return string

check_numeric_classification <- function(phyloseq_obj, ...) {
  classifications <- list(...)
  if (any(unlist(lapply(classifications, is.null))) |
      any(unlist(lapply(classifications, is.na)))) {
    return(NULL)
  } else {
    return(tryCatch({
      unlist(lapply(
        classifications,
        FUN = function(classification) {
          colnames(access(phyloseq_obj, 'tax_table')[, classification])
        }
      ))
    }, error = function(e) {
      stop(
        "taxa_filter(): `classification` must be at least one column name, or index, from the tax_table()",
        call. = FALSE
      )
    }))
  }
}

#' Creates color palettes for figures.
#'
#' Creates color palettes for figures using the the RColorBrewer, or the default
#' color palette I made. The first 8 colors of the palette are from B. Wong's 2011
#' publication. The additional colors I added and may change as time goes by.
#' @useDynLib phylosmith
#' @usage create_palette(color_count, colors = 'default')
#' @param color_count Number of colors to choose for palette.
#' @param colors Name of a color set from the
#' \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R-accepted
#' colors.
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom grDevices colorRampPalette
#' @return palette
#' @examples
#' #create_palette(8, 'Dark2')

create_palette <- function(color_count, colors = 'default') {
  options(warn = -1)
  mycolors <- c(
    "#A8B1CC",
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#9EDA8F",
    "#CC79A7",
    "#757575",
    "#DE9861",
    "#A6CBE0",
    "#B275D8",
    "#82BB47",
    "#e0503a",
    "#F5E56C",
    "#949696",
    "#4989DE",
    "#E2E2E2",
    "#565656",
    "#F7B04C",
    "#696bb2"
  )
  # image(1:length(mycolors), 1, as.matrix(1:length(mycolors)), col=mycolors, xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
  if (any(colors == 'default')) {
    colors <- mycolors
    if (color_count <= length(colors)) {
      return(colors[seq(color_count)])
    }
  }
  if (any(colors %in% 'rev') | any(colors %in% 'reverse')) {
    colors <- rev(mycolors)
    if (color_count <= length(mycolors)) {
      return(rev(mycolors[seq(color_count)]))
    }
  }
  if (any(!(colors %in% colors()))) {
    if (any(colors %in% rownames(brewer.pal.info))) {
      getPalette <- colorRampPalette(brewer.pal(min(c(
        color_count,
        brewer.pal.info[rownames(brewer.pal.info) == colors, 1]
      )), colors))
    } else {
      getPalette <- colorRampPalette(colors)
    }
  } else {
    getPalette <- colorRampPalette(colors)
  }
  return(getPalette(color_count))
}

#' Creates confidence interval ellipse for scatterplot
#'
#' Creates confidence interval ellipse for scatterplot. Can be created
#' for each group or entire set.
#' @useDynLib phylosmith
#' @usage CI_ellipse(points, groups = NULL, level = 0.95)
#' @param points dataframe of the x,y coordinates of the datapoints, along with
#' any metadata such as group
#' @param groups Name of the column in points that has group membership
#' @param level the confidence level to be plotted (default 95\%)
#' @importFrom stats cov.wt qf
#' @return ellipse points

CI_ellipse <- function(points,
                       groups = NULL,
                       level = 0.95) {
  ellipse_df <- data.frame()
  theta <- (0:99) * 2 * pi / 100
  unit_circle <- cbind(cos(theta), sin(theta))
  if (!is.null(groups)) {
    for (group in unique(points[, groups])) {
      sub_points <- points[points[, groups] == group, c('x', 'y')]
      ellipse_info <- cov.wt(sub_points[, c('x', 'y')])
      shape <- ellipse_info$cov
      center <- ellipse_info$center
      radius <- sqrt(2 * qf(level, 2, length(sub_points$x) - 1))
      Q <- chol(shape, pivot = TRUE)
      order <- order(attr(Q, "pivot"))
      ellipse <- t(center + radius * t(unit_circle %*% Q[, order]))
      ellipse_df <- rbind(ellipse_df, cbind(ellipse, group))
    }
  } else {
    ellipse_info <- cov.wt(points[, c('x', 'y')])
    shape <- ellipse_info$cov
    center <- ellipse_info$center
    radius <- sqrt(2 * qf(level, 2, length(points$x) - 1))
    Q <- chol(shape, pivot = TRUE)
    order <- order(attr(Q, "pivot"))
    ellipse <- t(center + radius * t(unit_circle %*% Q[, order]))
    ellipse_df <- cbind(ellipse)
  }
  return(ellipse_df)
}
