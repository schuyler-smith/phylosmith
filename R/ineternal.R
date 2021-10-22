#' Change names of factors in graph_data to change labels in plot objects
#'
#' Change names of factors in graph_data to change labels in plot objects
#' @useDynLib phylosmith
#' @usage change_labels(graph_data, treatment = NULL, treatment_labels = NULL,
#' sample_labels = NULL, classification = NULL, classification_labels = NULL)
#' @param graph_data A \code{\link[data.table:data.table]{data.table()}}) object
#' containing the melted phyloseq data.
#' @param treatment name of the column defined as the treatment.
#' @param treatment_labels a vector of names to be used as labels for
#' treatments/facets.
#' @param sample_labels a vector of names to be used as labels for Samples.
#' @param classification name of the column defined as the classification.
#' @param classification_labels a vector of names to be used as labels for the
#' taxonomic classifications.
#' @import data.table
#' @return data.table

change_labels <- function(graph_data,
                          treatment = NULL,
                          treatment_labels = NULL,
                          sample_labels = NULL,
                          classification = NULL,
                          classification_labels = NULL
                          ) {
  if(!is.null(treatment_labels)){
    set(graph_data, j = treatment,
        value = factor(treatment_labels[as.numeric(graph_data[[treatment]])],
                       levels = treatment_labels)
    )
  }
  if(!is.null(sample_labels)){
    set(graph_data, j = 'Sample',
        value = factor(sample_labels[as.numeric(graph_data[['Sample']])],
                       levels = sample_labels)
    )
  }
  if(!is.null(classification_labels)){
    set(graph_data, j = classification,
        value = factor(classification_labels[as.numeric(graph_data[[classification]])],
                       levels = classification_labels)
    )
  }
  return(graph_data)
}

#' Converts numeric values to column names in sample_data.
#'
#' Converts numeric values to column names in sample_data.
#' @useDynLib phylosmith
#' @usage check_index_treatment(phyloseq_obj, columns)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param columns Column name as a \code{string} or \code{numeric} in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be any number of
#' multiple columns and they will be combined into a new column.
#' @return string

check_index_treatment <- function(phyloseq_obj, columns) {
  argument <- deparse(substitute(columns))
  treatments <- list(columns)
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
        paste("{",argument,"}","must be at least one column name, or index, from the sample_data()"),
        call. = FALSE
      )
    }))
  }
}

#' Converts numeric values to column names in tax_table
#'
#' Converts numeric values to column names in tax_table.
#' @useDynLib phylosmith
#' @usage check_index_classification(phyloseq_obj, columns)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param columns Column name as a \code{string} or \code{numeric} in the
#' \code{\link[phyloseq:tax_table]{tax_table}} for the factor to conglomerate
#' by.
#' @return string

check_index_classification <- function(phyloseq_obj, columns) {
  argument <- deparse(substitute(columns))
  classifications <- list(columns)
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
        paste("{",argument,"}","must be at least one column name, or index, from the tax_table()"),
        call. = FALSE
      )
    }))
  }
}

#' phyloseq object so taxa are rows.
#'
#' phyloseq object so taxa are rows.
#' @useDynLib phylosmith
#' @usage check_TaR(phyloseq_obj)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @return phyloseq_obj

check_TaR <- function(phyloseq_obj){
  if(!(attributes(phyloseq_obj@otu_table)$taxa_are_rows)){
    phyloseq_obj@otu_table <- t(phyloseq_obj@otu_table)
  }
  return(phyloseq_obj)
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
    "#757575", "#E69F00", "#56B4E9", "#009E73",
    "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
    "#9EDA8F", "#DE9861","#565656",  "#A6CBE0",
    "#B275D8", "#82BB47", "#e0503a", "#F5E56C",
    "#949696", "#4989DE", "#E2E2E2",
    "#F7B04C", "#696bb2")
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
    group_members <- table(points[, groups])
    group_members <- as.numeric(names(group_members[group_members >= 3]))
    for (group in group_members) {
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


#' Bin data
#'
#' Assigns bins to a vector of values
#' @useDynLib phylosmith
#' @usage bin(data, nbins, labels = NULL)
#' @param data vector of data to bin
#' @param nbins number of bins to produce
#' @param labels names for the bins
#' @return bins of values

bin <- function(data, nbins = 9, labels = NULL) {
  vec <- FALSE
  if (is.atomic(data) == TRUE & is.null(dim(data)) == TRUE) {
    vec <- TRUE
    data <- data.frame(data)
  }
  if (is.list(data) == FALSE)
    data <- data.frame(data)
  data <- na.omit(data)
  if (!is.null(labels))
    if (nbins != length(labels))
      stop("number of 'nbins' and 'labels' differ")
  if (nbins <= 1)
    stop("nbins must be greater than 1")
  data <- lapply(data, function(x)
    if (is.numeric(x)) {
      if (length(unique(x)) <= nbins){
        as.factor(x)
      } else {
        CUT(x, breaks = unique(nbins), labels = labels)
      }
    } else as.factor(x))
  if (vec) {
    data <- unname(unlist(data))
  }
  return(data)
}

CUT <- function(x, breaks, ...) {
  if (length(breaks) == 1L) {
    nb <- as.integer(breaks + 1)
    rx <- range(x, na.rm = TRUE)
    dx <- diff(rx)
    if (dx == 0) {
      dx <- abs(rx[1L])
      breaks <- seq.int(rx[1L] - dx/1000, rx[2L] + dx/1000, length.out = nb)
    } else {
      breaks <- seq.int(rx[1L], rx[2L], length.out = nb)
      breaks[c(1L, nb)] <- c(floor(rx[1L]), ceiling(rx[2L]))
    }
  }
  breaks_f <- c(breaks[1], as.numeric(formatC(0 + breaks[2:(length(breaks)-1)], digits = 1, width = 1L)), breaks[length(breaks)])
  cut_x <- cut(x, breaks = unique(breaks_f), ...)
  return(cut_x)
}

ADDNA <- function(x) {
  if (is.factor(x) & !("NA" %in% levels(x))) x <- factor(x, levels = c(levels(x), "NA"))
  x[is.na(x)] <- "NA"
  x
}
