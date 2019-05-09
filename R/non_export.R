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
#' @examples
#' check_numeric_treatment <- function(phyloseq_obj, ...){
#' treatments <- list(...)
#' treatments <- unlist(sapply(treatments, FUN = function(treatment){
#'   if(is.numeric(treatment)){
#'     return(colnames(access(phyloseq_obj, 'sam_data')[,treatment]))
#'   } else {return(treatment)}
#' }))
#' return(treatments)
#' }
#' check_numeric_treatment(soil_column, 2)

check_numeric_treatment <- function(phyloseq_obj, ...){
    treatments <- list(...)
    if(any(unlist(sapply(treatments, is.null))) | any(unlist(sapply(treatments, is.na)))){
      return(NULL)
    } else {
      return(unlist(sapply(treatments, FUN = function(treatment){
        colnames(access(phyloseq_obj, 'sam_data')[,treatment])
      })))
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
#' @examples
#'
#' check_numeric_classification <- function(phyloseq_obj, ...){
#' classifications <- list(...)
#' classifications <- unlist(sapply(classifications,
#'     FUN = function(classification){
#'         if(is.numeric(classification)){
#'             return(colnames(access(phyloseq_obj,
#'             'tax_table')[,classification]))
#'         } else {return(classification)}
#'         }))
#' return(classifications)
#' }
#' check_numeric_classification(soil_column, 2)

check_numeric_classification <- function(phyloseq_obj, ...){
    classifications <- list(...)
    if(any(unlist(sapply(classifications, is.null))) | any(unlist(sapply(classifications, is.na)))){
      return(NULL)
    } else {
      return(unlist(sapply(classifications, FUN = function(classification){
            colnames(access(phyloseq_obj, 'tax_table')[,classification])
      })))
    }
}

#' Creates color palettes for figures.
#' 
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
#' @import RColorBrewer
#' @import grDevices
#' @return palette
#' @examples
#' #create_palette(8, 'Dark2')

create_palette <- function(color_count, colors = 'default'){
    options(warn = -1)
    mycolors <- c(
        "#A8B1CC", "#E69F00", "#56B4E9", "#009E73",
        "#F0E442", "#0072B2", "#D55E00", "#9EDA8F",
        "#CC79A7", "#757575", "#DE9861", "#A6CBE0",
        "#B275D8", "#82BB47", "#e0503a", "#F5E56C",
        "#949696", "#4989DE", "#E2E2E2", "#565656",
        "#F7B04C", "#696bb2")
    # image(1:length(mycolors), 1, as.matrix(1:length(mycolors)), col=mycolors, xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
    if(any(colors == 'default')){
        colors <- mycolors
        if(color_count <= length(colors)){
            return(colors[seq(color_count)])
        }
    }
    if(any(colors %in% 'rev') | any(colors %in% 'reverse')){
        colors <- rev(mycolors)
        if(color_count <= length(mycolors)){
            return(rev(mycolors[seq(color_count)]))
        }
    }
    if(any(!(colors %in% colors()))){
        if(any(colors %in% rownames(brewer.pal.info))){
            getPalette <- colorRampPalette(brewer.pal(min(c(color_count,
            brewer.pal.info[rownames(brewer.pal.info) == colors, 1])), colors))
        } else { getPalette <- colorRampPalette(colors)}
    } else { getPalette <- colorRampPalette(colors)}
    return(getPalette(color_count))
}
