#' Filter taxa based on proportion of samples they are observed in. Function from the phylosmith-package.
#'
#' This function takes a phyloseq object and finds which taxa are seen in a given proportion of samples, either in the entire dataset, by treatment, or a particular treatment of interest.
#' @useDynLib phylosmith
#' @usage taxa_filter(phyloseq_obj, treatment = NULL, subset = NULL,
#' frequency = 0, below = FALSE, drop_samples = FALSE)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with information about each sample, and it must contain \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each taxa/gene.
#' @param treatment Column name as a string or number in the \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any samples that to not contain this factor. This can be a vector of multiple factors to subset on.
#' @param frequency The proportion of samples the taxa is found in.
#' @param below Does frequency define the minimum (\code{FALSE}) or maximum (\code{TRUE}) proportion of samples the taxa is found in.
#' @param drop_samples Should the function remove samples that that are empty after removing taxa filtered by frequency (\code{TRUE}).
#' @keywords manip
#' @import data.table
#' @export

taxa_filter <- function(phyloseq_obj, treatment = NULL, subset = NULL, frequency = 0, below = FALSE, drop_samples = FALSE){
  # data("mock_phyloseq")
  # phyloseq_obj = mock_phyloseq; frequency = 0; treatment = c("treatment", "day"); subset = "5"; below = FALSE; drop_samples = FALSE
  if(!(is.null(phyloseq_obj@phy_tree))){phylo_tree <- phyloseq_obj@phy_tree} else {phylo_tree <- FALSE}
  if(!(is.null(phyloseq_obj@refseq))){refseq <- phyloseq_obj@refseq} else {refseq <- FALSE}
  phyloseq_obj <- phyloseq(phyloseq_obj@otu_table, phyloseq_obj@tax_table, phyloseq_obj@sam_data)

  if(!(is.null(treatment))){
    if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
    phyloseq_obj <- merge_treatments(phyloseq_obj, treatment)
    treatment_name <- paste(treatment, collapse = sep)
    treatment_classes <- sort(unique(phyloseq_obj@sam_data[[treatment_name]]))

    if(below == TRUE){
      phyloseq_obj <- do.call(merge_phyloseq,
                              apply(array(treatment_classes), 1, FUN = function(group){
                                sub_phy <- eval(parse(text=paste0("subset_samples(phyloseq_obj, ",treatment_name," == '",group,"')")))
                                cutoff <- floor(ncol(sub_phy@otu_table) * frequency)
                                sub_phy <- filter_taxa(sub_phy, function(x){sum(x != 0, na.rm = TRUE) <= cutoff}, TRUE)
                                return(sub_phy)}))
    } else {
      phyloseq_obj <- do.call(merge_phyloseq,
                              apply(array(treatment_classes), 1, FUN = function(group){
                                sub_phy <- eval(parse(text=paste0("subset_samples(phyloseq_obj, ",treatment_name," == '",group,"')")))
                                cutoff <- floor(ncol(sub_phy@otu_table) * frequency)
                                sub_phy <- filter_taxa(sub_phy, function(x){sum(x != 0, na.rm = TRUE) >= cutoff}, TRUE)
                                if(sum(taxa_sums(sub_phy), na.rm = TRUE) != 0){sub_phy <- filter_taxa(sub_phy, function(x){sum(x, na.rm = TRUE) != 0}, TRUE)}
                                return(sub_phy)}))
    }
  } else {
    cutoff <- floor(ncol(phyloseq_obj@otu_table) * frequency)
    if(below == TRUE){
      phyloseq_obj <- filter_taxa(phyloseq_obj, function(x){sum(x != 0, na.rm = TRUE) <= cutoff}, TRUE)
    } else {
      phyloseq_obj <- filter_taxa(phyloseq_obj, function(x){sum(x != 0, na.rm = TRUE) >= cutoff}, TRUE)
      phyloseq_obj <- filter_taxa(phyloseq_obj, function(x){sum(x, na.rm = TRUE) != 0}, TRUE)
      }
  }
  if(!(is.null(subset))){
    phyloseq_obj <- prune_samples(sample_names(phyloseq_obj)[unname(apply(phyloseq_obj@sam_data[,c(treatment, treatment_name)], 1, function(x){any(x %in% subset)}))], phyloseq_obj)
    phyloseq_obj <- filter_taxa(phyloseq_obj, function(x){sum(x, na.rm = TRUE) > 0}, TRUE)
  }
  if(drop_samples == TRUE){
    phyloseq_obj <- prune_samples(sample_sums(phyloseq_obj) > 0, phyloseq_obj)
  }
  if(!(is.logical(phylo_tree))){phyloseq_obj@phy_tree <- phylo_tree}
  if(!(is.logical(refseq))){phyloseq_obj@refseq <- refseq}
  return(phyloseq_obj)
}


#' Combine meta-data columns. Function from the phylosmith-package.
#'
#' Combines multiple columns of a \code{\link[phyloseq]{phyloseq-class}} object \code{\link[phyloseq:sample_data]{sample_data}} into a single-variable column.
#' @useDynLib phylosmith
#' @usage merge_treatments(phyloseq_obj, ...)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with information about each sample.
#' @param ... any number of column names as strings or numbers in the \code{\link[phyloseq:sample_data]{sample_data}} that are to be combined.
#' @keywords manip
#' @import data.table
#' @export

merge_treatments <- function(phyloseq_obj, ...){
  treatments <- list(...)
  treatments <- sapply(treatments, FUN = function(treatment){
    if(is.numeric(treatment)){
      return(colnames(phyloseq_obj@sam_data[,treatment]))
      } else {return(treatment)}
  })
  treatment_classes <- setDT(as(phyloseq_obj@sam_data[,colnames(phyloseq_obj@sam_data) %in% treatments], "data.frame"))
  treatment_name <- paste(treatments, collapse = sep)
  order <- apply(eval(parse(text=paste0("expand.grid(", paste0(paste0("levels(factor(phyloseq_obj@sam_data[['", treatments, "']]))", collapse = ', ')), ")"))), 1, FUN = function(combination){paste0(combination, collapse = sep)})
  eval(parse(text=paste0("treatment_classes[, '",treatment_name,"' := as.character(paste(", paste(treatments, collapse = ', '), ', sep = sep), by = treatment_classes)]')))
  phyloseq_obj@sam_data[[treatment_name]] <- factor(treatment_classes[[treatment_name]], levels = order)
  return(phyloseq_obj)
}


#' Transform abundance data in an \code{otu_table} to relative abundance, sample-by-sample. Function from the phylosmith-package.
#'
#' Transform abundance data into relative abundance, i.e. proportional data. This is an alternative method of normalization and may not be appropriate for all datasets, particularly if your sequencing depth varies between samples.
#' @useDynLib phylosmith
#' @usage relative_abundance(phyloseq_obj)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with information about each sample, and it must contain \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each taxa/gene.
#' @keywords manip
#' @export

relative_abundance <- function(phyloseq_obj){
  options(warnings=-1)
  phyloseq_obj <- transform_sample_counts(phyloseq_obj, function(sample) sample/sum(sample))
  return(phyloseq_obj)
}


#' order_treatment
#'
#' Reorders the levels of a metadata column in a \code{\link[phyloseq]{phyloseq-class}} object \code{\link[phyloseq:sample_data]{sample_data}}.
#' @useDynLib phylosmith
#' @usage order_treatment(phyloseq_obj, treatment, order)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with information about each sample, and it must contain \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each taxa/gene.
#' @param treatment Column name as a string or number in the \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param order The order of factors in \code{treatment} column as a vector of strings.
#' @keywords manip
#' @import data.table
#' @export

order_treatment <- function(phyloseq_obj, treatment, order){
  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  phyloseq_obj@sam_data[[treatment]] <- factor(phyloseq_obj@sam_data[[treatment]], levels = order)
  return(phyloseq_obj)
}

#' Merge samples based on common factor within sample_data. Function from the phylosmith-package.
#'
#' This function takes a phyloseq object and merges the samples that meet the specified criteria into a single sample. This is meant for replicates, or samples statistically proven to not be significantly different and should be used with caution as it may be a misleading representation of the data.
#' @useDynLib phylosmith
#' @usage merge_samples(phyloseq_obj, treatment, subset = NULL, merge_on = treatment)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with information about each sample, and it must contain \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each taxa/gene.
#' @param treatment Column name as a string or number in the \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any samples that to not contain this factor. This can be a vector of multiple factors to subset on.
#' @param merge_on Defines which variable the data is merged according to. This needs to be a column name as a string or number in the \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of multiple columns and they will be combined into a new column.
#' @keywords manip
#' @export

merge_samples <- function(phyloseq_obj, treatment, subset = NULL, merge_on = treatment){
  options(warn = -1)
  if(!(is.null(phyloseq_obj@phy_tree))){phylo_tree <- phyloseq_obj@phy_tree} else {phylo_tree <- FALSE}
  if(!(is.null(phyloseq_obj@refseq))){refseq <- phyloseq_obj@refseq} else {refseq <- FALSE}
  phyloseq_obj <- phyloseq(phyloseq_obj@otu_table, phyloseq_obj@tax_table, phyloseq_obj@sam_data)

  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  if(is.numeric(merge_on)){merge_on <- colnames(phyloseq_obj@sam_data[,merge_on])}
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, subset)
  phyloseq_obj <- merge_treatments(phyloseq_obj, merge_on)
  treatment_name <- paste(treatment, collapse = sep)
  merge_on <- paste(merge_on, collapse = sep)
  treatment_classes <- sort(unique(phyloseq_obj@sam_data[[treatment_name]]))
  treatment_classes <- eval(parse(text=paste0('treatment_classes[grepl("', paste0(subset), '", treatment_classes)]')))
  merge_sample_levels <- as.character(unique(sort(unlist(phyloseq_obj@sam_data[[merge_on]]))))
  if(merge_on != treatment){merge_sample_levels <- paste(sapply(treatment_classes,rep,times=length(merge_sample_levels)), rep(merge_sample_levels, length(treatment_classes)), sep = sep)}

  phyloseq_table <- melt_phyloseq(phyloseq_obj)
  if(merge_on != treatment){phyloseq_table[, 'Merged_Name' := do.call(paste0, list(phyloseq_table[[treatment_name]], sep, phyloseq_table[[merge_on]]))]
  } else {phyloseq_table[, 'Merged_Name' := phyloseq_table[[merge_on]]]}
  otu_tab <- dcast(phyloseq_table[,c('OTU','Abundance','Merged_Name'), with=FALSE], Merged_Name ~ OTU, value.var = 'Abundance', fun.aggregate = mean)

  sub_phy <- do.call(merge_phyloseq,
     sapply(treatment_classes, FUN = function(group){
       group_phy <- eval(parse(text=paste0('subset_samples(taxa_filter(phyloseq_obj, treatment), ', treatment_name, ' == "', group, '")')))
       if(nsamples(group_phy) > 1){sub_phy <- phyloseq::merge_samples(group_phy, merge_on)
       merge_names <- rownames(sub_phy@sam_data)
       if(group != merge_names){sample_names(sub_phy) <- paste0(group, sep, merge_names)}
       sam <- as(sub_phy@sam_data, 'data.frame')
       sam[,treatment_name] <- group
       for(i in treatment){
         sam[,i] <- unique(group_phy@sam_data[,i])}
       sam[, merge_on] <- factor(merge_names, levels = levels(phyloseq_obj@sam_data[[merge_on]]))
       sub_phy@sam_data <- sample_data(sam)}
       return(sub_phy)
     })
  )
  phyloseq_obj <- tryCatch({phyloseq_obj <- eval(parse(text=paste0('subset_samples(phyloseq_obj, !(', treatment_name,' %in% treatment_classes))')))},
     error = function(e){phyloseq_obj <- sub_phy},
     finally = {merge_phyloseq(phyloseq_obj, sub_phy)})
  phyloseq_obj <- phyloseq(otu_table(t(as.matrix(otu_tab[order(factor(otu_tab$Merged_Name, levels = merge_sample_levels)),], rownames = 'Merged_Name')), taxa_are_rows = TRUE),
                 phyloseq_obj@tax_table,
                 phyloseq_obj@sam_data[order(factor(rownames(phyloseq_obj@sam_data), levels = merge_sample_levels)),])
  if(!(is.logical(phylo_tree))){phyloseq_obj@phy_tree <- phylo_tree}
  if(!(is.logical(refseq))){phyloseq_obj@refseq <- refseq}
  return(phyloseq_obj)
}

#' Melt a phyloseq object data.table.
#'
#' melt_phyloseq takes a phyloseq object and melts its ou_table, taxa_tale, and sample_Data into a single into a data.table.
#' @useDynLib phylosmith
#' @usage melt_phyloseq(phyloseq_obj)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It must contain \code{\link[phyloseq:sample_data]{sample_data()}} with information about each sample, and it must contain \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each taxa/gene.
#' @keywords manip
#' @seealso \code{\link[phyloseq:psmelt]{psmelt()}}
#' @import data.table
#' @export

melt_phyloseq <- function(phyloseq_obj){
  melted_phyloseq <- melt.data.table(data.table(as(phyloseq_obj@otu_table, "matrix"), keep.rownames = TRUE), id.vars = 1)
  colnames(melted_phyloseq) <- c("OTU", "Sample", "Abundance")
  taxa <- data.table(phyloseq_obj@tax_table, OTU = taxa_names(phyloseq_obj))
  sample_data <- data.table(data.frame(phyloseq_obj@sam_data, stringsAsFactors = FALSE))
  sample_data[, 'Sample' := sample_names(phyloseq_obj)]

  melted_phyloseq <- merge(melted_phyloseq, sample_data, by = "Sample")
  melted_phyloseq <- merge(melted_phyloseq, taxa, by = "OTU")
  melted_phyloseq <- melted_phyloseq[order(melted_phyloseq$Abundance, decreasing = TRUE), ]

  return(melted_phyloseq)
}


#' Conglomerate taxa by sample on a given classification level
#'
#' Conglomerate taxa by sample on a given classification level from the tax_table.
#' @useDynLib phylosmith
#' @usage conglomerate_taxa(phyloseq_obj, classification, taxa_are_ordered = TRUE)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It must contain \code{\link[phyloseq:sample_data]{sample_data()}} with information about each sample, and it must contain \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each taxa/gene.
#' @param classification Column name as a string or number in the \code{\link[phyloseq:tax_table]{tax_table}} for the factor to conglomerate by.
#' @param taxa_are_ordered Whether the order of factors in the tax_table represent a decreasing heirarchy (TRUE) or are independant (FALSE). If FALSE, will only return the factor given by \code{classification}.
#' @keywords manip
#' @seealso \code{\link[phyloseq:tax_glom]{tax_glom()}}
#' @import data.table
#' @export

conglomerate_taxa <- function(phyloseq_obj, classification, taxa_are_ordered = TRUE){
  if(is.numeric(classification)){classification <- colnames(phyloseq_obj@tax_table[,classification])}

  if(taxa_are_ordered){phyloseq_obj@tax_table <- phyloseq_obj@tax_table[,1:which(rank_names(phyloseq_obj) %in% classification)]
  } else {phyloseq_obj@tax_table <- phyloseq_obj@tax_table[,classification]}

  phyloseq_table <- melt_phyloseq(phyloseq_obj)
  otus <- eval(parse(text=paste0("dcast(phyloseq_table, with=FALSE , Sample ~ ", paste(colnames(phyloseq_obj@tax_table), collapse = '+'), ", value.var = 'Abundance', fun = sum)")))
  otus <- as.matrix(otus, rownames = 1)
  taxa <- eval(parse(text=paste0("setkey(unique(phyloseq_table[, c('", paste(colnames(phyloseq_obj@tax_table), collapse = "', '"), "')]), ", paste(colnames(phyloseq_obj@tax_table), collapse = ', '), ")")))
  taxa <- as.matrix(taxa, rownames = colnames(otus))
  
  phyloseq_obj <- phyloseq(otu_table(t(otus), taxa_are_rows = TRUE), tax_table(taxa), sample_data(phyloseq_obj@sam_data))
  
  return(phyloseq_obj)
}