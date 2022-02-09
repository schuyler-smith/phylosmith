#' Merge samples based on common factor within sample_data.
#' Function from the phylosmith-package.
#'
#' Inputs a phyloseq object and merges the samples that meet the
#' specified criteria into a single sample. This is meant for replicates, or
#' samples statistically proven to not be significantly different. This should
#' be used with caution as it may be a misleading representation of the data.
#' @useDynLib phylosmith
#' @usage conglomerate_samples(phyloseq_obj, treatment, subset = NULL)
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
#' @seealso \code{\link[phyloseq:merge_samples]{merge_samples()}}
#' @export
#' @return phyloseq-object
#' @examples conglomerate_samples(soil_column, treatment = c('Day', 'Matrix', 'Treatment'))

conglomerate_samples <-
  function(phyloseq_obj,
           treatment,
           subset = NULL) {
  if (!inherits(phyloseq_obj, "phyloseq")) {
    stop("`phyloseq_obj` must be a phyloseq-class\n        object",
         call. = FALSE)
  }
  if (is.null(access(phyloseq_obj, "sam_data"))) {
    stop("`phyloseq_obj` must contain\n        sample_data() information",
         call. = FALSE)
  }
  treatment <- phylosmith:::check_index_treatment(phyloseq_obj, treatment)
  if (any(!(treatment %in% colnames(access(phyloseq_obj, "sam_data"))))) {
    stop("`treatment` must be at least one column\n        name, or index, from the sample_data()",
         call. = FALSE)
  }
  if (!(is.null(access(phyloseq_obj, "phy_tree")))) {
    phylo_tree <- access(phyloseq_obj, "phy_tree")
  }  else {
    phylo_tree <- FALSE
  }
  if (!(is.null(access(phyloseq_obj, "refseq")))) {
    refseqs <- access(phyloseq_obj, "refseq")
  }  else {
    refseqs <- FALSE
  }
  if (!(is.null(access(phyloseq_obj, "tax_table")))) {
    tax_table <- access(phyloseq_obj, "tax_table")
  }  else {
    tax_table <- FALSE
  }

  original_levels <- lapply(access(phyloseq_obj, "sam_data"), levels)
  treatment_name <- paste(treatment, collapse = sep)
  phyloseq_obj <- merge_treatments(phyloseq_obj, treatment)
  sam_data <- data.table(as(access(phyloseq_obj, "sam_data"), 'data.frame'), keep.rownames = "Sample")
  set(sam_data, j = "Sample", value = sam_data[[treatment_name]])
  if(length(treatment) > 1){sam_data[, (treatment_name) := NULL]}
  otu_table <- data.table(t(as(access(phyloseq_obj, "otu_table"), 'matrix')), keep.rownames = "Sample")
  otu_table[, Sample := sam_data[['Sample']]]
  if(!is.null(subset)){
    sam_data <- sam_data[sam_data[, Reduce(`|`, lapply(.SD, `%in%`, subset)), .SDcols = c(treatment)]]
    otu_table <- otu_table[Sample %in% sam_data$Sample]
  }
  otu_table <- otu_table[, lapply(.SD, mean, na.rm=TRUE), by=Sample]
  sam_data <- sam_data[, .SD[1], Sample]

  phyloseq_obj <- phyloseq(
    otu_table(t(as.matrix(otu_table, rownames = "Sample")),
            taxa_are_rows = TRUE), tax_table,
    phyloseq::sample_data(data.frame(sam_data, row.names = 1)))

  if (!(is.logical(refseqs))) {
    phyloseq_obj <- phyloseq(phyloseq_obj@otu_table, phyloseq_obj@tax_table,
                             phyloseq_obj@sam_data, refseq(refseqs))
  }
  if (!(is.logical(phylo_tree))) {
    phy_tree(phyloseq_obj) <- phylo_tree
  }
  return(phyloseq_obj)
}


#' Conglomerate taxa by sample on a given classification level.
#' Function from the phylosmith-package.
#'
#' Conglomerate taxa by sample on a given classification level from the
#' tax_table.
#' @useDynLib phylosmith
#' @usage conglomerate_taxa(phyloseq_obj, classification, hierarchical = TRUE,
#' use_taxonomic_names = TRUE)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}} with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param classification Column name as a \code{string} or \code{numeric} in
#' the \code{\link[phyloseq:tax_table]{tax_table}} for the factor to
#' conglomerate by.
#' @param hierarchical Whether the order of factors in the tax_table represent
#' a decreasing hierarchy (TRUE) or are independent (FALSE). If FALSE, will
#' only return the factor given by \code{classification}.
#' @param use_taxonomic_names If (TRUE), will use the hierarchical taxonomic
#' name from Domain to classification level. If (FALSE) will uses the
#' classification level with a numerical assignment.
#' @seealso \code{\link[phyloseq:tax_glom]{tax_glom()}}
#' @export
#' @return phyloseq-object
#' @examples conglomerate_taxa(soil_column, 'Phylum')

conglomerate_taxa <- function(phyloseq_obj,
                              classification,
                              hierarchical = TRUE,
                              use_taxonomic_names = TRUE) {
  if (!inherits(phyloseq_obj, "phyloseq")) {
    stop("`phyloseq_obj` must be a phyloseq-class
        object", call. = FALSE)
  }
  if (is.null(access(phyloseq_obj, 'tax_table'))) {
    stop("`phyloseq_obj` must contain tax_table()
        information",
         call. = FALSE)
  }
  phyloseq_obj <- check_TaR(phyloseq_obj)
  classification <- check_index_classification(phyloseq_obj, classification)
  if (any(!(classification %in% colnames(access(
    phyloseq_obj, 'tax_table'
  )))))
  {
    stop("`classification` must be a column from the
        tax_table()",
         call. = FALSE)
  }
  if (!(is.logical(hierarchical))) {
    stop("`hierarchical` must be either `TRUE` or
        `FALSE`", call. = FALSE)
  }
  if (!(is.null(access(phyloseq_obj, 'phy_tree')))) {
    phylo_tree <- access(phyloseq_obj, 'phy_tree')
  } else {
    phylo_tree <- FALSE
  }
  if (!(is.null(access(phyloseq_obj, 'refseq')))) {
    refseqs <- access(phyloseq_obj, 'refseq')
  } else {
    refseqs <- FALSE
  }
  if (!(is.null(access(phyloseq_obj, 'sam_data')))) {
    sam <- access(phyloseq_obj, 'sam_data')
  } else {
    sam <- FALSE
  }
  sample_order <- sample_names(phyloseq_obj)
  taxa <- as(access(phyloseq_obj,'tax_table'), 'matrix')
  taxa <- as.data.table(taxa, keep.rownames='OTU')
  for (j in seq_len(ncol(taxa))) set(taxa,which(is.na(taxa[[j]])),j, "Unclassified")
  for (j in seq_len(ncol(taxa))) set(taxa,which(taxa[[j]] %in% "Incertae Sedis"),j, "Unclassified")
  for (j in seq_len(ncol(taxa))) set(taxa,which(taxa[[j]] %in% ""),j, "Unclassified")
  
  if (hierarchical) {
    if(which(colnames(taxa) %in% classification) != length(colnames(taxa))){
      taxa[, `:=`(seq(ncol(taxa))[-c(seq(which(colnames(taxa) %in% 
                                                 classification)))], NULL)]
    }
    set(taxa, i = taxa[, .I[get(classification) %in% "Unclassified"]], j = classification,
        value = paste0('Unclassified ',
                       unlist(taxa[ taxa[, .I[get(classification) %in% "Unclassified"]],
                                    colnames(taxa)[which(colnames(taxa) %in% classification) - 1],
                                    with=F])))
    for(i in seq_len(ncol(taxa)-1)){
      set(taxa, i = taxa[, .I[get(classification) %in% "Unclassified Unclassified"]], j = classification,
          value = paste0('Unclassified ',
                         unlist(taxa[ taxa[, .I[get(classification) %in% "Unclassified Unclassified"]],
                                      colnames(taxa)[which(colnames(taxa) %in% classification) - i],
                                      with=F])))
    }
  } else {
    taxa[, seq(ncol(taxa))[-c(which(colnames(taxa) %in% classification))][-1] := NULL]
  }
  
  otus <- as(access(phyloseq_obj,'otu_table'), 'matrix')
  otus <- as.data.table(otus, keep.rownames='OTU')
  otus[, OTU := taxa[[classification]]]
  otus <- otus[, lapply(.SD, sum), by = OTU]
  otus <- as.matrix(otus, rownames = 'OTU')
  taxa[, OTU := get(classification)]
  taxa <- as.matrix(unique(taxa, by = 'OTU'), rownames = 'OTU')
  
  phyloseq_obj <- phyloseq(otu_table(otus[, sample_order], taxa_are_rows = TRUE),
                           tax_table(taxa))
  if(!use_taxonomic_names){
    taxa_names(phyloseq_obj) <- paste0(classification,"_",
                                       seq(length(taxa_names(phyloseq_obj))))
  }

  if (!(is.logical(sam))) {
    sample_data(phyloseq_obj) <- sam
  }
  if (!(is.logical(phylo_tree))) {
    warning('trees cannot be preserved after taxa conglomeration')
  }
  if (!(is.logical(refseqs))) {
    warning('reference sequences cannot be preserved after taxa conglomeration')
  }
  rm(list = c('taxa', 'sam', 'phylo_tree', 'refseqs'))
  gc()
  return(phyloseq_obj)
}

#' Melt a phyloseq object into a data.table.
#' Function from the phylosmith-package.
#'
#' melt_phyloseq inputs a phyloseq object and melts its otu_table, taxa_tale,
#' and sample_Data into a single into a data.table.
#' @useDynLib phylosmith
#' @usage melt_phyloseq(phyloseq_obj)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}} with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @seealso \code{\link[phyloseq:psmelt]{psmelt()}}
#' @export
#' @return data.table
#' @examples melt_phyloseq(soil_column)

melt_phyloseq <- function(phyloseq_obj) {
  options(warn = -1)
  if (!(inherits(phyloseq_obj, "phyloseq") | inherits(phyloseq_obj, "otu_table"))) {
    stop("`phyloseq_obj` must be a phyloseq-class object or otu_table to melt",
         call. = FALSE)
  }
  phyloseq_obj <- check_TaR(phyloseq_obj)
  melted_phyloseq <- data.table(as(
    access(phyloseq_obj, 'otu_table'), "matrix"
  ), keep.rownames = 'OTU')
  if (!(is.null(access(phyloseq_obj, 'tax_table')))) {
    taxa <- data.table(as(access(phyloseq_obj, 'tax_table'), 'matrix'), keep.rownames = 'OTU')
  } else {taxa <- NULL}

  if (!(is.null(access(phyloseq_obj, 'sam_data')))) {
    sample_data <-
      data.table(data.frame(access(phyloseq_obj, 'sam_data'),
                            stringsAsFactors = FALSE), keep.rownames = TRUE)
    sample_data[, Sample := NULL]
    setnames(sample_data, 'rn', 'Sample')
  } else {
    sample_data <- data.table(Sample = sample_names(phyloseq_obj))
  }
  if (!(is.null(taxa))) {
    melted_phyloseq <- merge(melted_phyloseq, taxa, by = "OTU")
  }
  melted_phyloseq <- melt(melted_phyloseq, variable.name = 'Sample', value.name = 'Abundance', id.vars = colnames(taxa))
  melted_phyloseq <-
    merge(melted_phyloseq, sample_data, by = "Sample")

  setorder(melted_phyloseq, -Abundance)
  rm(list = c('taxa', 'sample_data'))
  gc()
  return(melted_phyloseq)
}

#' Combine meta-data columns. Function from the phylosmith-package.
#'
#' Combines multiple columns of a \code{\link[phyloseq]{phyloseq-class}}
#' object \code{\link[phyloseq:sample_data]{sample_data}} into a
#' single-variable column.
#' @useDynLib phylosmith
#' @usage merge_treatments(phyloseq_obj, treatments)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample.
#' @param treatments A vector of any number of column names as \code{string}s or \code{numeric}s
#' in the \code{\link[phyloseq:sample_data]{sample_data}} that are to be
#' combined.
#' @export
#' @return phyloseq-object
#' @examples merge_treatments(soil_column, c('Matrix', 'Treatment', 'Day'))

merge_treatments <- function(phyloseq_obj, treatments) {
  if (!inherits(phyloseq_obj, "phyloseq")) {
    stop("`phyloseq_obj` must be a phyloseq-class
        object", call. = FALSE)
  }
  if (is.null(access(phyloseq_obj, 'sam_data'))) {
    stop("`phyloseq_obj` must contain sample_data()
        information",
         call. = FALSE)
  }
  treatments <- check_index_treatment(phyloseq_obj, treatments)
  if (any(!(treatments %in% colnames(access(
    phyloseq_obj, 'sam_data'
  ))))) {
    stop(
      "`treatments` must be at least two column
        names, or indices, from the sample_data()",
      call. = FALSE
    )
  }
  treatment_classes <- setDT(as(access(phyloseq_obj, 'sam_data')[,
                                                                 colnames(access(phyloseq_obj, 'sam_data')) %in% treatments],
                                "data.frame"))
  treatment_name <- paste(treatments, collapse = sep)
  order <- apply(
    eval(parse(text = paste0(
      "expand.grid(",
      paste0(
        paste0(
          "levels(factor(access(phyloseq_obj, 'sam_data')[['",
          treatments,
          "']]))",
          collapse = ', '
        )
      ), ")"
    ))),
    1,
    FUN = function(combination) {
      paste0(combination, collapse = ' ')
    }
  )
  eval(parse(
    text = paste0(
      "treatment_classes[, '",
      treatment_name,
      "' := as.character(paste(",
      paste(treatments, collapse = ', '),
      ', sep = " "), by = treatment_classes)]'
    )
  ))
  sample_data(phyloseq_obj)[[treatment_name]] <- factor(treatment_classes[[treatment_name]], levels = order)
  return(phyloseq_obj)
}

#' Transform abundance data in an \code{otu_table} to relative abundance,
#' sample-by-sample. Function from the phylosmith-package.
#'
#' Transform abundance data into relative abundance, i.e. proportional data.
#' This is an alternative method of normalization and may not be appropriate
#' for all datasets, particularly if your sequencing depth varies between
#' samples.
#' @useDynLib phylosmith
#' @usage relative_abundance(phyloseq_obj, sig.fig=4)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param phyloseq_obj Number of significant figures to round to.
#' @export
#' @return phyloseq-object
#' @examples relative_abundance(soil_column, 3)

relative_abundance <- function(phyloseq_obj, sig.fig=4) {
  if (!inherits(phyloseq_obj, "phyloseq")) {
    stop("`phyloseq_obj` must be a phyloseq-class
            object", call. = FALSE)
  }
  phyloseq_obj <- check_TaR(phyloseq_obj)
  abundance_table <- access(phyloseq_obj, 'otu_table')
  abundance_table <-
    apply(
      abundance_table,
      2,
      FUN = function(c) {
        round(c / sum(c), sig.fig)
      }
    )
  abundance_table[is.na(abundance_table)] <- 0
  otu_table(phyloseq_obj) <-
    otu_table(abundance_table, taxa_are_rows = TRUE)
  return(phyloseq_obj)
}

#' Re-orders the samples of a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and changes the
#' order of sample index either based on the metadata, or a given order. If
#' metadata columns are used, they will take the factor order.
#' @useDynLib phylosmith
#' @usage set_sample_order(phyloseq_obj, sort_on)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param sort_on Column name as a \code{string} or \code{numeric} in the
#' \code{\link[phyloseq:sample_data]{sample_data}}, or vector of sample names
#' or indices in particular order.
#' @export
#' @return phyloseq object
#' @examples
#' set_sample_order(soil_column, c('Matrix', 'Treatment'))

set_sample_order <- function(phyloseq_obj, sort_on) {
  if (!inherits(phyloseq_obj, "phyloseq")) {
    stop("`phyloseq_obj` must be a phyloseq-class
         object", call. = FALSE)
  }
  if (!(is.null(access(phyloseq_obj, 'phy_tree')))) {
    phylo_tree <- access(phyloseq_obj, 'phy_tree')
  } else {
    phylo_tree <- FALSE
  }
  if (!(is.null(access(phyloseq_obj, 'refseq')))) {
    refseqs <- access(phyloseq_obj, 'refseq')
  } else {
    refseqs <- FALSE
  }
  if (!(is.null(access(phyloseq_obj, 'sam_data')))) {
    sam <- access(phyloseq_obj, 'sam_data')
  } else {
    sam <- FALSE
  }
  if (!(is.null(access(phyloseq_obj, 'tax_table')))) {
    tax <- access(phyloseq_obj, 'tax_table')
  } else {
    tax <- FALSE
  }
  if (!(is.null(access(phyloseq_obj, 'otu_table')))) {
    otu <- access(phyloseq_obj, 'otu_table')
  } else {
    otu <- FALSE
  }
  metadata <- as(access(phyloseq_obj, 'sam_data'), 'data.frame')
  metadata <- data.table(samples = rownames(metadata), metadata)
  if (length(sort_on) < nsamples(phyloseq_obj)) {
    sort_on <- check_index_treatment(phyloseq_obj, sort_on)
    if (any(!(sort_on %in% colnames(access(
      phyloseq_obj, 'sam_data'
    ))))) {
      stop(
        "`sort_on` must be at least one column
           name, or index, from the sample_data(), or a vector of a set
           order of sample names or indices",
        call. = FALSE
      )
    }
    eval(parse(text = paste0(
      'setkey(metadata, ', paste(sort_on, collapse = ', '), ')'
    )))
  } else {
    if (is.character(sort_on)) {
      if (!(any(sort_on %in% sample_names(phyloseq_obj)))) {
        stop("`sort_on` must be at least one column\n        name, or index, from the sample_data(), or a vector of a set\n        order of sample names or indices",
             call. = FALSE)
      }
      metadata <- metadata[match(sort_on, metadata$samples), ]
    } else {
      if (is.numeric(sort_on)){
        metadata <- metadata[sort_on, ]
      }
    }
  }
  otu <- otu[, metadata$samples]
  metadata <- sample_data(data.frame(metadata, row.names = 1))

  phyloseq_obj <- phyloseq(otu, sam)
  if (!(is.logical(tax))) {
    tax_table(phyloseq_obj) <- tax
  }
  if (!(is.logical(refseqs))) {
    phyloseq_obj <-
      phyloseq(
        phyloseq_obj@otu_table,
        phyloseq_obj@tax_table,
        phyloseq_obj@sam_data,
        refseq(refseqs)
      )
  }
  if (!(is.logical(phylo_tree))) {
    phy_tree(phyloseq_obj) <- phylo_tree
  }
  return(phyloseq_obj)
}

#' set_treatment_levels
#'
#' Reorders the levels of a metadata factor column in a
#' \code{\link[phyloseq]{phyloseq-class}} object
#' \code{\link[phyloseq:sample_data]{sample_data}}.
#' @useDynLib phylosmith
#' @usage set_treatment_levels(phyloseq_obj, treatment, order)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param treatment Column name as a \code{string} or \code{numeric} in the
#' \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param order The order of factors in \code{treatment} column as a vector of
#' strings. If assigned "numeric" will set ascending numerical order.
#' @export
#' @return phyloseq-object
#' @examples set_treatment_levels(soil_column, treatment = 'Matrix',
#' order = c('Manure', 'Soil', 'Effluent'))
#' set_treatment_levels(soil_column, 'Day', 'numeric')

set_treatment_levels <- function(phyloseq_obj, treatment, order) {
  if (!inherits(phyloseq_obj, "phyloseq")) {
    stop("`phyloseq_obj` must be a phyloseq-class
        object", call. = FALSE)
  }
  if (is.null(access(phyloseq_obj, 'sam_data'))) {
    stop("`phyloseq_obj` must contain sample_data()
        information",
         call. = FALSE)
  }
  treatment <- check_index_treatment(phyloseq_obj, treatment)
  if (length(treatment) > 1 |
      !(treatment %in% colnames(access(phyloseq_obj, 'sam_data')))) {
    stop("`treatment` must be a single column name,
            or index, from the sample_data()",
         call. = FALSE)
  }
  if (order[1] == 'numeric') {
    order <- as.character(sort(as.numeric(unique(
      as.character(access(phyloseq_obj, 'sam_data')[[treatment]])
    ))))
  }
  sample_data(phyloseq_obj)[[treatment]] <-
    factor(access(phyloseq_obj,
                  'sam_data')[[treatment]], levels = order)
  return(phyloseq_obj)
}

#' Filter taxa based on proportion of samples they are observed in.
#' Function from the phylosmith-package.
#'
#' Inputs a phyloseq object and finds which taxa are seen in a
#' given proportion of samples, either in the entire dataset, by treatment, or
#' a particular treatment of interest.
#' @useDynLib phylosmith
#' @usage taxa_filter(phyloseq_obj, treatment = NULL, subset = NULL,
#' frequency = 0, below = FALSE, drop_samples = FALSE)
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
#' @param frequency The proportion of samples the taxa is found in.
#' @param below Does frequency define the minimum (\code{FALSE}) or maximum
#' (\code{TRUE}) proportion of samples the taxa is found in.
#' @param drop_samples Should the function remove samples that that are empty
#' after removing taxa filtered by frequency (\code{TRUE}).
#' @export
#' @return phyloseq-object
#' @examples taxa_filter(soil_column, frequency = 0.8)
#' taxa_filter(soil_column, treatment = c('Matrix', 'Treatment'),
#' subset = 'Soil Amended', frequency = 0.8)

taxa_filter <-
  function(phyloseq_obj,
           treatment = NULL,
           subset = NULL,
           frequency = 0,
           below = FALSE,
           drop_samples = FALSE) {
  if (!inherits(phyloseq_obj, "phyloseq")) {
    stop("`phyloseq_obj` must be a phyloseq-class object",
         call. = FALSE)
  }
  if (is.null(access(phyloseq_obj, 'sam_data'))) {
    stop("`phyloseq_obj` must contain sample_data()
          information",
         call. = FALSE)
  }
  treatment <- check_index_treatment(phyloseq_obj, treatment)
  if (!(is.null(treatment)) &
      any(!(treatment %in% colnames(access(
        phyloseq_obj, 'sam_data'
      ))))) {
    stop(
      "`treatment` must be at least one column name, or
      index, from the sample_data()",
      call. = FALSE
    )
  }
  if (!(is.numeric(frequency)) |
      !(frequency >= 0 & frequency <= 1)) {
    stop("`frequency` must be a numeric value between
      0 and 1", call. = FALSE)
  }
  if (!(is.logical(below))) {
    stop("`below` must be either `TRUE` or `FALSE`",
         call. = FALSE)
  }
  if (!(is.logical(drop_samples))) {
    stop("`drop_samples` must be either `TRUE` or `FALSE`",
         call. = FALSE)
  }
  phyloseq_obj <- check_TaR(phyloseq_obj)
  phyloseq_table <- data.table(as(
    access(phyloseq_obj, 'otu_table'), "matrix"
  ), keep.rownames = 'OTU')
  phyloseq_table <- melt(phyloseq_table, id.vars = 'OTU',
                         variable.name = 'Sample', value.name = 'Abundance')
  phyloseq_table <- phyloseq_table[Abundance > 0]
  if (!(is.null(treatment))) {
    phyloseq_obj <- merge_treatments(phyloseq_obj, treatment)
    treatment_name <- paste(treatment, collapse = sep)
    sam_data <- data.table(as(access(phyloseq_obj, 'sam_data'), 'data.frame'), keep.rownames = "Sample")
    sam_data <- sam_data[,c("Sample", treatment, treatment_name), with= FALSE]
    if (!(is.null(subset))) {
      sam_data <- sam_data[sam_data[, Reduce(`|`, lapply(.SD, `%in%`, subset)),
                                    .SDcols = c(treatment, treatment_name)]]
      phyloseq_table <- phyloseq_table[Sample %in% sam_data$Sample]
      phyloseq_obj <- prune_samples((sample_names(phyloseq_obj) %in% sam_data$Sample), phyloseq_obj)
    }
    sample_counts <- phyloseq_table[, .(n_samples = uniqueN(Sample))]
    taxa_counts <- phyloseq_table[, .(count = .N), by = c('OTU')]
    taxa_counts[, proportion := count/unlist(sample_counts)]
    phyloseq_table <- phyloseq_table[OTU %in% taxa_counts[proportion > frequency]$OTU]

    treatment_classes <- unique(sam_data[[treatment_name]])
    taxa <- data.table()
    for(trt in sam_data){
      sub_table <- phyloseq_table[Sample %in% sam_data[get(treatment_name) %in% trt]$Sample]
      sample_counts <- sub_table[, .(n_samples = uniqueN(Sample))]
      taxa_counts <- sub_table[, .(count = .N), by = c('OTU')]
      taxa_counts[, proportion := count/unlist(sample_counts)]
      taxa <- rbind(taxa, taxa_counts)
    }
    rm('sub_table', 'sam_data')
  } else {
    sample_counts <- phyloseq_table[, .(n_samples = uniqueN(Sample))]
    taxa_counts <- phyloseq_table[, .(count = .N), by = c('OTU')]
    taxa_counts[, n_samples := sample_counts$n_samples]
    taxa <- taxa_counts[, .(proportion = count/n_samples), by = c('OTU')]
  }
  rm('phyloseq_table', 'sample_counts', 'taxa_counts')
  if(below){
    taxa <- taxa[proportion <= frequency]
  } else {
    taxa <- taxa[proportion >= frequency]
  }
  taxa <- unique(taxa[['OTU']])
  phyloseq_obj <- taxa_prune(phyloseq_obj,
                             taxa_names(phyloseq_obj)[!(taxa_names(phyloseq_obj) %in% taxa)])
  if (drop_samples == TRUE) {
    phyloseq_obj <- prune_samples(sample_sums(phyloseq_obj) > 0, phyloseq_obj)
  }
  rm('taxa')
  gc()
  return(phyloseq_obj)
}

#' Remove specific taxa from all samples based on a given classification level.
#' Function from the phylosmith-package.
#'
#' Prunes taxa from phyloseq objects based on taxanomic names.
#' @useDynLib phylosmith
#' @usage taxa_prune(phyloseq_obj, taxa_to_remove,
#' classification=NULL, na.rm=FALSE)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}} with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param taxa_to_remove A vector of the classification names of the taxon to
#' be removed. If you are using ASVs and want to remove specific sequences, this
#' should be set OTU or NULL. By default, it will search for the names in all
#' taxanomic ranks
#' @param classification Column name as a \code{string} or \code{numeric} in
#' the \code{\link[phyloseq:tax_table]{tax_table}} for the factor to
#' conglomerate by.
#' @param na.rm if TRUE, and classification is specified, will remove taxon at
#' the classificaiton level that have NA values.
#' @seealso \code{\link[phyloseq:tax_glom]{tax_glom()}}
#' @export
#' @return phyloseq-object
#' @examples taxa_prune(soil_column, 'Firmicutes', 'Phylum')

taxa_prune <- function(phyloseq_obj,
                       taxa_to_remove,
                       classification = NULL,
                       na.rm = FALSE) {
  if (!inherits(phyloseq_obj, "phyloseq")) {
    stop("`phyloseq_obj` must be a phyloseq-class
        object", call. = FALSE)
  }
  if (is.null(access(phyloseq_obj, 'tax_table'))) {
    stop("`phyloseq_obj` must contain tax_table()
        information",
         call. = FALSE)
  }
  classification <- check_index_classification(phyloseq_obj, classification)
  if (any(!(classification %in% colnames(access(
    phyloseq_obj, 'tax_table'
  )))))
  {
    stop("`classification` must be a column from the
        tax_table()",
         call. = FALSE)
  }

  taxa <- as(access(phyloseq_obj,'tax_table'), 'matrix')
  taxa <- as.data.table(taxa, keep.rownames='OTU')
  if(is.null(classification)){classification <- colnames(taxa)}
  for (rank in classification) {
    set(taxa, i = which(taxa[[rank]] %in% taxa_to_remove), j = rank, value = 'REMOVE')
  }
  if(length(classification) == 1 & na.rm){
    set(taxa, i = which(is.na(taxa[[classification]])), j = classification, value = 'REMOVE')
  }

  taxa <- taxa[!(taxa[, Reduce(`|`, lapply(.SD, `%in%`, 'REMOVE')), .SDcols = classification])]
  taxa <- as.matrix(taxa, rownames = 'OTU')
  tax_table(phyloseq_obj) <- tax_table(taxa)

  rm(list = c('taxa'))
  gc()
  return(phyloseq_obj)
}

#' Extracts specific taxa from all samples based on a given classification level.
#' Function from the phylosmith-package.
#' 
#' Extracts taxa from phyloseq objects based on taxanomic names.
#' @useDynLib phylosmith
#' @usage taxa_extract(phyloseq_obj, taxa_to_extract, classification=NULL)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}} with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param taxa_to_extract A vector of the classification names of the taxon to
#' be extracted. If you are using ASVs and want to extract specific sequences, this
#' should be set OTU or NULL. By default, it will search for the names in all
#' taxanomic ranks
#' @param classification Column name as a \code{string} or \code{numeric} in
#' the \code{\link[phyloseq:tax_table]{tax_table}} for the factor to
#' conglomerate by.
#' @seealso \code{\link[phyloseq:tax_glom]{tax_glom()}}
#' @export
#' @return phyloseq-object
#' @examples taxa_extract(soil_column, 'Firmicutes', 'Phylum')

taxa_extract <- function(phyloseq_obj,
                         taxa_to_extract,
                         classification = NULL) {
  if (!inherits(phyloseq_obj, "phyloseq")) {
    stop("`phyloseq_obj` must be a phyloseq-class
        object", call. = FALSE)
  }
  if (is.null(access(phyloseq_obj, 'tax_table'))) {
    stop("`phyloseq_obj` must contain tax_table()
        information",
         call. = FALSE)
  }
  classification <- check_index_classification(phyloseq_obj, classification)
  if (any(!(classification %in% colnames(access(
    phyloseq_obj, 'tax_table'
  )))))
  {
    stop("`classification` must be a column from the
        tax_table()",
         call. = FALSE)
  }
  taxa <- as(access(phyloseq_obj,'tax_table'), 'matrix')
  taxa <- as.data.table(taxa, keep.rownames='OTU')
  if(is.null(classification)){classification <- colnames(taxa)}
  taxa <- cbind(taxa, ID = taxa$OTU)
  for (rank in classification) {
    set(taxa, i = unlist(lapply(tolower(taxa_to_extract), grep, tolower(taxa[[rank]]))), j = 'ID', value = 'EXTRACT')
  }
  taxa <- taxa[ID == 'EXTRACT'][, ID := NULL]
  if(nrow(taxa) == 0){
    stop("no taxa found matching `taxa_to_extract`",
         call. = FALSE)
  }
  taxa <- as.matrix(taxa, rownames = 'OTU')
  tax_table(phyloseq_obj) <- tax_table(taxa)
  
  rm(list = c('taxa'))
  gc()
  return(phyloseq_obj)
}

