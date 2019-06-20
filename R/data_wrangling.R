#' Find taxa shared between treatments of a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and finds which taxa
#' are shared between all of the specified treatments from data in the
#' \code{\link[phyloseq:sample_data]{sample_data()}}), or every sample in the
#' dataset.
#' @useDynLib phylosmith
#' @usage common_taxa(phyloseq_obj, treatment = NULL, subset = NULL, n = 'all')
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
#' @param n Number of treatment groups that need to share the taxa to be
#' considered a common taxa.
#' @seealso \code{\link{unique_taxa}}
#' @export
#' @return vector
#' @examples common_taxa(soil_column, treatment = 'Treatment',
#' subset = 'Amended', n = 'all')

common_taxa <-
  function(phyloseq_obj,
           treatment = NULL,
           subset = NULL,
           n = 'all') {
    if (!inherits(phyloseq_obj, "phyloseq")) {
      stop("`phyloseq_obj` must be a phyloseq-class object",
           call. = FALSE)
    }
    if (is.null(access(phyloseq_obj, 'sam_data'))) {
      stop("`phyloseq_obj` must contain sample_data()
        information",
           call. = FALSE)
    }
    treatment <- check_numeric_treatment(phyloseq_obj, treatment)
    if (any(!(treatment %in% colnames(access(
      phyloseq_obj, 'sam_data'
    ))))) {
      stop(
        "`treatment` must be at least one column name, or
        index, from the sample_data()",
        call. = FALSE
      )
    }
    if (!(is.numeric(n)) & n != 'all') {
      stop(
        "`n` must be either 'all' or a numeric value less
        than the number of treatments being compared",
        call. = FALSE
      )
    }
    phyloseq_obj <- phyloseq(
      access(phyloseq_obj, 'otu_table'),
      access(phyloseq_obj, 'tax_table'),
      access(phyloseq_obj, 'sam_data')
    )
    phyloseq_obj <-
      taxa_filter(phyloseq_obj, treatment, subset = subset)
    treatment_name <- paste(treatment, collapse = sep)
    treatment_classes <- unique(access(phyloseq_obj,
                                       'sam_data')[[treatment_name]])
    if (n == 'all') {
      n <- length(treatment_classes)
    }
    if (n > length(treatment_classes)) {
      stop(
        "`n` must be either 'all' or a numeric value less
        than the number of treatments being compared",
        call. = FALSE
      )
    }
    seen_taxa <- lapply(
      treatment_classes,
      FUN = function(trt) {
        taxa_names(taxa_filter(
          phyloseq_obj,
          treatment = treatment,
          subset = trt,
          frequency = 0
        ))
      }
    )
    if (n == 0) {
      seen_taxa <-
        tryCatch(
          taxa_names(taxa_filter(phyloseq_obj, frequency = .99)),
          error = function(e) {
            stop(
              "
      No taxa seen in every sample.
      To get a list of taxa seen in a certain proportion of samples use:
      taxa_names(taxa_filter(phyloseq_obj, frequency = x))"
            )
          }
        )
      n <- 1
    }
    names(seen_taxa) <- treatment_classes
    taxa_counts <- table(unlist(seen_taxa))
    shared_taxa <- names(taxa_counts[taxa_counts == round(n)])
    return(shared_taxa)
  }

#' Merge samples based on common factor within sample_data.
#' Function from the phylosmith-package.
#'
#' Inputs a phyloseq object and merges the samples that meet the
#' specified criteria into a single sample. This is meant for replicates, or
#' samples statistically proven to not be significantly different. This should
#' be used with caution as it may be a misleading representation of the data.
#' @useDynLib phylosmith
#' @usage conglomerate_samples(phyloseq_obj, treatment, subset = NULL,
#' merge_on = treatment)
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
#' @param merge_on Defines which variable the data is merged according to.
#' This needs to be a column name as a \code{string} or \code{numeric} in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @seealso \code{\link[phyloseq:merge_samples]{merge_samples()}}
#' @export
#' @return phyloseq-object
#' @examples conglomerate_samples(soil_column, treatment = c('Matrix', 'Treatment'), merge_on = 'Day')

conglomerate_samples <-
  function(phyloseq_obj,
           treatment,
           subset = NULL,
           merge_on = treatment) {
    if (!inherits(phyloseq_obj, "phyloseq")) {
      stop("`phyloseq_obj` must be a phyloseq-class
        object", call. = FALSE)
    }
    if (is.null(access(phyloseq_obj, 'sam_data'))) {
      stop("`phyloseq_obj` must contain
        sample_data() information",
           call. = FALSE)
    }
    treatment <- check_numeric_treatment(phyloseq_obj, treatment)
    if (any(!(treatment %in% colnames(access(
      phyloseq_obj, 'sam_data'
    ))))) {
      stop(
        "`treatment` must be at least one column
        name, or index, from the sample_data()",
        call. = FALSE
      )
    }
    merge_on <- check_numeric_treatment(phyloseq_obj, merge_on)
    if (any(!(merge_on %in% colnames(access(
      phyloseq_obj, 'sam_data'
    ))))) {
      stop("`merge_on` must be at least one column
        name, or index, from the sample_data()",
           call. = FALSE)
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
    original_levels <-
      lapply(access(phyloseq_obj, 'sam_data'), levels)
    phyloseq_obj <- phyloseq(
      access(phyloseq_obj, 'otu_table'),
      access(phyloseq_obj, 'tax_table'),
      access(phyloseq_obj, 'sam_data')
    )
    merge_on_name <- paste(merge_on, collapse = sep)

    phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, subset)
    phyloseq_obj <- merge_treatments(phyloseq_obj, merge_on_name)
    treatment_name <- paste(treatment, collapse = sep)
    merge_sample_levels <-
      as.character(unique(sort(unlist(
        access(phyloseq_obj,
               'sam_data')[[merge_on_name]]
      ))))

    treatment_classes <- sort(unique(access(phyloseq_obj,
                                            'sam_data')[[treatment_name]]))
    treatment_classes <-
      eval(parse(
        text = paste0(
          'treatment_classes[grepl("',
          paste0(subset),
          '", treatment_classes)]'
        )
      ))
    phyloseq_table <- melt_phyloseq(phyloseq_obj)

    if (any(merge_on != treatment)) {
      merge_sample_levels <- paste(
        vapply(
          as.character(treatment_classes),
          rep,
          character(length(merge_sample_levels)),
          times = length(merge_sample_levels)
        ),
        rep(merge_sample_levels,
            length(treatment_classes)),
        sep = ' '
      )
      phyloseq_table[, 'Merged_Name' := do.call(paste0,
                                                list(phyloseq_table[[treatment_name]], ' ',
                                                     phyloseq_table[[merge_on_name]]))]
    } else {
      phyloseq_table[, 'Merged_Name' := phyloseq_table[[merge_on_name]]]
    }
    otu_tab <-
      dcast(
        phyloseq_table[, c('OTU', 'Abundance', 'Merged_Name'),
                       with = FALSE],
        Merged_Name ~ OTU,
        value.var = 'Abundance',
        fun.aggregate = mean
      )

    sub_phy <- do.call(merge_phyloseq,
                       lapply(
                         treatment_classes,
                         FUN = function(group) {
                           group_phy <- eval(parse(
                             text =
                               paste0(
                                 'subset_samples(taxa_filter(phyloseq_obj,
                    treatment), ',
                                 treatment_name,
                                 ' == "',
                                 group,
                                 '")'
                               )
                           ))
                           if (nsamples(group_phy) > 1) {
                             sub_phy <- merge_samples(group_phy, merge_on_name)
                             merge_names <-
                               rownames(access(sub_phy, 'sam_data'))
                             if (any(group != merge_names)) {
                               sample_names(sub_phy) <- paste0(group, ' ',
                                                               merge_names)
                             }
                             sam <-
                               as(access(sub_phy, 'sam_data'), 'data.frame')
                             sam[, treatment_name] <- group
                             for (i in treatment) {
                               sam[, i] <- unique(group_phy@sam_data[, i])
                             }
                             sam[, merge_on_name] <- factor(merge_names,
                                                            levels = levels(access(phyloseq_obj,
                                                                                   'sam_data')[[merge_on_name]]))
                             sample_data(sub_phy) <- sample_data(sam)
                           }
                           return(sub_phy)
                         }
                       ))
    phyloseq_obj <- tryCatch({
      eval(parse(
        text = paste0(
          'subset_samples(phyloseq_obj, !(',
          treatment_name,
          ' %in% treatment_classes))'
        )
      ))
    },
    error = function(e) {
      phyloseq_obj <- sub_phy
    },
    finally = {
      merge_phyloseq(phyloseq_obj, sub_phy)
    })
    phyloseq_obj <- phyloseq(
      otu_table(t(
        as.matrix(otu_tab[order(factor(otu_tab$Merged_Name,
                                       levels = merge_sample_levels)), ],
                  rownames = 'Merged_Name')
      ),
      taxa_are_rows = TRUE),
      access(phyloseq_obj, 'tax_table'),
      access(phyloseq_obj, 'sam_data')[order(factor(rownames(access(
        phyloseq_obj, 'sam_data'
      )),
      levels = merge_sample_levels)), ]
    )
    for (i in seq_along(original_levels)) {
      if (!(is.null(unname(unlist(
        original_levels[i]
      ))))) {
        phyloseq_obj <- set_treatment_levels(phyloseq_obj,
                                             names(original_levels)[i], unname(unlist(original_levels[i])))
      }
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


#' Conglomerate taxa by sample on a given classification level.
#' Function from the phylosmith-package.
#'
#' Conglomerate taxa by sample on a given classification level from the
#' tax_table.
#' @useDynLib phylosmith
#' @usage conglomerate_taxa(phyloseq_obj, classification, hierarchical = TRUE)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}} with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param classification Column name as a \code{string} or \code{numeric} in
#' the \code{\link[phyloseq:tax_table]{tax_table}} for the factor to
#' conglomerate by.
#' @param hierarchical Whether the order of factors in the tax_table represent
#' a decreasing heirarchy (TRUE) or are independant (FALSE). If FALSE, will
#' only return the factor given by \code{classification}.
#' @seealso \code{\link[phyloseq:tax_glom]{tax_glom()}}
#' @export
#' @return phyloseq-object
#' @examples conglomerate_taxa(soil_column, 'Phylum')

conglomerate_taxa <- function(phyloseq_obj,
                              classification,
                              hierarchical = TRUE) {
  if (!inherits(phyloseq_obj, "phyloseq")) {
    stop("`phyloseq_obj` must be a phyloseq-class
        object", call. = FALSE)
  }
  if (is.null(access(phyloseq_obj, 'tax_table'))) {
    stop("`phyloseq_obj` must contain tax_table()
        information",
         call. = FALSE)
  }
  classification <- check_numeric_classification(phyloseq_obj,
                                                 classification)
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
  phyloseq_obj <- phyloseq(access(phyloseq_obj, 'otu_table'),
                           access(phyloseq_obj, 'tax_table'))
  if (hierarchical) {
    tax_table(phyloseq_obj) <- access(phyloseq_obj,
                                      'tax_table')[, seq(which(rank_names(phyloseq_obj) %in%
                                                                 classification))]
  } else {
    tax_table(phyloseq_obj) <- access(phyloseq_obj,
                                      'tax_table')[, classification]
  }

  phyloseq_table <- melt_phyloseq(phyloseq_obj)
  otus <- eval(parse(
    text = paste0(
      "dcast(phyloseq_table, with = FALSE, Sample ~ ",
      paste(colnames(access(
        phyloseq_obj, 'tax_table'
      )), collapse = '+'),
      ", value.var = 'Abundance', fun = sum)"
    )
  ))
  otus <- as.matrix(otus, rownames = 1)
  taxa <-
    eval(parse(
      text = paste0(
        "setkey(unique(phyloseq_table[, c('",
        paste(colnames(access(
          phyloseq_obj, 'tax_table'
        )), collapse = "', '"),
        "')]), ",
        paste(colnames(access(
          phyloseq_obj, 'tax_table'
        )), collapse = ', '),
        ")"
      )
    ))
  taxa <-
    as.matrix(taxa, rownames = unlist(taxa[, .(col_test = do.call(paste, c(.SD, sep = "_")))]))

  phyloseq_obj <- phyloseq(otu_table(t(otus)[, sample_order], taxa_are_rows = TRUE),
                           tax_table(taxa))
  if (!(is.logical(sam))) {
    sample_data(phyloseq_obj) <- sam
  }
  if (!(is.logical(phylo_tree))) {
    warning('trees cannot be preserved after taxa conglomeration')
  }
  if (!(is.logical(refseqs))) {
    warning('reference sequences cannot be preserved after taxa conglomeration')
  }
  return(phyloseq_obj)
}

#' Melt a phyloseq object into a data.table.
#' Function from the phylosmith-package.
#'
#' melt_phyloseq inputs a phyloseq object and melts its ou_table, taxa_tale,
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
  if (!inherits(phyloseq_obj, "phyloseq")) {
    stop("`phyloseq_obj` must be a phyloseq-class object",
         call. = FALSE)
  }
  if (is.null(access(phyloseq_obj, 'tax_table')) &
      is.null(access(phyloseq_obj
                     , 'sam_data'))) {
    stop(
      "`phyloseq_obj` must contain either tax_table()
            and/or sample_data() information",
      call. = FALSE
    )
  }
  melted_phyloseq <-
    melt.data.table(data.table(as(
      access(phyloseq_obj,
             'otu_table'), "matrix"
    ), keep.rownames = TRUE), id.vars = 1)
  colnames(melted_phyloseq) <- c("OTU", "Sample", "Abundance")
  if (!(is.null(access(phyloseq_obj, 'tax_table')))) {
    taxa <- data.table(as(access(phyloseq_obj, 'tax_table'), 'matrix'),
                       OTU = taxa_names(phyloseq_obj))
  }
  sample_data <- NULL
  if (!(is.null(access(phyloseq_obj, 'sam_data')))) {
    sample_data <-
      data.table(data.frame(access(phyloseq_obj, 'sam_data'),
                            stringsAsFactors = FALSE))
  }
  if (is.null(sample_data)) {
    sample_data <- data.table(Sample = sample_names(phyloseq_obj))
  } else {
    sample_data[, 'Sample' := sample_names(phyloseq_obj)]
  }

  melted_phyloseq <-
    merge(melted_phyloseq, sample_data, by = "Sample")
  if (!(is.null(access(phyloseq_obj, 'tax_table')))) {
    melted_phyloseq <- merge(melted_phyloseq, taxa, by = "OTU")
  }
  melted_phyloseq <-
    melted_phyloseq[order(melted_phyloseq$Abundance
                          , decreasing = TRUE),]
  return(melted_phyloseq)
}

#' Combine meta-data columns. Function from the phylosmith-package.
#'
#' Combines multiple columns of a \code{\link[phyloseq]{phyloseq-class}}
#' object \code{\link[phyloseq:sample_data]{sample_data}} into a
#' single-variable column.
#' @useDynLib phylosmith
#' @usage merge_treatments(phyloseq_obj, ...)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample.
#' @param ... any number of column names as \code{string}s or \code{numeric}s
#' in the \code{\link[phyloseq:sample_data]{sample_data}} that are to be
#' combined.
#' @export
#' @return phyloseq-object
#' @examples merge_treatments(soil_column, 'Matrix', 'Treatment', 'Day')

merge_treatments <- function(phyloseq_obj, ...) {
  if (!inherits(phyloseq_obj, "phyloseq")) {
    stop("`phyloseq_obj` must be a phyloseq-class
        object", call. = FALSE)
  }
  if (is.null(access(phyloseq_obj, 'sam_data'))) {
    stop("`phyloseq_obj` must contain sample_data()
        information",
         call. = FALSE)
  }
  treatments <- check_numeric_treatment(phyloseq_obj, ...)
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
#' @usage relative_abundance(phyloseq_obj)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @export
#' @return phyloseq-object
#' @examples relative_abundance(soil_column)

relative_abundance <- function(phyloseq_obj) {
  if (!inherits(phyloseq_obj, "phyloseq")) {
    stop("`phyloseq_obj` must be a phyloseq-class
            object", call. = FALSE)
  }
  abundance_table <- access(phyloseq_obj, 'otu_table')
  abundance_table <-
    apply(
      abundance_table,
      2,
      FUN = function(c) {
        c / sum(c)
      }
    )
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
    sort_on <- check_numeric_treatment(phyloseq_obj, sort_on)
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
        stop(
          "`sort_on` must be at least one column
        name, or index, from the sample_data(), or a vector of a set
        order of sample names or indices",
          call. = FALSE
        )
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
#' strings. If assinged "numeric" will set ascending numerical order.
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
  treatment <- check_numeric_treatment(phyloseq_obj, treatment)
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
    treatment <- check_numeric_treatment(phyloseq_obj, treatment)
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
    original_levels <-
      lapply(access(phyloseq_obj, 'sam_data'), levels)

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
    phyloseq_obj <- phyloseq(
      access(phyloseq_obj, 'otu_table'),
      access(phyloseq_obj, 'tax_table'),
      access(phyloseq_obj, 'sam_data')
    )

    treatment_name <- paste(treatment, collapse = sep)
    if (!(is.null(treatment))) {
      phyloseq_obj <- merge_treatments(phyloseq_obj, treatment)
      if (!(is.null(subset)) &
          any(!(subset %in% unlist(access(
            phyloseq_obj, 'sam_data'
          )[, c(treatment, treatment_name)])))) {
        stop("`subset` must be at least one factor from the `treatment`",
             call. = FALSE)
      }
      treatment_classes <- sort(unique(access(phyloseq_obj,
                                              'sam_data')[[treatment_name]]))
      if (below == TRUE) {
        phyloseq_obj <- do.call(merge_phyloseq,
                                apply(
                                  array(treatment_classes),
                                  1,
                                  FUN = function(group) {
                                    sub_phy <- eval(parse(
                                      text = paste0(
                                        "subset_samples(phyloseq_obj, ",
                                        treatment_name,
                                        " == '",
                                        group,
                                        "')"
                                      )
                                    ))
                                    cutoff <-
                                      floor(ncol(sub_phy@otu_table) * frequency)
                                    sub_phy <- filter_taxa(sub_phy, function(x) {
                                      sum(x != 0, na.rm = TRUE) <= cutoff
                                    }, TRUE)
                                    return(sub_phy)
                                  }
                                ))
      } else {
        phyloseq_obj <- do.call(merge_phyloseq,
                                apply(
                                  array(treatment_classes),
                                  1,
                                  FUN = function(group) {
                                    sub_phy <- eval(parse(
                                      text = paste0(
                                        "subset_samples(phyloseq_obj, ",
                                        treatment_name,
                                        " == '",
                                        group,
                                        "')"
                                      )
                                    ))
                                    cutoff <-
                                      floor(ncol(sub_phy@otu_table) * frequency)
                                    sub_phy <- filter_taxa(sub_phy, function(x) {
                                      sum(x != 0, na.rm = TRUE) >= cutoff
                                    }, TRUE)
                                    if (sum(taxa_sums(sub_phy), na.rm = TRUE) != 0) {
                                      sub_phy <- filter_taxa(sub_phy, function(x) {
                                        sum(x, na.rm = TRUE) != 0
                                      }, TRUE)
                                    }
                                    return(sub_phy)
                                  }
                                ))
      }
    } else {
      cutoff <- floor(ncol(access(phyloseq_obj, 'otu_table')) * frequency)
      if (below == TRUE) {
        phyloseq_obj <- filter_taxa(phyloseq_obj, function(x) {
          sum(x != 0, na.rm = TRUE) <= cutoff
        }, TRUE)
      } else {
        phyloseq_obj <- filter_taxa(phyloseq_obj, function(x) {
          sum(x != 0, na.rm = TRUE) >= cutoff
        }, TRUE)
        phyloseq_obj <- filter_taxa(phyloseq_obj, function(x) {
          sum(x, na.rm = TRUE) != 0
        }, TRUE)
      }
    }
    if (!(is.null(subset))) {
      phyloseq_obj <-
        prune_samples(sample_names(phyloseq_obj)[unname(apply(access(phyloseq_obj, 'sam_data')[, c(treatment, treatment_name)], 1,
                                                              function(x) {
                                                                any(x %in% subset)
                                                              }))], phyloseq_obj)
      phyloseq_obj <- filter_taxa(phyloseq_obj, function(x) {
        sum(x, na.rm = TRUE) > 0
      }, TRUE)
    }
    if (drop_samples == TRUE) {
      phyloseq_obj <- prune_samples(sample_sums(phyloseq_obj) > 0,
                                    phyloseq_obj)
    }
    for (i in seq_along(original_levels[!(names(original_levels) %in% treatment_name)])) {
      if (!(is.null(unname(unlist(
        original_levels[i]
      ))))) {
        phyloseq_obj <- set_treatment_levels(phyloseq_obj,
                                             names(original_levels)[i], unname(unlist(original_levels[i])))
      }
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

#' Compute proportions for taxa.
#'
#' Computes the proportion of a taxa classification. Function from the
#' phylosmith-package.
#' @useDynLib phylosmith
#' @usage taxa_proportions(phyloseq_obj, classification, treatment = NA)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. it
#' must contain \code{\link[phyloseq:tax_table]{tax_table()}}) with
#' information about each taxa/gene.
#' @param classification Column name as a \code{string} or \code{numeric} in
#' the \code{\link[phyloseq:tax_table]{tax_table}} for the prportions to be
#' reported on.
#' @param treatment Column name as a \code{string} or \code{numeric} in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column. If \code{NA},
#' then proprtions will be reported for the entire dataset. If set to "sample"
#' it will repoort proportions by sample.
#' @export
#' @return data.table
#' @examples taxa_proportions(soil_column, 'Phylum', treatment = NA)
#' taxa_proportions(soil_column, 'Phylum', treatment = 'sample')
#' taxa_proportions(soil_column, 'Phylum', treatment = c('Matrix', 'Treatment'))

taxa_proportions <-
  function(phyloseq_obj, classification, treatment = NA) {
    if (!inherits(phyloseq_obj, "phyloseq")) {
      stop("`phyloseq_obj` must be a phyloseq-class
        object", call. = FALSE)
    }
    if (is.null(access(phyloseq_obj, 'tax_table'))) {
      stop("`phyloseq_obj` must contain tax_table()
        information",
           call. = FALSE)
    }
    classification <- check_numeric_classification(phyloseq_obj,
                                                   classification)
    if (any(!(classification %in% colnames(access(
      phyloseq_obj, 'tax_table'
    )))))
    {
      stop("`classification` must be a column from the
        the tax_table()",
           call. = FALSE)
    }
    if (any(!(is.na(treatment))) &
        is.null(access(phyloseq_obj, 'sam_data'))) {
      stop(
        "`phyloseq_obj` must contain sample_data()
        information if `treatment` argument is used",
        call. = FALSE
      )
    }
    if (any(!(is.na(treatment))) & any(!(treatment %in% c('sample',
                                                          colnames(
                                                            access(phyloseq_obj, 'sam_data')
                                                          ))))) {
      stop(
        "`treatment` must be either NA, 'sample', or
        at least one column name, or index, from the sample_data()",
        call. = FALSE
      )
    }
    tax_table(phyloseq_obj) <- access(phyloseq_obj,
                                      'tax_table')[, classification]

    if (any(!(is.na(treatment))) & !('sample' %in% treatment)) {
      phyloseq_obj <- taxa_filter(phyloseq_obj, treatment)
      treatment_name <- paste(treatment, collapse = sep)
      phyloseq_obj <- phyloseq(
        access(phyloseq_obj, 'otu_table'),
        access(phyloseq_obj, 'tax_table'),
        access(phyloseq_obj, 'sam_data')[, treatment_name]
      )
    } else {
      phyloseq_obj <- phyloseq(access(phyloseq_obj, 'otu_table'),
                               access(phyloseq_obj, 'tax_table'))
    }
    class_table <- melt_phyloseq(conglomerate_taxa(phyloseq_obj,
                                                   classification))
    if (any(!(is.na(treatment))) & !('sample' %in% treatment)) {
      class_table <- class_table[,-2]
      class_table[, Abundance := sum(Abundance),
                  by = c(treatment_name, classification)]
      class_table <- unique(class_table)
      class_table[, Proportion := round(Abundance / sum(Abundance), 3),
                  by = treatment_name]
      class_table <- class_table[, c(3, 4, 5)]
      eval(parse(text = paste0(
        'setkey(class_table, ', treatment_name, ')'
      )))
    } else if ('sample' %in% treatment) {
      class_table[, Proportion := round(Abundance / sum(Abundance), 3),
                  by = Sample]
      class_table <- class_table[, c(2, 4, 5)]
      setkey(class_table, Sample)
    } else if (any(is.na(treatment))) {
      class_table <-
        class_table[, c(3, 4)][, lapply(.SD, sum, na.rm = TRUE),
                               by = classification]
      class_table <- class_table[,
                                 Proportion := round(Abundance / sum(Abundance), 3)][, -2]
    }
    return(class_table)
  }

#' Find unique taxa between treatments of a phyloseq object.
#' Function from the phylosmith-package.
#'
#' Inputs a \code{\link[phyloseq]{phyloseq-class}} object and
#' finds which taxa are taxa that are unique to a specific subset of the data.
#' @useDynLib phylosmith
#' @usage unique_taxa(phyloseq_obj, treatment, subset = NULL)
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
#' @seealso \code{\link{common_taxa}}
#' @export
#' @return list
#' @examples unique_taxa(soil_column, c('Matrix', 'Treatment'))

unique_taxa <- function(phyloseq_obj, treatment, subset = NULL) {
  if (!inherits(phyloseq_obj, "phyloseq")) {
    stop("`phyloseq_obj` must be a phyloseq-class object",
         call. = FALSE)
  }
  if (is.null(access(phyloseq_obj, 'sam_data'))) {
    stop("`phyloseq_obj` must contain sample_data()
        information",
         call. = FALSE)
  }
  treatment <- check_numeric_treatment(phyloseq_obj, treatment)
  if (any(!(treatment %in% colnames(access(
    phyloseq_obj, 'sam_data'
  ))))) {
    stop(
      "`treatment` must be at least one column name, or
        index, from the sample_data()",
      call. = FALSE
    )
  }
  phyloseq_obj <- phyloseq(
    access(phyloseq_obj, 'otu_table'),
    access(phyloseq_obj, 'tax_table'),
    access(phyloseq_obj, 'sam_data')
  )
  phyloseq_obj <-
    taxa_filter(phyloseq_obj, treatment, subset = subset)
  treatment_name <- paste(treatment, collapse = sep)
  treatment_classes <- unique(access(phyloseq_obj,
                                     'sam_data')[[treatment_name]])

  seen_taxa <- lapply(
    treatment_classes,
    FUN = function(trt) {
      taxa_names(taxa_filter(
        phyloseq_obj,
        treatment = treatment_name,
        subset = trt,
        frequency = 0
      ))
    }
  )
  names(seen_taxa) <- treatment_classes
  taxa_counts <- table(unlist(seen_taxa))
  unique_taxa <- names(taxa_counts[taxa_counts == 1])
  unique_taxa <- lapply(
    seen_taxa,
    FUN = function(toi) {
      toi[toi %in% unique_taxa]
    }
  )
  return(unique_taxa)
}
