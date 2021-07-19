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
    treatment <- check_index_treatment(phyloseq_obj, treatment)
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
#' the \code{\link[phyloseq:tax_table]{tax_table}} for the proportions to be
#' reported on.
#' @param treatment Column name as a \code{string} or \code{numeric} in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column. If \code{NA},
#' then proportions will be reported for the entire dataset. If set to "sample"
#' it will report proportions by sample.
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
    classification <- check_index_classification(phyloseq_obj,
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
    phyloseq_obj <- conglomerate_taxa(phyloseq_obj, classification, FALSE)
    class_table <- melt_phyloseq(phyloseq_obj)
    if (any(!(is.na(treatment))) & !('sample' %in% treatment)) {
      class_table <- class_table[, .(Abundance = sum(Abundance)),
                                 by = c(treatment_name, classification)]
      class_table[, Proportion := round(Abundance / sum(Abundance), 4),
                  by = treatment_name]
      setkeyv(class_table, treatment_name)
    } else if ('sample' %in% treatment) {
      class_table <- class_table[, .(classification = get(classification), Proportion = round(Abundance / sum(Abundance), 4)), by = Sample]
      setnames(class_table, 'classification', classification)
      setkey(class_table, Sample)
    } else if (any(is.na(treatment))) {
      class_table <- class_table[, .(Abundance = sum(Abundance)),by = classification]
      class_table[, Proportion := round(Abundance / sum(Abundance), 4)]
      setorder(class_table, -Proportion)
    }
    class_table[, Abundance := NULL][]
    return(class_table)
  }


#' Filter taxa in phyloseq-object to only include core taxa.
#' Function from the phylosmith-package.
#'
#' Inputs a phyloseq object and finds which taxa are seen in a
#' given proportion of samples at a minimum relative abundance, either in the
#' entire dataset, by treatment, or a particular treatment of interest.
#' @useDynLib phylosmith
#' @usage taxa_core(phyloseq_obj, treatment = NULL, subset = NULL,
#' frequency = 0.5, abundance_threshold = 0.01)
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
#' @param abundance_threshold The minimum relative abundance the taxa is found
#' in for each sample.
#' @export
#' @return phyloseq-object
#' @examples taxa_core(soil_column, frequency = 0.2, abundance_threshold = 0.01)

taxa_core <-
  function(phyloseq_obj,
           treatment = NULL,
           subset = NULL,
           frequency = 0.5,
           abundance_threshold = 0.01
  ) {
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
    if (!(is.numeric(abundance_threshold)) |
        !(abundance_threshold >= 0 & abundance_threshold <= 1)) {
      stop("`abundance_threshold` must be a numeric value between
        0 and 1", call. = FALSE)
    }
    phyloseq_obj <- taxa_filter(phyloseq_obj, treatment = treatment, subset = subset)
    treatment_name <- paste(treatment, collapse = sep)
    core_taxa <- relative_abundance(phyloseq_obj)
    core_taxa <- melt_phyloseq(phyloseq(core_taxa@otu_table, core_taxa@sam_data))

    if(is.null(treatment)){
      N <- nsamples(phyloseq_obj)
      core_taxa <- core_taxa[Abundance >= abundance_threshold]
      core_taxa <- core_taxa[, .(count = .N), by = OTU][count >= floor(N*frequency)][['OTU']]
    } else {
      taxa <- vector()
      treatments <- levels(core_taxa[[treatment_name]])
      for(treatment in treatments){
        N <- sum(phyloseq_obj@sam_data[[treatment_name]] == treatment)
        sub_table <- core_taxa[core_taxa[[treatment_name]] == treatment, ]
        sub_table <- sub_table[Abundance >= abundance_threshold]
        taxa <- c(taxa,
                  sub_table[, .(count = .N), by = OTU][count >= floor(N*frequency)][['OTU']])
      }
      core_taxa <- taxa
    }
    phyloseq_obj <- prune_taxa(core_taxa, phyloseq_obj)
    return(phyloseq_obj)
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
  treatment <- check_index_treatment(phyloseq_obj, treatment)
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
