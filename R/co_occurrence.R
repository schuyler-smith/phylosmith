#' pair-wise Spearman rank co-occurrence, written in efficient c++ code.
#' Function from the phylosmith-package.
#'
#' A rewrite of the pair-wise Spearman rank co-occurrence routine written by
#' \href{https://github.com/germs-lab/FastCoOccur}{Jin Choi}. The routine has
#' been adapted to integrate with the \code{\link[Rcpp]{Rcpp-package}} API.
#' @useDynLib phylosmith
#' @usage co_occurrence(phyloseq_obj, treatment = NULL, subset = NULL,
#' rho = 0, p = 0.05, method = 'spearman', cores = 1)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object.
#' @param treatment Column name as a \code{string} or \code{numeric} in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any
#' samples that to not contain this factor. This can be a vector of multiple
#' factors to subset on.
#' @param rho \code{numeric} The rho-value cutoff. All returned co-occurrences
#' will have a rho-value less than or equal to \code{rho} or less than or equal to -\code{rho}.
#' @param p \code{numeric} The p-value cutoff. All returned co-occurrences
#' will have a p-value less than or equal to \code{p}.
#' @param method Which correlation method to calculate, "pearson", "spearman".
#' @param cores \code{numeric} Number of CPU cores to use for the pair-wise
#' permutations. Default (0) uses max cores available. Parallelization not
#' available for systems running MacOS without openMP configuration.
#' @importFrom parallel detectCores
#' @keywords nonparametric
#' @seealso \code{\link{permute_rho}} \code{\link{phylosmith}}
#' @export
#' @return data.table
#' @examples
#' co_occurrence(soil_column, treatment = c('Matrix', 'Treatment'),
#' subset = 'Amended', rho = 0.8, p = 0.05, cores = 1)

# sourceCpp("src/correlations_Rcpp.cpp")

co_occurrence <-
  function(phyloseq_obj,
           treatment = NULL,
           subset = NULL,
           rho = 0,
           p = 0.05,
           method = 'spearman',
           cores = 1) {
    if (!inherits(phyloseq_obj, "phyloseq")) {
      stop("`phyloseq_obj` must be a phyloseq-class object",
           call. = FALSE)
    }
    treatment <- check_index_treatment(phyloseq_obj, treatment)
    if (!(is.null(treatment)) &
        is.null(access(phyloseq_obj, 'sam_data'))) {
      stop(
        "`phyloseq_obj` must contain sample_data()
        information if `treatment` argument is used",
        call. = FALSE
      )
    }
    if (any(!(treatment %in% colnames(access(
      phyloseq_obj, 'sam_data'
    ))))) {
      stop(
        "`treatment` must be at least one column name,
        or index, from the sample_data()",
        call. = FALSE
      )
    }
    if (!(is.numeric(p)) | !(p >= 0 & p <= 1)) {
      stop("`p` must be a numeric value between 0 and 1",
           call. = FALSE)
    }
    if (!(is.numeric(rho)) | !(rho >= 0 & rho <= 1)) {
      stop("`rho` must be a numeric value between 0 and 1",
           call. = FALSE)
    }
    if (!(is.numeric(cores)) |
        !(cores >= 0 & cores <= (detectCores() - 1))) {
      stop(
        "`cores` must be a numeric value between 0 and ",
        (detectCores() - 1),
        "\n(upper limit is set by the cores available on machine used to
        execute code)" ,
        call. = FALSE
      )
    }
    if (cores == 0) {
      cores <- (detectCores() - 1)
    }
    match.arg(method, c("pearson", "kendall", "spearman"))
    options(warnings = -1)

    phyloseq_obj <-
      taxa_filter(phyloseq_obj, treatment = treatment, subset = subset)
    treatment_name <- paste(treatment, collapse = sep)
 #   phyloseq_obj <- relative_abundance(phyloseq_obj)

    treatment_classes <- as.character(unique(access(phyloseq_obj,
                                                    'sam_data')[[treatment_name]]))
    treatment_indices <- lapply(
      treatment_classes,
      FUN = function(trt) {
        which(as.character(access(phyloseq_obj,
                                  'sam_data')[[treatment_name]]) %in% trt)
      }
    )
    if (is.null(treatment)) {
      treatment_classes <- 'Experiment_Wide'
      treatment_indices <- list(seq(nsamples(phyloseq_obj)))
    }
    phyloseq_obj <- access(phyloseq_obj, 'otu_table')
    co_occurrence <- data.table()
    for(i in seq_along(treatment_indices)){
      treatment_co_occurrence <- Correlation(
        X = phyloseq_obj[,treatment_indices[[i]]],
        cor_coef_cutoff = rho,
        p_cutoff = p,
        method = method,
        ncores = cores
      )
      treatment_co_occurrence[['X']] <- rownames(phyloseq_obj[,treatment_indices[[i]]])[treatment_co_occurrence[['X']]]
      treatment_co_occurrence[['Y']] <- rownames(phyloseq_obj[,treatment_indices[[i]]])[treatment_co_occurrence[['Y']]]
      if(length(treatment_indices) > 0){
        treatment_co_occurrence <- cbind(Treatment = treatment_classes[i], treatment_co_occurrence)
      }
      co_occurrence <- rbind(co_occurrence, treatment_co_occurrence)
    }
    return(as.data.table(co_occurrence))
  }

#' Permutes the pair-wise Spearman rank co-occurrence, to determine a
#' significant rho-cutoff. Function from the phylosmith-package.
#'
#' Permutes the pair-wise Spearman rank co-occurrence, to determine a
#' significant rho-cutoff.
#' @useDynLib phylosmith
#' @usage permute_rho(phyloseq_obj, treatment = NULL, subset = NULL,
#' replicate_samples = NULL, permutations = 10, method = 'spearman',
#' cores = 1)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object.
#' @param treatment Column name as a \code{string} or \code{numeric} in the
#' \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
#' multiple columns and they will be combined into a new column.
#' @param subset A factor within the \code{treatment}. This will remove any
#' samples that to not contain this factor. This can be a vector of multiple
#' factors to subset on.
#' @param replicate_samples Column name as a \code{string} or \code{numeric}
#' in the \code{\link[phyloseq:sample_data]{sample_data}} that indicates which
#' samples are non-independent of each other.
#' @param permutations \code{numeric} Number of iterations to compute.
#' @param method Which correlation method to calculate, "pearson", "spearman".
#' @param cores \code{numeric} Number of CPU cores to use for the pair-wise
#' permutations. Default (0) uses max cores available. Parallelization not
#' available for systems running MacOS without openMP configuration.
#' @keywords nonparametric
#' @importFrom parallel detectCores
#' @seealso \code{\link{co_occurrence}}
#' @export
#' @return table
#' @examples
#' permute_rho(soil_column, treatment = c('Matrix', 'Treatment'),
#' subset = 'Amended', replicate_samples = 'Day', permutations = 1,  cores = 0)

# sourceCpp('src/correlations_Rcpp.cpp')

permute_rho <-
  function(phyloseq_obj,
           treatment = NULL,
           subset = NULL,
           replicate_samples = NULL,
           permutations = 10,
           method = 'spearman',
           cores = 1) {
    if (!inherits(phyloseq_obj, "phyloseq")) {
      stop("`phyloseq_obj` must be a phyloseq-class object",
           call. = FALSE)
    }
    treatment <- check_index_treatment(phyloseq_obj, treatment)
    if (!(is.null(treatment)) &
        is.null(access(phyloseq_obj, 'sam_data'))) {
      stop(
        "`phyloseq_obj` must contain sample_data()
        information if `treatment` argument is used",
        call. = FALSE
      )
    }
    if (!(is.null(replicate_samples))){
      replicate_samples <- check_index_treatment(phyloseq_obj,
                                                   replicate_samples)
    }
    if (!(is.null(replicate_samples)) &
        is.null(access(phyloseq_obj, 'sam_data'))) {
      stop(
        "`phyloseq_obj` must contain sample_data()
        information if `replicate_samples` argument is used",
        call. = FALSE
      )
    }
    if (!(is.null(replicate_samples)) &
        any(!(replicate_samples %in% colnames(access(
          phyloseq_obj,
          'sam_data'
        ))))) {
      stop(
        "`replicate_samples` must be at least one column
        name, or index, from the sample_data()",
        call. = FALSE
      )
    }
    if (!(is.numeric(permutations)) | !(permutations >= 0)) {
      stop("`permutations` must be a numeric value greater
        than 0", call. = FALSE)
    }
    if (!(is.numeric(cores)) |
        !(cores >= 0 & cores <= (detectCores() - 1))) {
      stop(
        "`cores` must be a numeric value between 0 and ",
        (detectCores() - 1),
        "\n(upper limit is set by the cores available on machine used to
        xecute code)",
        call. = FALSE
      )
    }
    options(warnings = -1)
    phyloseq_obj <- taxa_filter(
      phyloseq_obj,
      treatment = treatment,
      subset = subset,
      frequency = 0
    )
    phyloseq_obj <- relative_abundance(phyloseq_obj)
    treatment_name <- paste(treatment, collapse = sep)
    treatment_classes <- as.character(unique(access(phyloseq_obj,
                                                    'sam_data')[[treatment_name]]))
    treatment_indices <- lapply(
      treatment_classes,
      FUN = function(trt) {
        which(as.character(access(phyloseq_obj,
                                  'sam_data')[[treatment_name]]) %in% trt) - 1
      }
    )
    if (is.null(treatment)) {
      treatment_classes <- 'Experiment_Wide'
      treatment_indices <- list(seq(nsamples(phyloseq_obj)) - 1)
    }
    replicate_sample_classes <- vector()
    if (is.null(replicate_samples)){
      replicate_indices <- seq(ncol(access(phyloseq_obj, 'otu_table')))
    } else if (!(is.null(replicate_samples))){
      phyloseq_obj_reps <- merge_treatments(phyloseq_obj,
                                            c(treatment, replicate_samples))
      replicate_name <- paste(c(treatment, replicate_samples),
                              collapse = sep)
      replicate_sample_classes <-
        as.character(unique(access(phyloseq_obj_reps,
                                   'sam_data')[[replicate_name]]))
      replicate_indices <- lapply(
        replicate_sample_classes,
        FUN = function(trt) {
          which(as.character(access(phyloseq_obj_reps,
                                    'sam_data')[[replicate_name]]) %in% trt)
        }
      )
    }
    rhos <- data.table()
 #   phyloseq_obj <- relative_abundance(phyloseq_obj)
    permuted_counts <- as(phyloseq_obj@otu_table, 'matrix')
    tryCatch({
      for (j in seq(permutations)){
        co_occurrence_table <- data.table()
        if(is.null(replicate_samples)){
          n <- nrow(permuted_counts)
          permuted_counts<- apply(permuted_counts, 2, FUN = function(x){
            sample(x, n)
          })
        } else {for(indices in replicate_indices){
          n <- nrow(permuted_counts)
          permuted_counts[, indices] <- permuted_counts[sample(seq(n), n), indices]
        }}
        for(i in seq_along(treatment_classes)){
          rep_co_occurrence_table <- data.table(permute_rho_Rcpp(
            phyloseq_obj@otu_table[, treatment_indices[[i]]],
            permuted_counts[, treatment_indices[[i]]],
            method,
            cores))
          rep_co_occurrence_table[, rho := round(rho, 2)]
          rep_co_occurrence_table[, Count := .N, by = .(rho)]
          rep_co_occurrence_table <- unique(rep_co_occurrence_table)
          if(length(treatment_classes) > 1){
            rep_co_occurrence_table <- cbind(Treatment = treatment_classes[i], rep_co_occurrence_table)
          }
          co_occurrence_table <- rbind(co_occurrence_table, unique(rep_co_occurrence_table))
        }

        if(length(treatment_classes) > 1){
          rhos <- rbindlist(list(rhos, co_occurrence_table))[,
                                                             lapply(.SD, sum, na.rm = TRUE), by = .(Treatment, rho)]
        } else {
          rhos <- rbindlist(list(rhos, co_occurrence_table))[,
                                                             lapply(.SD, sum, na.rm = TRUE), by = .(rho)]
        }
      }
    },
    interrupt = function(interrupt) {
      message('Interrupted after ', i - 1, ' permutations completed.')
      setkey(rhos, Treatment, rho)
      return(rhos)
    })
    if(length(treatment_classes) > 1){
      setkey(rhos, Treatment, rho)
    } else {
      setkey(rhos, rho)
    }
    return(rhos)
  } #else {
#     return(stats::quantile(rhos, 1-p, na.rm = TRUE))}
# }

#' Calculate quantiles for the permuted rho values from the Spearman-rank
#' co-occurrence. Function from the phylosmith-package.
#'
#' Calculate quantiles for the permuted rho values from the Spearman-rank
#' co-occurrence.
#' @useDynLib phylosmith
#' @usage quantile_permuted_rhos(permuted_rhos, p = 0.05, by_treatment = TRUE)
#' @param permuted_rhos A \code{data.table} output from
#' \code{\link[=permute_rho]{permute_rho}}.
#' @param p The significance threshold for setting cutoffs.
#' @param by_treatment Whether to find the rho cutoffs for each treatment
#' individually or for the entire experiment. Suggested to do by treatment
#' first, to see if there is any treatments that are outliers.
#' @seealso \code{\link[=permute_rho]{permute_rho}}
#' @export
#' @return data.frame
#' @examples
#' permuted_rhos <- permute_rho(soil_column,
#' treatment = c('Matrix', 'Treatment'), replicate_samples = 'Day',
#' permutations = 1,  cores = 0)
#' quantile_permuted_rhos(permuted_rhos)
#' quantile_permuted_rhos(permuted_rhos, by_treatment = FALSE)

quantile_permuted_rhos <- function(permuted_rhos,
                                   p = 0.05,
                                   by_treatment = TRUE) {
  if (!(is.data.frame(permuted_rhos))) {
    stop("`permuted_rhos` must be at data.frame
        object", call. = FALSE)
  }
  if (!(is.numeric(p)) | !(p >= 0 & p <= 1)) {
    stop("`p` must be a numeric value between 0
        and 1", call. = FALSE)
  }
  if (!(is.logical(by_treatment))) {
    stop("`by_treatment` must must be either
        `TRUE`, or `FALSE`",
         call. = FALSE)
  }
  if (by_treatment) {
    permuted_rhos[, Proportion := Count / sum(Count), by = 'Treatment']
    quantiles <-
      permuted_rhos[, list(lower = rho[sum(cumsum(Proportion) <=
                                             (p / 2))], upper = rho[sum(cumsum(Proportion) <= (1 - (p / 2)))]),
                    by = 'Treatment']
  } else {
    permuted_rhos[, Proportion := Count / sum(Count)]
    quantiles <-
      permuted_rhos[, list(lower = rho[sum(cumsum(Proportion) <=
                                             (p / 2))], upper = rho[sum(cumsum(Proportion) <= (1 - (p / 2)))])]
  }
  return(quantiles)
}

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

histogram_permuted_rhos <- function(permuted_rhos,
                                    p = NULL,
                                    x_breaks = 0.25,
                                    colors = 'default') {
  if (!(is.data.frame(permuted_rhos))) {
    stop("`permuted_rhos` must be at data.frame
        object", call. = FALSE)
  }
  if (!(is.null(p))) {
    if (!(is.numeric(p)) | !(p >= -0 & p <= 1)) {
      stop("`p` must be a numeric value
            between 0 and 1", call. = FALSE)
    }
  }
  if (!(is.numeric(x_breaks)) |
      !(x_breaks >= -1 & x_breaks <= 1)) {
    stop("`x_breaks` must be a numeric value
        between -1 and 1", call. = FALSE)
  }
  color_count <- length(unique(permuted_rhos[['Treatment']]))
  graph_colors <- create_palette(color_count, colors)

  if(color_count > 0){permuted_rhos[, Proportion := Count / sum(Count),
                                    by = Treatment]
    quantiles <- permuted_rhos[, list(lower = rho[sum(cumsum(Proportion) <= (p /
                                                                               2))],
                                      upper = rho[sum(cumsum(Proportion) <= (1 - (p / 2)))]),
                               by = Treatment]

    permuted_rhos[, bin := findInterval(rho, seq(-1, 1, (max(permuted_rhos$rho) - min(permuted_rhos$rho))/100)), by = Treatment]
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
    geom_bar(stat = 'identity',
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
  g <- g + theme_light() +
    theme(
      axis.line.x = element_line(
        colour = 'black',
        size = 1,
        linetype = 'solid'
      ),
      axis.line.y = element_line(
        colour = 'black',
        size = 1,
        linetype = 'solid'
      ),
      axis.text.x = element_text(
        size = 8,
        vjust = 0.7,
        hjust = 0,
        angle = -30
      ),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10, face = "bold"),
      legend.background = element_rect(fill = (alpha = 0)),
      legend.spacing.x = unit(0.01, 'npc')
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.0025, 0.002)))
  return(g)
}

#' Curate co-occurrence data. Function from the phylosmith-package.
#'
#' Used to curate a co-occurrence table from the \code{\link{co_occurrence}}
#' function. Takes a list of taxa and finds all pairs containing any of those
#' taxa.
#' @useDynLib phylosmith
#' @usage curate_co_occurrence(co_occurrence_table, taxa_of_interest,
#' number_of_treatments = 1)
#' @param co_occurrence_table Table of the co-occurrence of taxa/genes in the
#' \code{phyloseq_obj}, computed using \code{\link{co_occurrence}}.
#' @param taxa_of_interest A list or vector of taxa names, like those
#' generated with \code{\link{unique_taxa}}.
#' @param number_of_treatments How many treatments should the taxa of interest
#' be seen in? Requires \code{integer} or 'all' (default = 1).
#' @keywords manip
#' @seealso \code{\link{co_occurrence}}
#' @return data.table

curate_co_occurrence <-
  function(co_occurrence_table,
           taxa_of_interest,
           number_of_treatments = 1) {
    sub_co_occurrence <-
      co_occurrence_table[(
        co_occurrence_table[[2]] %in%
          taxa_of_interest |
          co_occurrence_table[[3]] %in% taxa_of_interest
      ), ]
    toi_table <- unique(cbind(rbindlist(
      list(sub_co_occurrence[, 1],
           sub_co_occurrence[, 1])
    ), rbindlist(
      list(sub_co_occurrence[, 2],
           sub_co_occurrence[, 3])
    )))
    toi_table <- toi_table[toi_table[[2]] %in% taxa_of_interest]
    if (number_of_treatments == 'all') {
      number_of_treatments <- length(unique(sub_co_occurrence[[1]]))
    } else {
      number_of_treatments <- number_of_treatments
    }
    toi <- names(table(toi_table[[2]])[table(toi_table[[2]]) >=
                                         number_of_treatments])
    sub_co_occurrence <-
      co_occurrence_table[(co_occurrence_table[[2]] %in%
                             toi | co_occurrence_table[[3]] %in% toi), ]

    # sourceCpp("src/arrange_co_occurrence_table_tbb.cpp")
    arranged_co_ocurrence <-
      as.data.table(arrange_co_occurrence_table(sub_co_occurrence, toi))

    setorder(arranged_co_ocurrence, 'Treatment')
    return(arranged_co_ocurrence)
  }
