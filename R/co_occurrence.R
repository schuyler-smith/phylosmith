#' pair-wise Spearman rank co-occurrence, written in efficient c++ code. Function from the phylosmith-package.
#'
#' A rewrite of the pair-wise Spearman rank co-occurrence routine written by \href{https://github.com/germs-lab/FastCoOccur}{Jin Choi}. The routine has been adapted to integrate with the \code{\link[Rcpp]{Rcpp-package}} API.
#' @useDynLib phylosmith
#' @usage co_occurrence(phyloseq_obj, treatment = NULL, p = 0.05, cores = 0)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object.
#' @param treatment Column name as a \code{string} or \code{numeric} in the \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of multiple columns and they will be combined into a new column.
#' @param p \code{numeric} The p-value cutoff. All returned co-occurrences will have a p-value less than or equal to \code{p}.
#' @param cores \code{numeric} Number of CPU cores to use for the pair-wise permutations. Default (0) uses max cores available. Parallelization not available for systems running MacOS without openMP configuration.
#' @aliases FastCoOccur
#' @import data.table
#' @import RcppArmadillo
#' @import RcppParallel
#' @import RcppProgress
#' @keywords nonparametric
#' @seealso \code{\link{bootstrap_rho}} \code{\link{phylosmith}}
#' @export

# sourceCpp("src/co_occurrence_Rcpp.cpp")

co_occurrence <- function(phyloseq_obj, treatment = NULL, p = 0.05, cores = 0){
  # phyloseq_obj = mock_phyloseq; treatment = c("treatment", "day"); p = 0.05
  options(warnings=-1)

  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment = treatment)
  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  treatment_name <- paste(treatment, collapse = sep)

  treatment_classes <- as.character(unique(phyloseq_obj@sam_data[[treatment_name]]))
  treatment_indices <- lapply(treatment_classes, FUN = function(trt){which(as.character(phyloseq_obj@sam_data[[treatment_name]]) %in% trt)-1})
  if(is.null(treatment)){treatment_classes <- 'Experiment_Wide'
  treatment_indices <- list(1:nsamples(phyloseq_obj)-1)}

  if(cores == 0){cores <- parallel::detectCores()}
  if(cores == 0){cores <- parallel::detectCores()}
  co_occurrence <- co_occurrence_Rcpp(phyloseq_obj@otu_table, treatment_indices, treatment_classes, p, cores)
  return(as.data.table(co_occurrence))
}

#' Bootstraps the pair-wise Spearman rank co-occurrence, to determine a significant rho-cutoff. Function from the phylosmith-package.
#'
#' Bootstraps the pair-wise Spearman rank co-occurrence, to determine a significant rho-cutoff.
#' @useDynLib phylosmith
#' @usage bootstrap_rho(phyloseq_obj, treatment = NULL,
#' replicate_samples = 'independent', permutations = 100)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object.
#' @param treatment Column name as a \code{string} or \code{numeric} in the \code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of multiple columns and they will be combined into a new column.
#' @param replicate_samples Column name as a \code{string} or \code{numeric} in the \code{\link[phyloseq:sample_data]{sample_data}} that indicates which samples are non-independent of each other.
#' @param permutations \code{numeric} Number of iterations to compute.
#' @keywords nonparametric
#' @import data.table
#' @import RcppArmadillo
#' @import RcppParallel
#' @import RcppProgress
#' @seealso \code{\link{co_occurrence}}
#' @export

# sourceCpp('src/co_occurrence_Rcpp.cpp')

bootstrap_rho <- function(phyloseq_obj, treatment = NULL, replicate_samples = 'independent', permutations = 100){
  # phyloseq_obj = mock_phyloseq; treatment = c("treatment", "day"); replicate_samples = 'independent'; permutations = 10; p = 0; cores = 0;
  options(warnings=-1)

  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  if(is.numeric(replicate_samples)){replicate_samples <- colnames(phyloseq_obj@sam_data[,replicate_samples])}

  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment = treatment, frequency = 0)
  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  treatment_name <- paste(treatment, collapse = sep)
  treatment_classes <- as.character(unique(phyloseq_obj@sam_data[[treatment_name]]))
  treatment_indices <- lapply(treatment_classes, FUN = function(trt){which(as.character(phyloseq_obj@sam_data[[treatment_name]]) %in% trt)-1})

  if(replicate_samples == 'independent' & is.null(treatment)){
    replicate_indices <- 1:ncol(phyloseq_obj@otu_table)
  } else if(replicate_samples == 'independent' & !(is.null(treatment))){
    phyloseq_obj_reps <- merge_treatments(phyloseq_obj, c(treatment))
    replicate_name <- paste(c(treatment), collapse = sep)
    replicate_samples <- as.character(unique(phyloseq_obj_reps@sam_data[[replicate_name]]))
    replicate_indices <- lapply(replicate_samples, FUN = function(trt){which(as.character(phyloseq_obj_reps@sam_data[[replicate_name]]) %in% trt)})
  } else if(replicate_samples != 'independent' & !(is.null(treatment))){
    phyloseq_obj_reps <- merge_treatments(phyloseq_obj, c(treatment, replicate_samples))
    replicate_name <- paste(c(treatment, replicate_samples), collapse = sep)
    replicate_samples <- as.character(unique(phyloseq_obj_reps@sam_data[[replicate_name]]))
    replicate_indices <- lapply(replicate_samples, FUN = function(trt){which(as.character(phyloseq_obj_reps@sam_data[[replicate_name]]) %in% trt)})
  }

  rhos <- data.table(Treatment = factor(levels = treatment_classes), rho = numeric(), Count = numeric())
  n <- nrow(phyloseq_obj@otu_table)
  permuted_phyloseq_obj <- phyloseq_obj

  tryCatch({
    for(i in 1:permutations){
      for(indices in replicate_indices){
        permuted_phyloseq_obj@otu_table[,indices] <- phyloseq_obj@otu_table[sample(1:n, n),indices]
      }
      co_occurence_table <- data.table(co_occurrence_rho_Rcpp(permuted_phyloseq_obj@otu_table, treatment_indices, treatment_classes))
      co_occurence_table[, rho  := round(.SD, 3), .SDcols = 'rho']
      co_occurence_table[, Count := .N, by = .(Treatment, rho)]
      co_occurence_table <- unique(co_occurence_table)
      rhos <- rbindlist(list(rhos, co_occurence_table))[, lapply(.SD, sum, na.rm = TRUE), by = .(Treatment, rho)]
    }},
    interrupt = function(interrupt){rhos <- rhos[-length(rhos)]; message('Interrupted after ', i, ' permutations.'); return(rhos)})

  return(setkey(rhos, Treatment, rho))
} #else {
#   return(stats::quantile(rhos, 1-p, na.rm = TRUE))}
# }

#' Curate co-occurrence data. Function from the phylosmith-package.
#'
#' Used to curate a co-occurrence table from the \code{\link{co_occurrence}} function. Takes a list of taxa and finds all pairs containing any of those taxa.
#' @useDynLib phylosmith
#' @usage curate_co_occurrence(co_occurrence_table, taxa_of_interest, number_of_treatments = 1)
#' @param co_occurrence_table Table of the co-occurrence of taxa/genes in the \code{phyloseq_obj}, computed using \code{\link{co_occurrence}}.
#' @param taxa_of_interest A list or vector of taxa names, like those generated with \code{\link{find_unique_taxa}}.
#' @param number_of_treatments How many treatments should the taxa of interest be seen in? Requires \code{integer} or 'all' (default = 1).
#' @keywords manip
#' @import data.table
#' @import RcppArmadillo
#' @import RcppParallel
#' @seealso \code{\link{co_occurrence}}

curate_co_occurrence <- function(co_occurrence_table, taxa_of_interest, number_of_treatments = 1){
  sub_co_occurrence <- co_occurrence_table[(co_occurrence_table[[2]] %in% taxa_of_interest | co_occurrence_table[[3]] %in% taxa_of_interest),]
  toi_table <- unique(cbind(rbindlist(list(sub_co_occurrence[,1], sub_co_occurrence[,1])), rbindlist(list(sub_co_occurrence[,2], sub_co_occurrence[,3]))))
  toi_table <- toi_table[toi_table[[2]] %in% taxa_of_interest]
  if(number_of_treatments == 'all'){number_of_treatments <- length(unique(sub_co_occurrence[[1]]))
  } else {number_of_treatments <- number_of_treatments}
  toi <- names(table(toi_table[[2]])[table(toi_table[[2]]) >= number_of_treatments])
  sub_co_occurrence <- co_occurrence_table[(co_occurrence_table[[2]] %in% toi | co_occurrence_table[[3]] %in% toi),]

  # sourceCpp("src/arrange_co_occurrence_table_tbb.cpp")
  arranged_co_ocurrence <- as.data.table(arrange_co_occurrence_table(sub_co_occurrence, toi))

  setorder(arranged_co_ocurrence, 'Treatment')
  return(arranged_co_ocurrence)
}

#' Calculate quantiles for the bootstrapped rho values from the Spearman-rank co-occurrence. Function from the phylosmith-package.
#'
#' Calculate quantiles for the bootstrapped rho values from the Spearman-rank co-occurrence.
#' @useDynLib phylosmith
#' @usage quantile_bootstrapped_rhos(bootstrapped_rhos, p = 0.05, by_treatment = TRUE)
#' @param bootstrapped_rhos A \code{data.table} output from \code{\link[=bootstrap_rho]{bootstrap_rho}}.
#' @param p The significance threshold for setting cutoffs.
#' @param by_treatment Whether to find the rho cutoffs for each treatment individually or for the entire experiment. Suggested to do by treatment first, to see if there is any treatments that are outliers.
#' @import data.table
#' @seealso \code{\link[=bootstrap_rho]{bootstrap_rho}}
#' @export
#'

quantile_bootstrapped_rhos <- function(bootstrapped_rhos, p = 0.05, by_treatment = TRUE){
  if(by_treatment){
    bootstrapped_rhos[, Proportion := Count/sum(Count), by = Treatment]
    quantiles <- bootstrapped_rhos[, list(lower = rho[sum(cumsum(Proportion) <= (p/2))], upper = rho[sum(cumsum(Proportion) <= (1-(p/2)))]), by = Treatment]
  } else {
    bootstrapped_rhos <- bootstrapped_rhos[,-1][, lapply(.SD, sum, na.rm = TRUE), by = rho]
    bootstrapped_rhos[, Proportion := Count/sum(Count)]
    quantiles <- bootstrapped_rhos[, list(lower = rho[sum(cumsum(Proportion) <= (p/2))], upper = rho[sum(cumsum(Proportion) <= (1-(p/2)))])]
  }
  return(quantiles)
}

#' Create a ggplot object of the distribution of rho values from bootstrap_rho(). Function from the phylosmith-package.
#'
#' Plots the output of \code{\link[=bootstrap_rho]{bootstrap_rho}} into a histogram with the distributions shown by treatment. This is a visualization tool to help show how the bootstrapping worked, and to see where the cutoffs lie.
#' @useDynLib phylosmith
#' @usage histogram_bootstrapped_rhos(bootstrapped_rhos, p = 0.05,
#' zeros = FALSE, x_breaks = 0.25, colors = 'default')
#' @param bootstrapped_rhos A \code{data.table} output from \code{\link[=bootstrap_rho]{bootstrap_rho}}.
#' @param p The significance threshold for setting cutoffs.
#' @param zeros In some cases either where a treatment is underpowered (few samples) there may be an extreme number of 0 values for rho, which can affect the scale an make the distribution curves look poor for treatments that are okay. The 0-rhos are still included in the calculations, but not displayed on the graph (\code{FALSE}).
#' @param x_breaks What intervals to set the ticks on the x-axis.
#' @param colors Name of a color set from the \link[=RColorBrewer]{RColorBrewer} package or a vector palete of R-accepted colors.
#' @import ggplot2
#' @seealso \code{\link[=bootstrap_rho]{bootstrap_rho}}
#' @export
#'

histogram_bootstrapped_rhos <- function(bootstrapped_rhos, p = 0.05, zeros = FALSE, x_breaks = 0.25, colors = 'default'){
  color_count <- length(unique(bootstrapped_rhos[,Treatment]))
  graph_colors <- create_palette(color_count, colors)

  bootstrapped_rhos[, Proportion := Count/sum(Count), by = Treatment]
  quantiles <- bootstrapped_rhos[, list(lower = rho[sum(cumsum(Proportion) <= (p/2))], upper = rho[sum(cumsum(Proportion) <= (1-(p/2)))]), by = Treatment]

  bootstrapped_rhos[, bin := findInterval(rho, seq(-1, 1, .03)), by = Treatment]
  bootstrapped_rhos <- bootstrapped_rhos[, list(rho = mean(rho), Count = sum(Count), Proportion = sum(Proportion)), by = .(Treatment, bin)]

  if(!(zeros)){bootstrapped_rhos <- bootstrapped_rhos[rho != 0]}
  g <- ggplot(bootstrapped_rhos, aes(x = rho, y = Proportion, fill = Treatment)) +
    scale_x_continuous(breaks = seq(-1, 1, x_breaks)) +
    scale_fill_manual(values = graph_colors) +
    theme_light() +
    geom_vline(data = quantiles, xintercept = c(quantiles$lower, quantiles$upper), color = c(graph_colors, graph_colors), size = 1.5, alpha = 0.8) +
    geom_bar(stat = 'identity', position = position_dodge(width = .02), width = .05)
  for(i in 1:nrow(quantiles)){
    g <- g + geom_text(data = quantiles, x = quantiles$lower[i]-.03, label = paste0(round(quantiles$lower[i],2)), y = (.075 + (i*.015)), color = graph_colors[i], size = 5) +
      geom_text(data = quantiles, x = quantiles$upper[i]+.03, label = paste0(round(quantiles$upper[i],2)), y = (.075 + (i*.015)), color = graph_colors[i], size = 5)
  }
  return(g)
}

#  ###RAM check .. because apparently some people have 36,000 genes in their studies...
#  n <- nrow(phyloseq_obj@otu_table)
#  memory_function <- function(n, t, r){(6.780e-07 + (2.867e-08 * (n*(n-1)*t))) - r}
#  required_memory <- round((6.780e-07 + (2.867e-08 * (n*(n-1)*length(treatments))))*p, 2)
#  if(Sys.info()['sysname'] == "Windows"){
#    available_memory <- round(as.numeric(gsub("\r","",gsub("FreePhysicalMemory=","",system('wmic OS get FreePhysicalMemory /Value',intern=TRUE)[3])))/(1024^2), 2)
#    } else {
#    available_memory <- round(as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern = TRUE))/(1024^2), 2)}
#  if(required_memory >= available_memory){
#    message("Calculating ", (n*(n-1)*length(treatments)), " pairs-wise co-occurrences.
# This could take ", required_memory, " GB of RAM, or more depending on the number of significant interactions.
# This machine has ", available_memory, " GB of RAM available.
# Recommend subsetting from ", n, " taxa to ",
#   floor(stats::uniroot(memory_function, t=length(treatments), r=available_memory, lower=0, upper=1000000)$root), " if your systems gets stuck."
#   )}

# start_time <- Sys.time()
# for(indices in replicate_indices){
#   permuted_phyloseq_obj@otu_table[,indices] <- phyloseq_obj@otu_table[sample(1:n, n),indices]
# }
# rhos[1] <- mean(FastCoOccur(permuted_phyloseq_obj, treatment, cores = cores)$rho)
# end_time <- Sys.time()
# time_taken <- end_time - start_time
# message('Allow up to ', round((time_taken*permutations)/360, 1), ' hours to complete these calculations.')

