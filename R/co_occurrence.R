#' pair-wise Spearman rank co-occurrence, written in efficient c++ code. Function from the phylosmith-package.
#'
#' A rewrite of the pair-wise Spearman rank co-occurrence routine written by \href{https://github.com/germs-lab/FastCoOccur}{Jin Choi}. The routine has been adapted to integrate with the \code{\link[Rcpp]{Rcpp-package}} API.
#' @useDynLib phylosmith
#' @usage co_occurrence(phyloseq_obj, treatment, p = 0.05, cores = 0)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package.
#' @param treatment Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param p the p-value cutoff. all returned co-occurrences must have a p-value less than or equal to p.
#' @param cores Number of CPU cores to use for the pair-wise permutations. Default uses all cores available.
#' @aliases FastCoOccur
#' @import data.table
#' @import RcppArmadillo
#' @import RcppParallel
#' @import RcppProgress
#' @keywords nonparametric
#' @seealso \code{\link{bootstrap_rho}} \code{\link{phylosmith}}
#' @examples
#' data(mock_phyloseq)
#' co_occurrence(mock_phyloseq, "day", 0.05)
#' @export

# sourceCpp("src/co_occurrence_Rcpp.cpp")

co_occurrence <- function(phyloseq_obj, treatment, p = 0.05, cores = 0){
  # phyloseq_obj = mock_phyloseq; treatment = c("treatment", "day"); p = 0.05
  options(warnings=-1)

  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment = treatment)
  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  treatment_name <- paste(treatment, collapse = ".")

  treatments <- as.character(unique(phyloseq_obj@sam_data[[treatment_name]]))
  treatment_indices <- lapply(treatments, FUN = function(trt){which(as.character(phyloseq_obj@sam_data[[treatment_name]]) %in% trt)-1})

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

  if(cores == 0){cores <- parallel::detectCores()}
  cooccurrence <- co_occurrence_Rcpp(phyloseq_obj@otu_table, treatment_indices, treatments, p, cores)
  return(as.data.table(cooccurrence))
}



#' Bootstraps the pair-wise Spearman rank co-occurrence, to determine a significant rho-cutoff. phylosmith
#'
#' Bootstraps the pair-wise Spearman rank co-occurrence, to determine a significant rho-cutoff.
#' @useDynLib phylosmith
#' @usage bootstrap_rho(phyloseq_obj, treatment, replicates = 'independent', permutations = 100,
#' cores = 0)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package.
#' @param treatment Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param replicates Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}} that indicates which samples are non-independent of each other.
#' @param permutations Number of iterations to compute.
#' @param cores Number of CPU cores to use for the pair-wise permutations. Default uses all cores available.
#' @keywords nonparametric
#' @import data.table
#' @import phyloseq
#' @import RcppArmadillo
#' @import RcppParallel
#' @import RcppProgress
#' @seealso \code{\link{co_occurrence}}
#' @examples
#' data(mock_phyloseq)
#' bootstrap_rho(mock_phyloseq, treatment = "day", permutations = 10)
#' @export

# sourceCpp('src/FastCoOccur_rho_Rcpp.cpp')

bootstrap_rho <- function(phyloseq_obj, treatment, replicates = 'independent', permutations = 100, cores = 0){
  # phyloseq_obj = mock_phyloseq; treatment = c("treatment", "day"); replicates = 'independent'; permutations = 10; p = 0; cores = 0;
  options(warnings=-1)

  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  if(is.numeric(replicates)){replicates <- colnames(phyloseq_obj@sam_data[,replicates])}

  if(replicates == 'independent'){
    replicate_indices <- 1:ncol(phyloseq_obj@otu_table)
  } else {
    phyloseq_obj_reps <- combine_treatments(phyloseq_obj, c(treatment, replicates))
    replicate_name <- paste(c(treatment, replicates), collapse = ".")
    replicates <- as.character(unique(phyloseq_obj_reps@sam_data[[replicate_name]]))
    replicate_indices <- lapply(replicates, FUN = function(trt){which(as.character(phyloseq_obj_reps@sam_data[[replicate_name]]) %in% trt)})
  }

  rhos<-list()
  n <- nrow(phyloseq_obj@otu_table)
  permuted_phyloseq_obj <- phyloseq_obj

  tryCatch({
    for(i in 1:permutations){
      for(indices in replicate_indices){
        permuted_phyloseq_obj@otu_table[,indices] <- phyloseq_obj@otu_table[sample(1:n, n),indices]
      }
      rhos[[i]] <- co_occurrence_rho(permuted_phyloseq_obj, treatment, cores = cores)
    }},
    interrupt = function(interrupt){rhos <- rhos[-length(rhos)]; message('Interrupted after ', length(rhos), ' permutations.'); return(rhos)})

  rhos <- unlist(rhos)
  # if(p == 0){
  return(rhos)
} #else {
#   return(stats::quantile(rhos, 1-p, na.rm = TRUE))}
# }

# start_time <- Sys.time()
# for(indices in replicate_indices){
#   permuted_phyloseq_obj@otu_table[,indices] <- phyloseq_obj@otu_table[sample(1:n, n),indices]
# }
# rhos[1] <- mean(FastCoOccur(permuted_phyloseq_obj, treatment, cores = cores)$rho)
# end_time <- Sys.time()
# time_taken <- end_time - start_time
# message('Allow up to ', round((time_taken*permutations)/360, 1), ' hours to complete these calculations.')


#' pair-wise Spearman rank co-occurrence, written in efficient c++ code.
#'
#' A rewrite of the pair-wise Spearman rank co-occurrence routine written by \href{https://github.com/germs-lab/FastCoOccur}{Jin Choi}. The routine has been adapted to integrate with the \code{\link[Rcpp]{Rcpp-package}} API.
#' @useDynLib phylosmith
#' @usage co_occurrence_rho(phyloseq_obj, treatment, cores = 0)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package.
#' @param treatment Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param cores Number of CPU cores to use for the pair-wise permutations. Default uses all cores available.
#' @import phyloseq
#' @import RcppArmadillo
#' @import RcppParallel
#' @import RcppProgress
#' @keywords nonparametric
#' @seealso \code{\link{bootstrap_rho}} \code{\link{phylosmith}}
#' @examples
#' data(mock_phyloseq)
#' phylosmith:::co_occurrence_rho(mock_phyloseq, "day", 0.05)

# sourceCpp("src/co_occurrence_rho_Rcpp.cpp")

co_occurrence_rho <- function(phyloseq_obj, treatment, cores = 0){
  # phyloseq_obj = mock_phyloseq; treatment = c("treatment", "day"); p = 0.05; cores = 0
  options(warnings=-1)

  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment = treatment)
  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  treatment_name <- paste(treatment, collapse = ".")

  treatments <- as.character(unique(phyloseq_obj@sam_data[[treatment_name]]))
  treatment_indices <- lapply(treatments, FUN = function(trt){which(as.character(phyloseq_obj@sam_data[[treatment_name]]) %in% trt)-1})

  if(cores == 0){cores <- parallel::detectCores()}
  rhos <- co_occurrence_rho_Rcpp(phyloseq_obj@otu_table, treatment_indices, treatments, cores)
  return(rhos)
}

#' Curate co-occurrence data. Function from the phylosmith-package.
#'
#' Used to curate a co-occurrence table from the \code{\link{co_occurrence}} function. Takes a list of taxa and finds all pairs containing any of those taxa.
#' @useDynLib phylosmith
#' @usage curate_cooccurrence(cooccurrence_table, taxa_of_interest, number_of_treatments = 1)
#' @param cooccurrence_table co-occurrence table generated with \code{\link{co_occurrence}}, or formatted in the same way.
#' @param taxa_of_interest a list or vector of taxa names, like those generated with \code{\link{find_unique_taxa}}.
#' @param number_of_treatments how many treatments should the taxa of interest be seen in? require \code{integer} or 'all' (default = 1).
#' @keywords manip
#' @export
#' @import data.table
#' @import phyloseq
#' @import RcppArmadillo
#' @import RcppParallel
#' @seealso \code{\link{co_occurrence}}
#' @examples
#' data(mock_phyloseq)
#' curate_cooccurrence(co_occurrence(mock_phyloseq, 'day', p=1), c('tat', 'cct'))

curate_cooccurrence <- function(cooccurrence_table, taxa_of_interest, number_of_treatments = 1){
  sub_cooccurrence <- cooccurrence_table[(cooccurrence_table[[2]] %in% taxa_of_interest | cooccurrence_table[[3]] %in% taxa_of_interest),]
  toi_table <- unique(cbind(rbindlist(list(sub_cooccurrence[,1], sub_cooccurrence[,1])), rbindlist(list(sub_cooccurrence[,2], sub_cooccurrence[,3]))))
  toi_table <- toi_table[toi_table[[2]] %in% taxa_of_interest]
  if(number_of_treatments == 'all'){number_of_treatments <- length(unique(sub_cooccurrence[[1]]))
  } else {number_of_treatments <- number_of_treatments}
  toi <- names(table(toi_table[[2]])[table(toi_table[[2]]) >= number_of_treatments])
  sub_cooccurrence <- cooccurrence_table[(cooccurrence_table[[2]] %in% toi | cooccurrence_table[[3]] %in% toi),]

  # sourceCpp("src/arrange_cooccurrence_table_tbb.cpp")
  arranged_coocurrence <- as.data.table(arrange_cooccurr_table(sub_cooccurrence, toi))

  setorder(arranged_coocurrence, 'Treatment')
  return(arranged_coocurrence)
}

