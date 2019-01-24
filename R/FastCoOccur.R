#' pair-wise Spearman rank co-occurrence, written in efficient c++ code. Function from the phylosmith-package.
#'
#' A rewrite of the pair-wise Spearman rank co-occurrence routine written by \href{https://github.com/germs-lab/FastCoOccur}{Jin Choi}. The routine has been adapted to integrate with the \link[=Rcpp]{Rcpp} API.
#' @useDynLib phylosmith
#' @usage FastCoOccur(phyloseq_obj, treatment, p = 0.05, cores = 0)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package.
#' @param treatment Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param p the p-value cutoff. all returned co-occurrences must have a p-value less than or equal to p.
#' @param cores Number of CPU cores to use for the pair-wise permutations. Default uses all cores available.
#' @import data.table
#' @import phyloseq
#' @import RcppArmadillo
#' @import RcppParallel
#' @import RcppProgress
#' @keywords nonparametric
#' @seealso \code{\link{bootstrap_rho}} \code{\link{phylosmith}}
#' @examples
#' data(mock_phyloseq)
#' FastCoOccur(mock_phyloseq, "day", 0.05)
#' @export

# sourceCpp("src/FastCoOccur_Rcpp.cpp")

FastCoOccur <- function(phyloseq_obj, treatment, p = 0.05, cores = 0){
  # phyloseq_obj = mock_phyloseq; treatment = c("treatment", "day"); p = 0.05
  options(warnings=-1)

  phyloseq_obj <- combine_treatments(phyloseq_obj, treatment)
  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  treatment_name <- paste(treatment, collapse = ".")

  treatments <- as.character(unique(phyloseq_obj@sam_data[[treatment_name]]))
  treatment_indices <- lapply(treatments, FUN = function(trt){which(as.character(phyloseq_obj@sam_data[[treatment_name]]) %in% trt)-1})

  ###RAM check .. because apparently some people have 36,000 genes in their studies...
  n <- nrow(phyloseq_obj@otu_table)
  memory_function <- function(n, t, r){(6.780e-07 + (2.867e-08 * (n*(n-1)*t))) - r}
  required_memory <- round((6.780e-07 + (2.867e-08 * (n*(n-1)*length(treatments)))), 2)
  available_memory <- round(as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern = TRUE))/(1024^2), 2)
  if(required_memory >= available_memory){
    message("Calculating ", (n*(n-1)*length(treatments)), " pairs-wise co-occurences.
 This may take ", required_memory, " GB of RAM, or more.
 This machine has ", available_memory, " GB of RAM available.
 Recommend subsetting from ", n, " taxa to ",
   floor(uniroot(memory_function, t=length(treatments), r=available_memory, lower=0, upper=1000000)$root), ".
 ... Until I decide a practical way to fix this.."
    )
  }

  if(cores == 0){cores <- parallel::detectCores()}
  cooccurrence <- FastCoOccur_Rcpp(phyloseq_obj@otu_table, treatment_indices = treatment_indices, treatment_names = treatments, p_cutoff = p, ncores = cores)
  return(as.data.table(cooccurrence))
}
