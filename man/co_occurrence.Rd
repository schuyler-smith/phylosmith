% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/co_occurrence.R
\name{co_occurrence}
\alias{co_occurrence}
\title{Pair-wise Spearman rank co-occurrence.}
\usage{
co_occurrence(phyloseq_obj, treatment = NULL, subset = NULL,
rho = 0, p = 0.05, method = 'spearman', cores = 1)
}
\arguments{
\item{phyloseq_obj}{A \code{\link[phyloseq]{phyloseq-class}} object.}

\item{treatment}{Column name as a \code{string} or \code{numeric} in the
\code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
multiple columns and they will be combined into a new column.}

\item{subset}{A factor within the \code{treatment}. This will remove any
samples that to not contain this factor. This can be a vector of multiple
factors to subset on.}

\item{rho}{\code{numeric} The rho-value cutoff. All returned co-occurrences
will have a rho-value less than or equal to \code{rho} or less than or
equal to -\code{rho}.}

\item{p}{\code{numeric} The p-value cutoff. All returned co-occurrences
will have a p-value less than or equal to \code{p}.}

\item{method}{Which correlation method to calculate, "pearson", "spearman".}

\item{cores}{\code{numeric} Number of CPU cores to use for the pair-wise
permutations. Default (0) uses max cores available. Parallelization not
available for systems running MacOS without openMP configuration.}
}
\value{
data.table
}
\description{
A rewrite of the pair-wise Spearman rank co-occurrence routine written by
\href{https://github.com/germs-lab/FastCoOccur}{Jin Choi}. The routine has
been adapted to integrate with the \code{\link[Rcpp]{Rcpp-package}} API.
}
\examples{
co_occurrence(soil_column, treatment = c('Matrix', 'Treatment'),
subset = 'Amended', rho = 0.8, p = 0.05, cores = 1)
}
\seealso{
\code{\link{permute_rho}} \code{\link{phylosmith}}
}
\keyword{nonparametric}
