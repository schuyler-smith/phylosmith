% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graphs.R
\name{network_ps}
\alias{network_ps}
\title{Create a node network ggplot object of the co-occurrence from a phyloseq object.
Function from the phylosmith-package.}
\usage{
network_ps(phyloseq_obj, classification = NULL,
treatment = NULL, subset = NULL, co_occurrence_table = NULL, layout = NULL,
nodes_of_interest = NULL, node_colors = 'default',
cluster = FALSE, cluster_colors = 'default', buffer = 0.5)
}
\arguments{
\item{phyloseq_obj}{A \code{\link[phyloseq]{phyloseq-class}} object. It
must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
information about each sample, and it must contain
\code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
taxa/gene.}

\item{classification}{Column name as a string or number in the
\code{\link[phyloseq:tax_table]{tax_table}} for the factor to use for node
colors.}

\item{treatment}{Column name as a string or number in the
\code{\link[phyloseq:sample_data]{sample_data}}. This can be a vector of
multiple columns and they will be combined into a new column.}

\item{subset}{A factor within the \code{treatment}. This will remove any
samples that to not contain this factor. This can be a vector of multiple
factors to subset on.}

\item{co_occurrence_table}{Table of the co-occurrence of taxa/genes in the
\code{phyloseq_obj}, computed using \code{\link{co_occurrence}}. If no
table is given, it will be computed with the \code{phyloseq_obj}, using the
given \code{treatment} and \code{p} = 0.05.}

\item{layout}{(optional) an igraph layout of the network, for reproducibility.
Can be created with \code{\link{network_layout_ps}}.}

\item{nodes_of_interest}{A vector of names of classes within the
\code{classification} to be labeled.}

\item{node_colors}{Name of a color set from the
\link[=RColorBrewer]{RColorBrewer} package or a vector palete of R accepted
colors.}

\item{cluster}{if \code{TRUE}, will use igraph's
\code{\link[igraph:cluster_fast_greedy]{cluster_fast_greedy}} method.
Alternatively, you may pass a vector of cluster assignments with order
corresponding to the order of the \code{taxa_names} in the
\code{phyloseq_obj}.}

\item{cluster_colors}{Name of a color set from the
\link[=RColorBrewer]{RColorBrewer} package or a vector palete of R accepted
colors to use for the clusters.}

\item{buffer}{Amount of space beyond the points to extend the cluster (
aesthetic argument).}
}
\value{
ggplot-object
}
\description{
Inputs a \code{\link[phyloseq]{phyloseq-class}} object and
creates a network from the co-occurrence. The co-occurrence can either be
input, or it will be calculated with the Spearman-rank correlation. Also,
the layout of the graph can be given as an argument as well for reprodusibility.
}
\examples{
#network_ps(soil_column, treatment = c('Matrix', 'Treatment'),
#subset = 'Soil Manure', co_occurrence_table = NULL, layout = NULL,
#classification = 'phylum')
}