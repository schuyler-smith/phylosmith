
## `phylosmith`

A conglomeration of functions that I have written, that I find useful, for analyzing phyloseq objects. Phyloseq objects are a great data-standard for microbiome and gene-expression data.

## Installation

```
library(devtools)
install_github('schuyler-smith/phylosmith')
library(phylosmith)
```

*for WINDOWS you need to install <a href="https://cran.r-project.org/bin/windows/Rtools/" target="_blank" >Rtools</a>, when prompted, select `add rtools to system PATH`.*

## Functions

### Data Curation

Call			 | Use
---------------- | ------------------------------------------------
combine_treatments | combines multiple columns in meta-data into a single column
find_common_taxa | find taxa common to each treatment
find_unique_taxa | find taxa unique to each treatment
merge_asvs       | combine ASVs to lowest common biological sequence
order_phyloseq_metadata | sets the orders of the factors in a sample_data column (for ordering graphs)
relative_abundance | transform abundance data to relative abundance
taxa_filter | filter taxa by proportion of samples seen in

### Graphs

Call			 | Use
---------------- | ------------------------------------------------
abundance_heatmap_ggplot | create a ggplot object of the heatmaps of the abndance table
abundance_lines_ggplot | create a ggplot object of the abundance data as a line graph
network_phyloseq | creates a ggplot object of the co-occurrence network of taxa
nmds_phyloseq_ggplot  | create a ggplot object of the NMDS from a phyloseq object
phylogeny_bars_ggplot | create a ggplot barplot object of the compositons of each sample at a taxonomic level
taxa_abundance_bars_ggplot | create a ggplot object of the abundance of taxa in each sample
tsne_phyloseq_ggplot  | create a ggplot object of the t-SNE from a phyloseq object

### Calculations

Call			 | Use
---------------- | ------------------------------------------------
co_occurrence | calculate co-occurrence between taxa
curate_cooccurrence | subsets the co-occurence table to specific taxa
bootstrap_rho | runs permutations of the otu_table to calculate a significant $\rho$ value

## Data

`phylosmith` includes two mock datasets. They are made to meet the phyloseq requirements, but do not represent real data. Therefore, they are not always perfect examples, but generally demonstrate how the functions operate.
