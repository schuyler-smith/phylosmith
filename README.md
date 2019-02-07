
## `phylosmith`

This package is a compilation of functions that I have written, that I find useful, for analyzing phyloseq objects.

## Installation

```
library(devtools)
install_github('schuyler-smith/phylosmith')
library(phylosmith)
```

*for WINDOWS you need to install <a href="https://cran.r-project.org/bin/windows/Rtools/" target="_blank" >Rtools</a>, when prompted, select `add rtools to system PATH`.*

## Functions

Call			 | Use
---------------- | ------------------------------------------------
FastCoOccur      | calculate co-occurrence between taxa
bootstrap_rho | runs permutations of the otu_table to calculate a significant $\rho$ value
concatenate_treatments | combines multiple columns in meta-data into a single column
curate_cooccurrence | subsets the co-occurence table to specific taxa
find_common_taxa | find taxa common to each treatment
find_unique_taxa | find taxa unique to each treatment
merge_asvs       | combine ASVs to lowest common biological sequence
nmds_phyloseq_ggplot  | create a ggplot object of the NMDS from a phyloseq object
relative_abundance | transform abundance data to relative abundance
taxa_filter | filter taxa by proportion of samples seen in
tsne_phyloseq_ggplot  | create a ggplot object of the t-SNE from a phyloseq object

## Data

`phylosmith` includes two mock datasets. They are made to meet the phyloseq requirements, but do not represent real data; and therefore are not always perfect examples, but generally demonstrate how the functions operate.
