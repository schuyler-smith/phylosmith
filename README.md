
## `phylosmith`
[vignette](https://schuyler-smith.github.io/phylosmith/)

[![Travis Build
Status](https://travis-ci.org/schuyler-smith/phylosmith.svg?branch=master)](https://github.com/schuyler-smith/phylosmith)

A conglomeration of functions that I have written, that I find useful, for analyzing phyloseq objects. Phyloseq objects are a great data-standard for microbiome and gene-expression data.

## Installation

```
library(devtools)
install_github('schuyler-smith/phylosmith')
library(phylosmith)
```

*for WINDOWS you need to install <a href="https://cran.r-project.org/bin/windows/Rtools/" target="_blank" >Rtools</a>, when prompted, select `add rtools to system PATH`.*

## Functions

### [Data Parsing](https://schuyler-smith.github.io/phylosmith/data_parsing.html)

Call			     | Description
-------------------- | ------------------------------------------------------------
[common_taxa](https://schuyler-smith.github.io/phylosmith/data_parsing.html#common_taxa) | find taxa common to each treatment
[conglomerate_samples](https://schuyler-smith.github.io/phylosmith/data_parsing.html#conglomerate_samples)   |   Merge samples based on common factor within sample_data
[conglomerate_taxa](https://schuyler-smith.github.io/phylosmith/data_parsing.html#conglomerate_taxa)  |  conglomerate taxa by sample on a given classification level
[melt_phyloseq](https://schuyler-smith.github.io/phylosmith/data_parsing.html#melt_phyloseq)   |   Melt a phyloseq object into a data.table.
[merge_treatments](https://schuyler-smith.github.io/phylosmith/data_parsing.html#merge_treatments) | combines multiple columns in meta-data into a single column
[order_treatment](https://schuyler-smith.github.io/phylosmith/data_parsing.html#order_treatment) | sets the orders of the factors in a sample_data column (for ordering graphs)
[relative_abundance](https://schuyler-smith.github.io/phylosmith/data_parsing.html#relative_abundance) | transform abundance data to relative abundance
[taxa_proportions](https://schuyler-smith.github.io/phylosmith/data_parsing.html#taxa_proportions) | computes the proportion of a taxa classification
[taxa_filter](https://schuyler-smith.github.io/phylosmith/data_parsing.html#taxa_filter) | filter taxa by proportion of samples seen in
[unique_taxa](https://schuyler-smith.github.io/phylosmith/data_parsing.html#unique_taxa) | find taxa unique to each treatment

### [Graphs](https://schuyler-smith.github.io/phylosmith/graphs.html)

Call                 | Description
-------------------- | ------------------------------------------------------------
[abundance_heatmap_ggplot](https://schuyler-smith.github.io/phylosmith/graphs.html#abundance_heatmap_ggplot) | create a ggplot object of the heatmaps of the abndance table
[abundance_lines_ggplot](https://schuyler-smith.github.io/phylosmith/graphs.html#abundance_lines_ggplot) | create a ggplot object of the abundance data as a line graph
[network_phyloseq](https://schuyler-smith.github.io/phylosmith/graphs.html#network_phyloseq) | creates a ggplot object of the co-occurrence network of taxa
[nmds_phyloseq_ggplot](https://schuyler-smith.github.io/phylosmith/graphs.html#nmds_phyloseq_ggplot)  | create a ggplot object of the NMDS from a phyloseq object
[phylogeny_bars_ggplot](https://schuyler-smith.github.io/phylosmith/graphs.html#phylogeny_bars_ggplot) | create a ggplot barplot object of the compositons of each sample at a taxonomic level
[taxa_abundance_bars_ggplot](https://schuyler-smith.github.io/phylosmith/graphs.html#taxa_abundance_bars_ggplot) | create a ggplot object of the abundance of taxa in each sample
[tsne_phyloseq_ggplot](https://schuyler-smith.github.io/phylosmith/graphs.html#tsne_phyloseq_ggplot)  | create a ggplot object of the t-SNE from a phyloseq object

### [Calculations](https://schuyler-smith.github.io/phylosmith/calculations.html)

Call                 | Description
-------------------- | ------------------------------------------------------------
[co_occurrence](https://schuyler-smith.github.io/phylosmith/calculations.html#permute_rho) | calculate co-occurrence between taxa
[permute_rho](https://schuyler-smith.github.io/phylosmith/calculations.html#co_occurrence) | runs permutations of the otu_table to calculate a significant $\rho$ value
[histogram_permuted_rhos](https://schuyler-smith.github.io/phylosmith/calculations.html##histogram_permuted_rhos) | Create a ggplot object of the distribution of $\rho$ values.
[quantile_permuted_rhos](https://schuyler-smith.github.io/phylosmith/calculations.html##quantile_permuted_rhos) | calculate quantiles for the permuted $\rho$ values from the Spearman-rank co-occurrence

## Datasets

Originally I had created 2 mock phyloseq objects (`mock_phyloseq` and `mock_phyloseq2`) that had no real-world data but served to show simple examples of how the functions worked. 

Then I decided that I should include a real example of microbiome data (`soil_column`) becasue it's always nice to see real examples. `soil_column` is a <a href="https://www.frontiersin.org/articles/10.3389/fmicb.2018.03197/full" target="_blank" >published dataset</a>  from my lab-group. The data is from an experiment where they looked at the microbial composition of farmland soil before and after manure application, over time, using 16S-sequencing.
