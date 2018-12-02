
## `phylosmith`

This package is a compilation of functions that I have written, that I find useful, for analyzing phyloseq objects.

## Installation

```
library(devtools)
install_github('schuyler-smith/phylosmith')
```

*for WINDOWS you need to install <a href="https://cran.r-project.org/bin/windows/Rtools/" target="_blank" >Rtools</a>, when prompted, select `add rtools to system PATH`.*

## Functions

Call			 | Use
---------------- | ------------------------------------------------
FastCoOccur      | calculate co-occurrence between taxa
find_generalists | filter taxa by proportion of samples seen in
find_unique_taxa | find taxa unique to each treatment
merge_asvs       | combine ASVs to lowest common biological sequence

## Data

`phylosmith` includes two mock datasets. They are made to meet the phyloseq requirements, but do not represent real data; and therefore are not always perfect examples, but generally demonstrate how the functions operate.
