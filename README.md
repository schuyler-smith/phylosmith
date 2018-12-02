
## `phylosmith`

This package is a compilation of functions that I have written, that I find useful, for analyzing phyloseq objects.

## Installation

```
require(devtools)
install_github('schuyler-smith/phylosmith')
```

for WINDOWS you need to install <a href="https://cran.r-project.org/bin/windows/Rtools/" target="_blank" >Rtools</a>, when prompted, select `add rtools to system PATH`.

## Data

`phylosmith` includes two mock datasets. They are made to meet the phyloseq requirements, but do not represent real data; and therefore are not always perfect examples, but generally demonstrate how the functions operate.

## Functions

<!-- ```
library(phylosmith)
data()

FastCoOccur(phyloseq_obj, treatment, p = 0.05)
find_generalists(phyloseq_obj, frequency = 0, treatments = NULL, subset = NULL, below = FALSE, drop_samples = FALSE)
find_unique_taxa(phyloseq_obj, column, keyword = NULL)
merge_asvs(...)
``` -->

---------------- | ------------------------------------------------
FastCoOccur      | co-occurrence
find_generalists | filter taxa by proportion of samples seen in
find_unique_taxa | find taxa seen only in particular treatment
merge_asvs       | combine ASV to lowest common biological sequence
