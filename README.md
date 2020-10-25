
## `phylosmith`
[vignette](https://schuyler-smith.github.io/phylosmith/)

[![Travis Build Status](https://travis-ci.org/schuyler-smith/phylosmith.svg?branch=master)](https://github.com/schuyler-smith/phylosmith) [![status](http://joss.theoj.org/papers/4d4780ab16c487764ebe108fa6bcdc2c/status.svg)](http://joss.theoj.org/papers/4d4780ab16c487764ebe108fa6bcdc2c)
[![DOI](https://zenodo.org/badge/159598780.svg)](https://zenodo.org/badge/latestdoi/159598780)

A conglomeration of functions that I have written, that I find useful, for analyzing phyloseq objects. Phyloseq objects are a great data-standard for microbiome and gene-expression data.

# Installation
## Requirements
### Linux Systems 
For some Linux systems you may need to install the following two programs through your terminal.

Ubuntu example:
```
sudo apt install libmysqlclient-dev libgdal-dev libudunits2-dev
```
These programs are required by some dependencies and may not come in your default OS distribution.

### Windows Systems 

if you are working on WINDOWS you likely need to install the CRAN program <a href="https://cran.r-project.org/bin/windows/Rtools/" target="_blank" >Rtools</a>.When prompted, select `add rtools to system PATH`.

### R

phylosmith depends on the usage of the [phyloseq package](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061217) released by [Dr. Paul McMurdie](https://github.com/joey711). The package is maintained on BioConductor, and can be installed through R using the following commands:

```
if(!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
} 
BiocManager::install("phyloseq")
```

Additionally, the package imports a number of other packages to use their advanced functions. These packages may install with the phylosmith installation, but it is always best to install independently.

```
install.packages(c("devtools", "RcppEigen", "RcppParallel", "Rtsne", "ggforce", "units"))
```

## phylosmith

The package is hosted on Github, and can be installed through R with:
```
devtools::install_github('schuyler-smith/phylosmith')
library(phylosmith)
```


