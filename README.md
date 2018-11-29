---
title: "phyloschuyler"
author: "Schuyler D. smith"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## `phyloschuyler` description

`phyloschuyler` is a compilation of functions that I have written, that I find very useful when analyzing phyloseq objects.

## Installation

```
require(devtools)
install_github('schuyler-smith/phyloschuyler')
```

## Data

`phyloschuyler` includes two mock datasets. They are made to meet the phyloseq requirements, but do not represent real data; and therefore are not always perfect examples, but generally demonstrate how the functions operate.

## Usage

```
library(phyloschuyler)
data(mock_phyloseq); data(mock_phyloseq_2)

find_generalists(mock_phyloseq, frequency = 0.3)
find_generalists(mock_phyloseq, frequency = 0.3, treatments = "day")
find_generalists(mock_phyloseq, frequency = 0.3, treatments = 3)
find_generalists(mock_phyloseq, frequency = 0.3, treatments = c("day", "treatment"))
find_generalists(mock_phyloseq, frequency = 0.3, treatments = c("day", "treatment"), subset = "5-soil")

find_unique_taxa(mock_phyloseq, column = 2)
find_unique_taxa(mock_phyloseq, column = "day")

merge_asvs(mock_phyloseq, mock_phyloseq_2)
```
