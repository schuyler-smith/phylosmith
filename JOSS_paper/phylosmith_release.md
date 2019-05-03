---
title: 'phylosmith: an R-package for reproducible and efficient microbiome analysis with phyloseq-objects'
authors:
  - name: Schuyler D. Smith
    orcid: 0000-0001-5720-5611
    affiliation: 1
affiliations:
 - name: Department of Bioinformatics and Computational Biology, Iowa State University
   index: 1
tags:
  - R
  - microbiology
  - microbiome
  - metagenomics
  - microbial ecology
date: 1 May 2019
bibliography: phylosmith_release.bib
csl: biomed-central.csl
output:
  pdf_document: default
  html_document:
    df_print: paged
---

# Summary

Presented here is ``phylosmith``, an ``R``-package that enables reproducible and efficient analysis of microbiome data with ``phyloseq-class`` objects by providing robust and efficient functions. ``phylosmith`` utilizes the standardized data format of ``phyloseq`` and ``R`` object accession methods to provide functions with simple and intuitive input arguments. 

The functions provided in ``phylosmith`` have been divided into 3 categories.

## Data Wrangling

Functions that will either return a transformed version of the input ``phyloseq`` object, a subset version, or extracted data from it on set parameters. In some cases the functions are practical rewrites of ones already available from ``phyloseq`` with additional features or a more efficient implementation with ``data.table`` for large datasets. Various functions include finding shared or exclusive taxa by treatments, agglomeration, factor handling, and filtering.

## Graphs

The graphs are designed to serve as a quick and easy way to visualise data for 
analysis and to provide a foundation for figures for publishing. Graphics include ordinations, phylogeny profiles, and co-occurrence networks. All images are produced as ``ggplot`` objects (@ggplot2), allowing for the image to be altered and additional layers given to tailor the graphic as desired. Additionally, the code for producing the graphs is readily accessible, allowing for the code to be reused and tailored to fit needs, providing a foundation to start from. The most novel, for the field of microbiome research, is the implementation of a t-SNE ordination. Most studies have 
used PCA or NMDS, which can suffers from convergence getting stuck in local minima on large datasets, t-SNE is designed for large datasets and is not susceptible to these same limitations.

## Calculations

As of this announcement, the functions in this section all pertain to calculating and analyzing the Spearman rank co-occurrence. The routine was written in efficient C++ code and interfaced with R using the ``Rcpp`` API (@Rcpp). The resulting co-occurrence table matches that produced by the ``cor()`` function in the ``R stats`` package, but is calculated much faster on a single thread, with a multi-threading options implemented as well.

# Need

Adoption of data-standards enable data that are readily available for sharing and also the creation and implementation of tools for reproducible research. It is commonly said that in the age of big-data that biologists are required to have computational proficiency and literacy [@carey]. It seems reasonable that their should be a large onus on bioinformaticians to create accessible and practical tools that enable the biologists.

For the field of microbiome research, there has developed a formulaic approach to analysis. A generic study will incorporate some combination and implementation of the same resulting figures; ordination, profile bar-chart, heatmap, network, etc [@generic1], [@generic2], [@generic3]. Each new scientist, often from a biology, microbiology, ecology, or environmental science background, is required to learn how to produce these analyses and figures. A lot of time is spent learning how to generate these plots. Even more time is spent learning how to process data; which can easily be done incorrectly without being apparently obvious (i.e. incorrect logical subsetting, factor levels set incorrectly, or even reordering of samples due to string sorting methods), leading to incorrect results and conclusions.

For microbiome researchers using the ``R`` statistical programming language [@R], a data-standard has been available in ``phyloseq`` [@phyloseq]. ``phyloseq`` provides an S4-class object that contains a count table, taxa table, and associated metadata, along with a phylogenetic tree slot and reference sequence slot. For beginning and medium, users of ``R``, S4 objects can be a barrier as they require an additional layer of accession methods compared to the base S3 objects. ``phyloseq`` offers several functions for handling its objects, as well as functions for producing some common figures, but is by no means a complete toolset. Additionally, when the authors originally wrote ``phyloseq``, advanced tools such as the package ``data.table`` [@datatable] were not practically available and thus had not been implemented within the program. 

Providing tools for reproducible and efficient research can help microbiome 
researchers to focus more effort on answering biological. Providing simple implementations of tools, such as t-SNE [@tsne], can increase the acceptance and adoption of new techniques in a field that is hesitant to do so. The importance of these tools should not be overlooked for the importance of science and understanding as a whole.

# Acknowledgement

This work was supported by the National Science Foundation Directorate of Biological Sciences under awards DEB 1737758 and DEB 1737765.

# References

