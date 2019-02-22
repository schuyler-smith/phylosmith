#' classify_ARG_genes
#'
#' Classifies ARGs from the \href{https://card.mcmaster.ca/home}{CARD database} using either accession number or gene name.
#' @usage classify_ARG_genes(phyloseq_obj, genes, obo = NULL)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package.
#' @param genes Column name or number in the \code{\link[phyloseq:tax_table]{tax_table}} where ARGs are as gene names or ARO accession number.
#' @param obo an object processed with \code{\link{processOBO}}.
#' @import stringr

classify_ARG_genes <- function(phyloseq_obj, genes, obo = NULL){
  if(is.numeric(genes)){genes <- colnames(phyloseq_obj@tax_table[,genes])}
  if(!(is.null(obo))){CARD <- obo}
  genes <- data.frame(phyloseq_obj@tax_table[,genes], stringsAsFactors = FALSE)
  classified_genes <- unlist(unname(apply(genes,1,FUN = function(gene){
    if(gene %in% CARD$ID){
      if(is.na(CARD$Resistance[gene])){return('Unclassified')}
      return(CARD$Resistance[gene])
    } else if(gene %in% CARD$Name){gene <- CARD$ID[which(CARD$Name %in% gene)]
      if(is.na(CARD$Resistance[gene])){return('Unclassified')}
      return(CARD$Resistance[gene])
    } else {return('Unclassified')}
  })))
  phyloseq_obj@tax_table <- tax_table(cbind(phyloseq_obj@tax_table, classified_genes))
  return(phyloseq_obj)
}

#' classify_ARG_mechanisms
#'
#' Classifies ARGs mechanism of function from the \href{https://card.mcmaster.ca/home}{CARD database} using either accession number or gene name.
#' @usage classify_ARG_mechanisms(phyloseq_obj, genes, obo = NULL)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package.
#' @param genes Column name or number in the \code{\link[phyloseq:tax_table]{tax_table}} where ARGs are as gene names or ARO accession number.
#' @param obo an object processed with \code{\link{processOBO}}.
#' @import stringr

classify_ARG_mechanisms <- function(phyloseq_obj, genes, obo = NULL){
  if(is.numeric(genes)){genes <- colnames(phyloseq_obj@tax_table[,genes])}
  if(!(is.null(obo))){CARD <- obo}
  genes <- data.frame(phyloseq_obj@tax_table[,genes], stringsAsFactors = FALSE)
  classified_genes <- unlist(unname(apply(genes,1,FUN = function(gene){
    if(gene %in% CARD$ID){
      if(is.na(CARD$Mechanism[gene])){return('Unclassified')}
      return(CARD$Mechanism[gene])
    } else if(gene %in% CARD$Name){gene <- CARD$ID[which(CARD$Name %in% gene)]
    if(is.na(CARD$Mechanism[gene])){return('Unclassified')}
    return(CARD$Mechanism[gene])
    } else {return('Unclassified')}
  })))
  phyloseq_obj@tax_table <- tax_table(cbind(phyloseq_obj@tax_table, classified_genes))
  return(phyloseq_obj)
}

#' processOBO
#'
#' Converts the aro.obo file from the \href{https://card.mcmaster.ca/home}{CARD database} into a data.table.
#' @usage processOBO(OBO_filepath)
#' @param OBO_filepath filepath to the aro.obo file.
#' @import stringr

processOBO <- function(OBO_filepath){
  con <- file(OBO_filepath, "r")
  OBO <- list()
  while(TRUE){
    line <- readLines(con, n = 1)
    if (length(line) == 0){break}
    if(length(grep('id: ARO', line)) > 0){
      ID <- str_split(line, ': ')[[1]][2]
      OBO$ID[ID] <- ID
      OBO$Resistance[ID] <- NA
      OBO$is_a[ID] <- NA
      OBO$part_of[ID] <- NA
      OBO$Mechanism[ID] <- NA
      res <- NULL}
    if(length(grep('name: ', line)) > 0){
      name <- str_split(line, ': ')[[1]][2]
      OBO$Name[ID] <- name}

    if(length(grep('confers_resistance_to', line)) > 0){
      res <- c(res, str_split(line, '! ')[[1]][2])
      OBO$Resistance[ID] <- gsub(' Antibiotic', '', str_to_title(paste(res, collapse = ', ')))}

    if(length(grep('is_a:', line)) > 0){
      if(is.na(OBO$is_a[ID])){OBO$is_a[ID] <- str_split(str_split(line, ' ! ')[[1]][1], ': ')[[1]][2]}}
    if(length(grep(': part_of', line)) > 0){
      if(is.na(OBO$part_of[ID])){OBO$part_of[ID] <- str_split(str_split(line, ' ! ')[[1]][1], 'part_of ')[[1]][2]}}
    if(length(grep(': regulates', line)) > 0){
      if(is.na(OBO$part_of[ID])){OBO$part_of[ID] <- str_split(str_split(line, ' ! ')[[1]][1], 'regulates ')[[1]][2]}}
    if(length(grep(': participates_in', line)) > 0){
      if(is.na(OBO$Mechanism[ID])){OBO$Mechanism[ID] <- str_split(line, ' ! ')[[1]][2]}}
  };  close(con)
  for(i in 1:length(OBO$ID)){
    gene <- OBO$ID[i]
    res <- OBO$Resistance[gene]
    part_of <- OBO$part_of[gene]
    is_a <- OBO$is_a[gene]
    for(j in 1:6){
      if(!(is.na(OBO$Resistance[part_of])) & is.na(res)){
        OBO$Resistance[gene] <- OBO$Resistance[part_of]}
      if(!(is.na(OBO$Resistance[is_a])) & is.na(res)){
        OBO$Resistance[gene] <- OBO$Resistance[is_a]}
      if(is.na(res)){
        part_of <- OBO$part_of[part_of]
        is_a <- OBO$is_a[is_a]}
      if(!(is.na(OBO$Resistance[gene]))){break}
    }
    mech <- OBO$Mechanism[gene]
    part_of <- OBO$part_of[gene]
    is_a <- OBO$is_a[gene]
    for(j in 1:10){
      if(!(is.na(OBO$Mechanism[part_of])) & is.na(mech)){
        OBO$Mechanism[gene] <- OBO$Mechanism[part_of]}
      if(!(is.na(OBO$Mechanism[is_a])) & is.na(mech)){
        OBO$Mechanism[gene] <- OBO$Mechanism[is_a]}
      if(is.na(mech)){
        part_of <- OBO$part_of[part_of]
        is_a <- OBO$is_a[is_a]}
      if(!(is.na(OBO$Mechanism[gene]))){break}
    }
  }
  return(OBO)
}

# OBO_filepath <- '~/Downloads/card-ontology/aro.obo'
# CARD <- processOBO(OBO_filepath)
# test <- readRDS('~/Dropbox/Co-occur_ARGs/data/arg_phy_updated.RDS')@tax_table[,4]
# classify_ARG_mechanism(test)
#
# unique(classify_ARG_genes(test))
# save(CARD, file = 'R/sysdata.rda')

