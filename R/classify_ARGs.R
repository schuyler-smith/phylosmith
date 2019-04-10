#' classify_ARG_classes
#'
#' Classifies ARGs from the \href{https://card.mcmaster.ca/home}{CARD
#' database} using either accession number or gene name.
#' @usage classify_ARG_classes(phyloseq_obj, genes, combine = 0, obo = NULL)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}})
#' with information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param genes Column name or number in the
#' \code{\link[phyloseq:tax_table]{tax_table}} where ARGs are as gene names or
#' ARO accession number.
#' @param combine The number of Classes a gene must belong to to be set to
#' 'Multiple_Resistance'.
#' @param obo an object processed with \code{\link{processOBO}}.
#' @import stringr
#' @import data.table
#' @return phyloseq-object

classify_ARG_classes <- function(phyloseq_obj, genes, combine = 0, obo = NULL){
    if(is.numeric(genes)){genes <- colnames(phyloseq_obj@tax_table[,genes])}
    if(!(is.null(obo))){CARD <- obo}
    genes <- data.frame(phyloseq_obj@tax_table[,genes],
        stringsAsFactors = FALSE)
    ARG_Class <- unlist(unname(apply(genes,1,FUN = function(gene){
        if(gene %in% CARD$ID){
            if(is.na(CARD$Resistance[gene])){res <- 'Unclassified'
            } else {res <- CARD$Resistance[gene]}
        } else if(gene %in% CARD$Name){
            gene <- CARD$ID[which(CARD$Name %in% gene)]
            if(is.na(CARD$Resistance[gene])){res <- 'Unclassified'
            } else {res <- CARD$Resistance[gene]}
        } else {res <- 'Unclassified'}
        if(combine == 0 | combine == 1){return(res)}
        if(length(str_split(res, ', ')[[1]]) >= combine){
            res <- 'Multiple_Resistance'}
        return(res)
    })))
    tax_tab <- data.table(as(phyloseq_obj@tax_table, 'matrix'))
    tax_tab[['ARG_Class']] = ARG_Class
    phyloseq_obj@tax_table <- tax_table(as.matrix(tax_tab,
        rownames = taxa_names(phyloseq_obj)))
    return(phyloseq_obj)
}

#' classify_ARG_mechanisms
#'
#' Classifies ARGs mechanism of function from the
#' \href{https://card.mcmaster.ca/home}{CARD database} using either accession
#' number or gene name.
#' @usage classify_ARG_mechanisms(phyloseq_obj, genes, obo = NULL)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object. It
#' must contain \code{\link[phyloseq:sample_data]{sample_data()}}) with
#' information about each sample, and it must contain
#' \code{\link[phyloseq:tax_table]{tax_table()}}) with information about each
#' taxa/gene.
#' @param genes Column name or number in the
#' \code{\link[phyloseq:tax_table]{tax_table}} where ARGs are as gene names or
#' ARO accession number.
#' @param obo an object processed with \code{\link{processOBO}}.
#' @import stringr
#' @import data.table
#' @return phyloseq-object

classify_ARG_mechanisms <- function(phyloseq_obj, genes, obo = NULL){
    if(is.numeric(genes)){genes <- colnames(phyloseq_obj@tax_table[,genes])}
    if(!(is.null(obo))){CARD <- obo}
    genes <- data.frame(phyloseq_obj@tax_table[,genes],
        stringsAsFactors = FALSE)
    ARG_Mechanism <- unlist(unname(apply(genes,1,FUN = function(gene){
        if(gene %in% CARD$ID){
            if(is.na(CARD$Mechanism[gene])){return('Unclassified')}
            return(CARD$Mechanism[gene])
        } else if(gene %in% CARD$Name){
            gene <- CARD$ID[which(CARD$Name %in% gene)]
        if(is.na(CARD$Mechanism[gene])){return('Unclassified')}
        return(CARD$Mechanism[gene])
        } else {return('Unclassified')}
    })))
    tax_tab <- data.table(as(phyloseq_obj@tax_table, 'matrix'))
    tax_tab[['ARG_Mechanism']] = ARG_Mechanism
    phyloseq_obj@tax_table <- tax_table(as.matrix(tax_tab,
        rownames = taxa_names(phyloseq_obj)))
    return(phyloseq_obj)
}

#' processOBO
#'
#' Converts the aro.obo file from the
#' \href{https://card.mcmaster.ca/home}{CARD database} into a data.table.
#' @usage processOBO(OBO_filepath)
#' @param OBO_filepath filepath to the aro.obo file.
#' @import stringr
#' @return list

processOBO <- function(OBO_filepath){
    con <- file(OBO_filepath, "r")
    OBO <- list()
    while(TRUE){
        line <- readLines(con, n = 1)
        if(length(line) == 0){break}
        if(line == '[Typedef]'){break}
        if(length(grep('id: ARO', line)) > 0){
            ID <- str_split(line, ': ')[[1]][2]
            OBO$ID[ID] <- ID
            OBO$Resistance[ID] <- NA
            OBO$is_a[ID] <- NA
            OBO$part_of[ID] <- NA
            OBO$Mechanism[ID] <- NA
            res <- NULL}
        if(length(grep('name: ', line)) > 0){
            name <- str_to_title(gsub(' antibiotic', '',
                str_split(line, ': ')[[1]][2]))
            OBO$Name[ID] <- name}

        if(length(grep('confers_resistance_to', line)) > 0){
            res <- c(res, str_split(line, ' ')[[1]][3])
            OBO$Resistance[ID] <- paste(res, collapse = ', ')}

        if(length(grep('is_a:', line)) > 0){
            if(is.na(OBO$is_a[ID])){
                OBO$is_a[ID] <- str_split(str_split(line, ' ! ')[[1]][1],
                ': ')[[1]][2]}}
        if(length(grep(': part_of', line)) > 0){
            if(is.na(OBO$part_of[ID])){
                OBO$part_of[ID] <- str_split(str_split(line, ' ! ')[[1]][1],
                'part_of ')[[1]][2]}}
        if(length(grep(': regulates', line)) > 0){
            if(is.na(OBO$part_of[ID])){
                OBO$part_of[ID] <- str_split(str_split(line, ' ! ')[[1]][1],
                'regulates ')[[1]][2]}}
        if(length(grep(': participates_in', line)) > 0){
            if(is.na(OBO$Mechanism[ID])){
                OBO$Mechanism[ID] <- str_split(line, ' ! ')[[1]][2]}}
    }; close(con)
    for(i in seq_along(OBO$ID)){
        gene <- OBO$ID[i]
        res <- OBO$Resistance[gene]
        part_of <- OBO$part_of[gene]
        is_a <- OBO$is_a[gene]
        if(is.na(res)){
            for(j in seq(6)){
                if(!(is.na(OBO$Resistance[part_of]))){
                    res <- OBO$Resistance[part_of]}
                if(!(is.na(OBO$Resistance[is_a]))){
                    res <- OBO$Resistance[is_a]}
                if(is.na(res)){
                    part_of <- OBO$part_of[part_of]
                    is_a <- OBO$is_a[is_a]}
                if(!(is.na(res))){break}
            }
        }
    ##removing this return specific resistance, instead of highest macro class
        if(!(is.na(res))){
            if(res != 'Undefined'){
                res <- paste(unique(sapply(str_split(res, ', ')[[1]],
                    FUN = function(class){
                    while(TRUE){
                        is_a <- OBO$is_a[class]
                        if(is_a != 'ARO:1000003'){class <- is_a
                        } else {break}}
                    return(class)
                })), collapse = ', ')
            }
        }
        OBO$Resistance[gene] <- res
        mech <- OBO$Mechanism[gene]
        part_of <- OBO$part_of[gene]
        is_a <- OBO$is_a[gene]
        if(is.na(mech)){
            for(j in seq(10)){
                if(!(is.na(OBO$Mechanism[part_of]))){
                    mech <- OBO$Mechanism[part_of]}
                if(!(is.na(OBO$Mechanism[is_a]))){
                    mech <- OBO$Mechanism[is_a]}
                if(is.na(mech)){
                    part_of <- OBO$part_of[part_of]
                    is_a <- OBO$is_a[is_a]}
                if(!(is.na(mech))){break}
            }
            OBO$Mechanism[gene] <- mech
        }
    }
    for(i in seq_along(OBO$Resistance)){
        if(is.na(OBO$Resistance[i])){OBO$Resistance[i] <- 'Undefined'
        } else {
            OBO$Resistance[i] <- paste(sapply(str_split(OBO$Resistance[i],
                ', ')[[1]], FUN = function(id){
            return(OBO$Name[id])
        }), collapse = ', ')
        }
    }
    return(OBO)
}



