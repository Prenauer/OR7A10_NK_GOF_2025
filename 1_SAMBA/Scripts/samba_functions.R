#!/usr/bin/env Rscript

library(plyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(edgeR)


## Mageck dataset preparation
Prep_Mageck <- function(counts, design, paired = T, DIR_OUT, PREFIX){
    PATHS <- list(counts = file.path(DIR_OUT,paste0(PREFIX,'_counts.txt')),
                  ntc = file.path(DIR_OUT,paste0(PREFIX,'_ntc.txt')),
                  output = file.path(DIR_OUT,paste0(PREFIX,'_output')),
                  script = file.path(DIR_OUT,paste0(PREFIX,'_script.sh')))
    
    ntc <- counts[which(counts$Gene == 'NTC'),'sgRNA']
    screen  <- paste0(design$Screen, collapse = ',')
    ifelse(paired, script.paired <- '--paired ', script.paired <- '')
    ifelse(paired,
           ctrl <- paste0(design$Control, collapse = ','),
           ctrl <- paste0(unique(design$Control), collapse = ','))
    script <- c('#!/bin/bash',
                paste0('mageck test --remove-zero both ',script.paired,'-k "',
                       PATHS$counts, '" -t ',screen,' -c ',ctrl,' --control-sgrna "', 
                       PATHS$ntc, '" -n "',PATHS$output, '"'))
    write.table(counts, PATHS$counts, sep = '\t', quote = F, row.names = F)
    write.table(ntc, PATHS$ntc, sep = '\t', quote = F, row.names = F, col.names = F)
    write.table(script, PATHS$script, sep = '\t', quote = F, row.names = F, col.names = F)
}



## Riger dataset preparation
Prep_Riger <- function(counts, design, paired = T, agg.column = 'Covariate1', ntc.name = 'NTC', DIR_OUT, PREFIX){
    PATHS <- list(in.up = file.path(DIR_OUT,paste0(PREFIX,'_in_up.txt')),
                  in.dn = file.path(DIR_OUT,paste0(PREFIX,'_in_dn.txt')),
                  out.up.ks = file.path(DIR_OUT,paste0(PREFIX,'_out_up_ks.txt')),
                  out.dn.ks = file.path(DIR_OUT,paste0(PREFIX,'_out_dn_ks.txt')),
                  out.up.ws = file.path(DIR_OUT,paste0(PREFIX,'_out_up_ws.txt')),
                  out.dn.ws = file.path(DIR_OUT,paste0(PREFIX,'_out_dn_ws.txt')),
                  out.up.sb = file.path(DIR_OUT,paste0(PREFIX,'_out_up_sb.txt')),
                  out.dn.sb = file.path(DIR_OUT,paste0(PREFIX,'_out_dn_sb.txt')),
                  script = file.path(DIR_OUT,paste0(PREFIX,'_script.sh')))
    # remove single-guide genes and filter undetected sgRNA  
    sg.keep <- counts[which(rowSums(counts[,3:ncol(counts)]) > 0), 'sgRNA']
    counts <- counts[which(counts$sgRNA %in% sg.keep),]
    gene.keep <- counts[with(counts, order(Gene)),1:2]
    gene.keep = ddply(gene.keep,'Gene', summarise, length(Gene))
    gene.keep <- gene.keep[which(gene.keep$`length(Gene)` > 1),'Gene']
    counts <- counts[which(counts$Gene %in% gene.keep),]
    # log-normalize
    lcpm <- cbind(counts[,1:2], apply(counts[,3:ncol(counts)], 2, function(x) { log2(1 + (x*1e6)/sum(x)) }))
    # get lfc data
    lfc = apply(design,1, function(x) {(lcpm[,x[1]] - lcpm[,x[2]])})
    colnames(lfc) <- design$Screen
    # get mean LFC 
    s <- dlply(design, agg.column, summarise, Screen)
    mlfc <- lcpm[,1:2]
    for(group in names(s)) mlfc <- cbind(mlfc, rowMeans(lfc[,s[[group]]$Screen]))
    colnames(mlfc)[3:ncol(mlfc)] <- names(s)
    mlfc$Pooled <- rowMeans(lfc)
    mlfc <- mlfc[which(mlfc$Gene != ntc.name),]
    # write riger input data
    write.table(Format_Riger(mlfc, col.score = 'Pooled', 'pos'), PATHS$in.up, sep = '\t', row.names = F, col.names = T, quote = F)
    write.table(Format_Riger(mlfc, col.score = 'Pooled', 'neg'), PATHS$in.dn, sep = '\t', row.names = F, col.names = T, quote = F)
    # write riger script
    script <- c('#!/bin/bash',
                'PATH_RIGER="01_SAMBA/rigerj-2.0.2-assembly.jar"',
                paste0('java -jar ${PATH_RIGER} -scoringMethod KSbyScore -randomSeed 42 -inputFile "',PATHS$in.up, '" -outputFile "',PATHS$out.up.ks,'"'), 
                paste0('java -jar ${PATH_RIGER} -scoringMethod KSbyScore -randomSeed 42 -inputFile "',PATHS$in.dn, '" -outputFile "',PATHS$out.dn.ks,'"')) 
    write.table(script, PATHS$script, sep = '\t', row.names = F, col.names = F, quote = F)
}

## Riger input reformating
Format_Riger <- function(data, col.score = 'score', direction = 'pos', eff = NULL){
    ifelse(direction == 'pos', data$rank <- rank(-data[,col.score], ties.method = 'last'),data$rank <- rank(data[,col.score], ties.method = 'last'))
    data <- data.frame('Construct' = data$sgRNA, 'GeneSymbol' = data$Gene, 'NormalizedScore' = data[,col.score], 'Construct Rank' = data$rank)
    ifelse(!is.null(eff), data <- merge(data, eff, by.x = 'Construct', by.y = 'sgRNA', all.x = T, all.y = F), data$HairpinWeight <- 1)
    colnames(data)[ncol(data)] <- 'HairpinWeight'
    colnames(data)[4] <- 'Construct Rank'
    return(data)
}

## Riger output reformating
Combine_Riger_Results <- function(path.up, path.dn, path.combined){
    pos <- read.delim(path.up)
    neg <- read.delim(path.dn)
    out <- merge(pos[,c(2:4,6)], neg[,c(2:4,6)], by = 'Gene.Name')
    colnames(out) <- c('Gene','Score_pos','Pvalue_pos','Ranks_pos','Score_neg','Pvalue_neg','Ranks_neg')
    out$Ranks_pos <- str_replace_all(out$Ranks_pos, ' ', '\\,')
    out$Ranks_neg <- str_replace_all(out$Ranks_neg, ' ', '\\,')
    out$Pvalue_pos <- abs(out$Pvalue_pos)
    out$Pvalue_neg <- abs(out$Pvalue_neg)
    write.table(out, path.combined, sep = '\t', row.names = F, col.names = T, quote = F)
    return(out)
}


## z-score analysis
Get_Zscores <- function(counts, design){
    # filter undetected sgRNA  
    sg.keep <- counts[which(rowSums(counts[,3:ncol(counts)]) > 0), 'sgRNA']
    counts <- counts[which(counts$sgRNA %in% sg.keep),]
    # log-normalize
    lcpm <- cbind(counts[,1:2], apply(counts[,3:ncol(counts)], 2, function(x) { log2(1 + (x*1e6)/sum(x)) }))
    # get lfc data
    lfc = apply(design,1, function(x) {(lcpm[,x[1]] - lcpm[,x[2]])})
    colnames(lfc) <- design$Screen
    # get z-scores with NTC as center
    ntc <- apply(lfc[which(lcpm$Gene == 'NTC'),], 2, median)
    z <- scale(lfc, center = ntc)
    z <- cbind(lcpm[,1:2],z)
    # get mean z 
    s <- dlply(design, agg.column, summarise, Screen)
    mz <- z[,1:2]
    for(group in names(s)) mz <- cbind(mz, rowMeans(z[,s[[group]]$Screen]))
    colnames(mz)[3:ncol(mz)] <- names(s)
    mz$Pooled <- rowMeans(z[,3:ncol(z)])
    mz$rank <- rank(mz$Pooled)
    return(mz)
}


Run_PBNPA <- function(counts = rc, design = d, DIR_OUT = file.path(DIR,'PBNPA'), PREFIX = paste0('PBNPA_',DATASET)){
    data <- lapply(1:nrow(design), 
                   function(x){ counts[,c('sgRNA','Gene',as.character(design[2,2:1]))] })
    names(data) <- design$Screen
    res <- PBNPA(data, sim.no = 10, alpha.threshold = 0.2, fdr = 0.05)
    return(res)
}

Run_CB2 <- function(counts = rc, design = d, DIR_OUT = file.path(DIR,'Riger'), PREFIX = paste0('riger_',DATASET)){

}


