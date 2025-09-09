#!/usr/bin/env Rscript

library(plyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(edgeR)


DIR <- '1_SAMBA/Datasets_Cancer'
DATASETS <- c('Aguirre','Hart','Meyer','Wang')

for(DATASET in DATASETS){
    rc <- read.delim(file.path(DIR,paste0('rc_',DATASET,'.txt')))
    d <- read.delim(file.path(DIR,paste0('DesignPaired_',DATASET,'.txt')))
    eff <- read.delim(file.path(DIR,paste0('eff_',DATASET,'.txt')))
    m <- read.delim(file.path(DIR,paste0('metadata_',DATASET,'.txt')))
    
    ## Mageck
    Prep_Mageck(rc, d, paired = F, file.path(DIR,'Analysis_Mageck'), paste0('mageck_',DATASET)) # Aguirre
    
    ## Riger-LFC
    Prep_Riger(counts = rc, design = d, DIR_OUT = file.path(DIR,'Analysis_Riger'), PREFIX = paste0('riger_',DATASET))
    r <- Combine_Riger_Results(path.up = file.path(DIR,'Analysis_Riger',paste0('riger_',DATASET, '_out_up_ks.txt')), 
                               path.dn = file.path(DIR,'Analysis_Riger',paste0('riger_',DATASET, '_out_dn_ks.txt')), 
                               path.combined = file.path(DIR,'Analysis_Riger',paste0('riger_',DATASET, '_out_ks.txt')))
    
    ## z-score
    mz <- Get_Zscores(counts = rc, design = d)
    
    ## CB2
    res <- Run_CB2()
    
    ## PBNPA
    res <- Run_PBNPA(counts = rc, design = d, DIR_OUT = file.path(DIR,'Analysis_PBNPA'), PREFIX = paste0('pbnpa_',DATASET))
}


