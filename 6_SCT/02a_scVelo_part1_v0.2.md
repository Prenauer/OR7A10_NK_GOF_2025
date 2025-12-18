### Load libraries and prepare environment

    ############################
    ## Set working directory
    ############################

    setwd('06_SCT')

    ############################
    ## Load required libraries
    ############################

    # Common libraries
    library(dplyr)
    library(stringr)
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(ggrastr)
    library(Matrix)
    library(reticulate)
    library(Seurat)
    library(future)
    # Libraries to convert SO to loom format
    library(SeuratDisk)
    library(SeuratWrappers)


    ## Source file with accessory functions
    source('00_Additional_Functions.R')

    ############################
    ## Set options
    ############################

    options(future.globals.maxSize= 10000*1024^2)
    options(matrixStats.useNames.NA = 'deprecated')
    options("ggrastr.default.dpi" = 750)

    ############################
    ## Load data
    ############################

    so <- readRDS('Datasets/so_proc2.rds')

### Export data for scvelo

    ############################
    ## Rename cells to match velocyto conventions
    ############################

    ## Define mapping between sample names and velocyto IDs
    id <- c(OR7A10st_Tumor = 'U947Q', OR7A10_Tumor = 'QSQ79')

    ## Construct new cell names using velocyto-compatible format
    new.names <- paste0(
      'possorted_genome_bam_',
      id[so$sample],
      ':',
      str_split_i(colnames(so), '-', 1),
      'x'
    )

    ## Rename cells in Seurat object
    so2 <- RenameCells(so, new.names = new.names)


    ############################
    ## Export data for RNA velocity analysis
    ############################

    ## Save Seurat object as loom file for scVelo
    SaveLoom(
      so2,
      'Datasets/so_proc1.loom',
      overwrite = TRUE
    )

    ## Write UMAP cell embeddings to file
    write.table(
      so2$umap@cell.embeddings,
      'Datasets/umap_embedding.dat',
      sep = '\t',
      quote = FALSE,
      row.names = FALSE
    )

    ## Write iNMF-normalized embeddings to file
    write.table(
      so2$inmfNorm@cell.embeddings,
      'Datasets/nmf_embedding.dat',
      sep = '\t',
      quote = FALSE,
      row.names = FALSE
    )
