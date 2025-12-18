### Load libraries and prepare environment

    ############################
    ## Set working directory
    ############################

    setwd('3_ORF_Miniscreen')

    ############################
    ## Load required libraries
    ############################

    library(dplyr)
    library(stringr)
    library(ggplot2)
    library(ggrastr)
    library(stringr)

    ## Source samba scripts for version 1.1
    source('Scripts/Samba_official_V1.1.R')

### Preprocessing

    ############################
    ## Sample ID conversion
    ############################

    ## Read sample barcode mapping table
    map.id <- read.delim('Data/SampleBarcode_Map.txt')

    ## Modify sample barcode strings to match expected format
    map.id[,2] <- paste0(
      substr(map.id[,2], 0, 10),
      '0',
      substr(map.id[,2], 11, 100)
    )

    ## Create named vector mapping filenames to sample IDs
    map.id <- structure(
      map.id[,1],
      names = paste0('Sample_', map.id[,2], '_readalign.txt')
    )

    ############################
    ## Reorganize raw count files into long-format dataframe
    ############################

    ## List all raw UMI count files
    filenames <- list.files('RawCounts/umi')

    ## Keep only files present in the sample ID mapping
    filenames <- intersect(names(map.id), filenames)

    ## Read and aggregate counts across all samples in parallel
    d <- do.call(
      rbind,
      BiocParallel::bplapply(
        filenames,
        function(filename){

          ## Retrieve sample ID corresponding to current file
          sample_id <- map.id[filename]

          ## Read raw sgRNA count file
          ct <- read.delim(
            file.path('RawCounts/umi', filename),
            header = FALSE,
            strip.white = TRUE
          )

          ## Extract sgRNA base name
          ct[,1] <- str_split_i(ct[,1], '\\.', 1)

          ## Construct full sgRNA identifier
          ct[,2] <- paste0(ct[,1], '_', ct[,2])

          ## Aggregate UMI counts per sgRNA
          ct <- reframe(
            ct,
            .by = 'V2',
            V3 = sum(V3)
          )

          ## Assemble tidy count dataframe for current sample
          ct <- data.frame(
            id    = sample_id,
            sgRNA = ct$V2,
            ct    = as.integer(ct$V3),
            row.names = NULL
          )

          ## Return per-sample count dataframe
          return(ct)

        },
        ## Configure parallel execution parameters
        BPPARAM = BiocParallel::MulticoreParam(
          workers = 4,
          progressbar = TRUE
        )
      )
    )

    ############################
    ## Cast long-format counts into matrix-style table
    ############################

    ## Convert long-format counts into wide matrix (sgRNA x sample)
    d <- reshape2::dcast(
      d,
      sgRNA ~ id,
      value.var = 'ct',
      fill = 0
    )

    ## Add gene annotation column derived from sgRNA IDs
    d <- data.frame(
      sgRNA = d$sgRNA,
      Gene  = str_split_i(d$sgRNA, '_', 1),
      d[,-1]
    )

    ## Write final count matrix to file
    write.table(
      d,
      'Data/data_counts.txt',
      sep = '\t',
      quote = FALSE,
      row.names = FALSE
    )
