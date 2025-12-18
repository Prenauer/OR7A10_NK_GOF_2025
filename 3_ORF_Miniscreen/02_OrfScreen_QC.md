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
    source('Samba_official_V1.1.R')

### PCA to assess between-sample variance

    ############################
    ## Load data and perform PCA
    ############################

    ## Read sgRNA count matrix
    d <- read.delim('Data/data_counts.txt')

    ## Apply variance-stabilizing transformation to count data
    d2 <- DESeq2::varianceStabilizingTransformation(
      as.matrix(d[,3:ncol(d)]),
      fitType = 'local'
    )

    ## Perform PCA on transposed transformed data
    pca <- prcomp(t(d2), scale = TRUE)

    ## Assemble PCA results with sample metadata
    pca <- data.frame(
      pca$x[,1:2],
      id    = rownames(pca$x),
      donor = str_split_i(rownames(pca$x), '_', 2),
      time  = str_split_i(rownames(pca$x), '_', 3)
    )


    ############################
    ## Generate PCA plots
    ############################

    ## Define color palette for grouping variables
    cols <- c('#B2B2B2','#E9AC4C','#496AB4','#D85446') %>% rev()

    ## Create PCA plot colored by donor
    p1 <- ggplot(pca, aes(x = PC1, y = PC2, fill = donor), color = 'gray20') + 
      geom_point(size = 2.5, shape = 21, stroke = 0.01) +
      theme_classic() +
      xlim(1.3 * range(pca$PC1)) +
      theme(legend.position = 'bottom') +
      scale_fill_manual(values = cols) +
      scale_color_manual(values = cols) +
      ggrepel::geom_text_repel(
        aes(label = id, color = donor),
        size = 2.5,
        max.overlaps = 100,
        min.segment.length = 0,
        segment.size = 0.2
      )

    ## Create PCA plot colored by time point
    p2 <- ggplot(pca, aes(x = PC1, y = PC2, fill = time), color = 'gray20') + 
      geom_point(size = 2.5, shape = 21, stroke = 0.01) +
      theme_classic() +
      xlim(1.3 * range(pca$PC1)) +
      theme(legend.position = 'bottom') +
      scale_fill_manual(values = cols) +
      scale_color_manual(values = cols) +
      ggrepel::geom_text_repel(
        aes(label = id, color = time),
        size = 2.5,
        max.overlaps = 100,
        min.segment.length = 0,
        segment.size = 0.2
      )

    ## Combine PCA plots into a single figure
    p <- cowplot::plot_grid(
      p1,
      p2,
      align = 'hv',
      axis = 'lbt',
      nrow = 1
    )

    ## Save PCA figure to file
    ggsave(
      plot = p,
      'Figures/qc_pca_samples.pdf',
      height = 2.5,
      width = 5,
      scale = 1.25
    )

    ## Save PCA coordinates and metadata
    write.table(
      pca,
      'Data/data_pca.txt',
      sep = '\t',
      quote = FALSE,
      row.names = FALSE
    )

### Correlation heatmap

    ############################
    ## Compute correlation matrix
    ############################

    ## Read sgRNA count matrix
    d <- read.delim('Data/data_counts.txt')

    ## Compute log-normalized CPM values
    d2 <- log2(
      1 + apply(d[,3:ncol(d)], 2, function(x) 1e6 * x / sum(x))
    )

    ## Compute Spearman correlation matrix
    corr <- cor(d2, method = 's')

    ## Write correlation matrix to file
    write.table(
      corr,
      'Data/qc_spearmanCorrelation.pdf',
      sep = '\t',
      row.names = TRUE,
      quote = FALSE
    )

    ############################
    ## Generate correlation heatmap
    ############################

    ## Create row/column annotation dataframe
    anno_row <- data.frame(
      Time = as.integer(
        substr(
          str_split_i(colnames(corr), '_', 3),
          2,
          2
        )
      ),
      row.names = colnames(corr)
    )

    ## Plot and save correlation heatmap
    pheatmap::pheatmap(
      corr,
      scale = 'none',
      treeheight_row = 12,
      cluster_cols = TRUE,
      treeheight_col = 12,
      clustering_method = 'ward.D2',
      cutree_cols = 3,
      cutree_rows = 3,
      annotation_row = anno_row,
      annotation_col = anno_row,
      angle_col = 90,
      width = 5.75,
      height = 5,
      filename = 'Figures/qc_heatmap_spearmanCorrelation.pdf'
    )

### Box-Whisker of sample read-count distributions

    ############################
    ## Prepare normalized count data
    ############################

    ## Read sgRNA count matrix
    d <- read.delim('Data/data_counts.txt')

    ## Compute log-normalized CPM values
    d2 <- log2(
      1 + apply(d[,3:ncol(d)], 2, function(x) 1e6 * x / sum(x))
    )

    ## Convert matrix to long format
    d2 <- reshape2::melt(d2)

    ## Rename melted dataframe columns
    colnames(d2) <- c('gene', 'sample', 'value')

    ## Extract time point from sample name
    d2$time <- str_split_i(as.character(d2$sample), '_', 3)

    ## Assign control label to control samples
    d2$time[grepl('ctrl', as.character(d2$sample))] <- 'ctrl'

    ## Convert time variable to ordered factor
    d2$time <- factor(d2$time, levels = c('ctrl','t0','t3','t6'))

    ## Extract sample type from sample name
    d2$type <- str_split_i(as.character(d2$sample), '_', 1) %>% factor()

    ############################
    ## Order samples for plotting
    ############################

    ## Determine sample ordering based on time, type, and abundance
    o <- unique(d2$sample[with(d2, order(time, type, -value))])

    ## Relevel sample factor for plotting
    d2$sample <- factor(d2$sample, levels = o)

    ## Define color palette
    cols <- c('#B2B2B2','#E9AC4C','#496AB4')

    ############################
    ## Generate box-and-jitter plots
    ############################

    ## Create vertical boxplot with jittered points
    p <- ggplot(d2, aes(x = sample, y = value, color = time, fill = time)) + 
      ggrastr::geom_jitter_rast(size = 0.1, pch = 19, alpha = 0.3) + 
      scale_color_manual(values = c('#B2B2B2',"#FCF171F1","#F48841FF","firebrick")) +
      scale_fill_manual(values  = c('#B2B2B2',"#FCF171F1","#F48841FF","firebrick")) +
      geom_boxplot(color = 'gray20', outlier.size = 0.2, outlier.alpha = 0.7) +
      coord_flip() +
      labs(x = '', y = 'UMI counts (log-norm)', title = 'Library distribution') + 
      cowplot::theme_cowplot()

    ## Save vertical boxplot
    ggsave(
      plot = p,
      filename = file.path('Figures/qc_umi-distrib.pdf'),
      height = 3.75,
      width = 3.75,
      scale = 1.5
    )

    ## Create horizontal boxplot with rotated labels
    p <- ggplot(d2, aes(x = sample, y = value, color = time, fill = time)) + 
      ggrastr::geom_jitter_rast(size = 0.1, pch = 19, alpha = 0.3) + 
      scale_color_manual(values = c('#B2B2B2',"#FCF171F1","#F48841FF","firebrick")) +
      scale_fill_manual(values  = c('#B2B2B2',"#FCF171F1","#F48841FF","firebrick")) +
      geom_boxplot(color = 'gray20', outlier.size = 0.2, outlier.alpha = 0.7) +
      labs(x = '', y = 'UMI counts (log-norm)', title = 'Library distribution') + 
      cowplot::theme_cowplot() +
      theme(axis.text.x = element_text(angle = 90))

    ## Save horizontal boxplot
    ggsave(
      plot = p,
      filename = file.path('Figures/qc_umi-distrib_horiz.pdf'),
      height = 3.75,
      width = 3.75,
      scale = 1.5
    )

### Cumulative prob functions for all groups

    ############################
    ## Compute ECDFs for each group
    ############################

    ## Load patchwork for plot composition
    library(patchwork)

    ## Define color palette for time groups
    cols <- structure(
      c('#B2B2B2',"#FCF171F1","#F48841FF","firebrick"),
      names = c('ctrl','t0','t3','t6')
    )

    ## Compute ECDFs for each time group
    cdf <- lapply(
      c('ctrl','t0','t3','t6'),
      function(group) {
        ecdf(d2[which(d2$time == group), 'value'])
      }
    ) %>% structure(., names = c('ctrl','t0','t3','t6'))

    ############################
    ## Plot and save ECDF figure
    ############################

    ## Create ECDF plot across groups
    p <- ggplot() +
      xlim(range(d2$value)) +
      geom_function(aes(color = 'ctrl'), fun = cdf[[1]], linewidth = 1.2) +
      geom_function(aes(color = 't0'),   fun = cdf[[2]], linewidth = 1.2) +
      geom_function(aes(color = 't3'),   fun = cdf[[3]], linewidth = 1.2) +
      geom_function(aes(color = 't6'),   fun = cdf[[4]], linewidth = 1.2) +
      scale_color_manual(values = cols) +
      theme_classic() +
      theme(legend.position = 'inside') +
      labs(
        y = 'Cumulative probability',
        x = 'Gene counts (log-norm)',
        title = 'Group ECDF'
      )

    ## Save ECDF plot
    ggsave(
      plot = p,
      filename = file.path('Figures/qc_cdf_samples.pdf'),
      height = 3.5,
      width = 2.25
    )

    ############################
    ## Wasserstein tests between groups
    ############################

    ## Load two-sample testing library
    library(twosamples)

    ## Compute pairwise Wasserstein tests
    ws <- do.call(
      rbind,
      lapply(1:3, function(i) {
        do.call(
          rbind,
          lapply((i+1):4, function(j) {

            ## Identify group names
            t1 <- levels(d2$time)[i]
            t2 <- levels(d2$time)[j]

            ## Set seed for reproducibility
            set.seed(1)

            ## Run Monte Carlo Wasserstein test
            r <- wass_test(
              a = d2$value[which(d2$time == t1)],
              b = d2$value[which(d2$time == t2)],
              nboots = 2000,
              p = 1e-6
            )

            ## Assemble test results
            data.frame(
              Group1 = t1,
              Group2 = t2,
              Statistic = r[1],
              pval = r[2],
              method = 'Monte-Carlo two-sample Wasserstein test'
            )
          })
        )
      })
    )

    ## Remove self-comparisons
    ws <- ws[which(ws$Group1 != ws$Group2),]

    ## Adjust p-values for multiple testing
    ws$padj <- p.adjust(ws$pval)

    ## Save Wasserstein test results
    write.table(
      ws,
      'Data/qc_wass_distribution.txt',
      sep = '\t',
      quote = FALSE,
      row.names = FALSE
    )
