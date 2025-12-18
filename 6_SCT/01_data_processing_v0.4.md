### Load libraries and prepare environment

    ############################
    ## Set working directory
    ############################

    setwd('6_SCT')

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
    # Library for SCT data integration
    library(rliger)


    ## Source file with accessory functions
    source('00_Additional_Functions.R')

    ############################
    ## Set options
    ############################

    options(future.globals.maxSize= 10000*1024^2)
    options(matrixStats.useNames.NA = 'deprecated')
    options("ggrastr.default.dpi" = 750)

### Setup Seurat Object (SO)

    ############################
    ## Create merged count matrix
    ############################

    ## List CellRanger output directories
    rawDataPaths <- dir('Datasets/CellRanger', recursive = FALSE, full.names = TRUE)

    ## Append filtered matrix subdirectory to each path
    rawDataPaths <- paste0(rawDataPaths, '/filtered_feature_bc_matrix')

    ## Read first dataset
    m <- Read10X(rawDataPaths[1])

    ## Merge second dataset into first matrix
    m <- combineMatrices(m, Read10X(rawDataPaths[2]), paste0('-', 2))


    ############################
    ## Create Seurat object
    ############################

    ## Initialize Seurat object with filtering thresholds
    so <- CreateSeuratObject(counts = m, min.cells = 3, min.features = 200)

    ## Remove raw matrix from memory
    rm(m)

    ## Compute percent mitochondrial reads per cell
    so$Percent.mt <- PercentageFeatureSet(so, pattern = '^MT-')

    ## Compute percent low-quality transcripts per cell
    so$Low_Quality <- PercentageFeatureSet(
      so,
      features = c('MALAT1', 'KCNQ1OT1')
    )

    ## Filter cells based on quality thresholds
    so <- subset(
      so,
      subset = Percent.mt < 20 & Low_Quality < 5
    )

    ############################
    ## Annotate sample identity
    ############################

    ## Define sample names
    samples <- c('OR7A10_Tumor', 'OR7A10st_Tumor')

    ## Assign sample identity based on barcode suffix
    so$sample <- samples[as.integer(str_split_i(colnames(so), '-', 2))]

### Integrate sample datasets

    ############################
    ## LIGER iNMF integration
    ############################

    ## Trigger garbage collection
    gc()

    ## Use sequential execution for integration
    plan('sequential')

    ## Normalize data by sample
    so <- rliger::normalize(so)

    ## Select variable genes across datasets
    so <- selectGenes(so, datasetVar = 'sample', thresh = 0.1)

    ## Scale data without centering
    so <- scaleNotCenter(so, datasetVar = 'sample')

    ## Run integrative non-negative matrix factorization
    so <- runINMF(
      so,
      datasetVar = 'sample',
      k = 20,
      nCores = 20,
      lambda = 10
    )

    ## Align factors across datasets
    so <- alignFactors(so, method = "centroidAlign")


    ############################
    ## Dimensional reduction
    ############################

    ## Enable multicore execution
    plan('multicore')

    ## Run UMAP on integrated embedding
    so <- RunUMAP(
      so,
      reduction = "inmfNorm",
      dims = 1:20,
      min.dist = 0.01,
      seed.use = 42
    )


    ############################
    ## Store UMAP coordinates
    ############################

    ## Save UMAP x-coordinate to metadata
    so$x <- so$umap@cell.embeddings[,1]

    ## Save UMAP y-coordinate to metadata
    so$y <- so$umap@cell.embeddings[,2]

    ## Compute nearest-neighbor graph
    so <- FindNeighbors(
      so,
      reduction = 'inmfNorm',
      dims = 1:20,
      l2.norm = TRUE,
      graph.name = c('RNA_nn','RNA_snn')
    )

    ############################
    ## Save intermediate object
    ############################

    ## Save unprocessed integrated Seurat object
    saveRDS(so, 'Datasets/so_unproc.rds')

### Find and exclude tumor cells

    ############################
    ## Visualize tumor marker expression
    ############################

    ## Generate feature plots for tumor-associated genes
    p <- FeaturePlot(
      so,
      c('EGFR','EPCAM','KRT18','KRT19'),
      pt.size = 5,
      order = TRUE,
      raster = TRUE
    )

    ## Save tumor marker UMAP plot
    ggsave(
      plot = p,
      'Figures/01_tumorRemoval_umap.pdf',
      height = 2.5 * 2,
      width = 2.9 * 2,
      scale = 1.5
    )

    ############################
    ## Remove tumor cells
    ############################

    ## Subset cells based on UMAP y-coordinate threshold
    so <- subset(so, subset = (y > -9))

### Cluster cells

    ############################
    ## Construct SNN graph
    ############################

    ## Compute nearest-neighbor graph using UMAP coordinates
    so <- FindNeighbors(
      so,
      reduction = 'umap',
      dims = 1:2,
      graph.name = c('RNA_nn','RNA_snn')
    )

    ############################
    ## Determine optimal clustering resolution
    ############################

    ## Switch to sequential execution
    plan('sequential')

    ## Evaluate clustering resolutions across wide range
    oc <- DetermineOptimalClusters(
      so,
      graph.name = 'RNA_snn',
      resolution = seq(0.05, 0.3, 0.01)
    )

    ## Display resolution diagnostic plot
    oc$plot

    ## Refine resolution search
    oc <- DetermineOptimalClusters(
      so,
      graph.name = 'RNA_snn',
      resolution = seq(0.1, 0.15, 0.002)
    )

    ## Display refined diagnostic plot
    oc$plot

    ############################
    ## Perform final clustering
    ############################

    ## Run Leiden clustering at selected resolution
    so <- FindClusters(
      so,
      graph.name = 'RNA_snn',
      algorithm = 4,
      random.seed = 42,
      resolution = 0.116,
      cluster.name = 'clust_numbered'
    )

    ############################
    ## Save clustered object
    ############################

    ## Save clustered Seurat object
    saveRDS(so, 'Datasets/so_proc1.rds')

### Classify iNK/mNK cell populations

    ############################
    ## Marker visualization
    ############################

    ## Plot NCAM1 expression
    p1 <- FeaturePlot(
      so,
      raster = TRUE,
      order = TRUE,
      raster.dpi = c(750, 750),
      max.cutoff = 'q99',
      'NCAM1',
      pt.size = 4.5
    ) +
      scale_color_viridis()

    ## Plot FCGR3A expression
    p2 <- FeaturePlot(
      so,
      raster = TRUE,
      order = TRUE,
      raster.dpi = c(750, 750),
      max.cutoff = 'q90',
      'FCGR3A',
      pt.size = 4.5
    ) +
      scale_color_viridis()

    ## Save marker feature plots
    ggsave(
      plot = p1 | p2,
      'Figures/featPlot_markers.pdf',
      height = 2.5,
      width = 2.9 * 2,
      scale = 1.5
    )


    ############################
    ## Define marker gene sets
    ############################

    ## Set default assay to RNA
    DefaultAssay(so) <- 'RNA'

    ## Load curated NK marker list from Tang_Cell2023
    subset.markerlist <- read.delim('Data/01_cellmarker_list_Tang_Cell2023.txt')

    ## Define functional NK marker groups manually
    function.markerlist <- do.call(
      rbind,
      list(
        data.frame(
          subset = 'inhib.hla.dep',
          marker = c('KIR2DL1','KIR2DL3','KIR3DL1','KIR3DL2','LILRB1','LAG3')
        ),
        data.frame(
          subset = 'inhib.hla.indep',
          marker = c('PDCD1','SIGLEC7','CD300A','CD96','IL1RAPL1','TIGIT','HAVCR2')
        ),
        data.frame(
          subset = 'activ.hla.dep',
          marker = c('KIR2DL4','CD160','KLRC2')
        ),
        data.frame(
          subset = 'activ.hla.indep',
          marker = c('NCR3','NCR1','KLRK1','CRTAM','FCGR3A','HCST')
        ),
        data.frame(
          subset = 'coreceptors',
          marker = c('CD226','SLAMF6','SLAMF7','CD244','TNFRSF9','CD59')
        ),
        data.frame(
          subset = 'cytotoxicity',
          marker = c('GZMA','GZMB','GZMH','GZMM','GZMK','GNLY','PRF1','CTSW')
        ),
        data.frame(
          subset = 'inflammatory',
          marker = c('CCL2','CCL3','CCL4','CCL5','CXCL10','CXCL9','IL1B','IL6',
                     'IL7','IL15','IL18')
        ),
        data.frame(
          subset = 'stress',
          marker = c('BAG3','CALU','DNAJB1','DUSP1','EGR1','FOS','FOSB','HIF1A',
                     'HSP90AA1','HSP90AB1','HSP90B1','HSPA1A','HSPA1B','HSPA6',
                     'HSPB1','HSPH1','IER2','JUN','JUNB','NFKBIA','NFKBIZ',
                     'RGS2','SLC2A3','SOCS3','UBC','ZFAND2A','ZFP3','ZFP36L1')
        )
      )
    )

    ## Combine curated and functional marker lists
    mlist <- rbind(subset.markerlist, function.markerlist)

    ## Convert marker list into named list format
    mlist <- lapply(
      unique(mlist$subset),
      function(x){
        return(mlist$marker[which(mlist$subset == x)])
      }
    ) %>% structure(names = unique(mlist$subset))

    ## Replace dashes with underscores in marker set names
    names(mlist) <- str_replace(names(mlist), '-', '_')


    ############################
    ## Compute AUCell signature scores
    ############################

    ## Identify highly variable genes overlapping marker lists
    hvg <- intersect(VariableFeatures(so$RNA), unique(unlist(mlist)))

    ## Extract normalized expression matrix
    m <- so$RNA$data[hvg,]

    ## Load AUCell package
    library(AUCell)

    ## Build gene expression rankings
    m2 <- AUCell_buildRankings(
      exprMat = m,
      featureType = "genes",
      splitByBlocks = TRUE,
      BPPARAM = BiocParallel::MulticoreParam(
        tasks = 100,
        force.GC = TRUE,
        progressbar = TRUE
      ),
      verbose = TRUE
    )

    ## Calculate AUC scores for marker gene sets
    m <- AUCell_calcAUC(
      mlist,
      m2,
      nCores = 20,
      verbose = TRUE
    )

    ## Save AUCell results
    saveRDS(m, 'Datasets/aucell_subsetMarkerSignatures.rds')

    ## Remove intermediate ranking object
    rm(m2)


    ############################
    ## Add signature scores to Seurat metadata
    ############################

    ## Transpose AUC matrix to cell-by-signature format
    m <- t(m@assays@data$AUC[, colnames(so)])

    ## Append signature scores to metadata
    so@meta.data <- cbind(so@meta.data, m)

### Make dot plots of markers and signatures to characterize clusters

    ############################
    ## Generate cluster-level dot plots
    ############################

    ## Assemble list of dot plots
    p <- list(

      ## Plot general NK markers
      customDotPlot(
        so,
        c('NCAM1','FCGR3A','KLRF1','RGS1'),
        group.by = 'clust_numbered',
        scale.max = 50
      ),

      ## Plot major NK functional programs
      customDotPlot(
        so,
        grep('^CD56', names(mlist), invert = TRUE, value = TRUE),
        group.by = 'clust_numbered',
        scale.max = 50
      ),

      ## Plot subtype-discriminating genes
      customDotPlot(
        so,
        c('KLRC2','IL32','CX3CR1','BATF3','MKI67','EGR2','PRF1','GZMB','CREM'),
        group.by = 'clust_numbered',
        scale.max = 50
      )
    )

    ## Combine dot plots into single figure
    p <- cowplot::plot_grid(
      plotlist = p,
      align = 'h',
      axis = 'lbt',
      nrow = 1,
      rel_widths = c(0.4, 1, 0.6, 0.7)
    )

    ## Save dot plot figure
    ggsave(
      plot = p,
      'Figures/dotplot_celltypes_clusters.pdf',
      height = 3.7,
      width = 1.7 * 4,
      scale = 1.5
    )

### Name clusters

    ############################
    ## Assign biologically meaningful cluster names
    ############################

    ## Map numeric clusters to annotated NK subsets
    so$ct <- c(
      'mNK_c5-KLRC2',
      'mNK_c3-BATF3',
      'mNK_c1-MKI67',
      'mNK_c2-MKI67-BATF3',
      'mNK_c6-EGR2',
      'mNK_c4-IL32',
      'iNK_c1-Cytotoxic',
      'iNK_c2-CREM'
    )[as.integer(so$clust_numbered)]

    ## Derive NK type (iNK vs mNK)
    so$nk_type <- str_split_i(so$ct, '_', 1)

    ## Set Seurat identities to cell type labels
    Idents(so) <- so$ct


    ############################
    ## Save annotated object
    ############################

    ## Save processed Seurat object
    saveRDS(so, 'Datasets/so_proc2.rds')

### Visualize marker genes and signatures across cell subsets

    ############################
    ## Dot plots across cell subsets
    ############################

    ## Define ordering of cell types for plotting
    group.order <- rev(sort(unique(so$ct)))

    ## Assemble list of dot plots across cell types
    p <- list(

      ## Plot general NK markers and signatures
      customDotPlot(
        so,
        c(
          'NCAM1','FCGR3A','KLRF1','RGS1',
          grep('^CD56', names(mlist), invert = TRUE, value = TRUE)
        ),
        group.order = group.order,
        group.by = 'ct',
        scale.max = 50
      ),

      ## Plot iNK-associated markers
      customDotPlot(
        so,
        c('cytotoxicity', 'PRF1','CREM','stress'),
        group.order = group.order,
        group.by = 'ct',
        scale.max = 50
      ),

      ## Plot mNK-associated markers
      customDotPlot(
        so,
        c('MKI67', 'BATF3', 'IL32','KLRC2','EGR2'),
        group.order = group.order,
        group.by = 'ct',
        scale.max = 50
      )
    )

    ## Combine dot plots into a single figure
    p <- cowplot::plot_grid(
      plotlist = p,
      align = 'h',
      axis = 'lbt',
      nrow = 1,
      rel_widths = c(1, 0.7, 0.7)
    )

    ## Save dot plot figure by cell subset
    ggsave(
      plot = p,
      'Figures/dotplot_celltypes_ct.pdf',
      height = 2.75,
      width = 2.5 * 3,
      scale = 1.5
    )


    ############################
    ## Dot plots across subsets and genotypes
    ############################

    ## Create combined subset-genotype identifier
    so$id <- paste0(
      so$ct,
      ':',
      (
        str_replace(so$sample, 'st', 'stop') %>%
          str_remove('_Tumor')
      )
    )

    ## Assemble list of dot plots across subset-genotype groups
    p <- list(

      ## Plot general NK markers
      customDotPlot(
        so,
        c(
          'NCAM1','FCGR3A','KLRF1','RGS1',
          grep('^CD56', names(mlist), invert = TRUE, value = TRUE)
        ),
        group.order = rev(sort(unique(so$id))),
        group.by = 'id',
        scale.max = 50
      ),

      ## Plot iNK-associated markers
      customDotPlot(
        so,
        c('cytotoxicity', 'PRF1','CREM','stress'),
        group.order = rev(sort(unique(so$id))),
        group.by = 'id',
        scale.max = 50
      ),

      ## Plot mNK-associated markers
      customDotPlot(
        so,
        c('MKI67', 'BATF3', 'IL32','KLRC2','EGR2'),
        group.order = rev(sort(unique(so$id))),
        group.by = 'id',
        scale.max = 50
      )
    )

    ## Combine subset-genotype dot plots
    p <- cowplot::plot_grid(
      plotlist = p,
      align = 'h',
      axis = 'lbt',
      nrow = 1,
      rel_widths = c(1, 0.75, 0.75)
    )

    ## Save dot plot figure by subset and genotype
    ggsave(
      plot = p,
      'Figures/dotplot_celltypes_ct-geno.pdf',
      height = 4,
      width = 3.1 * 3,
      scale = 1.5
    )


    ############################
    ## Dot plots of signature genes
    ############################

    ## Generate dot plot of all signature genes across cell types
    p <- customDotPlot(
      so,
      as.character(
        unlist(
          mlist[
            grep('^CD56', names(mlist), invert = TRUE, value = TRUE)
          ]
        )
      ),
      group.order = group.order,
      group.by = 'ct',
      scale.max = 20
    )

    ## Save dot plot of all signature genes by cell type
    ggsave(
      plot = p,
      'Figures/dotplot_celltypes_allGenes.pdf',
      height = 2.75,
      width = 12,
      scale = 1.5
    )

    ## Generate dot plot of all signature genes by subset-genotype
    p <- customDotPlot(
      so,
      as.character(
        unlist(
          mlist[
            grep('^CD56', names(mlist), invert = TRUE, value = TRUE)
          ]
        )
      ),
      group.order = rev(sort(unique(so$id))),
      group.by = 'id',
      scale.max = 20
    )

    ## Save dot plot of all signature genes by subset-genotype
    ggsave(
      plot = p,
      'Figures/dotplot_celltypes-geno_allGenes.pdf',
      height = 3.5,
      width = 13,
      scale = 1.5
    )

### Make subset-labeled UMAPs

    ############################
    ## Compute label coordinates for cluster annotations
    ############################

    ## Compute centroid coordinates for each cell type label
    centroids <- reframe(so@meta.data, .by = c('ct'), label = ct[1], x = mean(x), y = mean(y))

    ############################
    ## Build full-dataset UMAP with subset labels
    ############################

    ## Create UMAP scatter plot for all cells with cell-type coloring
    p1 <- ggplot(so@meta.data, aes(x = x, y = y)) +
        ggrastr::geom_point_rast(size = 2, color = 'gray25') +
        ggrastr::geom_point_rast(size = 0.75, aes(color = ct)) +
        theme_void() +
        ggrepel::geom_text_repel(
            data = centroids,
            aes(x = x, y = y, label = label),
            min.segment.length = 0
        ) +
        theme(legend.position = 'none')

    ############################
    ## Define consistent color palette for cell types
    ############################

    ## Define colors for each cell type
    cols <- structure(ggcolor(8), names = sort(unique(so$ct)))

    ############################
    ## Build per-sample UMAPs with subset labels
    ############################

    ## Generate per-sample UMAP plots with consistent coloring
    p2 <- lapply(unique(so$sample), function(samples){

        ## Subset metadata to current sample
        d <- so@meta.data[which(so$sample == samples),]

        ## Compute centroid coordinates for current sample
        centroids <- reframe(d, .by = c('ct'), label = ct[1], x = mean(x), y = mean(y))

        ## Build UMAP plot for current sample
        ggplot(d, aes(x = x, y = y)) +
            theme_void() +
            ggrastr::geom_point_rast(data = so@meta.data, size = 1, color = 'gray25') +
            ggrastr::geom_point_rast(data = so@meta.data, size = 0.25, color = 'gray98') +
            ggrastr::geom_point_rast(size = 1, color = 'gray25') +
            ggrastr::geom_point_rast(size = 0.25, aes(color = ct)) +
            labs(title = samples) +
            scale_color_manual(values = cols) +
            theme(legend.position = 'none', title = element_text(vjust = 0.25))
    })

    ############################
    ## Arrange panels and save figure
    ############################

    ## Stack per-sample plots vertically
    p2 <- cowplot::plot_grid(plotlist = p2, align = 'hv', axis = 'lb', ncol = 1)

    ## Combine full-dataset and per-sample panels
    p <- cowplot::plot_grid(p1, p2, align = 'hv', axis = 'lb', ncol = 2, rel_widths = c(1, 0.6))

    ## Save combined UMAP figure
    ggsave(plot = p, 'Figures/umap_celltypes_samples.pdf', height = 3 * 1, width = 2.5 * 2)

### Compare distributions of cell subsets

    ############################
    ## Compute subset proportions for stacked barplot
    ############################

    ## Extract relevant metadata columns
    d <- so@meta.data[, c('ct', 'nk_type', 'sample')]

    ## Compute per-sample totals
    d <- mutate(d, .by = c('sample'), total = length(sample))

    ## Compute per-celltype proportions per sample
    d <- reframe(d, .by = c('ct', 'sample'), pct = length(ct) / total[1])

    ## Simplify sample labels for plotting
    d$sample <- str_split_i(d$sample, '_', 1) %>% str_replace(., 'st', 'stop')

    ## Set factor ordering for sample labels
    d$sample <- factor(d$sample, levels = sort(unique(d$sample)))

    ## Define ordered cell type list
    celltypes <- sort(unique(so$ct))

    ## Set factor ordering for cell types
    d$ct <- factor(d$ct, levels = (celltypes))

    ## Define colors for cell types
    cols <- structure(ggcolor(8), names = levels(d$ct))

    ############################
    ## Make stacked bar plot across samples
    ############################

    ## Build stacked barplot of subset proportions by sample
    p1 <- ggplot(d, aes(y = sample, x = pct, fill = ct)) +
        geom_bar(stat = 'identity', position = 'stack', color = 'gray20') +
        theme_classic() +
        scale_fill_manual(values = cols) +
        labs(y = NULL, x = ('Proportion of NK cells')) +
        theme(
            legend.key.spacing = unit(0.5, 'lines'),
            legend.key.size = unit(0.9, 'lines'),
            legend.key.height = unit(0.5, 'lines'),
            legend.key.width = unit(0.5, 'lines'),
            legend.text = element_text(size = 7),
            legend.key.spacing.y = unit(0.01, 'lines'),
            legend.title = element_text(size = 10)
        )

    ############################
    ## Compute subset proportions for grouped barplot
    ############################

    ## Extract relevant metadata columns
    d <- so@meta.data[, c('ct', 'nk_type', 'sample')]

    ## Compute per-sample totals
    d <- mutate(d, .by = c('sample'), total = length(sample))

    ## Compute per-celltype proportions per sample
    d <- reframe(d, .by = c('ct', 'sample'), pct = length(ct) / total[1])

    ## Set factor ordering for cell types
    d$ct <- factor(d$ct, levels = sort(celltypes))

    ## Set factor ordering for samples
    d$sample <- factor(d$sample, levels = rev(samples))

    ## Define colors for sample groups
    cols <- structure(c('#B2B2B2', '#E9AC4C'), names = rev(samples))

    ## Create wrapped x-axis labels for cell types
    xlabs <- paste0(substr(celltypes, 0, 6), '\n', substr(celltypes, 8, 100))

    ## Assign names to x-axis label mapping
    names(xlabs) <- celltypes

    ############################
    ## Make grouped bar plot across cell types
    ############################

    ## Build grouped barplot of subset proportions by genotype
    p2 <- ggplot(d, aes(x = ct, y = pct, fill = sample)) +
        geom_bar(stat = 'identity', position = 'dodge', color = 'gray20') +
        theme_classic() +
        scale_fill_manual(values = cols) +
        scale_x_discrete(labels = xlabs) +
        labs(x = NULL, y = str_wrap('Proportion of NK cells')) +
        theme(
            legend.key.spacing = unit(0.5, 'lines'),
            legend.key.size = unit(0.9, 'lines'),
            axis.text.x = element_text(lineheight = 0.75),
            legend.key.height = unit(0.5, 'lines'),
            legend.key.width = unit(0.5, 'lines'),
            legend.text = element_text(size = 7),
            legend.key.spacing.y = unit(0.01, 'lines'),
            legend.title = element_text(size = 10)
        )

    ############################
    ## Combine barplots and save
    ############################

    ## Combine stacked and grouped barplots vertically
    p <- cowplot::plot_grid(p1, p2, ncol = 1, align = 'hv', axis = 'lbt', rel_heights = c(0.75, 1))

    ## Save barplot figure
    ggsave(plot = p, 'Figures/barplot_subset_pct.pdf', height = 1.2 * 2, width = 8)


    ############################
    ## Statistical tests for subset proportions between genotypes
    ############################

    ## Load rstatix for pairwise tests
    library(rstatix)

    ## Extract subset and sample metadata
    d <- so@meta.data[, c('ct', 'sample')]

    ## Count cells per subset and sample
    d <- reframe(d, .by = c('ct', 'sample'), x = length(ct))

    ## Compute total cells per sample
    d <- mutate(d, .by = 'sample', total = sum(x))

    ## Compute complement counts for contingency testing
    d$y <- d$total - d$x

    ## Compute proportions for reference
    d$p <- d$x / d$total

    ## Set rownames for convenience
    rownames(d) <- paste0(d$sample, '_', d$ct)

    ## Run pairwise Fisher tests for each cell type
    r <- do.call(rbind, lapply(celltypes, function(ct){

        ## Skip excluded cell type
        if(ct == 'CD4_TCM') return(NULL)

        ## Compute pairwise Fisher tests for current cell type
        r <- pairwise_fisher_test(d[which(d$ct == ct), c('x', 'y')])

        ## Return annotated result table
        return(data.frame(ct = ct, r[,1:4]))
    }))

    ## Adjust p-values across tests
    r$p.adj <- p.adjust(r$p)

    ## Reorder rows and drop redundant column
    r <- r[order(r$ct), -4]

    ## Assign significance stars based on adjusted p-values
    r$signif <- case_when(
        r$p.adj < 1e-4 ~ '****',
        r$p.adj < 1e-3 ~ '***',
        r$p.adj < 1e-2 ~ '**',
        r$p.adj < 0.05 ~ '*',
        r$p.adj >= 0.05 ~ 'ns'
    )

    ## Save statistical test results
    write.table(r, 'Data/subset_pct_comp.txt', sep = '\t', quote = FALSE, row.names = FALSE)

### Cell cycle analysis

    ############################
    ## Cell cycle scoring
    ############################

    ## Compute cell cycle scores without changing identities
    so <- CellCycleScoring(
        so,
        s.features = cc.genes$s.genes,
        g2m.features = cc.genes$g2m.genes,
        set.ident = FALSE
    )

    ############################
    ## Visualize cell cycle phase on UMAP
    ############################

    ## Compute centroid coordinates for each cell cycle phase
    centroids <- reframe(
        so@meta.data,
        .by = 'Phase',
        label = Phase[1],
        x = mean(x),
        y = mean(y)
    )

    ## Build UMAP plot colored by cell cycle phase
    p <- ggplot(so@meta.data, aes(x = x, y = y)) +
        ggrastr::geom_point_rast(size = 2, color = 'gray25') +
        ggrastr::geom_point_rast(size = 0.75, aes(color = Phase)) +
        theme_void() +
        ggrepel::geom_text_repel(
            data = centroids,
            aes(x = x, y = y, label = label),
            min.segment.length = 0
        ) +
        theme(legend.position = 'bottom')

    ## Save cell cycle UMAP plot
    ggsave(plot = p, 'Figures/umap_cellcycle.pdf', height = 2.5 * 1, width = 2.5 * 1)
