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
    # Libraries for signature analyses
    library(AUCell)
    library(fgsea)
    # libraries for network plots
    library(ggraph)
    library(igraph)
    library(ggnetwork)
    # libraries for upset plots
    library(ggdendroplot)
    library(ggvenn)
    library(ggVennDiagram)
    # library for GAMs
    library(mgcv)


    ## Source file with accessory functions
    source('00_Additional_Functions.R')

    ############################
    ## Set options
    ############################

    options(future.globals.maxSize= 10000*1024^2)
    options(matrixStats.useNames.NA = 'deprecated')
    options("ggrastr.default.dpi" = 750)


    ############################
    ## Load and prepare data
    ############################

    so <- readRDS('Datasets/so_proc2.rds')
    m <- readRDS('Datasets/aucell_PidSignatures.rds')
    m <- (m@assays@data$AUC[, colnames(so)])
    so$pw <- CreateAssay5Object(count=m, data=m)

### DSR Analyses

    ############################
    ## Select pathway predictors based on DE results
    ############################

    ## Load differential pathway results
    de <- read.delim("Data/de_pw_pid.txt")

    ## Select top upregulated pathways per cell type
    up <- slice_head(
        de[which(de$logFC > 0.2 & de$padj < 0.01),],
        by = 'celltype',
        n = 2
    )$feature %>% unique()

    ## Select top downregulated pathways per cell type
    dn <- slice_head(
        de[which(de$logFC < -0.2 & de$padj < 0.01),],
        by = 'celltype',
        n = 2
    )$feature %>% unique()

    ## Define commonly regulated upregulated pathways
    pws.common.up <- c(
        'PID-ATF2-PATHWAY','PID-CDC42-REG-PATHWAY','PID-HNF3A-PATHWAY',
        'PID-IL1-PATHWAY','PID-IL23-PATHWAY','PID-IL3-PATHWAY',
        'PID-NFAT-TFPATHWAY','PID-RB-1PATHWAY','PID-TCR-CALCIUM-PATHWAY',
        'PID-TCR-JNK-PATHWAY'
    )

    ## Define commonly regulated downregulated pathways
    pws.common.dn <- c(
        'PID-CXCR3-PATHWAY','PID-ERBB1-RECEPTOR-PROXIMAL-PATHWAY',
        'PID-GLYPICAN-1PATHWAY','PID-INTEGRIN-CS-PATHWAY',
        'PID-LYMPH-ANGIOGENESIS-PATHWAY'
    )

    ## Combine common pathways
    pws.common <- c(pws.common.up, pws.common.dn)

    ## Combine top DE pathways
    pws.top <- c(up, dn)

    ## Define final pathway predictor set
    pws <- c(pws.common, pws.top) %>% unique()


    ############################
    ## Select DE genes for DSR response variables
    ############################

    ## Load differential gene expression results
    de <- read.delim('Data/de_genes.txt')

    ## Filter genes by effect size and significance
    deg <- de[which(abs(de$logFC) > 0.5 & de$padj < 0.01), 'feature'] %>%
        unique()


    ############################
    ## Run DSR using pathways as predictors and genes as responses
    ############################

    ## Initialize result table
    r <- data.frame()

    ## Iterate over pathway predictors
    for (x in pws) {

        ## Run DSR in parallel over DE genes
        tmp <- bplapply(
            deg,
            function(y){

                ## Prepare DSR input data
                d <- setupDsrData(x, y)

                ## Fit GAM-based DSR model
                fit <- dsrFitGam(d, select = FALSE)

                ## Skip failed fits
                if (is.null(fit)) return(NULL)

                ## Extract DSR statistics
                return(extractDsrResults(fit))
            },
            BPPARAM = MulticoreParam(workers = 18, force.GC = FALSE,
                                     progressbar = TRUE, RNGseed = 1)
        )

        ## Remove NULL results
        tmp <- tmp[!is.null(tmp)]

        ## Combine results into a data frame
        tmp <- do.call(rbind, tmp)

        ## Add predictor and response identifiers
        tmp <- data.frame(predictor = x, response = deg, tmp)

        ## Trigger garbage collection
        gc()

        ## Print progress
        cat('\n', x, ' ', which(pws == x), '\n')

        ## Append to result table
        r <- rbind(r, tmp)

        ## Remove temporary object
        rm(tmp)
    }

    ############################
    ## Adjust p-values for multiple testing
    ############################

    ## Apply multiple testing correction per predictor
    r <- mutate(
        r,
        .by = 'predictor',
        cor_xy_padj = p.adjust(cor_xy_p, n = length(rownames(so$RNA))),
        fit_padj = p.adjust(fit_p, n = length(rownames(so$RNA))),
        geno_padj = p.adjust(geno_p, n = length(rownames(so$RNA)))
    )

    ## Save DSR results
    write.table(
        r,
        'Data/dsr_sig-gene.txt',
        sep = '\t',
        row.names = FALSE,
        quote = FALSE
    )

### Summarize DSR results in a heatmap

    ############################
    ## Define gene-level DE direction for annotations
    ############################

    ## Load differential gene expression results
    de <- read.delim('Data/de_genes.txt')

    ## Filter genes by DE criteria
    de <- de[which(
        abs(de$logFC) > 0.5 &
        de$padj < 0.01 &
        (de$pct.1 > 0.25) &
        (de$pct.2 > 0.25)
    ), ]

    ## Average logFC per gene
    de <- reframe(de, .by = 'feature', lfc = mean(logFC))

    ## Encode gene direction
    de.mean.lfc <- structure(
        c('dn','up')[1 + as.integer(de$lfc > 0)],
        names = de$feature
    )


    ############################
    ## Define pathway-level DE direction for annotations
    ############################

    ## Load differential pathway results
    pw <- read.delim("Data/de_pw_pid.txt")

    ## Filter pathways by DE criteria
    pw <- pw[which(abs(pw$logFC) > 0.25 & pw$padj < 0.01), ]

    ## Average logFC per pathway
    pw <- reframe(pw, .by = 'feature', lfc = mean(logFC))

    ## Encode pathway direction
    pw.mean.lfc <- structure(
        c('dn','up')[1 + as.integer(pw$lfc > 0)],
        names = pw$feature
    )


    ############################
    ## Prepare DSR results for heatmap visualization
    ############################

    ## Load DSR results
    r <- read.delim('Data/dsr_sig-gene.txt')

    ## Annotate response gene direction
    r$dir.resp <- de.mean.lfc[r$response]

    ## Annotate predictor pathway direction
    r$dir.pred <- pw.mean.lfc[r$predictor]


    ############################
    ## Filter DSR results by model quality and biological relevance
    ############################

    ## Apply quality and relevance filters
    f <- r[which(
        r$fit_padj < 0.01 &
        r$fit_dev > 0.25 &
        r$fit_r2 > 0.1 &
        abs(r$cor_xy_r2) > 0.25 &
        r$predictor %in% unique(pw$feature) &
        r$response %in% unique(de$feature) &
        !grepl('^ENSG0', r$response)
    ), ]

    ## Subset full results to filtered predictors and responses
    r <- r[which(
        r$predictor %in% unique(f$predictor) &
        r$response %in% unique(f$response)
    ), ]


    ############################
    ## Cast DSR correlations into matrix form
    ############################

    ## Create response-by-predictor correlation matrix
    d <- reshape2::dcast(r, response ~ predictor,
                         value.var = 'cor_xy_r2', fill = 0)

    ## Assign gene names as row names
    d <- data.frame(d[, -1], row.names = d[, 1])


    ############################
    ## Create annotation data frames for heatmap
    ############################

    ## Create gene direction annotation
    anno_gene <- unique(r[, c('response','dir.resp')])
    anno_gene <- data.frame(
        'Gene_dir.' = anno_gene$dir.resp,
        row.names = anno_gene$response
    )

    ## Create pathway direction annotation
    anno_pw <- unique(r[, c('predictor','dir.pred')])
    anno_pw <- data.frame(
        'PW_dir.' = anno_pw$dir.pred,
        row.names = str_replace_all(anno_pw$predictor, '-', '.')
    )

    ## Generate row labels for pathways
    labels_row <- str_remove(colnames(d), 'PID.') %>%
        str_replace_all(., '\\.', '-') %>%
        str_replace_all(., '-PATHWAY|PATHWAY', ' pathway')

    ## Define annotation color schemes
    anno_cols <- list(
        'PW_dir.' = c(dn = '#62D7DF', up = '#F1988D'),
        'Gene_dir.' = c(dn = '#62D7DF', up = '#F1988D')
    )


    ############################
    ## Generate and save heatmaps
    ############################

    ## Draw horizontal heatmap
    p <- cowplot::ggdraw() + cowplot::draw_grob(
        pheatmap::pheatmap(
            t(d),
            clustering_method = 'ward.D2',
            cutree_rows = 6,
            cutree_cols = 7,
            treeheight_col = 12,
            treeheight_row = 12,
            annotation_colors = anno_cols,
            annotation_col = anno_gene,
            annotation_row = anno_pw,
            labels_row = labels_row,
            angle_col = '90'
        )$gtable
    )

    ## Save horizontal heatmap
    ggsave(
        plot = p,
        'Figures/hm_dsr.pdf',
        height = 2.4,
        width = 5.75,
        scale = 1.4
    )

    ## Draw vertical heatmap
    p <- cowplot::ggdraw() + cowplot::draw_grob(
        pheatmap::pheatmap(
            d,
            clustering_method = 'ward.D2',
            cutree_rows = 7,
            cutree_cols = 6,
            treeheight_col = 12,
            treeheight_row = 12,
            annotation_colors = anno_cols,
            annotation_row = anno_gene,
            annotation_col = anno_pw,
            labels_col = labels_row,
            angle_col = '90'
        )$gtable
    )

    ## Save vertical heatmap
    ggsave(
        plot = p,
        'Figures/hm_dsr_vert.pdf',
        height = 4.5,
        width = 3.4,
        scale = 1.4
    )

### Make scatter plots of top DSR relationships

    ############################
    ## Define top DSR predictor–response pairs
    ############################

    ## Specify pathway–gene and pathway–pathway relationships
    top.dsm <- list(
        c('PID-ATF2-PATHWAY', 'JUN'),
        c('PID-IL2-STAT5-PATHWAY', 'LTA'),
        c('PID-NFAT-TFPATHWAY', 'IFNG'),
        c('PID-NFAT-TFPATHWAY', 'TNF'),
        c('PID-IL1-PATHWAY', 'NFKB1'),
        c('PID-TCR-PATHWAY', 'PID-GLYPICAN-1PATHWAY'),
        c('PID-GLYPICAN-1PATHWAY', 'TGFBR2'),
        c('PID-GLYPICAN-1PATHWAY', 'TGFB1'),
        c('PID-GLYPICAN-1PATHWAY', 'FGR'),
        c('PID-NFAT-TFPATHWAY', 'TBX21'),
        c('PID-MYC-ACTIV-PATHWAY', 'GZMB'),
        c('PID-IL23-PATHWAY', 'NFKB1')
    )

    ############################
    ## Generate scatter plots for top DSR hits
    ############################

    ## Create scatter plots using DSR testing function
    plist <- lapply(top.dsm, function(x){
        testDsr(x[1], x[2], include_geno = TRUE)
    })

    ## Combine scatter plots into grid
    p <- cowplot::plot_grid(
        plotlist = plist,
        align = 'hv',
        axis = 'lbt',
        ncol = 2
    )

    ## Save scatter plot figure
    ggsave(
        plot = p,
        'Figures/dsm_tophits.pdf',
        height = 0.75 * length(top.dsm),
        width = 3.5 * 2,
        scale = 2
    )

### PCA and dot plots of key pathway genes

    ############################
    ## Define pathways and load gene sets
    ############################

    ## Specify pathways of interest
    pws <- c(
        "PID-NFAT-TFPATHWAY","PID-HNF3A-PATHWAY","PID-IL2-STAT5-PATHWAY",
        "PID-ATF2-PATHWAY","PID-GLYPICAN-1PATHWAY","PID-IL23-PATHWAY"
    )

    ## Load canonical pathway gene sets
    gs <- fgsea::gmtPathways('ref/c2.canonical_pathways.v7.4.symbols.gmt.txt')

    ## Restrict to PID pathways
    gs <- gs[grepl('^PID_', names(gs))]

    ## Harmonize pathway naming convention
    names(gs) <- str_replace_all(names(gs), '_', '-')


    ############################
    ## Select highly variable and well-detected genes
    ############################

    ## Retrieve highly variable genes
    hvg <- VariableFeatures(so$RNA)

    ## Compute detection rate per gene
    det_rate <- apply(so$RNA$data[hvg,], 1, PercentAbove, threshold = 0)

    ## Filter genes by detection rate
    hvg <- hvg[which(det_rate > 0.2)]

    ## Intersect HVGs with pathway gene sets
    glist <- lapply(gs[pws], function(x) intersect(hvg, x))


    ############################
    ## Construct combined cell identity
    ############################

    ## Create combined cell type and genotype identifier
    so$id <- paste0(
        str_split_i(so$ct, '-', 1), ':',
        str_replace(str_split_i(so$sample, '_', 1), 'st', 'stop')
    )


    ############################
    ## Generate mNK dot plots and PCA plots
    ############################

    ## Set RNA as default assay
    DefaultAssay(so) <- 'RNA'

    ## Save full Seurat object
    so2 <- so

    ## Subset to mNK cells
    so <- subset(so2, subset = (nk_type == 'mNK'))

    ## Generate dot plots and PCA plots for mNK
    plist.mnk <- lapply(names(glist), function(n) {

        ## Generate plot title
        ptitle <- str_remove(n, 'PID-') %>%
            str_replace(., '-PATHWAY|PATHWAY', ' path.')

        ## Create dot plot and PCA plot
        MakeDotPlot(glist[[n]], title = ptitle, thresh = 0.1)
    }) %>% unlist(recursive = FALSE)

    ## Combine mNK dot plots
    p <- cowplot::plot_grid(
        plotlist = plist.mnk[seq(1, length(plist.mnk), 2)],
        align = 'hv',
        axis = 'lbt',
        ncol = 3
    )

    ## Save mNK dot plot figure
    ggsave(
        plot = p,
        'Figures/heatdot_signatures_mNK_sizeB.pdf',
        limitsize = FALSE,
        height = 1.2 * 2,
        width = 7.25,
        scale = 3.5
    )


    ############################
    ## Generate iNK dot plots and PCA plots
    ############################

    ## Subset to iNK cells
    so <- subset(so2, subset = (nk_type == 'iNK'))

    ## Generate dot plots and PCA plots for iNK
    plist.ink <- lapply(names(glist), function(n) {

        ## Generate plot title
        ptitle <- str_remove(n, 'PID-') %>%
            str_replace(., '-PATHWAY|PATHWAY', ' path.')

        ## Create dot plot and PCA plot
        MakeDotPlot(glist[[n]], title = ptitle, thresh = 0.1)
    }) %>% unlist(recursive = FALSE)

    ## Combine iNK dot plots
    p <- cowplot::plot_grid(
        plotlist = plist.ink[seq(1, length(plist.ink), 2)],
        align = 'hv',
        axis = 'lbt',
        ncol = 3
    )

    ## Save iNK dot plot figures
    ggsave(
        plot = p,
        'Figures/heatdot_signatures_iNK_sizeA.pdf',
        limitsize = FALSE,
        height = 0.8 * 2,
        width = 5,
        scale = 3.5
    )

    ggsave(
        plot = p,
        'Figures/heatdot_signatures_iNK_sizeB.pdf',
        limitsize = FALSE,
        height = 0.8 * 2,
        width = 7.25,
        scale = 3.5
    )

    ## Restore full Seurat object
    so <- so2

    ## Remove temporary object
    rm(so2)


    ############################
    ## Generate PCA plots across mNK and iNK
    ############################

    ## Combine PCA plots from mNK and iNK
    plist <- append(
        plist.mnk[seq(2, length(plist.mnk), 2)],
        plist.ink[seq(2, length(plist.ink), 2)]
    )

    ## Combine PCA plots into grid
    p <- cowplot::plot_grid(
        plotlist = plist,
        align = 'hv',
        axis = 'lbt',
        ncol = 6
    )

    ## Save PCA plot figure
    ggsave(
        plot = p,
        'Figures/pca_signatures.pdf',
        limitsize = FALSE,
        height = 0.9 * 2,
        width = 1.05 * 6,
        scale = 2
    )
