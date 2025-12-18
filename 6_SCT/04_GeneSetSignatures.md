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

### Pathway signature analysis

    ############################
    ## Build AUCell pathway signatures
    ############################

    ## Extract raw UMI counts matrix
    m <- so$RNA$counts

    ## Build gene expression rankings for AUCell
    m2 <- AUCell_buildRankings(
        exprMat = m,
        featureType = "genes",
        splitByBlocks = TRUE,
        BPPARAM = BiocParallel::MulticoreParam(tasks = 100, force.GC = TRUE, progressbar = TRUE),
        verbose = TRUE
    )

    ############################
    ## Load and filter pathway gene sets
    ############################

    ## Load canonical pathway gene sets from GMT file
    gs <- fgsea::gmtPathways('ref/c2.canonical_pathways.v7.4.symbols.gmt.txt')

    ## Filter gene sets to PID pathways only
    gs <- gs[grepl('^PID_', names(gs))]

    ############################
    ## Compute AUCell AUC scores for pathways
    ############################

    ## Calculate AUC scores for pathway signatures
    m <- AUCell_calcAUC(gs, m2, nCores = 20, verbose = TRUE)

    ## Save AUCell pathway signature object
    saveRDS(m, 'Datasets/aucell_PidSignatures.rds')

    ## Remove rankings object from memory
    rm(m2)

    ############################
    ## Add pathway signature assay to Seurat object
    ############################

    ## Extract AUC matrix aligned to Seurat cell order
    m <- (m@assays@data$AUC[, colnames(so)])

    ## Create assay object containing pathway signatures
    so$pw <- CreateAssay5Object(count = m, data = m)

### Differential signature analysis (pathways)

    ############################
    ## Set up differential signature comparison
    ############################

    ## Define unique cell type groups
    groups <- sort(unique(so$ct))

    ## Define group 1 genotype label
    gp1 <- 'OR7A10_Tumor'

    ## Define group 2 genotype label
    gp2 <- 'OR7A10st_Tumor'

    ## Set default assay to pathway signature assay
    DefaultAssay(so) <- 'pw'

    ############################
    ## Run differential signature analysis per cell type
    ############################

    ## Compute differential pathway signatures for each cell type
    pw <- do.call(rbind, lapply(groups, function(celltype){

        ## Run differential signature test between genotypes within a cell type
        de <- FindMarkers(
            subset(so, subset = (ct == celltype)),
            group.by = 'sample',
            logfc.threshold = 0,
            densify = TRUE,
            features = rownames(so),
            ident.1 = gp1,
            ident.2 = gp2
        )

        ## Standardize column names
        colnames(de) <- c('pval','logFC','pct.1','pct.2','padj')

        ## Add metadata columns to differential results
        de <- data.frame(
            celltype = celltype,
            feature = rownames(de),
            comp = paste0(gp1, '-', gp2),
            group1 = gp1,
            group2 = gp2,
            de
        )

        ## Order results by significance and effect size
        de <- de[with(de, order(-log10(pval), abs(logFC), decreasing = TRUE)),]

        ## Return result table for current cell type
        return(de)
    }))

    ############################
    ## Add gene set sizes and save results
    ############################

    ## Sort pathway results for clean output
    pw <- pw[with(pw, order(celltype, pval)),]

    ## Compute gene set lengths intersecting measured genes
    gs.length <- lapply(gs, function(x) length(intersect(rownames(so$RNA), x))) %>%
        unlist()

    ## Set names for gene set size mapping
    names(gs.length) <- str_replace_all(names(gs.length), '_', '-')

    ## Add gene set size column to pathway DE table
    pw$set_size <- gs.length[pw$feature]

    ## Save differential pathway table
    write.table(
        pw,
        'Data/de_pw_pid.txt',
        sep = '\t',
        quote = FALSE,
        row.names = FALSE
    )

### Pathway volcano plots

    ############################
    ## Load and format pathway DE results
    ############################

    ## Load differential pathway signature results
    de <- read.delim('Data/de_pw_pid.txt')

    ## Simplify pathway labels for plotting
    de$feature <- str_remove_all(de$feature, 'PID-|PATHWAY|-PATHWAY')

    ## Create group column required by volcano plotting function
    de$group <- de$celltype

    ############################
    ## Generate volcano plots per cell type
    ############################

    ## Build volcano plot list across cell types
    plist <- lapply(sort(unique(de$celltype)), function(celltype){

            ## Define plot title
            plot.title <- celltype

            ## Print plot title to console
            cat(plot.title, '\n')

            ## Generate volcano plot for current cell type
            plotVolcano(
                de[which(de$celltype == celltype),],
                text.size = 2,
                xlim.expansion = 3,
                pt.size = 1,
                rastr = TRUE,
                lfc.thresh = 0.2,
                title = plot.title
            )
        })

    ############################
    ## Combine and save volcano plot grid
    ############################

    ## Combine pathway volcano plots into a grid
    p <- cowplot::plot_grid(
        plotlist = plist,
        align = 'hv',
        axis = 'lbt',
        ncol = 4
    )

    ## Save volcano plot figure
    ggsave(
        plot = p,
        'Figures/vol_de_pw.pdf',
        height = 1.25 * 2,
        width = 1.25 * 4,
        scale = 2
    )

### Pathway upset plots

    ############################
    ## Prepare data for pathway upset plotting
    ############################

    ## Load differential pathway signature results
    de <- read.delim("Data/de_pw_pid.txt")

    ## Extract unique comparisons
    comps <- unique(de$comp)

    ## Assign regulation direction based on log fold change
    de$dir <- c('up','dn')[1 + as.integer(de$logFC < 0)]

    ## Filter pathways by significance and effect size
    de <- de[which(de$padj < 0.01 & abs(de$logFC) > 0.2),]

    ## Create set identifiers combining direction and cell type
    de$id <- paste0(de$dir, ':', de$celltype)

    ## Define universe of pathways
    g <- unique(de$feature)

    ## Build logical membership matrix for each set
    d <- do.call(cbind, lapply(unique(de$id), function(id){
        g %in% de$feature[which(de$id == id)]
    }))

    ## Assign set names to membership matrix columns
    colnames(d) <- unique(de$id)

    ## Convert membership matrix to list format for venn/upset
    d <- ggvenn::data_frame_to_list(data.frame(d))

    ## Clean set labels for plotting
    names(d) <- str_replace(names(d), 'up.', 'UP:') %>%
        str_replace(., 'dn.', 'DN:') %>%
        str_split_i(., '\\.', 1)

    ############################
    ## Generate upset plots for UP and DN pathways
    ############################

    ## Define colors for cell types
    cols <- structure(ggcolor(8), names = str_split_i(sort(unique(so$ct)), '-', 1))

    ## Build upset plot list for upregulated and downregulated sets
    plist <- lapply(c('UP','DN'), function(direction){

        ## Construct venn object for direction-specific sets
        v <- Venn(d[grep(direction, names(d))])

        ## Compute number of non-empty intersections
        n.intersects <- sum(process_region_data(v)$count > 0)

        ## Extract region-level intersection data
        pr <- process_region_data(v)

        ## Filter to non-empty regions
        pr <- pr[which(pr$count > 0),]

        ## Order intersections by size and identifier
        pr <- pr[order(-pr$count, pr$id),]

        ## Filter out very small intersection labels
        pr <- pr[which(nchar(pr$id) > 5),]

        ## Extract pathway lists per intersection
        glist <- lapply(pr$item, function(i) g[i]) %>% structure(., names = pr$id)

        ## Print intersection pathway lists
        print(glist)
        
        ## Define bar colors based on cell type group
        bar.col <- cols[str_split_i(v@names, ':', 2)]

        ## Generate upset plot
        p <- plot_upset(
            v,
            nintersects = n.intersects,
            sets.bar.color = bar.col,
            order.set.by = 'name',
            relative_height = 2,
            sets.bar.x.label = '# DE Pathways'
        )

        ## Wrap plot for cowplot compatibility
        p <- cowplot::ggdraw() + cowplot::draw_grob(aplot::gglistGrob(p))
        
        ## Return upset plot
        return(p)
    })

    ############################
    ## Combine and save upset plots
    ############################

    ## Combine UP and DN upset plots vertically
    p <- cowplot::plot_grid(
        plotlist = plist,
        align = 'hv',
        axis = 'lbt',
        ncol = 1,
        rel_widths = c(0.4, 1)
    )

    ## Save upset plot figure
    ggsave(
        plot = p,
        'Figures/upset_de_pathways.pdf',
        height = 2 * 2,
        width = 1 * 5,
        scale = 2.2
    )

### Pathway network analyses

    ############################
    ## Load and format pathway results for network visualization
    ############################

    ## Load differential pathway signature results
    de <- read.delim("Data/de_pw_pid.txt")

    ## Order results by cell type and effect size
    de <- de[with(de, order(celltype, -abs(logFC))),]

    ## Construct group label for network plotting
    de$group <- paste0(de$celltype, ': ', de$group1, '-', de$group2)

    ## Store cell type column expected by network function
    de$ct <- de$celltype

    ## Reformat pathway names for network plotting
    de$pathway <- str_replace_all(de$feature, '-', '_')

    ## Store enrichment score as logFC
    de$NES <- de$logFC

    ############################
    ## Split into upregulated and downregulated pathways
    ############################

    ## Subset significant upregulated pathways
    up <- de[which(de$logFC > 0 & de$padj < 0.01),]

    ## Subset significant downregulated pathways
    dn <- de[which(de$logFC < 0 & de$padj < 0.01),]

    ############################
    ## Generate network plots and save
    ############################

    ## Combine upregulated and downregulated network plots
    p <- cowplot::plot_grid(
        RunNetPW(up, 'Upreg. pathways', 2, 1),
        RunNetPW(dn, 'Downreg. pathways', 2, 1),
        align = 'hv',
        axis = 'lbt',
        ncol = 2,
        rel_widths = c(1, 1)
    )

    ## Save pathway network figure
    ggsave(
        plot = p,
        paste0('Figures/net_pw_pid.pdf'),
        height = 3 * 1,
        width = 5.5 * 2,
        scale = 1.5
    )

### Phenotype signatures

    ############################
    ## Define pathways and genes for visualization
    ############################

    ## Define pathway signatures to plot
    pw2plot <- c(
        'PID-TCR-JNK-PATHWAY',
        'PID-NFAT-TFPATHWAY',
        'PID-IL2-STAT5-PATHWAY',
        'PID-GLYPICAN-1PATHWAY',
        'PID-IL1-PATHWAY'
    )

    ## Define phenotype genes to plot
    genes2plot <- c('GZMB','IFNG','LTA','BCL2','BCL2L1','BAX')

    ## Define fill colors by sample
    cols <- structure(
        c('#B2B2B2', '#E9AC4C'),
        names = rev(unique(so$sample))
    )

    ############################
    ## Prepare labels and factor ordering
    ############################

    ## Set cell type factor ordering for plots
    so$ct <- factor(so$ct, levels = rev(sort(unique(so$ct))))

    ## Create cell type label map for y-axis
    labs.ct <- structure(
        str_split_i(rev(sort(unique(so$ct))), '-', 1),
        names = rev(levels(so$ct))
    )

    ## Create pathway label map for plotting
    labs.pw <- structure(
        str_remove(pw2plot, 'PID-'),
        names = pw2plot
    )

    ############################
    ## Plot pathway signature violins
    ############################

    ## Switch default assay to pathway signatures
    DefaultAssay(so) <- 'pw'

    ## Create stacked violin plot for pathway signatures
    p1 <- VlnPlot(
        so,
        pw2plot,
        group.by = 'ct',
        pt.size = 0,
        stack = TRUE,
        fill.by = 'ident',
        adjust = 0.6,
        split.by = 'sample',
        same.y.lims = FALSE,
        flip = FALSE
    ) +
        geom_vline(xintercept = 0.1, linetype = 'dashed', linewidth = 0.1) +
        scale_fill_manual(values = cols) +
        scale_y_discrete(labels = labs.ct) +
        labs(x = 'Pathway signature', y = NULL) +
        theme(
            legend.position = 'none',
            legend.title = element_text(size = 10),
        )

    ############################
    ## Plot gene expression violins
    ############################

    ## Switch default assay to RNA expression
    DefaultAssay(so) <- 'RNA'

    ## Create stacked violin plot for phenotype genes
    p2 <- VlnPlot(
        so,
        genes2plot,
        group.by = 'ct',
        pt.size = 0,
        stack = TRUE,
        fill.by = 'ident',
        adjust = 0.6,
        split.by = 'sample',
        same.y.lims = FALSE,
        flip = FALSE
    ) +
        geom_vline(xintercept = 0.1, linetype = 'dashed', linewidth = 0.1) +
        scale_fill_manual(values = cols) +
        scale_y_discrete(labels = labs.ct) +
        labs(x = 'Gene expr.', y = NULL) +
        theme(
            legend.position = 'none',
            axis.text.x.top = element_text(angle = '90')
        )

    ############################
    ## Combine and save violin plot figure
    ############################

    ## Combine pathway and gene violin plots
    p <- cowplot::plot_grid(
        p1,
        p2,
        align = 'hv',
        axis = 'lbt',
        ncol = 2,
        rel_widths = c(0.75, 1)
    )

    ## Save combined violin plot figure
    ggsave(
        plot = p,
        'Figures/topPW_Vln_Signaling.pdf',
        height = 3.5,
        width = 2.7 * 2,
        scale = 1.5
    )


    ############################
    ## Extract significance labels for pathway signatures
    ############################

    ## Load differential pathway signature results
    de <- read.delim("Data/de_pw_pid.txt")

    ## Subset to pathways of interest
    de <- de[which(de$feature %in% pw2plot),]

    ## Set pathway factor ordering
    de$feature <- factor(de$feature, levels = pw2plot)

    ## Order rows and select output columns
    de <- de[with(de, order(feature, celltype)), c(1, 2, 7, 10)]

    ## Convert adjusted p-values into significance stars
    de$sig <- sapply(de$padj, function(x) {
        case_when(
            x < 1e-4 ~ '****',
            x < 1e-3 ~ '***',
            x < 1e-2 ~ '**',
            x < 0.05 ~ '*',
            x >= 0.05 ~ 'ns'
        )
    })

    ############################
    ## Extract significance labels for phenotype genes
    ############################

    ## Load differential gene expression results
    de <- read.delim("Data/de_genes.txt")

    ## Subset to genes of interest
    de <- de[which(de$feature %in% genes2plot),]

    ## Set gene factor ordering
    de$feature <- factor(de$feature, levels = genes2plot)

    ## Order rows and select output columns
    de <- de[with(de, order(feature, celltype)), c(1, 2, 7, 10)]

    ## Convert adjusted p-values into significance stars
    de$sig <- sapply(de$padj, function(x) {
        case_when(
            x < 1e-4 ~ '****',
            x < 1e-3 ~ '***',
            x < 1e-2 ~ '**',
            x < 0.05 ~ '*',
            x >= 0.05 ~ 'ns'
        )
    })
