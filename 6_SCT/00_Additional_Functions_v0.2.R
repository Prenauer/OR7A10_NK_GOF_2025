#!/usr/bin/env Rscript

############################
## Required libraries
############################

library(stringr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(cluster)
library(Seurat)
library(future)
library(ggrastr)
library(ggrepel)
library(viridis)
library(grid)

############################
## Combine two count matrices
############################

## Integrate genes and cells from two matrices
combineMatrices <- function(m1, m2, suffix = "-1") {
    
    ## Identify genes unique to each matrix
    g1 <- setdiff(rownames(m1), rownames(m2))
    g2 <- setdiff(rownames(m2), rownames(m1))
    
    ## Pad missing genes with zeros in m1
    m1 <- rbind(
        m1,
        matrix(
            0,
            nrow = length(g2),
            ncol = ncol(m1),
            dimnames = list(g2, colnames(m1))
        )
    )
    
    ## Pad missing genes with zeros in m2
    m2 <- rbind(
        m2,
        matrix(
            0,
            nrow = length(g1),
            ncol = ncol(m2),
            dimnames = list(g1, colnames(m2))
        )
    )
    
    ## Ensure same gene order
    m2 <- m2[rownames(m1), , drop = FALSE]
    
    ## Rename second matrix columns
    colnames(m2) <- stringr::str_replace(colnames(m2), "-1", suffix)
    
    ## Merge matrices
    m1 <- cbind(m1, m2)
    
    return(m1)
}

############################
## Determine optimal clustering resolution
############################

DetermineOptimalClusters <- function(
        so,
        resolution,
        plot.only = FALSE,
        graph.name = "SCT_snn",
        plot.res = NULL,
        optimalResolution = NULL,
        seed = 42
) {
    
    if (!plot.only) {
        
        ## Ensure resolution values are provided
        if (is.null(resolution)) {
            stop("`resolution` must be provided when plot.only = FALSE")
        }
        
        ## Sequential execution for reproducibility
        future::plan("sequential")
        
        ## Extract UMAP coordinates
        xy <- as.data.frame(so@reductions$umap@cell.embeddings)
        
        ## Distance matrix for silhouette calculation
        d <- dist(xy)
        
        ## Within-cluster sum-of-squares
        CalcWSS <- function(x, y) {
            (length(x) - 1) * (var(x) + var(y))
        }
        
        ## Evaluate clustering across resolutions
        optimalResolution <- do.call(
            rbind,
            lapply(resolution, function(res) {
                
                tmp <- FindClusters(
                    so,
                    graph.name = graph.name,
                    algorithm = 4,
                    resolution = res,
                    random.seed = seed
                )
                
                wss <- reframe(
                    cbind(xy, cluster = Idents(tmp)),
                    .by = "cluster",
                    wss = CalcWSS(umap_1, umap_2)
                )$wss %>% sum()
                
                asw <- mean(
                    cluster::silhouette(
                        as.integer(Idents(tmp)),
                        dist = d
                    )[, 3]
                )
                
                cat(res, "\n")
                
                data.frame(
                    res = res,
                    wss = wss,
                    asw = asw,
                    n_clust = length(levels(Idents(tmp)))
                )
            })
        )
    }
    
    ############################
    ## Select optimal resolution
    ############################
    
    if (is.null(plot.res)) {
        
        res <- optimalResolution[
            with(optimalResolution, order(n_clust, -asw, wss)),
            ,
            drop = FALSE
        ]
        
        r <- reframe(
            res,
            .by = "n_clust",
            asw.pick = head(res[which(asw == max(asw)), ], 1)
        )[, 2] %>% unlist()
        
    } else {
        r <- plot.res
    }
    
    ############################
    ## Diagnostic plots
    ############################
    
    p.resCheck <- cowplot::plot_grid(
        
        ggplot(optimalResolution, aes(x = res, y = n_clust)) +
            geom_point(color = "tomato2") +
            geom_path(color = "tomato2") +
            geom_vline(xintercept = r, linetype = "dashed",
                       color = "gray25", alpha = 0.75) +
            labs(
                x = "Clustering resolution",
                y = stringr::str_wrap("# of clusters", 20)
            ) +
            theme_test(),
        
        ggplot(optimalResolution, aes(x = res, y = wss)) +
            geom_point(color = "cyan2") +
            geom_path(color = "cyan2") +
            geom_vline(xintercept = r, linetype = "dashed",
                       color = "gray25", alpha = 0.75) +
            labs(
                x = "Clustering resolution",
                y = stringr::str_wrap(
                    "Within-cluster sum-of-squares (WSS)", 20
                )
            ) +
            theme_test(),
        
        ggplot(optimalResolution, aes(x = res, y = asw)) +
            geom_point(color = "orange2") +
            geom_path(color = "orange2") +
            geom_vline(xintercept = r, linetype = "dashed",
                       color = "gray25", alpha = 0.75) +
            labs(
                x = "Clustering resolution",
                y = stringr::str_wrap(
                    "Average silhouette width", 20
                )
            ) +
            theme_test(),
        
        ncol = 1,
        align = "hv"
    )
    
    cat("Optimal res = ", paste0(sort(r), collapse = ","), "\n")
    
    return(list(res = optimalResolution, plot = p.resCheck))
}

############################
## Customized Seurat DotPlot
############################

customDotPlot <- function(
        so,
        genes,
        group.order = NULL,
        group.by = "ct",
        scale.max = 60
) {
    
    ## Reorder groups if provided
    if (!is.null(group.order)) {
        so@meta.data[, group.by] <- factor(
            so@meta.data[, group.by],
            levels = group.order
        )
    }
    
    ## Generate dot plot
    DotPlot(
        so,
        genes,
        group.by = group.by,
        scale.max = scale.max
    ) +
        labs(x = NULL, y = NULL) +
        theme(
            legend.position = "bottom",
            legend.key.spacing = unit(0.5, "lines"),
            legend.key.height = unit(0.5, "lines"),
            legend.key.width = unit(0.5, "lines"),
            legend.text = element_text(size = 7),
            legend.key.spacing.y = unit(0.01, "lines"),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.title = element_text(size = 10)
        ) +
        scale_color_viridis_c(
            option = "inferno",
            direction = -1,
            end = 0.9
        ) +
        guides(
            size = guide_legend(nrow = 2),
            alpha = guide_legend(nrow = 2)
        )
}

############################
## Volcano plot function
############################

plotVolcano <- function(
        de,
        lfc.thresh = 1,
        p.thresh = 0.01,
        pt.size = 2,
        title = NULL,
        text.size = 2.2,
        rastr = TRUE,
        xlim.expansion = 1.3
) {
    
    ## Ensure required columns exist
    stopifnot(all(c("padj", "logFC", "feature") %in% colnames(de)))
    
    ## Assign plotting variables
    de$p <- de$padj
    de$lfc <- de$logFC
    
    ## Order by raw p-value if present
    if ("pval" %in% colnames(de)) {
        de <- de[order(de$pval), ]
    }
    
    ## Select genes for labeling
    genes2label <- c(
        head(de$feature[de$lfc > lfc.thresh & de$p < p.thresh], 5),
        head(de$feature[de$lfc < -lfc.thresh & de$p < p.thresh], 5)
    ) %>% unique()
    
    ## Replace zero p-values
    plot.min.p <- min(de$p[de$p != 0], na.rm = TRUE) * 0.1
    de$p[de$p == 0] <- plot.min.p
    
    ## Axis limits
    plot.xlim <- xlim.expansion * range(de$lfc, na.rm = TRUE)
    plot.ylim <- max(-log10(de$p), na.rm = TRUE) * c(-0.02, 1.2)
    
    ## Base plot
    p <- ggplot(de, aes(x = lfc, y = -log10(p)))
    
    ## Points
    if (rastr) {
        p <- p +
            geom_point_rast(color = "gray", size = pt.size) +
            geom_point_rast(
                data = de[de$lfc > lfc.thresh & de$p < p.thresh, ],
                color = "firebrick",
                size = pt.size
            ) +
            geom_point_rast(
                data = de[de$lfc < -lfc.thresh & de$p < p.thresh, ],
                color = "steelblue",
                size = pt.size
            )
    } else {
        p <- p +
            geom_point(color = "gray", alpha = 0.5, size = pt.size)
    }
    
    ## Labels
    p <- p +
        ggrepel::geom_text_repel(
            data = de[de$feature %in% genes2label, ],
            aes(label = feature),
            size = text.size,
            min.segment.length = 0,
            max.iter = 1e5,
            segment.alpha = 0.7,
            fontface = "italic"
        ) +
        xlim(plot.xlim) +
        ylim(plot.ylim) +
        theme_classic() +
        labs(
            x = "Log2 fold-change",
            y = "Adj. p (-log10)",
            title = title
        )
    
    return(p)
}

############################
## Network plot function
############################

## Function makes a network plot from pathway analysis results
RunNetPW <- function(pw, plot.title = NULL, n.pw.per.group = 5, seed = 42) {
    
    ## Ensure required packages are loaded
    require(ggraph)
    require(igraph)
    require(ggnetwork)
    require(ggrepel)
    
    ## Setup graph (select top pathways per group)
    top.pw <- slice_head(pw, by = 'group', n = n.pw.per.group)$pathway
    pw <- pw[which(pw$pathway %in% top.pw), ]
    
    ## Create edge list for (celltype -> comparison) and (comparison -> pathway)
    e <- rbind(
        data.frame(from = pw$ct, to = pw$group, size = 0.1),
        data.frame(from = pw$group, to = pw$pathway, size = pw$NES / max(pw$NES, TRUE))
    )
    
    ## Create vertex list
    v <- rbind(
        data.frame(name = pw$pathway, type = 'pathway'),
        data.frame(name = pw$ct, type = 'celltype'),
        data.frame(name = pw$group, type = 'comparison')
    )
    v <- unique(v)
    
    ## Build igraph object
    ig <- graph_from_data_frame(e, vertices = v)
    
    ## Compute layout (Fruchterman–Reingold)
    set.seed(seed)
    ig2 <- (ig %>% add_layout_(with_fr()))
    
    ## Convert graph to ggnetwork object
    net <- ggnetwork(ig2)
    
    ## Initialize class labels
    net$class <- 'edge'
    
    ## Mark nodes (where start and end coordinates match)
    net$class[which(net$x == net$xend & net$y == net$yend)] <- 'node'
    
    ## Reformat pathway names
    net$name[which(net$type == 'pathway')] <- net$name[which(net$type == 'pathway')] %>%
        str_remove(., 'GOBP_') %>% str_remove(., 'REACTOME_') %>%
        str_remove(., 'KEGG_') %>% str_remove(., 'WP_') %>%
        str_remove(., 'BIOCARTA_') %>% str_replace_all(., '_', ' ')
    
    ## Exclude celltype nodes and edges from the visualization
    net <- net[which(net$type != 'celltype'), ]
    
    ## Initialize node labels
    net$node_label <- NA
    
    ## Set node labels for nodes
    net$node_label[which(net$class == 'node')] <-
        str_split_i(net$name[which(net$class == 'node')], ': ', 1)
    
    ## Define colors for subsets
    cols <- structure(
        c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC"),
        names = c("iNK_c1-Cytotoxic", "iNK_c2-CREM", "mNK_c1-MKI67", "mNK_c2-MKI67-BATF3",
                  "mNK_c3-BATF3", "mNK_c4-IL32", "mNK_c5-KLRC2", "mNK_c6-EGR2")
    )
    
    ## Extract group label prefix for coloring
    net$gp2 <- str_split_i(net$name, ':', 1)
    
    ## Compute node sizes (per node label group)
    net <- mutate(net, .by = 'node_label', nodesize = length(node_label))
    
    ## Calculate pathway node sizes relative to number of connecting edges
    node_size <- reframe(pw, .by = 'feature', n = length(feature))
    node_size <- structure(node_size$n,
                           names = str_replace_all(node_size$feature, '-', ' '))
    
    ## Map pathway sizes to network nodes
    net$node_size <- node_size[net$gp2]
    
    ## Fix pathway label names
    net[which(net$class == 'node' & net$type == 'pathway'), 'name'] <-
        str_remove_all(
            net[which(net$class == 'node' & net$type == 'pathway'), 'name'],
            'PID '
        )
    
    ## Fix comparison label names
    net[which(net$class == 'node' & net$type == 'comparison'), 'node_label'] <-
        str_split_i(
            net[which(net$class == 'node' & net$type == 'comparison'), 'node_label'],
            '-',
            1
        )
    
    ## Generate plot
    ggplot(net, aes(x = x, y = y)) +
        geom_edges(
            color = "#DDDDDD",
            aes(linewidth = size, xend = xend, yend = yend)
        ) +
        geom_edges(
            data = net[which(net$class == 'edge' & net$type == 'comparison'), ],
            arrow = arrow(length = unit(0.25, 'lines'), type = "closed"),
            aes(linewidth = size, color = gp2, xend = xend, yend = yend)
        ) +
        scale_linewidth(range = c(0.1, 1)) +
        theme_blank() +
        geom_nodes(
            data = net[which(net$class == 'node' & net$type == 'comparison'), ],
            shape = 22, color = 'gray30', size = 6,
            aes(fill = gp2)
        ) +
        geom_nodes(
            data = net[which(net$class == 'node' & net$type == 'pathway'), ],
            shape = 21, color = 'gray30', fill = '#20854EFF',
            aes(size = node_size)
        ) +
        scale_size_identity() +
        geom_nodetext_repel(
            data = net[which(net$class == 'node' & net$type == 'comparison'), ],
            aes(label = str_wrap(node_label, 20)),
            size = 3, min.segment.length = 0, nudge_y = 0.5,
            segment.size = 0.1, segment.alpha = 0.5, max.overlaps = 50
        ) +
        geom_nodetext_repel(
            data = net[which(net$class == 'node' & net$type == 'pathway'), ],
            aes(label = str_wrap(name, 20)),
            size = 3, min.segment.length = 0,
            nudge_y = 0.5, lineheight = 0.8, fontface = 'bold',
            segment.size = 0.1, segment.alpha = 0.5, max.overlaps = 50
        ) +
        scale_fill_manual(values = cols) +
        scale_color_manual(values = cols) +
        labs(
            title = plot.title,
            fill = 'Subset',
            linewidth = str_wrap('Pathway enrichment score', 8)
        ) +
        guides(
            fill = guide_legend(order = 1),
            linewidth = guide_legend(order = 2),
            color = FALSE
        ) +
        theme(
            legend.key.height = unit(0.2, 'lines'),
            legend.spacing.x = unit(0.2, 'lines')
        )
}

############################
## Extended dot plot function
############################

## Function makes a hclust dot plot with a dendrogram
MakeDotPlot <- function(g, title = '', thresh = 1/5, method = 'ward.D2',
                        cluster.axis = 'xy') {
    
    ## Ensure required packages are loaded
    require(ggdendroplot)
    
    ## Restrict genes to HVGs present in input set
    g <- intersect(hvg, g)
    
    ## Build matrix for PCA from average expression by id
    m <- t(so$RNA$data[g, ])
    m <- S4Vectors::splitAsList(m, factor(so$id, levels = unique(so$id)))
    m <- do.call(rbind, lapply(m, colMeans))
    
    ## Harmonize rownames for downstream labels
    rownames(m) <- str_replace(rownames(m), '_', '-')
    
    ## Run PCA on averaged expression
    res.pca <- prcomp(m, scale = TRUE)
    
    ## Filter genes by detection rate threshold
    g <- g[(apply(so$RNA$data[g, ], 1, PercentAbove, threshold = 0) > thresh)]
    
    ## Build matrix for dotplot from averaged expression by id
    m <- t(so$RNA$data[g, ])
    m <- S4Vectors::splitAsList(m, factor(so$id, levels = unique(so$id)))
    m <- do.call(rbind, lapply(m, colMeans)) %>% t()
    
    ## Initialize clustering defaults
    cl.x <- g
    cl.y <- c(
        paste0(as.character(sort(unique(so$ct))), ':OR7A10'),
        paste0(as.character(sort(unique(so$ct))), ':OR7A10st')
    )
    
    ## Cluster along y-axis (samples / ids) if requested
    if (cluster.axis %in% c('y', 'xy')) {
        cl.y <- hclust(dist(t(m)), method = method)
    }
    
    ## Cluster along x-axis (genes) if requested
    if (cluster.axis %in% c('x', 'xy')) {
        cl.x <- hclust(dist(m), method = method)
    }
    
    ## Apply y-axis ordering to Seurat metadata id factor
    so$id <- factor(so$id, levels = cl.y$labels[cl.y$order])
    
    ## Define axis color (computed but not used in plotting here)
    axis.color <- c('#D98500', 'gray30')[1 + as.integer(grepl('OR7A10st', levels(so$id)))]
    
    ## Apply x-axis ordering to gene list
    g <- g[cl.x$order]
    
    ## Make dot plot with ordered ids and genes
    p1 <- DotPlot(so, g, group.by = 'id') +
        labs(x = NULL, y = NULL, title = title) +
        scale_color_distiller(palette = "RdYlBu") +
        guides(size = guide_legend(nrow = 1), alpha = guide_legend(nrow = 2)) +
        theme(
            legend.position = 'bottom',
            legend.key.spacing = unit(0.5, 'lines'),
            legend.key.height = unit(0.5, 'lines'),
            legend.key.width = unit(0.5, 'lines'),
            legend.text = element_text(size = 7),
            legend.key.spacing.y = unit(0.01, 'lines'),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = 'italic'),
            legend.title = element_text(size = 10),
            axis.text.y = ggtext::element_markdown(colour = '#D98500')
        )
    
    ## Create dendrogram plot for y-axis clustering
    p1.dend <- ggplot() + geom_dendro(cl.y, pointing = 'side') + theme_void()
    
    ## Combine dotplot and dendrogram
    p1 <- cowplot::plot_grid(
        p1,
        p1.dend,
        align = 'h',
        axis = 'lbt',
        ncol = 2,
        rel_widths = c(0.6, 0.2)
    )
    
    ## Build PCA plot data frame
    d <- data.frame(
        res.pca$x[, 1:2],
        id = rownames(res.pca$x),
        ct = substr(rownames(res.pca$x), 0, 6),
        geno = str_split_i(rownames(res.pca$x), ':', 2)
    )
    
    ## Generate PCA plot colored by genotype
    p2 <- ggplot(d, aes(x = PC1, y = PC2, fill = geno), color = 'gray20') +
        geom_point(size = 2.5, shape = 21, stroke = 0.01) +
        xlim(1.3 * range(d$PC1)) +
        theme_classic() +
        theme(legend.position = 'none') +
        scale_fill_manual(values = structure(c('#B2B2B2', '#E9AC4C'),
                                             names = c('OR7A10st', 'OR7A10'))) +
        scale_color_manual(values = structure(c('gray30', '#D98500'),
                                              names = c('OR7A10st', 'OR7A10'))) +
        ggrepel::geom_text_repel(
            aes(label = ct, color = geno),
            size = 2.5,
            min.segment.length = 0,
            segment.size = 0.2
        ) +
        labs(title = title)
    
    ## Return dotplot+dendrogram and PCA plot
    p <- list(p1, p2)
    
    return(p)
}

############################
## Matrix preprocessing utilities
############################

## Apply symmetric winsorization to each row of a matrix
winsorize <- function(m, thresh) {
    
    ## Apply winsorization row-wise
    m.win <- apply(m, 1, function(x) {
        
        ## Cap positive values
        x[x > thresh] <- thresh
        
        ## Cap negative values
        x[x < -thresh] <- -thresh
        
        ## Return modified vector
        return(x)
    })
    
    ## Restore original orientation and return as data frame
    return(as.data.frame(t(m.win)))
}


############################
## DSR model fitting
############################

## Fit a GAM for signature–signature or signature–gene relationships
dsrFitGam <- function(d, use.nb = NULL, use.umi = NULL, f = NULL,
                      filter.data = TRUE, select = TRUE) {
    
    ## Load required package
    require(mgcv)
    
    ## Optionally filter zero values
    if (filter.data) {
        d <- d[d$x > 0 & d$y > 0, ]
    }
    
    ## Infer data type from attributes if not explicitly provided
    if (all(c("x.type", "y.type") %in% names(attributes(d))) &&
        is.null(use.nb) && is.null(use.umi)) {
        
        ## Predictor is counts
        if (attr(d, "x.type") == "counts") {
            use.umi <- TRUE
        }
        
        ## Response is counts
        if (attr(d, "y.type") == "counts") {
            use.umi <- TRUE
            use.nb  <- TRUE
        }
    }
    
    ## Default to non-count models if unspecified
    if (is.null(use.nb) && is.null(use.umi)) {
        use.nb  <- FALSE
        use.umi <- FALSE
    }
    
    ############################
    ## Negative binomial GAM
    ############################
    
    if (use.nb) {
        
        ## Define model formula if missing
        if (is.null(f)) {
            f <- y ~ geno + ct + s(x, bs = "ps", k = 10) + offset(umi)
        }
        
        ## Fit model
        fit <- try(
            mgcv::bam(
                f,
                data    = d,
                family  = "nb",
                select  = select,
                control = list(keepData = TRUE),
                method  = "fREML"
            ),
            silent = TRUE
        )
    }
    
    ############################
    ## Gamma GAM
    ############################
    
    if (!use.nb) {
        
        ## Define formula with or without UMI offset
        if (use.umi) {
            f <- y ~ geno + ct + s(x, bs = "ps", k = 10) + offset(umi)
        } else {
            f <- y ~ geno + ct + s(x, bs = "ps", k = 10)
        }
        
        ## Fit model
        fit <- try(
            mgcv::bam(
                f,
                data    = d,
                family  = Gamma(link = "log"),
                select  = select,
                control = list(keepData = FALSE),
                method  = "fREML"
            ),
            silent = TRUE
        )
    }
    
    ## Exit gracefully on failure
    if (attempt::is_try_error(fit)) {
        return(NULL)
    }
    
    ## Preserve original row order for downstream indexing
    fit$order <- rownames(d)
    
    return(fit)
}


############################
## Extract DSR model statistics
############################

## Extract model statistics and correlation metrics
extractDsrResults <- function(fit) {
    
    ## Load required package
    require(mgcv)
    
    ## Extract model frame
    d <- fit$model
    
    ## Extract term-level predictions
    pr <- predict(fit, type = "terms")
    
    ## Abort if insufficient data
    if (nrow(pr) < 20) {
        return(NULL)
    }
    
    ############################
    ## Compute partial residuals
    ############################
    
    ## Predictor-only partial residuals
    d$p_x <- (pr[, "s(x)"] + residuals(fit))[fit$order]
    
    ## Predictor + genotype partial residuals
    d$p_gx <- (rowSums(pr[, c("geno", "s(x)")]) + residuals(fit))[fit$order]
    
    ############################
    ## Correlation statistics
    ############################
    
    ## Correlation between predictor and partial residuals
    ctx  <- cor.test(d$x, d$p_x)[c("p.value", "estimate")]
    
    ## Correlation including genotype effect
    ctxg <- cor.test(d$x, d$p_gx)[c("p.value", "estimate")]
    
    ############################
    ## Genotype coefficient
    ############################
    
    genoEst <- summary(fit)$p.table["genoOR7A10_Tumor", c(4, 1)]
    
    ct <- structure(
        c(ctx, ctxg, genoEst),
        names = c(
            "cor_xy_p", "cor_xy_r2",
            "cor_xgy_p", "cor_xgy_r2",
            "geno_p", "geno_est"
        )
    )
    
    ############################
    ## Model-level statistics
    ############################
    
    r <- structure(
        c(
            summary(fit)$s.table[1, 3:4],
            summary(fit)$r.sq,
            summary(fit)$dev.expl
        ),
        names = c("fit_F", "fit_p", "fit_r2", "fit_dev")
    )
    
    ## Combine results into a single-row data frame
    res <- as.data.frame(matrix(unlist(c(ct, r)), nrow = 1))
    colnames(res) <- c(names(ct), names(r))
    
    return(res)
}


############################
## Extract fitted effects for visualization
############################

## Extract predictor and genotype-adjusted effects
getDsrEffects <- function(fit) {
    
    ## Load required package
    require(mgcv)
    
    ## Extract model data
    d <- fit$model
    
    ## Extract terms and residuals
    pr <- cbind(
        predict(fit, type = "terms"),
        resid = residuals(fit)
    )[fit$order, ]
    
    ## Construct output table
    output <- data.frame(
        d,
        x_effect     = rowSums(pr[, c("s(x)", "resid")]),
        geno_effect  = rowSums(pr[, c("geno", "resid")]),
        xGeno_effect = rowSums(pr[, c("s(x)", "geno", "resid")])
    )
    
    return(output)
}


############################
## Prepare DSR input data
############################

## Create DSR input table from Seurat object
setupDsrData <- function(x, y, force.x.counts = FALSE, force.y.norm = FALSE) {
    
    ## Identify assays containing predictor and response
    assays <- Assays(so)
    
    a.x <- assays[sapply(assays, function(a) x %in% rownames(so[[a]]))][1]
    a.y <- assays[sapply(assays, function(a) y %in% rownames(so[[a]]))][1]
    
    
    ## Resolve metadata sources
    if (is.na(a.x) && x %in% colnames(so@meta.data)) a.x <- "metadata"
    if (is.na(a.y) && y %in% colnames(so@meta.data)) a.y <- "metadata"
    
    ## Extract data vectors
    if (a.y == "metadata") y.data <- so@meta.data[, y]
    if (a.x == "metadata") x.data <- so@meta.data[, x]
    
    if (a.y == "RNA") y.data <- so[[a.y]]$counts[y, ]
    if (!(a.y %in% c("RNA", "metadata"))) y.data <- so[[a.y]]$data[y, ]
    
    if (force.x.counts && a.x != "metadata") x.data <- so[[a.x]]$counts[x, ]
    if (!force.x.counts && a.x != "metadata") x.data <- so[[a.x]]$data[x, ]
    
    if (a.y != "metadata" && force.y.norm) y.data <- so[[a.y]]$data[y, ]
    
    ## Assemble DSR data frame
    d <- data.frame(
        x    = x.data,
        y    = y.data,
        ct   = factor(so$ct),
        umi  = log(so$nCount_RNA),
        geno = factor(
            so$sample,
            levels = c("OR7A10st_Tumor", "OR7A10_Tumor")
        )
    )
    
    ## Annotate data types
    attr(d, "x.type") <- c("data", "counts")[1 + as.integer(force.x.counts)]
    attr(d, "y.type") <- c("data", "counts")[1 + as.integer(a.y == "RNA")]
    
    return(d)
}


############################
## DSR visualization and diagnostics
############################

## Analyze DSR relationship and generate diagnostic scatter plots
testDsr <- function(x, y, raster = TRUE, include_geno = FALSE) {
    
    ############################
    ## Fit DSR model
    ############################
    
    ## Prepare data for model fitting
    d <- setupDsrData(x, y)
    
    ## Fit GAM model
    fit <- dsrFitGam(d, select = FALSE)
    
    ## Extract model statistics
    r <- extractDsrResults(fit)
    
    ## Extract fitted effects
    eff <- getDsrEffects(fit)
    
    ## Report model being tested
    cat("Model: ", x, " as a predictor of ", y, "\n")
    
    
    ############################
    ## Prepare data for visualization
    ############################
    
    ## Reload data using normalized response
    d <- setupDsrData(x, y, force.y.norm = TRUE)
    
    ## Filter zero values
    d <- d[d$x > 0 & d$y > 0, ]
    
    ## Generate smoothing grid
    smooth.x <- seq(
        quantile(d$x, 0.01),
        quantile(d$x, 0.99),
        diff(quantile(d$x, c(0.01, 0.99))) / 100
    )
    
    
    ############################
    ## Generate statistical annotations
    ############################
    
    ## Raw correlation statistics
    p1.stats <- cor.test(d$x, d$y)[c("estimate", "p.value")] %>% unlist()
    
    ## Annotation for raw scatter
    p1.label <- make_scatter_anno(
        d$x,
        d$y,
        p1.stats[1],
        p1.stats[2]
    )
    
    ## Annotation for adjusted-effect scatter
    p2.label <- make_scatter_anno(
        eff$x,
        eff$x_effect,
        r$cor_xy_r2,
        r$cor_xy_p
    )
    
    
    ############################
    ## Initialize scatter plots
    ############################
    
    ## Randomize plotting order so that one color doesn't overlap all others
    set.seed(1)
    p1 <- ggplot(
        d[sample(seq_len(nrow(d))), ],
        aes(x = x, y = y)
    )
    
    ## Select adjusted effect to visualize
    if (include_geno) {
        set.seed(1)
        p2 <- ggplot(
            eff[sample(seq_len(nrow(eff))), ],
            aes(x = x, y = xGeno_effect)
        )
    } else {
        set.seed(1)
        p2 <- ggplot(
            eff[sample(seq_len(nrow(eff))), ],
            aes(x = x, y = x_effect)
        )
    }
    
    
    ############################
    ## Add points (rasterized or vector)
    ############################
    
    if (raster) {
        p1 <- p1 + geom_point_rast(aes(color = geno), size = 0.01)
        p2 <- p2 + geom_point_rast(aes(color = geno), size = 0.01)
    } else {
        p1 <- p1 + geom_point(aes(color = geno), size = 0.01)
        p2 <- p2 + geom_point(aes(color = geno), size = 0.01)
    }
    
    
    ############################
    ## Add smoothing curves and styling
    ############################
    
    ## Raw scatter plot
    p1 <- p1 +
        geom_smooth(
            alpha = 0.1,
            method = "gam",
            formula = y ~ x + s(x, bs = "ps"),
            fullrange = FALSE,
            xseq = smooth.x
        ) +
        theme_classic() +
        scale_color_manual(
            values = c(
                OR7A10st_Tumor = "gray40",
                OR7A10_Tumor   = "#E9AC4C"
            ),
            labels = c(
                OR7A10st_Tumor = "OR7A10stop",
                OR7A10_Tumor   = "OR7A10"
            ),
            name = "Genotype"
        ) +
        labs(x = x, y = y) +
        p1.label +
        theme(
            legend.key.height = unit(1, "mm"),
            legend.key.width  = unit(1, "mm"),
            legend.position   = "bottom"
        ) +
        guides(color = guide_legend(override.aes = list(size = 2)))
    
    ## Adjusted-effect scatter plot
    p2 <- p2 +
        geom_smooth(
            alpha = 0.1,
            method = "gam",
            formula = y ~ x + s(x, bs = "ps"),
            fullrange = FALSE,
            xseq = smooth.x
        ) +
        theme_classic() +
        scale_color_manual(
            values = c(
                OR7A10st_Tumor = "gray40",
                OR7A10_Tumor   = "#E9AC4C"
            ),
            labels = c(
                OR7A10st_Tumor = "OR7A10stop",
                OR7A10_Tumor   = "OR7A10"
            ),
            name = "Genotype"
        ) +
        labs(x = x, y = paste0(y, " (adj. effect)")) +
        p2.label +
        theme(
            legend.key.height = unit(1, "mm"),
            legend.key.width  = unit(1, "mm"),
            legend.position   = "bottom"
        ) +
        guides(color = guide_legend(override.aes = list(size = 2)))
    
    
    ############################
    ## Add marginal distributions
    ############################
    
    p1 <- ggExtra::ggMarginal(
        p1,
        groupColour = TRUE,
        groupFill   = TRUE,
        margins     = "both",
        size        = 4,
        alpha       = 0.2
    ) %>% cowplot::as_grob()
    
    p2 <- ggExtra::ggMarginal(
        p2,
        groupColour = TRUE,
        groupFill   = TRUE,
        margins     = "both",
        size        = 4,
        alpha       = 0.2
    ) %>% cowplot::as_grob()
    
    ## Combine plots and return
    p <- cowplot::plot_grid(
        p1,
        p2,
        align = "hv",
        axis  = "lbt",
        ncol  = 2,
        rel_widths = c(1, 1)
    )
    
    ## Print model statistics
    print(r)
    
    return(p)
}


############################
## Annotation helper for DSR scatter plots
############################

## Create correlation annotation for scatter plots
make_scatter_anno <- function(x, y, estimate, pval) {
    
    ## Determine label position
    if (sign(estimate) < 0) {
        p.df <- data.frame(
            x = max(x) - diff(range(x)) * 0.25,
            y = max(y) - diff(range(y)) * 0.05,
            label = paste0(
                "Correlation:\n",
                "rho = ", formatC(estimate, digits = 2, format = "f"),
                "\np = ", formatC(pval, digits = 2, format = "e")
            )
        )
    }
    
    if (sign(estimate) > 0) {
        p.df <- data.frame(
            x = min(x) + diff(range(x)) * 0.01,
            y = max(y) - diff(range(y)) * 0.05,
            label = paste0(
                "Correlation:\n",
                "rho = ", formatC(estimate, digits = 2, format = "f"),
                "\np = ", formatC(pval, digits = 2, format = "e")
            )
        )
    }
    
    ## Return ggplot annotation
    return(
        geom_text(
            data = p.df,
            aes(x = x, y = y, label = label),
            size = 3,
            lineheight = 0.7,
            hjust = 0
        )
    )
}

