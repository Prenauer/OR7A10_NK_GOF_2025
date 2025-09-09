#!/usr/bin/env Rscript

## Function integrates genes and cells from two different matrices into one.
combineMatrices <- function(m1, m2, suffix='-1'){
    g1 <- setdiff(rownames(m1), rownames(m2))
    g2 <- setdiff(rownames(m2), rownames(m1))
    
    m1 <- rbind(m1, matrix(0, nrow=length(g2),ncol=ncol(m1), 
                           dimnames=list(g2,colnames(m1))))
    m2 <- rbind(m2, matrix(0, nrow=length(g1),ncol=ncol(m2), 
                           dimnames=list(g1,colnames(m2))))
    m2 <- m2[rownames(m1),]
    colnames(m2) <- stringr::str_replace(colnames(m2), '-1', suffix)
    m1 <- cbind(m1,m2)
    return(m1)
}

## Function calculates leiden clustering over a range of resolutions and reports
##   the WSS and ASW for each resolution value.
DetermineOptimalClusters <- function(so, resolution=NULL, plot.only=F, graph.name='SCT_snn',
                                     plot.res=NULL, optimalResolution=NULL, seed=42){
    if(!plot.only){
        ## Optimize clustering (calculating WSS and silhouette-width)
        plan('sequential')
        xy <- so@reductions$umap@cell.embeddings %>% data.frame()
        d <- dist(xy)#Matrix::Matrix(so[[graph.name]])
        CalcWSS <- function(x,y) (length(x)-1)*(var(x)+var(y))
        optimalResolution <- do.call(rbind, lapply(resolution, function(res){
            tmp <- FindClusters(so, graph.name = graph.name, algorithm = 4, 
                                resolution = res, random.seed = seed)
            wss <- reframe(cbind(xy,cluster=Idents(tmp)), .by='cluster', 
                           wss=CalcWSS(umap_1,umap_2))$wss %>% sum()
            asw <- mean(cluster::silhouette(as.integer(Idents(tmp)), dist=d)[,3])
            cat(res,'\n')
            return(data.frame(res=res, wss=wss, asw=asw, n_clust=length(levels(Idents(tmp)))))
        }))
    }
    
    ## plot cluster-n, WSS, and ASW as a function of resolution
    if(is.null(plot.res)){
        res <- optimalResolution[with(optimalResolution, 
                                      order(n_clust,-asw,wss)),]
        r <- reframe(res, .by='n_clust', 
                     asw.pick=head(res[which(asw==max(asw))],1))[,2] %>% unlist()
    }
    else{
        r <- plot.res
    }
    
    p.resCheck <- cowplot::plot_grid(ggplot(optimalResolution, aes(x=res, y=n_clust)) + 
                                         geom_point(color='tomato2') + 
                                         geom_path(color='tomato2') + 
                                         geom_vline(xintercept=r, linetype='dashed', 
                                                    color='gray25', alpha=0.75) +
                                         labs(x='Clustering resolution', 
                                              y=stringr::str_wrap('# of clusters',20)) + 
                                         theme_test(),
                                     ggplot(optimalResolution, aes(x=res, y=wss)) + 
                                         geom_point(color='cyan2') + 
                                         geom_path(color='cyan2') + 
                                         geom_vline(xintercept=r, linetype='dashed', 
                                                    color='gray25', alpha=0.75) +
                                         labs(x='Clustering resolution', 
                                              y=stringr::str_wrap('Within-cluster sum-of-squares (wss)',20)) +
                                         theme_test(),
                                     ggplot(optimalResolution, aes(x=res, y=asw)) + 
                                         geom_point(color='orange2') + 
                                         geom_path(color='orange2') + 
                                         geom_vline(xintercept=r, linetype='dashed', 
                                                    color='gray25', alpha=0.75) +
                                         labs(x='Clustering resolution', 
                                              y=stringr::str_wrap('Average silhouette width',20)) + 
                                         theme_test(),
                                     ncol=1, align='hv')
    cat('Optimal res = ', paste0(sort(r),collapse=','),'\n')
    
    return(list(res=optimalResolution, plot=p.resCheck))
}


## Function generates a customized Dot Plot
customDotPlot <- function(so, genes, group.order=NULL, 
                          group.by='ct', scale.max=60){
    ## Reorder the cells according to the provided order
    if(!is.null(group.order)){
        so@meta.data[,group.by] <- factor(so@meta.data[,group.by], 
                                          levels=group.order)
    }
    ## Make Dotplot of the provided genes
    DotPlot(so,genes, 
            group.by=group.by, scale.max=scale.max) + 
        # remove labels
        labs(x=NULL, y=NULL) +
        # customize legend presentation
        theme(legend.position='bottom', legend.key.spacing=unit(0.5,'lines'), 
              legend.key.height=unit(0.5,'lines'), 
              legend.key.width=unit(0.5,'lines'), 
              legend.text=element_text(size=7), 
              legend.key.spacing.y=unit(0.01,'lines'),
              axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
              legend.title=element_text(size=10)) +
        # change dot color-scheme
        scale_color_viridis_c(option='inferno', direction=-1, end=0.9) +
        # set order of which lengend is shown first
        guides(size = guide_legend(nrow = 2), alpha = guide_legend(nrow = 2))
}

## Function makes a detailed volcano plot from DE results
plotVolcano <- function(de, lfc.thresh=1, p.thresh=0.01, pt.size=2, title=NULL,
                        text.size=2.2, rastr=T, xlim.expansion=1.3){
    de$p <- de[,'padj']
    de$lfc <- de[,'logFC']
    de$feature <- de[,'feature']
    de <- de[order(de$pval),]
    genes2label <- c(head(de$feature[which((de$lfc > lfc.thresh)&(de$p < p.thresh))],5),
                     head(de$feature[which((de$lfc < -lfc.thresh)&(de$p < p.thresh))],5))
    de <- de[order(abs(de$lfc), decreasing=T),]
    genes2label <- c(genes2label,
                     head(de$feature[which((de$lfc > lfc.thresh)&(de$p < p.thresh))],5),
                     head(de$feature[which((de$lfc < -lfc.thresh)&(de$p < p.thresh))],5))
    genes2label <- unique(genes2label)
    
    # set min pval limit for plotting
    plot.min.p <- min(de$p[which(de$p != 0)], na.rm=T)*0.1
    de$p[which(de$p == 0)] <- plot.min.p
    
    # get plot size limits
    plot.xlim <- xlim.expansion*range(de$lfc, na.rm=T)
    plot.ylim <- max(-log10(de$p), na.rm=T)*c(-0.02,1.2)
    
    p <- ggplot(de, aes(x=lfc, y=-log10(p))) 
    if(rastr){
        p <- p + geom_point_rast(color='gray20', alpha=1, stroke=0, size=pt.size*1.5) +
            geom_point_rast(color='gray', alpha=1, stroke=0, size=pt.size) +
            geom_point_rast(data=de[which(de$lfc > lfc.thresh & de$p < p.thresh),],
                            color='firebrick', alpha=1, stroke=0, size=pt.size) +
            geom_point_rast(data=de[which(de$lfc < -lfc.thresh & de$p < p.thresh),],
                            color='steelblue', alpha=1, stroke=0, size=pt.size)
    }
    if(!rastr){
        p <- p + geom_point(color='gray', alpha=0.5, stroke=0, size=pt.size) +
            geom_point(data=de[which(de$lfc > lfc.thresh & de$p < p.thresh),],
                       color='firebrick', alpha=0.8, stroke=0, size=pt.size) +
            geom_point(data=de[which(de$lfc < -lfc.thresh & de$p < p.thresh),],
                       color='steelblue', alpha=0.8, stroke=0, size=pt.size)
    }
    p <- p + ggrepel::geom_text_repel(data=de[which(de$lfc < 0 & de$feature %in% genes2label),],
                                      aes(label=feature), min.segment.length=0,
                                      nudge_x=-0.5, nudge_y=0.5, max.iter = 1e5,
                                      segment.alpha=0.7, segment.size=0.01,
                                      size=text.size, max.overlaps=50, fontface='italic') +
        ggrepel::geom_text_repel(data=de[which(de$lfc > 0 & de$feature %in% genes2label),],
                                 aes(label=feature), min.segment.length=0,
                                 nudge_x=0.5, nudge_y=0.5, max.iter = 1e5,
                                 segment.alpha=0.7, segment.size=0.01,
                                 size=text.size, max.overlaps=50, fontface='italic') +
        xlim(plot.xlim) + ylim(plot.ylim) +
        theme_classic() + labs(x='Log2 fold-change', y='Adj. p (-log10)', title=title)
    return(p)
}

## Function makes a network plot from pathway analysis results
RunNetPW <- function(pw, plot.title=NULL, n.pw.per.group=5, seed=42){
    # Ensure required packages are loaded 
    require(ggraph)
    require(igraph)
    require(ggnetwork)
    require(ggrepel)
    
    # setup graph
    top.pw <- slice_head(pw, by='group', n=n.pw.per.group)$pathway
    pw <- pw[which(pw$pathway %in% top.pw),]
    e <- rbind(
        data.frame(from=pw$ct, to=pw$group, size=0.1),
        data.frame(from=pw$group, to=pw$pathway, size=pw$NES/max(pw$NES, T)))
    v <- rbind(data.frame(name=pw$pathway,type='pathway'),
               data.frame(name=pw$ct,type='celltype'),
               data.frame(name=pw$group,type='comparison'))
    v <- unique(v)
    
    ig <- graph_from_data_frame( e, vertices=v )
    set.seed(seed)
    ig2 <- (ig %>% add_layout_(with_fr()))
    net <- ggnetwork(ig2)
    net$class <- 'edge'
    net$class[which(net$x==net$xend & net$y==net$yend)] <- 'node'
    
    # Reformat pathway names
    net$name[which(net$type=='pathway')] <- net$name[which(net$type=='pathway')] %>%
        str_remove(.,'GOBP_') %>% str_remove(.,'REACTOME_') %>%
        str_remove(.,'KEGG_') %>% str_remove(.,'WP_') %>% 
        str_remove(.,'BIOCARTA_') %>% str_replace_all(., '_', ' ')
    
    ## exclude celltype nodes and edges from the visualization
    net <- net[which(net$type != 'celltype'),]
    
    # adjust names of 'comparison' nodes
    net$node_label=NA
    net$node_label[which(net$class=='node')] <- 
        str_split_i(net$name[which(net$class=='node')], ': ',1) 
    
    ## Make plot
    cols <- structure(c("#F8766D","#CD9600","#7CAE00","#00BE67","#00BFC4","#00A9FF","#C77CFF","#FF61CC"),
                      names=c("iNK_c1-Cytotoxic","iNK_c2-CREM","mNK_c1-MKI67","mNK_c2-MKI67-BATF3","mNK_c3-BATF3","mNK_c4-IL32","mNK_c5-KLRC2","mNK_c6-EGR2"))
    net$gp2 <- str_split_i(net$name, ':', 1)   
    net <- mutate(net, .by='node_label', nodesize=length(node_label))
    
    ## calculate pathway node size relative to the number of connecting edges
    node_size <- reframe(pw, .by='feature', n=length(feature))
    node_size <- structure(node_size$n, 
                           names=str_replace_all(node_size$feature,'-',' '))
    net$node_size <- node_size[net$gp2]
    
    ## Fix label names
    net[which(net$class=='node' & net$type=='pathway'), 'name'] <- 
        str_remove_all(net[which(net$class=='node' & 
                                     net$type=='pathway'), 'name'], 
                       'PID ')
    net[which(net$class=='node' & net$type=='comparison'),'node_label'] <- 
        str_split_i(net[which(net$class=='node' & 
                                  net$type=='comparison'),'node_label'], '-',1)
    
    ggplot(net, aes(x = x, y = y)) + 
        geom_edges(color="#DDDDDD",
                   aes(linewidth=size, xend = xend, yend = yend)) +
        geom_edges(data=net[which(net$class=='edge' & net$type=='comparison'),], 
                   arrow = arrow(length = unit(0.25, 'lines'), 
                                 type = "closed"), #curvature=0.4, 
                   aes(linewidth=size, color=gp2, xend = xend, yend = yend)) +
        scale_linewidth(range=c(0.1,1)) +
        theme_blank() +
        geom_nodes(data=net[which(net$class=='node' & net$type=='comparison'),], 
                   shape=22, color='gray30', size=6,
                   aes(fill=gp2)) + 
        geom_nodes(data=net[which(net$class=='node' & net$type=='pathway'),], 
                   shape=21, color='gray30',fill='#20854EFF', 
                   aes(size=node_size)) + 
        scale_size_identity() + 
        geom_nodetext_repel(data=net[which(net$class=='node' & net$type=='comparison'),], 
                            aes(label=str_wrap(node_label,20)), size=3, 
                            min.segment.length=0, nudge_y=0.5, 
                            segment.size=0.1, segment.alpha=0.5, max.overlaps=50) + 
        geom_nodetext_repel(data=net[which(net$class=='node' & net$type=='pathway'),], 
                            aes(label=str_wrap(name,20)), size=3, min.segment.length=0, 
                            nudge_y=0.5, lineheight=0.8,fontface='bold',
                            segment.size=0.1, segment.alpha=0.5, max.overlaps=50) + 
        scale_fill_manual(values=cols) + 
        scale_color_manual(values=cols) + 
        labs(title=plot.title, fill='Subset', 
             linewidth=str_wrap('Pathway enrichment score',8)) + 
        guides(fill=guide_legend(order = 1),linewidth=guide_legend(order=2), 
               color=F) +
        theme(legend.key.height=unit(0.2,'lines'), legend.spacing.x=unit(0.2,'lines')) 
}

## Function makes a hclust dot plot with a dendrogram 
MakeDotPlot <- function(g, title='', thresh=1/5, method='ward.D2', 
                        cluster.axis='xy'){
    require(ggdendroplot)
    g <- intersect(hvg, g)
    m <- t(so$RNA$data[g,])
    m <- S4Vectors::splitAsList(m, factor(so$id, levels=unique(so$id)))
    m <- do.call(rbind, lapply(m, colMeans))
    rownames(m) <- str_replace(rownames(m), '_','-')
    res.pca <- prcomp(m, scale=T)
    g <- g[(apply(so$RNA$data[g,], 1, PercentAbove, threshold=0) > thresh)]
    m <- t(so$RNA$data[g,])
    m <- S4Vectors::splitAsList(m, factor(so$id, levels=unique(so$id)))
    #m <- do.call(rbind, lapply(m, function(x){apply(x, 2, PercentAbove, threshold=0)}))
    m <- do.call(rbind, lapply(m, colMeans)) %>% t()
    #rownames(m) <- unique(so$id)
    # Set default order
    cl.x <- g
    cl.y <- c(paste0(as.character(sort(unique(so$ct))),':OR7A10'), 
              paste0(as.character(sort(unique(so$ct))),':OR7A10st'))
    if(cluster.axis %in% c('y','xy')) {
        cl.y <- hclust(dist(t(m)), method=method)
    }
    if(cluster.axis %in% c('x','xy')) cl.x <- hclust(dist(m), method=method)
    so$id <- factor(so$id, levels=cl.y$labels[cl.y$order])
    axis.color <- c('#D98500','gray30')[1+as.integer(grepl('OR7A10st',levels(so$id)))]
    g <- g[cl.x$order]
    p1 <- DotPlot(so,g,group.by='id') + #coord_flip() + 
        labs(x=NULL, y=NULL, title=title) +scale_color_distiller(palette = "RdYlBu") +
        guides(size = guide_legend(nrow = 1), alpha = guide_legend(nrow = 2)) +
        theme(legend.position='bottom', legend.key.spacing=unit(0.5,'lines'), 
              legend.key.height=unit(0.5,'lines'), legend.key.width=unit(0.5,'lines'), 
              legend.text=element_text(size=7), legend.key.spacing.y=unit(0.01,'lines'),
              axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face='italic'),
              legend.title=element_text(size=10),
              theme(axis.text.y = ggtext::element_markdown(colour ='#D98500')))
    p1.dend <- ggplot() + geom_dendro(cl.y, pointing='side') + theme_void()
    p1 <- cowplot::plot_grid(p1,p1.dend, align='h', axis='lbt',ncol=2,
                             rel_widths=c(0.6,0.2))
    d <- data.frame(res.pca$x[,1:2], id=rownames(res.pca$x),
                    ct=substr(rownames(res.pca$x),0,6),
                    geno=str_split_i(rownames(res.pca$x),':',2))
    p2 <- ggplot(d, aes(x=PC1, y=PC2, fill=geno), color='gray20') + 
        geom_point(size=2.5, shape=21, stroke=0.01) + 
        xlim(1.3*range(d$PC1)) + 
        theme_classic() + theme(legend.position='none') + 
        scale_fill_manual(values=structure(c('#B2B2B2','#E9AC4C'),
                                           names=c('OR7A10st','OR7A10'))) +
        scale_color_manual(values=structure(c('gray30','#D98500'),
                                            names=c('OR7A10st','OR7A10'))) +
        ggrepel::geom_text_repel(aes(label=ct, color=geno), size=2.5, min.segment.length=0, segment.size=0.2) + 
        labs(title=title)
    p <- list(p1,p2)
    return(p)
}

## Function processes a matrix so that there is a min and max cutoff
winsorize <- function(m, thresh){
    apply(m, 1, function(x){
        x[which(x > thresh)] <- thresh
        x[which(x < -thresh)] <- -thresh
        return(x)
    }) %>% t() %>% data.frame()
}


###########
###########
## DSR functions


## Function fits a GAM for signature-signature or signature-gene combinations.
dsrFitGam <- function(d, use.nb=NULL, use.umi=NULL, f=NULL, 
                      filter.data=T, select=T){
    require(mgcv)
    if(filter.data) d <- d[which(d$x > 0 & d$y > 0),]
    
    if(sum(c('x.type', 'y.type') %in% names(attributes(d))) == 2 &
       is.null(use.nb) & is.null(use.umi)){
        if(attr(d, 'x.type')=='counts') use.umi <- T
        if(attr(d, 'y.type')=='counts'){
            use.umi <- T
            use.nb <- T
        } 
    }
    if(is.null(use.nb) & is.null(use.umi)){
        use.nb=F
        use.umi=F
    }
    
    if(use.nb){
        if(is.null(f)) f <- formula(y ~ geno + ct + s(x, bs="ps",k=10) + offset(umi))
        fit <- try(mgcv::bam(f,  select=select,  control=list(keepData=T),
                             data=d, family='nb', method='fREML'), silent=T)
    }
    if(!use.nb){
        if(use.umi){
            f <- formula(y ~ geno + ct + s(x, bs="ps",k=10) + offset(umi))
        }
        if(!use.umi){
            f <- formula(y ~ geno + ct + s(x, bs="ps",k=10))
        }
        if(is.null(f)) f <- formula(y ~ geno + ct + s(x, bs="ps",k=10))
        fit <- try(mgcv::bam(f,  select=select,  control=list(keepData=F),
                             data=d, Gamma(link='log'), method='fREML'), silent=T)
        
    }
    if(attempt::is_try_error(fit)) return(NULL)
    
    ## Add order of cell barcodes to fit object
    fit$order <- rownames(d)
    return(fit)
    
}

## Function takes dsrFitGam fit and extracts model parameters. This also
##  calculates correlation between the predictor and partial residuals.
extractDsrResults <- function(fit){
    require(mgcv)
    d <- fit$model
    #pr <- (predict(fit, newdata=fit$model, type='terms'))
    pr <- (predict(fit, type='terms'))
    if(nrow(pr) < 20) return(NULL)
    d$p_x <- (pr[,'s(x)'] + residuals(fit))[fit$order]
    d$p_gx <- (rowSums(pr[,c('geno', 's(x)')]) + residuals(fit))[fit$order]
    ctx <- cor.test(d$x,d$p_x)[c('p.value','estimate')]
    ctxg <- cor.test(d$x,d$p_gx)[c('p.value','estimate')]
    genoEst <- summary(fit)$p.table['genoOR7A10_Tumor',c(4,1)]
    ct <- structure(c(ctx,ctxg, genoEst), names=c('cor_xy_p','cor_xy_r2','cor_xgy_p','cor_xgy_r2','geno_p','geno_est'))
    r <- structure(c(summary(fit)$s.table[1,3:4],summary(fit)$r.sq, summary(fit)$dev.expl) , names=c('fit_F','fit_p','fit_r2','fit_dev'))
    res <- matrix(unlist(c(ct, r)), nrow=1) %>% data.frame()
    colnames(res) <- c(names(ct), names(r))
    #res[,3:ncol(res)] <- apply(res[,3:ncol(res)], 2, as.numeric)
    
    return(res)
    
}

## Function takes dsrFitGam fit and extracts predictor and partial residuals.
getDsrEffects <- function(fit){
    require(mgcv)
    d <- fit$model
    pr <- cbind(predict(fit, type='terms'),resid=residuals(fit))[fit$order,]
    output <- data.frame(d, x_effect=rowSums(pr[,c('s(x)','resid')]), 
                         geno_effect=rowSums(pr[,c('geno','resid')]),
                         xGeno_effect=rowSums(pr[,c('s(x)','geno','resid')]))
    return(output)
}


## Function creates the input for dsrFitGam using the names of the predictor 
##  and response. 
setupDsrData <- function(x, y, force.x.counts=F, force.y.norm=F){
    assays <- Assays(so)
    a.x <- assays[sapply(assays, function(a) (x %in% rownames(so[[a]])))][1]
    a.y <- assays[sapply(assays, function(a) (y %in% rownames(so[[a]])))][1]
    
    # check if in metadata
    if(is.na(a.x)) if(x %in% colnames(so@meta.data)) a.x <- 'metadata'
    if(is.na(a.y)) if(y %in% colnames(so@meta.data)) a.y <- 'metadata'
    
    if(a.y=='metadata') y.data <- so@meta.data[,y]
    if(a.x=='metadata') x.data <- so@meta.data[,x]
    if(a.y=='RNA') y.data <- so[[a.y]]$counts[y,]
    if(!(a.y %in% c('RNA','metadata'))) y.data <- so[[a.y]]$data[y,]
    if(force.x.counts & a.x!='metadata') x.data <- so[[a.x]]$counts[x,]
    if(!force.x.counts & a.x!='metadata') x.data <- so[[a.x]]$data[x,]
    if(a.y != 'metadata' & force.y.norm) y.data <- so[[a.y]]$data[y,]
    
    d <- data.frame(x=x.data, y=y.data, 
                    ct=factor(so$ct), umi=log(so$nCount_RNA),
                    geno=factor(so$sample, 
                                levels=c('OR7A10st_Tumor', 'OR7A10_Tumor')))
    attr(d, 'x.type') <- c('data','counts')[1+as.integer(force.x.counts)]
    attr(d, 'y.type') <- c('data','counts')[1+as.integer(a.y=='RNA')]
    return(d)
    
}


## Function takes names of the predictor and response, analyzes DSR, then 
##    makes a scatter plot of predictor vs adjusted effect with stats printed
##    in the upper left side. 
testDsr <- function(x,y, raster=T, include_geno=F){
    d=setupDsrData(x,y)
    
    fit=dsrFitGam(d,select=F)
    r=extractDsrResults(fit)
    eff <- getDsrEffects(fit)
    
    cat('Model: ',x, ' as a predictor of ', y,'\n')

    ## setup basic scatter plot
    d=setupDsrData(x,y, force.y.norm=T)
    d <- d[which(d$x > 0 & d$y > 0),]
    smooth.x <- seq(quantile(d$x, 0.01),quantile(d$x, 0.99), 
                    diff(quantile(d$x, c(0.01,0.99)))/100)
    
    
    ## Stat labels
    p1.stats <- cor.test(d$x, d$y)[c('estimate', 'p.value')] %>% unlist()
    p1.label <- make_scatter_anno(d$x, d$y, p1.stats[1], p1.stats[2])
    p2.label <- make_scatter_anno(eff$x, eff$x_effect, r$cor_xy_r2, r$cor_xy_p)
    
    
    
    set.seed(1)
    p1 <- ggplot(d[sample(1:nrow(d), nrow(d)),], aes(x=x,y=y)) 
    
    if(include_geno) {
        set.seed(1)
        p2 <- ggplot(eff[sample(1:nrow(eff), nrow(eff)),], 
                     aes(x=x,y=xGeno_effect)) 
    }
    if(!include_geno) {
        set.seed(1)
        p2 <- ggplot(eff[sample(1:nrow(eff), nrow(eff)),], aes(x=x,y=x_effect)) 
    }
    
    
    if(raster){
        p1 <- p1 + geom_point_rast(aes(color=geno), size=0.01)
        p2 <- p2 + geom_point_rast(aes(color=geno), size=0.01)
    }
    if(!raster){
        p1 <- p1 + geom_point(aes(color=geno), size=0.01)
        p2 <- p2 + geom_point(aes(color=geno), size=0.01)
    }
    p1 <- p1 +  geom_smooth(alpha=0.1, 
                            method='gam', formula=y~ x + s(x, bs="ps"), 
                            fullrange=F,xseq=smooth.x) + theme_classic() + 
        scale_color_manual(values=structure(c('gray40','#E9AC4C'),
                                           names=c("OR7A10st_Tumor",
                                                   "OR7A10_Tumor")),
                           labels=structure(c("OR7A10stop",
                                              "OR7A10"),
                                            names=c("OR7A10st_Tumor",
                                                    "OR7A10_Tumor")),
                           name='Genotype') +
        labs(x=x,y=y) + p1.label + 
        theme(legend.key.height=unit(1, 'mm'), legend.key.width=unit(1,'mm'),
              legend.position='bottom') + 
        guides(color = guide_legend(override.aes = list(size = 2)))
        
    p2 <- p2 + geom_smooth(alpha=0.1, 
                           method='gam', formula=y~ x + s(x, bs="ps"), 
                           fullrange=F,xseq=smooth.x) + theme_classic() + 
        scale_color_manual(values=structure(c('gray40','#E9AC4C'),
                                           names=c("OR7A10st_Tumor",
                                                   "OR7A10_Tumor")),
                           labels=structure(c("OR7A10stop",
                                              "OR7A10"),
                                            names=c("OR7A10st_Tumor",
                                                    "OR7A10_Tumor")),
                           name='Genotype') +
        labs(x=x,y=paste0(y, ' (adj. effect)')) + p2.label + 
        theme(legend.key.height=unit(1, 'mm'), legend.key.width=unit(1,'mm'),
              legend.position='bottom') + 
        guides(color = guide_legend(override.aes = list(size = 2)))
    
    p1 <- ggExtra::ggMarginal(p1, groupColour=T, margins='both', groupFill=T, 
                              size=4, alpha=0.2) %>% cowplot::as_grob()
    p2 <- ggExtra::ggMarginal(p2, groupColour=T, margins='both', groupFill=T, 
                              size=4, alpha=0.2) %>% cowplot::as_grob()
    
    p <- cowplot::plot_grid(p1, p2,
                       align='hv', axis='lbt',ncol=2, rel_widths=c(1,1))
    print(r)
    
    return(p)
}

## Function is an accessory for 'testDsr' plotting that preps labels for stats
make_scatter_anno <- function(x, y, estimate, pval){
    if(sign(estimate) < 0){
        p.label <- data.frame(
            x=max(x)-(diff(range(x))*0.25),
            y=max(y)-(diff(range(y))*0.05),
            label=paste0('Correlation:\nrho = ',
                         formatC(estimate, digits=2, format='f'), 
                         '\np = ',formatC(pval, digits=2, format='e'))
        )
    }    
    if(sign(estimate) > 0){
        p.label <- data.frame(
            x=min(x)+(diff(range(x))*0.01),
            y=max(y)-(diff(range(y))*0.05),
            label=paste0('Correlation:\nrho = ',
                         formatC(estimate, digits=2, format='f'), 
                         '\np = ',formatC(pval, digits=2, format='e'))
        )
    }
    p.label <- geom_text(data=p.label, aes(x=x,y=y, label=label), size=3,
                         lineheight=0.7, hjust = 0)
    return(p.label)
    
}    

