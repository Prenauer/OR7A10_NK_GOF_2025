###This script analyze the RNA-seq data with the DESeq2 packages###
library(dplyr)
library(tidyverse)
library(purrr)
library(stringr)
library(data.table)
library(ggrepel)
library(pheatmap)

#BiocManager::install("DESeq2")
library(DESeq2)
setwd("~/Yale University Dropbox/Kaiyuan Tang/Sidi_Chen_Lab/Project/NKGOF/RNAseq")
count_raw <- fread(file = "NKGOF_RNAseq_Count_Matrix.txt") %>%
  column_to_rownames(var = "gene")#set the gene to rownames so that it can be read by edgeR

#filter out unexpressed genes across all samples
count <- count_raw[rowSums(count_raw[, 1:36] > 20) > 0, ] 
count_filter <- count %>% #a subset version of only the stimulated population
  select("77-24_Rep1", "77-24_Rep2", "77-24_Rep3", "77b-24_Rep1", "77b-24_Rep2","77b-24_Rep3", "79-24_Rep1","79-24_Rep2", "79-24_Rep3", "79b-24_Rep1", "79b-24_Rep2", "79b-24_Rep3", "90-24_Rep1", "90-24_Rep2", "90-24_Rep3", "90b-24_Rep1", "90b-24_Rep2", "90b-24_Rep3")

#1. Readin with DESeq2
#generate metadata for Metastasis data matrix 
sample_name <- colnames(count)
sample_name1 <- colnames(count_filter) #stimulated subset

condition <- c(rep("HER2CAR_OR7A10_0hr", 3), rep("HER2CAR_OR7A10_24hr", 3), rep("HER2CAR_OR7A10Stop_0hr", 3), rep("HER2CAR_OR7A10Stop_24hr", 3), rep("OR7A10_0hr", 3), rep("OR7A10_24hr", 3), rep("OR7A10Stop_0hr", 3), rep("OR7A10Stop_24hr", 3), rep("tHER2CAR_OR7A10_0hr", 3), rep("tHER2CAR_OR7A10_24hr", 3), rep("tHER2CAR_OR7A10Stop_0hr", 3), rep("tHER2CAR_OR7A10Stop_24hr", 3))

condition1 <- c(rep("HER2CAR_OR7A10_24hr", 3), rep("HER2CAR_OR7A10Stop_24hr", 3), rep("OR7A10_24hr", 3), rep("OR7A10Stop_24hr", 3), rep("tHER2CAR_OR7A10_24hr", 3), rep("tHER2CAR_OR7A10Stop_24hr", 3)) #stimulated subset

stim <- c(rep("0hr", 3), rep("24hr", 3), rep("0hr", 3), rep("24hr", 3), rep("0hr", 3), rep("24hr", 3), rep("0hr", 3), rep("24hr", 3), rep("0hr", 3), rep("24hr", 3), rep("0hr", 3), rep("24hr", 3))
CAR <- c(rep("HER2CAR", 12), rep("NoCAR", 12), rep("tHER2CAR", 12))
OE <- c(rep("OR7A10", 6), rep("OR7A10Stop", 6), rep("OR7A10", 6), rep("OR7A10Stop", 6), rep("OR7A10", 6), rep("OR7A10Stop", 6))

#stimulated subset
stim1 <- c(rep("24hr", 3), rep("24hr", 3), rep("24hr", 3), rep("24hr", 3), rep("24hr", 3), rep("24hr", 3))
CAR1 <- c(rep("HER2CAR", 6), rep("NoCAR", 6), rep("tHER2CAR", 6))
OE1 <- c(rep("OR7A10", 3), rep("OR7A10Stop", 3), rep("OR7A10", 3), rep("OR7A10Stop", 3), rep("OR7A10", 3), rep("OR7A10Stop", 3))


meta <- data.frame(cbind(sample_name, condition, stim, CAR, OE)) %>%
  column_to_rownames(var = "sample_name")

meta1 <- data.frame(cbind(sample_name1, condition1, stim1, CAR1, OE1)) %>%
  column_to_rownames(var = "sample_name1")

meta$condition <- as.factor(meta$condition)
meta$stim <- as.factor(meta$stim)
meta$CAR <- as.factor(meta$CAR)
meta$OE <- as.factor(meta$OE)

meta1$condition1 <- as.factor(meta1$condition1)
meta1$stim1 <- as.factor(meta1$stim1)
meta1$CAR1 <- as.factor(meta1$CAR1)
meta1$OE1 <- as.factor(meta1$OE1)

#create DESeq2 object
dds_inter <- DESeqDataSetFromMatrix(countData = count,
                                    colData = meta,
                                    design = ~ OE * CAR * stim) #Full model with interaction between variables

dds_filter <- DESeqDataSetFromMatrix(countData = count_filter,
                                       colData = meta1,
                                       design = ~ OE1*CAR1) #Full model with just the stimulation population


#2. Normalization 
dds_inter <- estimateSizeFactors(dds_inter)
sizeFactors(dds_inter)

dds_filter <- estimateSizeFactors(dds_filter)
sizeFactors(dds_filter)
NKGOF_normalized_count1 <- counts(dds_filter, normalized = T)

#3. QC with correlation heatmap and PCA
#correlation heatmap 
vsd <- vst(dds_inter, blind = T) #log transformation 
vsd_mat <- assay(vsd) #extract log transformed matrix
vsd_cor <- cor(vsd_mat) #pairwise correlation between samples 

vsd1 <- vst(dds_filter, blind = T) #log transformation 
vsd_mat1 <- assay(vsd1) #extract log transformed matrix
vsd_cor1 <- cor(vsd_mat1) #pairwise correlation between samples 

#install.packages("pheatmap")
library(pheatmap)
pheatmap(vsd_cor, annotation = dplyr::select(meta, condition))
pheatmap(vsd_cor1, annotation = dplyr::select(meta, condition))

library(paletteer)
library(ggsci)
#PCA plot
pca_stim <- plotPCA(vsd, intgroup="stim") + 
  theme_classic()+
  #scale_color_paletteer_d("nationalparkcolors::Acadia")+
  scale_color_npg()+
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + 
  labs(title = "PCA Plot by Stimulation")

pca_OE <- plotPCA(vsd, intgroup="OE") + 
  theme_classic()+
  scale_color_npg()+
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + 
  labs(title = "PCA Plot by OE")

pca_CAR <- plotPCA(vsd, intgroup="CAR") + 
  theme_classic()+
  scale_color_npg()+
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + 
  labs(title = "PCA Plot by CAR Construct")

subset <- plotPCA(vsd1, intgroup="OE1") + 
  theme_classic()+
  scale_color_npg()+
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + 
  labs(title = "PCA Plot by OE")

ggsave("PCA_stim.pdf", plot = pca_stim, width = 4, height = 4)
ggsave("PCA_OE.pdf", plot = pca_OE, width = 4, height = 4)
ggsave("PCA_CAR.pdf", plot = pca_CAR, width = 4, height = 4)
ggsave("PCA_OE_subset.pdf", plot = subset, width = 4, height = 4)

#4. Differential expression analysis 
#setting the controls
dds_inter$CAR <- relevel(dds_inter$CAR, ref = "tHER2CAR")
dds_inter$stim <- relevel(dds_inter$stim, ref = "0hr")
dds_inter$OE <- relevel(dds_inter$OE, ref = "OR7A10Stop")
dds_inter <- DESeq(dds_inter)
resultsNames(dds_inter)

#setting the controls
dds_filter$CAR1 <- relevel(dds_filter$CAR1, ref = "tHER2CAR")
dds_filter$OE1 <- relevel(dds_filter$OE1, ref = "OR7A10Stop")
dds_filter <- DESeq(dds_filter)
resultsNames(dds_filter)

#custom function to curate and filter the data 
result_curation <- function(res, lfc.thresh=0.5, p.thresh=0.05) {
  #order
  res_ordered <- res[order(res$pvalue),]
  #convert to df
  res_df <- as.data.frame(res_ordered) %>%
    filter(complete.cases(.)) #remove NA in the dataframe
  #filter
  res_df$significant <- 'Not Significant'
  res_df[res_df$padj < p.thresh & res_df$log2FoldChange > lfc.thresh,]$significant <- "Up"
  res_df[res_df$padj < p.thresh & res_df$log2FoldChange < -lfc.thresh,]$significant <- "Down"
  
  return(res_df)
}

#Get the lfc shrinked results 
#main effect
resLFC_stim <- lfcShrink(dds_inter, coef="stim_24hr_vs_0hr", type="apeglm")
resLFC_CAR <- lfcShrink(dds_inter, coef="CAR_HER2CAR_vs_tHER2CAR", type="apeglm")
resLFC_OE <- lfcShrink(dds_inter, coef="OE_OR7A10_vs_OR7A10Stop", type="apeglm")

#Stim subset
resLFC_OE1 <- lfcShrink(dds_filter, coef="OE1_OR7A10_vs_OR7A10Stop", type="apeglm")
resLFC_OE1OR7A10.CAR1HER2CAR <-lfcShrink(dds_filter, coef="OE1OR7A10.CAR1HER2CAR", type="apeglm")

#data curation
resLFC_stim_df <- result_curation(resLFC_stim)

resLFC_OE1_df <- result_curation(resLFC_OE1)
resLFC_OE1OR7A10.CAR1HER2CAR_df <- result_curation(resLFC_OE1OR7A10.CAR1HER2CAR) 

#Custom function for creating volcano plot
library(ggrastr)
# Create the volcano plot
DESeq2_volcano <- function(df, pt.size=2, title=NULL,lfc.thresh=0.5, p.thresh=0.05,
                           text.size=2.2, xlim.expansion=1.3, head5=T) {
  # --- 1. Data Preparation ---
  # Add gene names from rownames to a new 'feature' column
  df$feature <- rownames(df)
  
  # Standardize column names for plotting
  de <- df %>%
    dplyr::rename(p = padj, lfc = log2FoldChange)

  # Choose between rasterized and vector points
  if (head5 == T) {
    # This logic prioritizes genes by both p-value and fold-change
    de <- de[order(de$p), ]
    genes2label <- c(head(de$feature[which((de$lfc > lfc.thresh) & (de$p < p.thresh))], 5),
                     head(de$feature[which((de$lfc < -lfc.thresh) & (de$p < p.thresh))], 5))
    
    de <- de[order(abs(de$lfc), decreasing=T), ]
    genes2label <- c(genes2label,
                     head(de$feature[which((de$lfc > lfc.thresh) & (de$p < p.thresh))], 5),
                     head(de$feature[which((de$lfc < -lfc.thresh) & (de$p < p.thresh))], 5))
    genes2label <- unique(genes2label)
  } 

  # --- 3. Handle P-values of Zero ---
  # Replace p=0 with a tiny value for log transformation
  if (any(de$p == 0, na.rm=T)) {
    plot.min.p <- min(de$p[which(de$p != 0)], na.rm=T) * 0.1
    de$p[which(de$p == 0)] <- plot.min.p
  }
  
  # --- 4. Plotting ---
  # Get plot size limits
  plot.xlim <- xlim.expansion * range(de$lfc, na.rm=T)
  plot.ylim <- max(-log10(de$p), na.rm=T) * c(-0.02, 1.2)
  
  p <- ggplot(de, aes(x=lfc, y=-log10(p), color=significant))
  p <- p + geom_point_rast(alpha=0.8, stroke=0, size=pt.size)
  
  p <- p +
    # Define custom colors and legend labels
    scale_color_manual(name = "Significance",
                       values = c("Up" = "firebrick", "Down" = "steelblue", "Not Significant" = "grey"),
                       labels = c("Downregulated", "Not Significant", "Upregulated")) +
    
    # Add labels using ggrepel to avoid overlap
    ggrepel::geom_text_repel(data = subset(de, feature %in% genes2label),
                             aes(label = feature),
                             min.segment.length = 0,
                             nudge_x = 0.2, nudge_y = 0.2,
                             max.iter = 1e5,
                             segment.alpha = 0.7, segment.size = 0.01,
                             size = text.size,
                             max.overlaps = 50,
                             fontface = 'italic',
                             color = 'black') + # Label color
    
    geom_vline(xintercept = c(-lfc.thresh, lfc.thresh), linetype = "dashed", color = "darkgrey") +
    geom_hline(yintercept = -log10(p.thresh), linetype = "dashed", color = "darkgrey") +
    
    # Set plot limits and theme
    xlim(plot.xlim) + ylim(plot.ylim) +
    theme_classic() +
    labs(x = 'Log2 Fold-Change', y = '-Log10(Adjusted P-value)', title = title) +
    theme(legend.position = "top", plot.title = element_text(hjust = 0.6, size = 14))
  
  return(p)
}

#plot main effect
p_stim <- DESeq2_volcano(resLFC_stim_df, title = "Stim vs Unstim Main Effect")
p_CAR <- DESeq2_volcano(resLFC_CAR_df, title = "HER2CAR vs tHER2CAR Main Effect")
p_OE <- DESeq2_volcano(resLFC_OE_df, title = "OR7A10 vs OR7A10stop Main Effect")

#plot stimulation subset
p_OE1 <- DESeq2_volcano(resLFC_OE1_df, title = "Subset OR7A10 vs OR7A10stop Main Effect")
p_OE1OR7A10.CAR1HER2CAR <- DESeq2_volcano(resLFC_OE1OR7A10.CAR1HER2CAR_df, title = "Subset Interaction OE1OR7A10.CAR1HER2CAR")


#save ggplot
ggsave("Volcano Stim vs Unstim Main Effect.pdf", plot = p_stim, width = 5, height = 5)
ggsave("Volcano OR7A10 vs OR7A10stop Main Effect subset.pdf", plot = p_OE1, width = 5, height = 5)
ggsave("Volcano Interaction OE1OR7A10.CAR1HER2CAR.pdf", plot = p_OE1OR7A10.CAR1HER2CAR, width = 5, height = 5)

#Plot interaction plot 
# --- Step 2: Filter for Your Gene and Prepare for Plotting ---
interaction_metadata <- meta1 %>%
  rownames_to_column(var = "sample_id")
interaction_long <- NKGOF_normalized_count1 %>%
  as.data.frame() %>%
  # Filter for only the gene of interest
  filter(row.names(.) == "OR7A10") %>%
  dplyr::select(interaction_metadata$sample_id) %>%
  # Convert from wide format to long format
  tibble::rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample_id",
    values_to = "expression"
  ) %>% 
  # Add the experimental conditions by joining with the metadata
  left_join(interaction_metadata, by = "sample_id") %>%
  mutate(CAR1 = factor(CAR1, levels = c("NoCAR", "tHER2CAR", "HER2CAR"))) %>%
  group_by(OE1, CAR1) %>%
  mutate(avg = mean(expression))

# --- Step 3: Calculate Summary Statistics ---
# This part is the same as before.
summary_for_plot <- interaction_long %>%
  group_by(CAR1, OE1) %>%
  summarise(
    mean_expr = mean(expression),
    se = sd(expression) / sqrt(n())
  ) %>%
  ungroup() %>%
  mutate(CAR = factor(CAR1, levels = c("HER2CAR", "tHER2CAR", "NoCAR")))


# --- Step 4: Create the Interaction Plot ---
inter <- ggplot(interaction_long, aes(x = OE1, y = expression, group = CAR1, color = CAR1)) +
  #geom_line(linewidth = 1) +
  #paletteer_d("nationalparkcolors::Acadia") +
  geom_point(size = 4, position = position_jitterdodge(jitter.width = 0.08, jitter.height = 0.1, dodge.width = 0.4), shape =21) +
  #geom_hline(yintercept = avg) +
  scale_color_npg() +
  #geom_errorbar(
    #aes(ymin = mean_expr - se, ymax = mean_expr + se),
    #width = 0.1,
    #linewidth = 0.8
 # ) +
  #scale_y_log10() +
  #scale_color_paletteer_d(`"awtools::mpalette"`) +
  labs(
    title = "Interaction between OE and CAR Type",
    subtitle = "Expression of Gene OR7A10",
    y = "Normalized Expression",
    color = "CAR Type"
  ) +
  # Set the y-axis to start at 0
  #scale_y_continuous(
  #expand = c(0, 0), 
  #limits = c(0, NA)
  #) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )
ggsave("Interaction plot between OE and CAR.pdf", plot = inter, width = 6, height = 5)

#6. Save to files

write.csv(resLFC_stim_df, file="Differential_expression_result/DEG_stim_24hr_vs_0hr.csv")
write.csv(resLFC_stim_df[resLFC_stim_df$significant != "Not Significant", ], file="Differential_expression_result/DEG_stim_24hr_vs_0hr_padj005_LFC05.csv")

write.csv(resLFC_OE1_df, file="Differential_expression_result/DEG_subset_OE_OR7A10_vs_OR7A10Stop.csv")
write.csv(resLFC_OE1_df[resLFC_OE1_df$significant != "Not Significant", ], file="Differential_expression_result/DEG_subset_OE_OR7A10_vs_OR7A10Stop_padj005_LFC05.csv")

write.csv(resLFC_OE1OR7A10.CAR1HER2CAR_df, file="Differential_expression_result/DEG_subset_Interaction_OE1OR7A10.CAR1HER2CAR.csv")
