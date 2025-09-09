library("dplyr")
library("fgsea")
library("GSEABase")
library("tibble")
library("clusterProfiler")
library("enrichplot")
library("GO.db")
library("ggplot2")
library(data.table)

###stim main effect
#read in file
genes <- read.csv("Differential_expression_result/DEG_stim_24hr_vs_0hr_padj005_LFC05.csv")

gene_list <- genes$log2FoldChange #the input of sort has to be a vector 
names(gene_list) <- genes$X
gene_list <- sort(gene_list, decreasing = TRUE)
library("org.Hs.eg.db", character.only = TRUE)
#BiocManager::install("org.Hs.eg.db", character.only = TRUE)
library("org.Hs.eg.db", character.only = TRUE)
keytypes(org.Hs.eg.db)
library(GSEABase)
library(clusterProfiler)
gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db)
require(DOSE)
stim <- dotplot(gse, showCategory=6, split=".sign") + facet_grid(.~.sign) + 
  ggtitle("Stim 24hr vs 0hr Main Effect Pathway Analysis") +
  scale_color_npg()+
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave("Pathway_Stim_Main_Effect.pdf", plot = stim, width = 9, height = 8)

#OE main effect for stimulation subset
genes <- read.csv("Differential_expression_result/DEG_subset_OE_OR7A10_vs_OR7A10Stop_padj005_LFC05.csv")

gene_list <- genes$log2FoldChange #the input of sort has to be a vector 
names(gene_list) <- genes$X
gene_list <- sort(gene_list, decreasing = TRUE)
library("org.Hs.eg.db", character.only = TRUE)
#BiocManager::install("org.Hs.eg.db", character.only = TRUE)
library("org.Hs.eg.db", character.only = TRUE)
keytypes(org.Hs.eg.db)
library(GSEABase)
library(clusterProfiler)
gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db)
require(DOSE)
OE <- dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign) + 
  ggtitle("OE OR7A10OE_vs_OR7A10Stop Main Effect Pathway Analysis") +
  #theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
ggsave("Pathway_OE_Main_Effect.pdf", plot = OE, width = 9, height = 8)

#GSEA Plot for the enriched pathways
gseaplot(gse, by = "all",title = gse$Description[15], geneSetID = 15)