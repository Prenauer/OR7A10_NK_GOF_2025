# NKGOF Bulk RNA-seq Analysis

![R](https://img.shields.io/badge/R-%3E%3D4.1-blue)
![DESeq2](https://img.shields.io/badge/DESeq2-RNA--seq-green)
![GSEA](https://img.shields.io/badge/GSEA-GO%20Biological%20Process-orange)
![Status](https://img.shields.io/badge/Status-Manuscript%20Analysis-success)

This directory contains R scripts used for **bulk RNA-seq differential expression and pathway analysis** in the NKGOF study, supporting figures and results reported in the associated manuscript.

---

## Overview of Analysis

The workflow consists of two major components:

1. **Differential expression analysis** using **DESeq2**
2. **Pathway enrichment analysis** using **GSEA / clusterProfiler**

These scripts analyze transcriptional changes associated with:
- OR7A10 vs OR7A10Stop
- CAR vs no-CAR conditions
- Stimulation (24 hr vs 0 hr)
- Higher-order interactions between CAR, OR7A10 status, and stimulation

---

## Workflow Summary

### 1. Differential expression analysis
[01_RNAseq_analysis.Rmd](https://github.com/Prenauer/OR7A10_NK_GOF_2025/blob/main/4_BulkRNAseq/01_RNAseq_analysis.md)
- Import of raw count matrices
- Filtering of low-expression genes
- Construction of sample metadata
- DESeq2 model fitting with interaction terms  
  (`OE * CAR * stim`)
- Variance-stabilizing transformation (VST)
- Quality control:
  - Sample correlation heatmaps
  - PCA colored by stimulation, CAR, and OR7A10 status
- Differential expression testing with **LFC shrinkage (apeglm)**
- Volcano plots for:
  - Main effects
  - Subset analyses
  - Interaction terms
- Interaction plots for selected genes (e.g. *OR7A10*)
- Export of full and filtered DEG tables

**Key outputs**
- PCA plots (`PCA_*.pdf`)
- Volcano plots (`Volcano_*.pdf`)
- Interaction plots
- DEG tables (`.csv`)

---

### 2. Gene set enrichment and pathway analysis
[02_Pathway_analysis.Rmd](https://github.com/Prenauer/OR7A10_NK_GOF_2025/blob/main/4_BulkRNAseq/02_Pathway_analysis.md)

- Import of curated DEG tables
- Construction of ranked gene lists by log2 fold change
- GSEA using `clusterProfiler::gseGO`
- Visualization with dot plots and enrichment plots

Analyses performed for:
- Stimulation main effect (24 hr vs 0 hr)
- OR7A10 vs OR7A10Stop (stimulated subset)

**Key outputs**
- GO enrichment dot plots
- GSEA enrichment plots
- Publication-ready pathway figures

---

## Software Requirements

- R (≥ 4.1)
- Bioconductor packages:
  - DESeq2
  - apeglm
  - clusterProfiler
  - fgsea
  - GSEABase
  - org.Hs.eg.db
- CRAN packages:
  - tidyverse
  - data.table
  - ggplot2
  - ggrepel
  - pheatmap
  - ggrastr
  - paletteer
  - ggsci

---

## Reproducibility Notes

- Scripts assume gene symbols as rownames
- Factor reference levels are explicitly set before DE testing
- Interaction terms are explicitly modeled
- Output directories (e.g. `Differential_expression_result/`) must exist prior to saving

---

## 📌 Citation

If you use this code, analysis framework, or adapt the Dynamic Signature Relationship (DSR) methodology, please cite:

**OR7A10 GPCR engineering boosts CAR-NK therapy against solid tumors**  
Luojia Yang*, Paul A. Renauer*, Kaiyuan Tang, Josh Saskin, Liqun Zhou, Charles Zou,  
Seok-Hoon Lee, Madison Fox, Samuel Johnson-Noya, Benedict Weiss, Stephanie Deng,  
Paris Fang, Binfan Chen, Giacomo Sferruzza, Saba Fooladi, Kai Zhao, Daniel Park,  
Feifei Zhang, Jiayi Tu, Jing Chen, Jennifer Moliterno, Murat Gunel,  
Lei Peng#, and Sidi Chen#. *Accepted-in-principle at **Nature**, December 12, 2025.*

\* Co–first authors  
\# Corresponding authors

---

## License

This project is released under the **MIT License**. See the `LICENSE` file for details.

---

## ✉️ Contact

For questions or collaboration:
- **Name:** Kaiyuan Tang
- **Email:** kaiyuan.tang@yale.edu

