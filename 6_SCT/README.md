# 🧬 Dynamic Gene and Pathway Analysis of NK Cell States

[![R](https://img.shields.io/badge/R-%E2%89%A54.2-blue.svg)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-%E2%89%A53.9-yellow.svg)](https://www.python.org/)
[![Seurat](https://img.shields.io/badge/Seurat-single--cell-red.svg)](https://satijalab.org/seurat/)
[![scVelo](https://img.shields.io/badge/scVelo-RNA%20velocity-green.svg)](https://scvelo.readthedocs.io/)
[![License](https://img.shields.io/badge/License-MIT-lightgrey.svg)](LICENSE)

---

## Overview

This repository contains the full computational pipeline for analyzing **single-cell RNA-seq data from NK cell populations**, integrating:

- Seurat-based preprocessing and clustering  
- RNA velocity and trajectory inference (velocyto + scVelo)  
- Differential gene and pathway analysis  
- **Dynamic Signature Relationship (DSR)** modeling using generalized additive models  

The workflow enables **gene-level, pathway-level, and dynamic modeling** of NK cell state transitions.

---

## 📁 Repository Structure

```text
.
├── Data
│ ├─ 01_cellmarker_list_Tang_Cell2023.txt
│ ├─ de_genes_v0.1.txt
│ ├─ de_pw_pid_v0.1.txt
│ ├─ dsr_sig-gene_v0.7.txt
│ └─ subset_pct_comp.txt
│
├── ref
│ ├─ c2.canonical_pathways.v7.4.symbols.gmt.txt
│ ├─ c2.reactome_pathways.v7.4.symbols.gmt.txt
│ ├─ c3.TF_targets.gtrd.v7.4.symbols.gmt.txt
│ ├─ c3.TF_targets.v7.4.symbols.gmt.txt
│ ├─ c4.3ca.v2024.1.Hs.symbols.gmt.txt
│ ├─ c5.go.bp.v7.4.symbols.gmt.txt
│ ├─ c7.immunesigdb.v2024.1.Hs.symbols.gmt.txt
│ ├─ c8.all.v2024.1.Hs.symbols.gmt.txt
│ ├─ MSigDB_Hallmark_2020.txt
│ └─ scsig.all.v1.0.1.symbols.gmt.txt
│
├── 00_Additional_Functions.md
├── 00_Additional_Functions.R
│ └── Utility Functions
│
├── 01_data_processing.md
├── 01_data_processing.Rmd
│ └── Data processing and clustering
│
├── 02a_scVelo_part1.md
├── 02a_scVelo_part1.Rmd
│ └── Export processed data 
│
├── 02b_velocyto_script.txt
│ └── Generate spliced/unspliced counts
│
├── 02c_Trajectory_scVelo.txt
│ └── RNA velocity and trajectory analysis
│
├── 03_DE-Genes.md
├── 03_DE-Genes.Rmd
│ └── DE analyses and plots
│
├── 04_GeneSetSignatures.md
├── 04_GeneSetSignatures.Rmd
│ └── Pathway Signature Analysis and plots
│
├── 05_DynamicSignatureRelationships.md
├── 05_DynamicSignatureRelationships.Rmd
│ └── Modeling pathway-gene relationships
│
├── LICENSE
│
└── README.md
```

---

## Workflow Summary

### 1. Data Processing and Clustering
[01_data_processing.Rmd](https://github.com/Prenauer/OR7A10_NK_GOF_2025/blob/main/6_SCT/01_data_processing.md)

- Import CellRanger outputs
- Quality control and filtering
- Sample integration (iNMF)
- Dimensionality reduction (UMAP)
- Graph-based clustering
- NK cell subset annotation (iNK / mNK)

**Key outputs**
- Processed Seurat objects (`so_proc.rds`)
- UMAP embeddings
- Cluster and cell-type labels

### 2. RNA Velocity and Trajectory Inference
[02a_scVelo_part1.Rmd](https://github.com/Prenauer/OR7A10_NK_GOF_2025/blob/main/6_SCT/02a_scVelo_part1.md)

[02b_velocyto_script.sh](https://github.com/Prenauer/OR7A10_NK_GOF_2025/blob/main/6_SCT/02b_velocyto_script.sh)

[02c_Trajectory_scVelo.py](https://github.com/Prenauer/OR7A10_NK_GOF_2025/blob/main/6_SCT/02c_Trajectory_scVelo.py)

**Steps**
- Generate spliced/unspliced counts with `velocyto`
- Merge velocity data with Seurat embeddings
- Recover transcriptional dynamics
- Compute RNA velocity, pseudotime, and latent time
- Visualize lineage trajectories

**Key outputs**
- `scvelo_analysis.h5ad`
- Velocity stream plots
- Latent time and pseudotime maps

### 3. Differential Gene Expression
[03_DE-Genes.Rmd](https://github.com/Prenauer/OR7A10_NK_GOF_2025/blob/main/6_SCT/03_DE-Genes.md)

- Subset-specific differential expression
- Volcano and upset plots
- Identification of genotype-dependent DE genes

### 4. Pathway Signature Analysis
[04_GeneSetSignatures.Rmd](https://github.com/Prenauer/OR7A10_NK_GOF_2025/blob/main/6_SCT/04_GeneSetSignatures.md)

- AUCell-based pathway scoring
- Differential pathway activity
- Pathway volcano, upset, and network plots

### 5. Dynamic Signature Relationships (DSR)
[05_DynamicSignatureRelationships.Rmd](https://github.com/Prenauer/OR7A10_NK_GOF_2025/blob/main/6_SCT/05_DynamicSignatureRelationships.md)

- GAM-based modeling of pathway–gene dynamics
- Identification of predictive regulators
- Heatmaps and scatter plots of dynamic effects

---

## 🔧 Utility Functions
[00_Additional_Functions.R](https://github.com/Prenauer/OR7A10_NK_GOF_2025/blob/main/6_SCT/00_Additional_Functions.md)

Reusable functions for:
- Matrix integration
- Clustering optimization
- Volcano and network plots
- DSR model fitting and diagnostics

These functions are sourced across multiple analysis notebooks.

---

## Requirements

**R (≥ 4.2)**
- Seurat
- AUCell
- mgcv
- cowplot
- ggrastr

**Python (≥ 3.9)**
- scvelo
- scanpy
- loompy
- numpy
- pandas

**System tools**
- velocyto
- samtools

---

## Workflow Order

1. `01_data_processing.Rmd`
2. `02a_scVelo_part1.Rmd`
3. `02b_velocyto_script.sh`
4. `02c_Trajectory_scVelo.py`
5. `03_DE-Genes.Rmd`
6. `04_GeneSetSignatures.Rmd`
7. `05_DynamicSignatureRelationships.Rmd`

---

## Reproducibility

All analyses are fully scripted and can be rerun end-to-end by executing the R Markdown files in numerical order.  
Paths, thresholds, and analysis parameters are explicitly defined within each script to ensure reproducibility.

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

For questions:
- **Name:** Paul Renauer 
- **Email:** paul.renauer@yale.edu
