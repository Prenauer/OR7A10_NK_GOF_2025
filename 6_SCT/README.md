# 🧬 Dynamic Gene and Pathway Analysis of NK Cell States

[![R](https://img.shields.io/badge/R-%E2%89%A54.2-blue.svg)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-%E2%89%A53.9-yellow.svg)](https://www.python.org/)
[![Seurat](https://img.shields.io/badge/Seurat-single--cell-red.svg)](https://satijalab.org/seurat/)
[![scVelo](https://img.shields.io/badge/scVelo-RNA%20velocity-green.svg)](https://scvelo.readthedocs.io/)
[![License](https://img.shields.io/badge/License-MIT-lightgrey.svg)](LICENSE)

This repository contains the full computational pipeline for analyzing **single-cell RNA-seq data from NK cell populations**, integrating:

- Seurat-based preprocessing and clustering  
- RNA velocity and trajectory inference (velocyto + scVelo)  
- Differential gene and pathway analysis  
- **Dynamic Signature Relationship (DSR)** modeling using generalized additive models  

The workflow enables **gene-level, pathway-level, and dynamic modeling** of NK cell state transitions.

---

## 📁 Repository Structure


---

## 🧪 Analysis Overview

### 1. Data Processing and Clustering
**Script:** `01_data_processing.Rmd`

- Import CellRanger outputs
- Quality control and filtering
- Sample integration (iNMF)
- Dimensionality reduction (UMAP)
- Graph-based clustering
- NK cell subset annotation (iNK / mNK)

**Key outputs**
- Processed Seurat objects (`so_proc*.rds`)
- UMAP embeddings
- Cluster and cell-type labels

---

### 2. RNA Velocity and Trajectory Inference
**Scripts**
- `02b_velocyto_script.sh`
- `02a_scVelo_part1.Rmd`
- `02c_Trajectory_scVelo.py`

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

---

### 3. Differential Gene Expression
**Script:** `03_DE-Genes.Rmd`

- Subset-specific differential expression
- Volcano and upset plots
- Identification of genotype-dependent DE genes

---

### 4. Pathway Signature Analysis
**Script:** `04_GeneSetSignatures.Rmd`

- AUCell-based pathway scoring
- Differential pathway activity
- Pathway volcano, upset, and network plots

---

### 5. Dynamic Signature Relationships (DSR)
**Script:** `05_DynamicSignatureRelationships.Rmd`

- GAM-based modeling of pathway–gene dynamics
- Identification of predictive regulators
- Heatmaps and scatter plots of dynamic effects

**Core idea**
> Model continuous transcriptional relationships while accounting for genotype, cell state, and nonlinearity.

---

## 🔧 Core Utility Functions

**Script:** `00_Additional_Functions.R`

Reusable functions for:
- Matrix integration
- Clustering optimization
- Volcano and network plots
- DSR model fitting and diagnostics

These functions are sourced across multiple analysis notebooks.

---

## ▶️ Running the Pipeline

### Requirements

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

### Recommended Execution Order

1. `01_data_processing.Rmd`
2. `02b_velocyto_script.sh`
3. `02a_scVelo_part1.Rmd`
4. `02c_Trajectory_scVelo.py`
5. `03_DE-Genes.Rmd`
6. `04_GeneSetSignatures.Rmd`
7. `05_DynamicSignatureRelationships.Rmd`

---

## 📊 Outputs

- **Figures/**  
  Publication-quality figures (UMAPs, volcano plots, velocity streams, networks)

- **Data/**  
  Intermediate and final result tables

---

## 🧠 Conceptual Notes

- Pathway-level modeling improves robustness over gene-only analyses
- RNA velocity provides a continuous temporal axis for dynamic inference
- DSR explicitly separates predictor effects, genotype effects, and residual structure

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
- **Name:** Paul Renauer 
- **Email:** paul.renauer@yale.edu
