# ORF Screen Analysis Pipeline

![Language](https://img.shields.io/badge/language-R-blue)
![Workflow](https://img.shields.io/badge/workflow-R%20Markdown-informational)
![Status](https://img.shields.io/badge/status-accepted--in--principle-green)
![Journal](https://img.shields.io/badge/journal-Nature-lightgrey)
![License](https://img.shields.io/badge/license-MIT-brightgreen)

---

## Overview

This repository contains the complete analysis pipeline for **ORF (Open Reading Frame) screening experiments**, including preprocessing, quality control, and downstream statistical analysis. The workflow is implemented in **R Markdown** to ensure transparency, reproducibility, and ease of extension.

The pipeline is designed for pooled ORF perturbation screens and supports:

- Barcode and sample ID harmonization  
- Read count aggregation and normalization  
- Quality control diagnostics  
- Hit identification and statistical testing  
- Visualization of ORF-level effects  

---

## Repository Structure

.
├── 01_Preprocess_OrfScreen.Rmd
│ └── Sample ID mapping, count aggregation, and preprocessing
│
├── 02_OrfScreen_QC.Rmd
│ └── Quality control, normalization checks, and replicate assessment
│
├── 03_OrfScreen_analysis.Rmd
│ └── Differential analysis, hit calling, and visualization
│
└── README.md


---

## Workflow Summary

### 1. ORF Screen Preprocessing  
**`01_Preprocess_OrfScreen.Rmd`**

- Converts raw ORF count files into a unified matrix  
- Harmonizes sample identifiers across batches  
- Aggregates ORF-level counts  
- Outputs analysis-ready count tables  

---

### 2. Quality Control  
**`02_OrfScreen_QC.Rmd`**

- Assesses library complexity and sequencing depth  
- Identifies low-quality samples and ORFs  
- Evaluates replicate concordance  
- Generates diagnostic and QC plots  

---

### 3. ORF Screen Analysis  
**`03_OrfScreen_analysis.Rmd`**

- Performs differential ORF enrichment/depletion analysis  
- Identifies statistically significant ORF hits  
- Generates volcano plots and summary visualizations  
- Exports ranked ORF-level results  

---

## Requirements

- **R ≥ 4.2**
- Key R packages (non-exhaustive):
  - `tidyverse`
  - `data.table`
  - `BiocParallel`
  - `ggplot2`
  - `cowplot`
  - `rstatix`
  - `reshape2`

Exact package versions and parameters are documented within each R Markdown file.

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

For questions or collaboration:
- **Name:** Paul Renauer 
- **Email:** paul.renauer@yale.edu

