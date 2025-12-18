# CRISPRa Screen Analysis Pipeline

![Language](https://img.shields.io/badge/language-R-blue)
![Workflow](https://img.shields.io/badge/workflow-R%20Markdown-informational)
![Status](https://img.shields.io/badge/status-accepted--in--principle-green)
![Journal](https://img.shields.io/badge/journal-Nature-lightgrey)
![License](https://img.shields.io/badge/license-MIT-brightgreen)

---

## Overview

This repository contains the analysis workflow for a **CRISPRa (CRISPR activation) pooled screening experiment**, implemented in **R Markdown** for transparency and reproducibility.

The attached script (**`01_CRISPRa_ScreenAnalysis.Rmd`**) performs end-to-end analysis using a **SAMBA-based pipeline**, including:

- Input ingestion and preprocessing into a screen-aware object  
- Guide-level differential analysis (sgRNA statistics)  
- Gene-level aggregation and hit calling  
- Generation of ranked hit plots (gene score vs adjusted significance)  
- Export of gene-level results and publication-ready figures  

---

## Repository Structure

```text
.
├── Data
│   └─ NK_Screen_LY3_CountTable.txt
│   └─ NK_Screen_LY3_GeneLevelResults.txt
│   └─ NK_Screen_LY3_GuideLevelResults.txt
│   └─ NK_Screen_LY3_lcpm.txt
│   └─ NK_Screen_LY3_mds_plotdata.txt
│
├── 01_CRISPRa_ScreenAnalysis.md
├── 01_CRISPRa_ScreenAnalysis.Rmd
│   └─── CRISPRa screen preprocessing, guide/gene analysis, and rank-plot visualization
│
├── LICENSE
│
├── README.md
│
└── Samba_official_V1.1.R
```

---

## Workflow Summary

### 1. Environment Setup and Inputs
[01_Preprocess_OrfScreen.Rmd](https://github.com/Prenauer/OR7A10_NK_GOF_2025/blob/main/3_ORF_Miniscreen/01_Preprocess_OrfScreen.md)
- Loads required R packages and analysis utilities
- Reads in CRISPRa screen count/design inputs
- Defines the experimental design used for differential testing

### 2. Screen Preprocessing (SAMBA)
- Converts raw sgRNA counts into an analysis-ready object via Preprocess_Samba()
- Encodes design information needed for downstream modeling

### 3. Differential Analysis and Hit Calling
- Performs guide-level testing using Analyze_Samba_Guides()
- Aggregates sgRNA evidence into gene-level statistics using Analyze_Samba_Genes()
- Writes gene-level outputs to disk (e.g., *_GeneLevelResults.txt)

### 4. Visualization and Output
- Constructs a ranked plot (gene score / z-score vs -log10(FDR) or analogous adjusted metric)
- Overlays null distributions and thresholds where applicable
- Saves publication-quality plots (e.g., NK_Screen_LY_volcano.pdf)

---

## Requirements

- R ≥ 4.2
- Key R packages (non-exhaustive):
   - tidyverse, dplyr, data.table
   - ggplot2
   - ggrastr
   - ggpointdensity
   - viridis
   - cowplot
- SAMBA pipeline functions referenced in the script:
   - Preprocess_Samba()
   - Analyze_Samba_Guides()
   - Analyze_Samba_Genes()

*Exact parameters, thresholds, and output filenames are defined inside the R Markdown file.

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
