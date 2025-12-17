# SAMBA (Manuscript Version)

[![R](https://img.shields.io/badge/R-%E2%89%A54.2-blue.svg)](https://cran.r-project.org/)
[![CRISPR Screen Analysis](https://img.shields.io/badge/Application-CRISPR%20Screens-purple)](#)
[![Manuscript Code](https://img.shields.io/badge/Code-Manuscript%20Snapshot-orange)](#)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

---

## Overview

This directory contains the **exact version of the SAMBA (Statistical Analysis of Multiplexed Barcode Assays) code used in the manuscript**:

> **OR7A10 GPCR engineering boosts CAR-NK therapy against solid tumors**

This version is provided **for reproducibility of the published analyses only**.  
It reflects the implementation and defaults used at the time of manuscript submission and **is not actively maintained**.

👉 **For the actively developed and maintained version of SAMBA, please visit:**

**https://github.com/Prenauer/SAMBA**

---

## Description

**SAMBA** is an end-to-end statistical framework for analyzing pooled CRISPR screening data, supporting:

- Guide-level enrichment analysis using edgeR (QLF or LRT)
- Gene-level aggregation via weighted gene-score models
- Robust null-distribution construction
- Flexible design matrices and contrasts
- Built-in preprocessing, normalization, and dispersion estimation

This implementation is optimized for **CRISPRa/i and knockout screens** with replicated designs.

---

## Contents

### Core Functions

| Function | Description |
|--------|-------------|
| `SAMBA()` | Full pipeline wrapper (preprocess → guide-level → gene-level) |
| `Preprocess_Samba()` | Filtering, normalization, dispersion estimation |
| `Analyze_Samba_Guides()` | Guide-level differential analysis |
| `Analyze_Samba_Genes()` | Gene-level aggregation |
| `MetaAnalysisSamba()` | Meta-analysis gene scoring |
| `GeneScoreSamba()` | Weighted gene-score aggregation |
| `FilterCountData()` | Guide filtering |
| `WeighGuides()` | Guide weighting by detectability |
| `RandomIndexGenerator()` | Null distribution construction |
| `WtSumScore()` | Weighted gene scoring |

---

## Intended Use

This code is intended to:

- **Reproduce figures and tables in the manuscript**
- Serve as an **archival snapshot** of the analysis pipeline
- Enable reviewers and readers to audit statistical methods

🚫 **Not recommended** for new analyses or pipeline development.

---

## Maintained Version

For bug fixes, new features, documentation, and support, use:

🔗 **https://github.com/Prenauer/SAMBA**

The maintained version may differ in:
- Default parameters
- Function signatures
- Statistical options
- Performance optimizations

---

## Requirements

- R ≥ 4.2
- edgeR
- limma
- plyr
- dplyr
- stringr
- metap

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

For software questions or feature requests, use the maintained repository:
👉 https://github.com/Prenauer/SAMBA
