# Structural Variant Analysis

[![Language](https://img.shields.io/badge/language-Python-blue)]  
[![Format](https://img.shields.io/badge/format-Jupyter_Notebook-orange)]  
[![Analysis](https://img.shields.io/badge/analysis-Structural_Variants-purple)]  
[![Status](https://img.shields.io/badge/status-Manuscript_Version-lightgrey)]  

This repository contains the **structural variant (SV) analysis notebooks** used in the manuscript  
**_OR7A10 GPCR engineering boosts CAR-NK therapy against solid tumors_**.

The code here reflects the **exact analysis state used for the manuscript** and is provided for **reproducibility and transparency**, harmonized with the RNA-seq and scRNA-seq analysis repositories.

---

## Overview

Structural variants were analyzed to quantify genomic alterations associated with OR7A10-engineered CAR-NK cells and matched controls. The workflow focuses on:

- Quantification of structural variant burden  
- Identification of unique and genotype-specific variants  
- Generation of curated SV tables for downstream interpretation  

This repository does **not** perform structural variant calling itself and assumes that SV calls were generated upstream using external tools.

---

## Repository Contents

### `Structural_Variant_quantification.ipynb`

**Purpose:**  
Quantify structural variants across samples and variant classes.

**Key steps:**
- Import structural variant callsets  
- Classify variants by type (e.g., deletions, duplications, inversions)  
- Aggregate variant counts per sample  
- Export summary tables used in downstream analyses  

---

### `Unique_Variant_Filtration.ipynb`

**Purpose:**  
Identify **unique or condition-specific structural variants**.

**Key steps:**
- Remove shared or common variants between sample groups  
- Perform set-based filtering across genotypes  
- Export filtered variant lists for interpretation and visualization  

---

## Dependencies

Analyses were performed using standard scientific Python tooling:

- Python ≥ 3.9  
- pandas  
- numpy  
- matplotlib / seaborn  
- jupyter  

Exact package versions were not pinned; results should be reproducible with standard Python environments.

---

## Usage

Run notebooks in the following order:

1. `Structural_Variant_quantification.ipynb`
2. `Unique_Variant_Filtration.ipynb`

Input structural variant files must be placed in the directories specified within each notebook.

---

## 📌 Citation

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

