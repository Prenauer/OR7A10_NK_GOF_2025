# Structural Variant Analysis Pipeline

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

## Repository Structure
```text
.
├── WGS_Data
│ └─ Donor0958_OR7A10OE_Unique_Indel_Variant.vcf
│ └─ Donor0958_OR7A10OE_Unique_SV1000.vcf
│ └─ Donor0958_OR7A10Stop_Unique_Indel_Variant.vcf
│ └─ Donor0958_OR7A10Stop_Unique_SV1000.vcf
│
├── 01_Structural_Variant_quantification.ipynb
│ └─── Import, classify, aggregate, export
│
├── 02_Unique_Variant_Filtration.ipynb
│ └─── Format, filter, export
│
├── 03_StructuralVar_Analysis.md
├── 03_StructuralVar_Analysis.Rmd
│ └─── Analyze, annotate, plot
│
├── LICENSE
└── README.md
```

---


### `1. Structural VariantQuantification`
[01_Structural_Variant_quantification.ipynb](https://github.com/Prenauer/OR7A10_NK_GOF_2025/blob/main/5_WGS/01_Structural_Variant_quantification.ipynb)

**Purpose:**  
Quantify structural variants across samples and variant classes.

**Key steps:**
- Import structural variant callsets  
- Classify variants by type (e.g., deletions, duplications, inversions)  
- Aggregate variant counts per sample  
- Export summary tables used in downstream analyses  

### `2. Unique Variant Filtering`
[02_Unique_Variant_Filtration.ipynb](https://github.com/Prenauer/OR7A10_NK_GOF_2025/blob/main/5_WGS/02_Unique_Variant_Filtration.ipynb)

**Purpose:**  
Identify **unique or condition-specific structural variants**.

**Key steps:**
- Remove shared or common variants between sample groups  
- Perform set-based filtering across genotypes  
- Export filtered variant lists for interpretation and visualization  

### `3. Structural Variant Analysis`
[03_StructuralVar_Analysis.Rmd](https://github.com/Prenauer/OR7A10_NK_GOF_2025/blob/main/5_WGS/03_StructuralVar_Analysis.md)

**Purpose:**
- Characterize and visualize condition-specific genomic variation—including indels and structural variants (SVs)—from whole-genome sequencing (WGS) data, and compare their genomic distribution between experimental conditions.

**Key steps:**
- Import VCF files for condition-specific indels and structural variants
- Remove VCF metadata and harmonize chromosome annotations
- Bin the genome into fixed-size windows to compute variant density
- Count indels per genomic bin for each condition
- Extract SV breakpoint pairs and filter out short-distance, likely artifactual events
- Convert breakpoint data into BEDPE format for downstream visualization
- Load GRCh38 cytoband ideograms and configure Circos plotting parameters
- Generate Circos plots showing genome-wide indel density and SV links for each condition
- Export publication-ready Circos PDFs for comparative interpretation

---
## Dependencies

Analyses were performed using standard scientific Python tooling:

- Python ≥ 3.9  
- pandas  
- numpy  
- matplotlib / seaborn  
- jupyter  

Exact package versions were not pinned; results should be reproducible with standard Python environments.

- R (≥ 4.1)
- vcfR
- ape
- dplyr
- tidyverse
- data.table
- CMplot
- VariantAnnotation
- StructuralVariantAnnotation
- RCircos

---

## Usage

Run notebooks in the following order:

1. `01_Structural_Variant_quantification.ipynb`
2. `02_Unique_Variant_Filtration.ipynb`
3. `03_StructuralVar_Analysis.Rmd`

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

