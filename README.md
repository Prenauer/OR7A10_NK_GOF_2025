# 🧬 OR7A10 CAR-NK Computational Analyses

[![Language](https://img.shields.io/badge/language-R%20%7C%20Python-blue)]  
[![Analysis](https://img.shields.io/badge/analysis-single--cell%20%7C%20RNA--seq%20%7C%20genomics-purple)]  
[![Status](https://img.shields.io/badge/status-Manuscript_Version-lightgrey)]  
[![License](https://img.shields.io/badge/license-MIT-lightgrey)]  

This umbrella repository documents and links **all computational analysis code** used in the manuscript:

> **_OR7A10 GPCR engineering boosts CAR-NK therapy against solid tumors_**  
> *Accepted-in-principle at **Nature**, December 2025*

Each linked repository reflects the **exact analysis state used for the manuscript** and 
is provided for **reproducibility, transparency, and review**.  
Code is **not maintained** beyond the manuscript version.

---

## 📁 Analysis Repositories

### 1️⃣ SAMBA (Manuscript Version)

**Scope**
- Screen Analysis Method with emperical Bayes estimation for Aggregated gene scoring (SAMBA)
- Manuscript implementation only

⚠️ **Important**
This repository contains **only the version of SAMBA used in the manuscript**.  
The **actively maintained version** is available at:

➡️ https://github.com/Prenauer/SAMBA

➡️ **Repository:** [Samba_screen_analysis](https://github.com/Prenauer/OR7A10_NK_GOF_2025/tree/main/1_SAMBA)

---

### 2️⃣ CRISPRa Screen Analysis

**Scope**
- CRISPRa screen data processing  
- Gene-level and pathway-level enrichment  

➡️ **Repository:** [CRISPRa_screen](https://github.com/Prenauer/OR7A10_NK_GOF_2025/tree/main/2_CRISPRa_screen)

---

### 3️⃣ ORF Screen Analysis

**Scope**
- ORF screen preprocessing  
- Quality control and normalization  
- Hit identification and downstream analysis  

**Execution order**
1. `01_Preprocess_OrfScreen.Rmd`
2. `02_OrfScreen_QC.Rmd`
3. `03_OrfScreen_analysis.Rmd`

➡️ **Repository:** [ORF_screen](https://github.com/Prenauer/OR7A10_NK_GOF_2025/tree/main/3_ORF_Miniscreen)

---

### 4️⃣ Bulk RNA-seq (NK GOF / LOF)

**Scope**
- Differential expression with DESeq2  
- Pathway enrichment and visualization  

**Key scripts**
- `NKGOF_RNAseq_analysis_DESeq2_Git.R`
- `Pathway_analysis_Git.R`

➡️ **Repository:** [RNA-seq_analysis](https://github.com/Prenauer/OR7A10_NK_GOF_2025/tree/main/4_BulkRNAseq)

---

### 5️⃣ Structural Variant Analysis

**Scope**
- Quantification of structural variant burden  
- Identification of unique and genotype-specific SVs  

**Notebooks**
- `Structural_Variant_quantification.ipynb`
- `Unique_Variant_Filtration.ipynb`

➡️ **Repository:** [Genetic_Variant_Analysis](https://github.com/Prenauer/OR7A10_NK_GOF_2025/tree/main/5_WGS)  

---

### 6️⃣ Single-cell RNA-seq & Dynamic Modeling (Core Analysis)

**Scope**
- scRNA-seq preprocessing and clustering  
- RNA velocity and trajectory inference  
- Differential gene and pathway analysis  
- Dynamic Signature Relationship (DSR) modeling  

➡️ **Repository:** [Single-cell RNA-seq & Dynamic Modeling](https://github.com/Prenauer/OR7A10_NK_GOF_2025/tree/main/6_SCT)

---

## ▶️ General Notes

- Scripts and notebooks must be run **in the order specified** within each repository README.
- Input data paths are defined **inside each script/notebook**.
- External preprocessing (e.g., CellRanger, variant calling) is assumed complete.

---

## 📌 Citation

If you use or reference any code from these repositories, please cite:

**OR7A10 GPCR engineering boosts CAR-NK therapy against solid tumors**  
Luojia Yang*, Paul A. Renauer*, Kaiyuan Tang, Josh Saskin, Liqun Zhou, Charles Zou,  
Seok-Hoon Lee, Madison Fox, Samuel Johnson-Noya, Benedict Weiss, Stephanie Deng,  
Paris Fang, Binfan Chen, Giacomo Sferruzza, Saba Fooladi, Kai Zhao, Daniel Park,  
Feifei Zhang, Jiayi Tu, Jing Chen, Jennifer Moliterno, Murat Gunel,  
Lei Peng#, and Sidi Chen#.  
*Accepted-in-principle at **Nature**, December 2025.*

\* Co–first authors  
\# Corresponding authors  

---

## 📄 License

All repositories are released under the **MIT License**.

---

## ✉️ Contact

For scientific or technical questions:

- **Paul A. Renauer** — paul.renauer@yale.edu  
- **Kaiyuan Tang** — kaiyuan.tang@yale.edu  

