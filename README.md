# Integrated-Transcriptomics-and-Survival-Analysis-of-Breast-Cancer-using-TCGA-RNA-seq-data-in-Rstudio
Comprehensive bioinformatics pipeline for TCGA breast cancer RNA-seq data analysis involving data preprocessing, differential expression analysis using DESeq2, pathway enrichment (GO and KEGG), and Kaplan–Meier survival analysis for biomarker discovery.
# TCGA-BRCA RNA-Seq Analysis Pipeline

## Overview

This project analyzes RNA-seq data from TCGA-BRCA (breast cancer) to identify differentially expressed genes, explore biological pathways, and evaluate the survival impact of selected genes. The workflow is implemented in R using standard Bioconductor packages.

---

## Data Acquisition

RNA-seq count data for TCGA-BRCA was obtained using the `TCGAbiolinks` package. The dataset was saved locally as an RDS file to avoid repeated downloads.

---

## Differential Expression Analysis

- Raw count data was extracted from the SummarizedExperiment object  
- Samples were grouped into **Tumor** and **Normal** based on sample type  
- Differential expression analysis was performed using `DESeq2`  
- Genes were ranked based on adjusted p-values  
- Significant genes were defined as those with `padj < 0.05`  

---

## Visualization

A volcano plot was generated using `ggplot2` to visualize:

- Log2 fold changes  
- Statistical significance of genes  
- Highlighted significant vs non-significant genes  

---

## Functional Enrichment Analysis

Significant genes were used for pathway analysis:

- ENSEMBL gene IDs were converted to ENTREZ IDs using `org.Hs.eg.db`  
- Gene Ontology (GO) enrichment analysis was performed (Biological Process)  
- KEGG pathway analysis was conducted to identify affected signaling pathways  
- Results were visualized using dot plots  

---

## Survival Analysis

- Clinical data was extracted from TCGA metadata  
- Expression of a selected gene (**TP53**) was analyzed  
- Patients were divided into **high** and **low expression groups** based on median expression  
- Survival analysis was performed using the Kaplan–Meier method  
- Survival curves were generated using `survminer`  

---

## Tools and Packages

- TCGAbiolinks  
- DESeq2  
- clusterProfiler  
- org.Hs.eg.db  
- enrichplot  
- ggplot2  
- survival  
- survminer  

---

## Outputs

- Differential expression results  
- Volcano plot  
- GO enrichment results  
- KEGG pathway analysis  
- Kaplan–Meier survival plot for TP53  

---

## Notes

This pipeline can be adapted to other TCGA datasets by changing the project ID. It provides a basic framework for RNA-seq analysis including expression profiling, pathway analysis, and survival evaluation.
