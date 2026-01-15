# Single-Cell RNA-seq (scRNA-seq) Analysis Pipeline

**Author:** Divya Mishra, Ph.D.  
**Expertise:** Molecular Genetics, NGS Data Analysis, Clinical Genomics

---

## Overview

This repository provides a **full end-to-end scRNA-seq pipeline**, from raw FASTQ files to clustering and cell type annotation:

- Quality control (FastQC / MultiQC)
- Adapter trimming (Trim Galore)
- Alignment and UMI counting (STARsolo)
- Normalization, scaling, dimensionality reduction
- Clustering and UMAP visualization
- Marker gene detection and cell type annotation
- R (Seurat) and Python (Scanpy) workflows included

---

##  Repository Structure
```
scRNAseq_pipeline/
├── README.md
├── data/raw/ # Raw FASTQ files
├── results/qc/
├── results/trimmed/
├── results/alignments/
├── results/counts/
├── results/clusters/
├── results/figures/
├── scripts/scRNAseq_commands.sh
├── R/scrnaseq_analysis.R
└── Python/scrnaseq_analysis.py

```
---
##  Running the Pipeline

#### Bash QC, trimming, alignment
```bash
bash scripts/scRNAseq_commands.sh

R Seurat analysis
Rscript R/scrnaseq_analysis.R

Python Scanpy analysis
python Python/scrnaseq_analysis.py
```
### Requirements
```
Command-line: FastQC, MultiQC, Trim Galore, STAR

R: Seurat, dplyr, ggplot2, patchwork

Python: scanpy, matplotlib
```
### Outputs
```
QC reports
Trimmed FASTQ
BAM alignment files
Gene-cell count matrix
Clustering results (UMAP / marker genes)
Publication-ready figures
