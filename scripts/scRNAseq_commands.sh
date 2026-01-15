#!/bin/bash
# ================================================================
# scRNA-seq Pipeline: QC, trimming, alignment
# Author: Divya Mishra, Ph.D.
# ================================================================

# 1. Quality check
fastqc data/raw/sample1_R1.fastq.gz data/raw/sample1_R2.fastq.gz -o results/qc/
multiqc results/qc/ -o results/qc/

# 2. Adapter trimming (if needed)
trim_galore --paired data/raw/sample1_R1.fastq.gz data/raw/sample1_R2.fastq.gz -o results/trimmed/

# 3. Alignment & gene quantification using STARsolo
STAR --runThreadN 8 \
     --genomeDir /path/to/star_genome_index \
     --readFilesIn results/trimmed/sample1_R1_val_1.fq.gz results/trimmed/sample1_R2_val_2.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix results/alignments/sample1_ \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 --soloCBlen 16 \
     --soloUMIstart 17 --soloUMIlen 10 \
     --soloFeatures Gene
### Note: STARsolo is used for droplet-based 10X-like scRNA-seq for UMI counting and gene-cell matrix.
