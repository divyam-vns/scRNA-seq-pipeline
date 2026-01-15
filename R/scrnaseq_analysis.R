# ================================================================
# scRNA-seq analysis using Seurat
# Author: Divya Mishra, Ph.D.
# ================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# 1. Load gene-cell count matrix (from STARsolo)
counts <- Read10X(data.dir = "results/counts/")
sc <- CreateSeuratObject(counts = counts, project = "scRNAseq", min.cells = 3, min.features = 200)

# 2. QC and filtering
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
VlnPlot(sc, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol=3)

sc <- subset(sc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# 3. Normalization and variable features
sc <- NormalizeData(sc)
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)

# 4. Scaling and PCA
all.genes <- rownames(sc)
sc <- ScaleData(sc, features = all.genes)
sc <- RunPCA(sc, features = VariableFeatures(sc))

# 5. Clustering and UMAP
sc <- FindNeighbors(sc, dims = 1:20)
sc <- FindClusters(sc, resolution = 0.5)
sc <- RunUMAP(sc, dims = 1:20)

# 6. Identify markers
markers <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "results/clusters/marker_genes.csv")

# 7. Cell type annotation (manual or automated)
# Example: based on known markers
FeaturePlot(sc, features = c("CD3D","MS4A1","LYZ"))

# 8. Save Seurat object
saveRDS(sc, file="results/clusters/scRNAseq_seurat.rds")

# 9. Visualization
pdf("results/figures/UMAP_clusters.pdf")
DimPlot(sc, reduction = "umap", label = TRUE)
dev.off()
