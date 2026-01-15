"""
scRNA-seq analysis using Scanpy
Author: Divya Mishra, Ph.D.
"""

import scanpy as sc
import matplotlib.pyplot as plt

# 1. Load gene-cell matrix
adata = sc.read_10x_mtx("results/counts/", var_names='gene_symbols', cache=True)

# 2. QC filtering
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.n_genes_by_counts > 200, :]
adata = adata[adata.obs.pct_counts_mt < 10, :]

# 3. Normalization and log transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable]

# 4. Scaling and PCA
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

# 5. Neighbors and clustering
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
sc.tl.leiden(adata, resolution=0.5)
sc.tl.umap(adata)

# 6. Marker genes and visualization
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=10, save='_markers.png')
sc.pl.umap(adata, color=['leiden'], save='_clusters.png')

# Save AnnData object
adata.write("results/clusters/scRNAseq_scanpy.h5ad")
