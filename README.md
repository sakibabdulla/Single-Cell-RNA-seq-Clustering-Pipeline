**Single-Cell RNA-seq Clustering Pipeline**

Welcome to the official  repository that contains a Python script for performing end-to-end single-cell RNA-seq (scRNA-seq) 
clustering analysis using the Scanpy library, with an optional Seurat workflow via rpy2.

----------------------------------------------------------------------------------------------------------------------------------------

![image alt](https://github.com/sakibabdulla/Single-Cell-RNA-seq-Clustering-Pipeline/blob/413fb93a7b4d4c5e2884869616d2332d1e85d018/workflow.png)

------------------------------------------------------------------------------------------------------------------------------------------
Overview

The script processes scRNA-seq data through a standard pipeline:
- Loads data (10x format or built-in PBMC3k example)
- Filters and performs quality control (QC)
- Normalizes data and identifies highly variable genes (HVGs)
- Reduces dimensions using PCA and UMAP
- Cluster cells using Louvain or Leiden algorithms
- Identifies marker genes per cluster
- Generates visualizations (UMAP plots, violin plots) and saves results

Optionally, it integrates with Seurat (via rpy2) for additional analysis. 

----------------------------------------------------------

Features
- Supports 10x mtx directories or h5 files
- Customizable filtering, normalization, and clustering parameters
- Outputs include cluster assignments, marker genes, and plots
- Optional Seurat integration for cross-validation

  ------------------------------------------------------------------------------------

Requirements

Python Dependencies
scanpy, anndata , pandas , numpy ,matplotlib , seaborn , scikit-learn , louvain   ,leidenalg 
, rpy2 (optional for Seurat)

---------------------------------------------------------------------------------------------------------
Install via pip:
pip install scanpy anndata pandas numpy matplotlib seaborn scikit-learn louvain leidenalg rpy2

-----------------------------------------------------------------------------------------------------------


1 Run the script with default settings (using PBMC3k example):

python sc_rnaseq_clustering_full_script.py --use_example --outdir results_pbmc --use_leiden

--------------------------------------------------------------------------------------------------

2 Or, provide your own 10x data:

python scRNAseq_clustering_full_script.py --input_10x_dir /path/to/10x/data --outdir /path/to/output

-----------------------------------------------------------------------------------------------------------

Arguments

--use_example: Use the built-in PBMC3k dataset (default: False)
--input_10x_dir: Path to 10x mtx directory
--input_h5: Path to 10x h5 file
--outdir: Output directory (default: scRNA_results)
--min_genes: Minimum genes per cell (default: 200)
--min_cells: Minimum cells per gene (default: 3)
--n_pcs: Number of principal components (default: 30)
--cluster_resolution: Clustering resolution (default: 1.0)
--use_leiden: Use Leiden instead of Louvain clustering (default: False)

-------------------------------------------------------------------------------------------------------------


Output
- scRNA_results/ directory containing:
- cell_clusters.csv: Cluster assignments for each cell
- markers_all_clusters.csv: Marker genes for each cluster
- adata_processed.h5ad: Processed AnnData object
- UMAP plots (e.g., _clusters.png, _gene_name.png)
- Violin plots and heatmaps
- Optional Seurat outputs (if rpy2 is available):
- seurat_clusters.csv
- seurat_obj.rds

---------------------------------------------------------------------------------------------------


Notes
Ensure input data is in 10x format. Other formats may require preprocessing.
Seurat integration requires R and rpy2 setup; skip if not configured.
Check the script for detailed comments on each step.

------------------------------------------------------------------------------

Contributing

Feel free to fork this repository, submit issues, or send pull requests for improvements!

Author: M Sakib Abdullah 
mohdsakib219@gmail.com


