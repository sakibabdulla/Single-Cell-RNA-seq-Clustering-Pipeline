"""
Full end-to-end single-cell RNA-seq clustering pipeline (Scanpy + optional Seurat via rpy2)

Filename: scRNAseq_clustering_full_script.py


"""

import os
import argparse
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def parse_args():
    parser = argparse.ArgumentParser(description='Single-cell RNA-seq clustering pipeline (Scanpy)')
    parser.add_argument('--use_example', action='store_true', help='Use built-in PBMC3k example from scanpy.datasets')
    parser.add_argument('--input_10x_dir', type=str, default=None, help='Path to 10x mtx/filtered_feature_bc_matrix directory')
    parser.add_argument('--input_h5', type=str, default=None, help='Path to 10x h5 file (optional)')
    parser.add_argument('--outdir', type=str, default='scRNA_results', help='Output directory')
    parser.add_argument('--min_genes', type=int, default=200, help='Minimum genes per cell')
    parser.add_argument('--min_cells', type=int, default=3, help='Minimum cells per gene')
    parser.add_argument('--n_pcs', type=int, default=30, help='Number of PCs for PCA')
    parser.add_argument('--cluster_resolution', type=float, default=1.0, help='Clustering resolution (Louvain/Leiden)')
    parser.add_argument('--use_leiden', action='store_true', help='Use Leiden instead of Louvain')
    return parser.parse_args()


def save_umap_and_violin(adata, outdir, marker_genes):
    sc.pl.umap(adata, color=['louvain'] if 'louvain' in adata.obs.columns else ["leiden"], save='_clusters.png', show=False)
    plt.close('all')
    # plot a few markers
    for g in marker_genes[:6]:
        sc.pl.umap(adata, color=g, save=f'_{g}.png', show=False)
        plt.close('all')
    # violin for top markers
    sc.pl.violin(adata, keys=marker_genes[:6], groupby='louvain' if 'louvain' in adata.obs.columns else 'leiden', save='_markers_violin.png', show=False)
    plt.close('all')



def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

   
    if args.use_example:
        print('Loading PBMC3k example from scanpy.datasets...')
        adata = sc.datasets.pbmc3k()
    elif args.input_h5 is not None:
        print(f'Loading 10x h5 file: {args.input_h5}')
        adata = sc.read_10x_h5(args.input_h5)
    elif args.input_10x_dir is not None:
        print(f'Loading 10x mtx directory: {args.input_10x_dir}')
        adata = sc.read_10x_mtx(args.input_10x_dir, var_names='gene_symbols', cache=True)
    else:
        raise ValueError('No input provided. Use --use_example or provide --input_10x_dir or --input_h5')

    print('Initial data shape:', adata.shape)

   
    adata.raw = adata

  
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_genes(adata, min_cells=args.min_cells)

  
    adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    adata = adata[adata.obs['pct_counts_mt'] < 20, :]
    adata = adata[adata.obs['n_genes_by_counts'] < 7500, :]

    print('After filtering shape:', adata.shape)

   
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True, flavor='seurat_v3')

  
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack', n_comps=args.n_pcs)
    sc.pl.pca_variance_ratio(adata, log=True, save='_pca_variance.png', show=False)

 
    sc.pp.neighbors(adata, n_pcs=args.n_pcs)
    sc.tl.umap(adata)

   
    if args.use_leiden:
        print('Running Leiden clustering...')
        sc.tl.leiden(adata, resolution=args.cluster_resolution)
        adata.obs['leiden'] = adata.obs['leiden'].astype('category')
    else:
        print('Running Louvain clustering...')
        sc.tl.louvain(adata, resolution=args.cluster_resolution)
        adata.obs['louvain'] = adata.obs['louvain'].astype('category')

  
    groupby_key = 'louvain' if 'louvain' in adata.obs.columns else 'leiden'
    sc.tl.rank_genes_groups(adata, groupby=groupby_key, n_genes=25, method='wilcoxon')
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    markers.to_csv(os.path.join(args.outdir, 'markers_all_clusters.csv'), index=False)

  
    cluster_assign = adata.obs[[groupby_key]]
    cluster_assign.to_csv(os.path.join(args.outdir, 'cell_clusters.csv'))

  
    adata.write_h5ad(os.path.join(args.outdir, 'adata_processed.h5ad'))

   
    sc.pl.umap(adata, color=groupby_key, title=f'UMAP_{groupby_key}', save=f'_{groupby_key}.png', show=False)
    plt.close('all')

   
    top_markers = []
    rg = adata.uns['rank_genes_groups']
    groups = rg['names'].dtype.names
    for g in groups:
        top = [x for x in rg['names'][g][:5]]
        top_markers.extend(top)
    top_markers = list(dict.fromkeys(top_markers))

    for gene in top_markers:
        if gene in adata.var_names:
            sc.pl.umap(adata, color=gene, save=f'_{gene}.png', show=False)
            plt.close('all')

    sc.pl.rank_genes_groups_heatmap(adata, groupby=groupby_key, show=False, save='_ranked_heatmap.png')
    plt.close('all')

    # -------------------------
    # Optional: Run Seurat via rpy2 (example using exported counts CSV)
    # -------------------------
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        pandas2ri.activate()
        has_rpy2 = True
    except Exception as e:
        print('rpy2 not available or failed to import. Skipping Seurat section. Error:', e)
        has_rpy2 = False

    if has_rpy2:
        print('Preparing counts for Seurat...')
        # Save raw counts (genes x cells) to CSV for R/Seurat example
        counts = adata.raw.X if adata.raw is not None else adata.X
        # Convert sparse matrix to array if needed
        try:
            counts_arr = counts.toarray()
        except Exception:
            counts_arr = counts
        counts_df = pd.DataFrame(counts_arr.T, index=adata.obs_names, columns=adata.var_names)
        counts_csv = os.path.join(args.outdir, 'counts_for_seurat.csv')
        counts_df.to_csv(counts_csv)

        # Minimal R code executed via rpy2 to create a Seurat object, normalize, find variable features and cluster
        r_code = f"""
        library(Seurat)
        counts <- read.csv('{counts_csv}', row.names=1)
        # counts currently: cells in rows, genes in columns. Seurat expects genes x cells
        counts_t <- t(as.matrix(counts))
        seu <- CreateSeuratObject(counts = counts_t)
        seu <- NormalizeData(seu)
        seu <- FindVariableFeatures(seu, selection.method = 'vst', nfeatures = 2000)
        seu <- ScaleData(seu)
        seu <- RunPCA(seu, npcs = {args.n_pcs})
        seu <- FindNeighbors(seu, dims = 1:{min(20, args.n_pcs)})
        seu <- FindClusters(seu, resolution = {args.cluster_resolution})
        seu <- RunUMAP(seu, dims = 1:{min(20, args.n_pcs)})
        # save cluster assignments
        clusters <- Idents(seu)
        write.csv(as.data.frame(clusters), file = '{os.path.join(args.outdir, 'seurat_clusters.csv')}')
        # save seurat object
        saveRDS(seu, file = '{os.path.join(args.outdir, 'seurat_obj.rds')}')
        """
        try:
            print('Running Seurat (R) code via rpy2...')
            ro.r(r_code)
            print('Seurat section finished. Output: seurat_clusters.csv and seurat_obj.rds')
        except Exception as e:
            print('Running Seurat via rpy2 failed. Error:', e)

    print('All done. Results saved to', args.outdir)


if __name__ == '__main__':
    main()
