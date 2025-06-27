"""
This script performs preprocessing of single-cell RNA-seq data stored in an `.h5ad` AnnData file.
It follows the standard pipeline for scDRS-compatible preprocessing:
    - Quality control filtering of cells and genes
    - Normalization and log-transformation
    - Selection of highly variable genes (HVGs)
    - Dimensionality reduction (PCA)
    - Neighborhood graph construction
    - UMAP embedding
It also generates basic UMAP plots for key metadata categories.

Usage:
------
python preprocess_adata.py -i <input.h5ad> -o <output.h5ad> -u <umap_output_dir>

Arguments:
----------
-i / --input       : Path to the input AnnData `.h5ad` file.
-o / --output      : Path to save the processed `.h5ad` file.
-u / --umap_dir    : Directory to save UMAP visualization figures.

Output:
-------
- Processed `.h5ad` file written to the specified output path.
- UMAP PNGs stored under:
    <umap_dir>/<output_stem>/general_umaps/*.png

Notes:
------
- The script intentionally skips gene scaling (no `sc.pp.scale`) to preserve log-normal structure 
  required by scDRS gene score modeling.
- Cell filtering uses min_genes=200 and gene filtering uses min_cells=50.
- Highly variable gene selection is performed with `batch_key='batch'` to account for batch variability.
- UMAPs are colored by several metadata fields including `celltype_imputed_lowerres` and `batch`.

Author: Hessel Dijkstra

"""
import anndata as an
import scanpy as sc
import pandas as pd
from pathlib import Path
import argparse

sc.settings.set_figure_params(dpi=300)

def preprocess_adata(adata):
    print("Processing data...")
    adata = filter_adata(adata)
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # sc.pp.highly_variable_genes(adata, min_mean=0.001, max_mean=3, min_disp=0.5, batch_key="batch")
    # adata = adata[:, adata.var.highly_variable]

    # NOTE: No mean-variance scaling for scDRS cell scoring because the scDRS preprocessing step models the gene-level log mean-variance relationship
    #       and should therefore not contain negative values. (See scdrs.pp.compute_stats documentation for details). Also

    # sc.pp.scale(adata)
    
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    print("Done")
    return adata


def filter_adata(adata):
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=50)
    return adata


def save_adata(adata, output_file):
    adata.write_h5ad(filename=output_file)
    print(f"Data written to {output_file}")


def save_general_umap_plots(
    adata,
    figdir,
    colors=["celltype_imputed_lowerres", "celltype_imputed", "batch", "sample_final"],
    color_map="RdBu_r",
    vmin=-5,
    vmax=5,
    wspace=0.7,
    titles = ["Cell Type", "Cell Type", "Batch", "Samples"]
):
    figdir.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = figdir
    
    for idx, color in enumerate(colors):
        sc.pl.umap(
            adata,
            color=[color],
            ncols=1,
            color_map=color_map,
            vmin=vmin,
            vmax=vmax,
            wspace=wspace,
            show=False,
            save=f"_{color}.png",
            title=titles[idx]  
        )
    print(f"UMAPs saved to {figdir}")



def main():
    parser = argparse.ArgumentParser(description="Preprocess and save an AnnData object with scDRS preprocessing.")

    parser.add_argument("-i", "--input", type=str, required=True, help="Path to input .h5ad file")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output .h5ad file name")
    parser.add_argument("-u", "--umap_dir", type=str, required=True, help="Output .h5ad file name")

    args = parser.parse_args()

    input_file = Path(args.input)
    output_file = Path(args.output)
    gen_umap_dir = Path(args.umap_dir) / output_file.stem / "general_umaps"

    # Load input
    adata = an.read_h5ad(input_file)

    # Process
    adata = preprocess_adata(adata)
    
    # Save data
    save_adata(adata, output_file)
    save_general_umap_plots(adata, gen_umap_dir)
    

if __name__ == "__main__":
    main()