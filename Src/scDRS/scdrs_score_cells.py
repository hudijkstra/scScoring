"""
This script performs single-modality scDRS(single-cell Disease Relevance Score) scoring using an input `.h5ad` file and a GWAS trait-associated gene set. This script computes
individual-cell level scDRS scoress, performs downstream group analysis, and generates UMAP plots of the normalized scDRS scores.

For more information, please visit the official scDRS documentation: https://martinjzhang.github.io/scDRS/

Usage:
    python scdrs_score_single.py \
        -d /path/to/data.h5ad \
        -gs /path/to/geneset.gs \
        -o /path/to/output_dir \
        -t "Trait Name"


Author: Hessel Dijkstra
"""

import scdrs
import anndata as an
import scanpy as sc
import pandas as pd
from pathlib import Path
import argparse
import warnings

warnings.filterwarnings("ignore")
sc.settings.set_figure_params(dpi=300)


def preprocess(adata, covariates = ["sex_code", "age", "n_genes", "batch_code"], adj_prop = "celltype_imputed_lowerres"):
    """
    Preprocesses an AnnData object using scDRS preprocessing.

    1. Correct covariates by regressing out the covariates (including a constant term) and adding back the original mean for each gene.

    2. Compute gene-level and cell-level statistics for the covariate-corrected data.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object with `.obs` containing metadata columns.
    covariates : list of str, optional
        List of covariate column names in `adata.obs`
    adj_prop : str, optional
        The adj_prop option is used for adjusting for cell group proportions, 
        where each cell is inversely weighted proportional to its corresponding cell group size for 
        computing expression mean and variance for genes.

    Returns
    -------
    AnnData
        Modified AnnData object with preprocessing information attached.
    """
    try:
        df_cov = adata.obs[covariates].copy()
    except KeyError as e:
        missing = [key for key in covariates if key not in adata.obs.columns]
        raise KeyError(f"Missing covariates in adata.obs: {missing}") from e

    scdrs.preprocess(adata, cov=df_cov, adj_prop=adj_prop, n_mean_bin=20, n_var_bin=20, copy=False)
    return adata


def load_geneset(geneset_path, adata, src_species="human", dst_species="human"):
    """
    Loads a gene set from file and intersects it with the genes in the AnnData object.

    Parameters
    ----------
    geneset_path : str or Path
        Path to a gene set file in `.gmt` format.
    adata : AnnData
        AnnData object whose gene set should be intersected.
    src_species : str, optional
        Source species of the gene set (default: "human").
    dst_species : str, optional
        Target species for mapping genes (default: "human").

    Returns
    -------
    dict
        Dictionary of intersected gene sets.
    """
    dict_gs = scdrs.util.load_gs(
        geneset_path,
        src_species=src_species,
        dst_species=dst_species,
        to_intersect=adata.var_names,
    )
    return dict_gs


def score_cells(adata, gene_list, gene_weights, trait):
    """
    Run scdrs score_cell function

    Parameters
    ----------
    adata : AnnData
        Preprocessed AnnData object.
    gene_list : list
        List of gene names in the gene set.
    gene_weights : list or np.ndarray
        Corresponding gene weights.
    trait : str
        Trait or gene set name used for labeling outputs.

    Returns
    -------
    pd.DataFrame
        DataFrame containing normalized cell-level enrichment scores.
    """
    print(f"Scoring cells for trait: {trait}...")
    score_df = scdrs.score_cell(
        data=adata,
        gene_list=gene_list,
        gene_weight=gene_weights,
        ctrl_match_key="mean_var",
        n_ctrl=1000,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=True,
        verbose=True,
    )
    return score_df


def run_group_analysis(adata, score_df, group_cols=["celltype_imputed_lowerres", "celltype_imputed"]):
    """
    Performs downstream group-level analysis on normalized cell scores.

    Parameters
    ----------
    adata : AnnData
        AnnData object used for grouping metadata.
    score_df : pd.DataFrame
        DataFrame containing cell-level scores.
    group_cols : list of str, optional
        List of columns in `adata.obs` used to define cell groups.

    Returns
    -------
    dict
        Dictionary with keys as group names and values as group-level statistics.
    """
    print("Running group analysis...")
    stats_dict = scdrs.method.downstream_group_analysis(
        adata=adata,
        df_full_score=score_df,
        group_cols=group_cols,
    )
    return stats_dict


def save_scores(stats_dir, trait, score_df, stats_dict):
    """
    Saves cell scores and downstream group analysis results to disk.

    Parameters
    ----------
    stats_dir : Path
        Directory where results will be saved.
    trait : str
        Trait or gene set name used in file naming.
    score_df : pd.DataFrame
        DataFrame containing normalized cell scores.
    stats_dict : dict
        Dictionary containing group-level statistics.
    """
    stats_dir.mkdir(parents=True, exist_ok=True)
    # Normalized cell scores
    score_df["norm_score"].to_csv(stats_dir / f"{trait}_scores.csv", index=True)
    # Downstream group analysis results
    for group, stats in stats_dict.items():
        stats.to_csv(stats_dir/ f"{trait}_{group}_stats.csv", index=True)
    print(f"Saved results for trait: {trait} to {stats_dir}")


def save_trait_umap_plot(
    adata,
    score_df,
    trait,
    figdir,
    title=None,
    color_map="RdBu_r",
    vmin=-5,
    vmax=5,
):
    """
    Saves UMAP plots of normalized scDRS scores for a specific trait.

    Parameters
    ----------
    adata : AnnData
        AnnData object with UMAP coordinates.
    score_df : pd.DataFrame
        DataFrame with `norm_score` for each cell.
    trait : str
        Trait or gene set name used for coloring and file naming.
    figdir : Path
        Output directory to save plots.
    title : str, optional
        Custom title for the plot. Defaults to trait name.
    color_map : str, optional
        Color map used for the UMAP (default: "RdBu_r").
    vmin : float, optional
        Minimum value for color scale.
    vmax : float, optional
        Maximum value for color scale.
    """

    figdir.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = figdir
    
     # Check if title is provided
    title = title if title else trait
    
    # Add score to adata
    adata.obs[trait] = score_df["norm_score"]
    
    # Save umap plots
    for s in [2, 2.5, 3, 3.5]:
        sc.pl.umap(
            adata,
            color=trait,
            color_map="RdBu_r",
            vmin=-5,
            vmax=5,
            s=s,
            show=False,
            save=f"_{trait}_{s}.png",
            title=title
        )


def parse_args():
    parser = argparse.ArgumentParser(description="Run scDRS scoring on AnnData.")
    parser.add_argument("-d", "--adata", type=str, required=True, help="Path to .h5ad file")
    parser.add_argument("-gs", "--geneset", type=str, required=True, help="Path to gene set file")
    parser.add_argument("-o", "--outdir", type=str, default="scdrs_results", help="Output directory")
    parser.add_argument("-t", "--title", type=str, default=None, help="Title for umap plot. Uses trait name in geneset if none specified. (Default=None)")
    parser.add_argument(
        "-cov",
        "--covariates",
        nargs="+",
        default=["sex_code", "age", "n_genes", "batch_code"],
        help="List of covariates to include (default: sex_code age n_genes batch_code)"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Extract the basename of the input file (without extension) for output folder
    adata_file_basename = Path(args.adata).stem
    gs_file_basename = Path(args.geneset).stem
    output_dir = Path(args.outdir) / adata_file_basename / gs_file_basename
    fig_output_dir = output_dir / "figures"
    stats_dir = output_dir / "stats"

    # Load data
    print(f"Loading data from {args.adata}...")
    adata = sc.read_h5ad(args.adata)

    # Preprocess
    preprocess(adata = adata, covariates = args.covariates)

    # Load gene sets
    print(f"Loading gene sets from {args.geneset}...")
    genesets = load_geneset(args.geneset, adata)

    # NOTE: probably better to force one geneset per run; easier for downstream processing.
    for trait, (gene_list, gene_weights) in genesets.items():
        # Score cells for each trait
        score_df = score_cells(adata, gene_list, gene_weights, trait)
        # Run group analysis
        stats_dict = run_group_analysis(adata, score_df)
        # Save the results
        save_trait_umap_plot(adata, score_df, trait, fig_output_dir, args.title)
        save_scores(stats_dir, trait, score_df, stats_dict)


if __name__ == "__main__":
    main()
