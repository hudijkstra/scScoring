"""
Performs multi-modal scDRS(single-cell Disease Relevance Score) analysis using gene expression and ATAC-derived gene scores from single-cell data. This script computes
individual-cell level scDRS scores for scRNA-seq and scATAC Gene Score data, and the combination of both, using a modified verion of scDRS score_cell function.
Afterwards, performs downstream group-level enrichment analysis, and generates UMAPs of normalized scDRS scores.

Main functionality:
- Loads RNA and ATAC `.h5ad` files
- Loads a single gene set file
- Computes scDRS scores for RNA and ATAC and for a combination of both using a modified version of scDRS `score_cell` function
- Performs group analysis for RNA, ATAC and the combination
- Saves normalized scDRS scores and group analysis results
- Generates UMAP plots of normalized scDRS scores for RNA and ATAC. Also generates a UMAP for the combination if concat-adata is provided.

Usage:
    python scdrs_score_modalities.py \
        -rna /path/to/rna.h5ad \
        -atac /path/to/atac.h5ad \
        -gs /path/to/geneset.gs \
        -cov /path/to/covariates.csv \
        -concat /path/to/concat.h5ad \
        -o /path/to/output_dir \
        -t "Trait Name" \
        -c 3

Author: Hessel Dijkstra

"""

import scdrs
import anndata as an
import scanpy as sc
import pandas as pd
from pathlib import Path
import argparse
import warnings
from score_cell_per_modality import score_cell_per_modality
from score_group_per_modality import score_group_per_modality


warnings.filterwarnings("ignore")
sc.settings.set_figure_params(dpi=300)


def preprocess(data_dict, covariates_file=None, adj_prop="celltype_imputed_lowerres"):
    """
    Preprocesses each AnnData object in the data_dict using scDRS preprocessing.

    If a covariates file is provided, the specified covariates are extracted and used for adjustment.

    Parameters:
        data_dict (dict): Dictionary with modality names as keys and AnnData objects as values.
        covariates_file (str or Path, optional): Path to covariates CSV file. Expected format:
            - No header
            - First column: dataset key matching data_dict keys
            - Remaining columns: covariate names present in `adata.obs`
        adj_prop (str): Column in `adata.obs` used for adjusting for cell group proportion.
    """
    if covariates_file is not None:
        cov_df = pd.read_csv(covariates_file, header=None, index_col=0)

    for key, adata in data_dict.items():
        cov = None
        if covariates_file is not None:
            if key not in cov_df.index:
                raise ValueError(f"{key} not found in covariates file.")
            covariates = cov_df.loc[key].dropna().tolist()

            try:
                cov = adata.obs[covariates].copy()
            except KeyError:
                missing = [c for c in covariates if c not in adata.obs.columns]
                raise KeyError(f"Missing covariates in adata.obs for {key}: {missing}")

        scdrs.preprocess(adata, cov=cov, adj_prop=adj_prop, n_mean_bin=20, n_var_bin=20, copy=False)

  
def load_geneset(data_dict, geneset_path, src_species="human", dst_species="human"):
    """
    Loads a gene set and intersects it with the genes in each AnnData object.

    Parameters:
        data_dict (dict): Dictionary with modality names as keys and AnnData objects as values.
        geneset_path (str or Path): Path to a gene set file in scDRS format (.gs).
        src_species (str): Species of the gene set.
        dst_species (str): Species of the data.

    Returns:
        trait (str): Trait name parsed from the gene set.
    """
    print(f"Loading gene set from {geneset_path}...")
    for modality, adata in data_dict.items():
        gs_dict = scdrs.util.load_gs(
            geneset_path,
            src_species=src_species,
            dst_species=dst_species,
            to_intersect=adata.var_names,
        )

        assert len(gs_dict) == 1, f"Expected exactly one gene set, but found {len(gs_dict)} in {args.geneset}."
        
        for trait, (gene_list, gene_weight) in gs_dict.items():
            data_dict[modality] = [
                adata,
                gene_list,
                gene_weight
            ]
            
    return trait


def save_scores(stats_dir, trait, cell_score_dict, group_stats_dict):
    """
    Saves computed scDRS scores and group-level statistics for each modality.

    Parameters:
        stats_dir (Path): Directory where results will be saved.
        trait (str): Trait name used for labeling.
        cell_score_dict (dict): Cell-level scDRS score results by modality.
        group_stats_dict (dict): Group-level (e.g., cell type) enrichment stats by modality.
    """
    modals = cell_score_dict.keys()
    
    group_stats_dir = stats_dir / "group_stats"
    cell_score_dir = stats_dir / "cell_scores"
    group_stats_dir.mkdir(parents=True, exist_ok=True)
    cell_score_dir.mkdir(parents=True, exist_ok=True)

    for modality in modals:
         # Normalized cell scores
        cell_score_dict[modality]["norm_score"].to_csv(cell_score_dir / f"{modality}_cell_scores.csv", index=True)
        
        # Downstream group analysis results
        for group, stats in group_stats_dict[modality].items():
            stats.to_csv(group_stats_dir / f"{modality}_{group}_stats.csv", index=True)
    
    print(f"Saved stats for trait: {trait} to {stats_dir}")


def generate_modality_umap(
    adata,
    trait,
    modality,
    title=None,
    color_map="RdBu_r",
    vmin=-5,
    vmax=5,
):
    """
    Generates and saves UMAP plots of scDRS scores for a given modality.

    Parameters:
        adata (AnnData): Annotated data object.
        trait (str): Trait used to label color axis.
        modality (str): Name of the modality ("RNA", "ATAC", etc.).
        title (str): Title for the plot.
        color_map (str): Colormap to use.
        vmin (float): Min value for color scaling.
        vmax (float): Max value for color scaling.
    """
    
    title = title if title else trait

    for s in [2, 3]:
        sc.pl.umap(
            adata,
            color=trait,
            color_map=color_map,
            vmin=vmin,
            vmax=vmax,
            s=s,
            show=False,
            save=f"_{modality}_{s}.png",
            title=f"{title} - {modality}"
        )


def parse_args():
    parser = argparse.ArgumentParser(description="Run multiomics scDRS.")
    
    parser.add_argument("-rna", "--rna_file", type=Path, required=True, help="Path to gene expression .h5ad file")
    parser.add_argument("-atac", "--atac_file", type=Path, required=True, help="Path to ATAC gene score .h5ad file")
    parser.add_argument("-concat", "--concat_adata", type=Path, default=None,
        help="Optional path to a precomputed .h5ad AnnData file. If provided, this will be used to generate aggregated modality UMAP of scDRS scores.")
    parser.add_argument("-gs", "--geneset", type=Path, required=True, help="Path to gene set file")
    parser.add_argument("-o", "--outdir", type=Path, default="scdrs_results", help="Output directory")
    parser.add_argument("-cov", "--covariates", type=Path, default=None,
        help="Optional path to a covariates file (no header). First column: dataset key; remaining columns: covariates")

    parser.add_argument("-t", "--title", type=str, default=None, help="Title for umap plot. Uses trait name in geneset if none specified. (Default=None)")
    parser.add_argument("-c", "--cores", type=int, default=3, help="Cores for parallel processing (default=3)")

    return parser.parse_args()


def print_args(args):
    """
    Prints the user-specified arguments in a clean format.
    Automatically adjusts column width for better alignment.
    """
    arg_dict = vars(args)
    max_key_length = max(len(k) for k in arg_dict.keys())  # Find longest argument name

    print("\n### Specified Arguments ###")
    for arg, value in arg_dict.items():
        print(f"{arg.ljust(max_key_length)} : {value}")  
    print("###########################\n")


def main():
    args = parse_args()
    print_args(args)

    # Create output directories
    gs_file_basename = args.geneset.stem
    output_dir = args.outdir / "multiomics" / gs_file_basename
    fig_dir = output_dir / "figures"
    stats_dir = output_dir / "stats"

    print("Loading data...")
    data_dict = {
        "ATAC" : sc.read_h5ad(args.atac_file),
        "RNA" : sc.read_h5ad(args.rna_file)
    } 
    
    preprocess(data_dict, args.covariates)

    # Load geneset and retrieve trait
    trait = load_geneset(data_dict, args.geneset)

    # Compute cell scores
    cell_score_dict = score_cell_per_modality(data_dict, return_ctrl_norm_score=True, n_cores=args.cores)

    # Cell type association analysis
    group_stats_dict = score_group_per_modality(data_dict, cell_score_dict, n_cores=args.cores)

    # Save results
    save_scores(stats_dir, trait, cell_score_dict, group_stats_dict)

    fig_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = fig_dir

    # NOTE: maybe change to just the coordinates
    if args.concat_adata:
        print(f"Loading precomputed data from {args.concat_adata}")
        adata = sc.read_h5ad(args.concat_adata)
        data_dict[list(cell_score_dict.keys())[-1]] = [adata, None, None]


    # Generate UMAP plots
    for modality, (adata, _, _) in data_dict.items():
        adata.obs[trait] = cell_score_dict[modality]["norm_score"]
        generate_modality_umap(
            adata=adata,
            trait=trait,
            modality=modality,
            title=args.title
        )

    print("All modalities processed successfully.\n")


if __name__ == "__main__":
    main()
