"""
This module extends scDRS downstream group analysis to support multiple data modalities.

Main Functions:
- score_group_per_modality: Runs group-level association analysis in parallel across modalities.
- _score_group_for_modality: Helper function to compute FDR counts and empirical statistics
  (Monte Carlo p-value and z-score) for each group.

Returns:
- Dictionary mapping each modality name to a dictionary of pandas DataFrames,
  keyed by grouping column name (e.g., "celltype_imputed"), containing:
  - Number of cells (`n_cell`)
  - Number of control samples (`n_ctrl`)
  - Monte Carlo p-value (`assoc_mcp`)
  - Monte Carlo z-score (`assoc_mcz`)
  - Cell count passing FDR thresholds (`n_fdr_0.05`, etc.)
"""

import numpy as np
import pandas as pd
import time
from statsmodels.stats.multitest import multipletests
from multiprocessing import Pool
import anndata


def score_group_per_modality(
    data_dict,
    modal_score_dict,
    n_cores=3,
    group_cols=["celltype_imputed_lowerres", "celltype_imputed"],
    fdr_thresholds=[0.05, 0.1, 0.2],
):
    """
    Perform group-level disease association analysis per modality.

    Parameters:
        data_dict (dict): Dictionary mapping modality names to tuples of:
                          (AnnData, gene list, gene weights).
                          Used only to get shared cell indices and obs annotations.
        modal_score_dict (dict): Dictionary mapping modality names to DataFrames
                                 with per-cell scores (raw, norm, pval, ctrl_*).
        n_cores (int): Number of parallel processes to use (default: 3).
        group_cols (list of str): Columns in `adata.obs` used for group annotation.
        fdr_thresholds (list of float): FDR significance thresholds.

    Returns:
        dict: Dictionary mapping modality name to a dictionary of DataFrames.
              Each inner dictionary maps `group_col` → results DataFrame with:
              - n_cell: number of cells in group
              - n_ctrl: number of Monte Carlo control samples
              - assoc_mcp: empirical p-value for association
              - assoc_mcz: z-score based on empirical null
              - n_fdr_<threshold>: number of cells with FDR below threshold
    """

    adata = next(iter(data_dict.values()))[0]
    obs_df = adata.obs.copy()
    
    assert (
        len(set(group_cols) - set(obs_df.columns)) == 0
    ), "Missing `group_cols` variables from `adata.obs.columns`."

    args_list = [
        (modality, df_score, obs_df, group_cols, fdr_thresholds)
        for modality, df_score in modal_score_dict.items()
    ]
    
    start_time = time.time()
    print("Analysing groups per modality...")
    with Pool(processes=n_cores) as pool:
        results = pool.map(_score_group_for_modality, args_list)

    end_time = time.time()
    print(f"Done in {end_time - start_time:.2f} seconds.")
    
    return dict(results)


def _score_group_for_modality(args):
    """
    Compute group-level disease association for a single modality.

    Parameters:
        args (tuple): A tuple of:
            - modality (str): Name of the modality.
            - df_full_score (pd.DataFrame): Cell-level score DataFrame with columns:
                - 'norm_score', 'pval', and many 'ctrl_norm_score_*' columns.
            - obs_df (pd.DataFrame): Metadata (adata.obs) used for grouping.
            - group_cols (list of str): obs columns to group by.
            - fdr_thresholds (list of float): Significance cutoffs.

    Returns:
        tuple: (modality_name, result_dict), where result_dict maps
               each group_col to a result DataFrame indexed by group name.
               Each result DataFrame includes:
               - n_cell: total number of cells in group
               - n_ctrl: number of Monte Carlo controls
               - assoc_mcp: empirical p-value (quantile test)
               - assoc_mcz: Monte Carlo z-score
               - n_fdr_<threshold>: number of cells with FDR < threshold
    """

    modality, df_full_score, obs_df, group_cols, fdr_thresholds = args

    cell_list = sorted(set(obs_df.index) & set(df_full_score.index))
    
    control_list = [x for x in df_full_score.columns if x.startswith("ctrl_norm_score")]
    n_ctrl = len(control_list)
    df_reg = obs_df.loc[cell_list, group_cols].copy()
    df_reg = df_reg.join(
        df_full_score.loc[cell_list, ["norm_score"] + control_list + ["pval"]]
    )

    dict_df_res = {}
    for group_col in group_cols:
        group_list = sorted(set(obs_df[group_col]))
        res_cols = ["n_cell", "n_ctrl", "assoc_mcp", "assoc_mcz"]
     
        for fdr_threshold in fdr_thresholds:
            res_cols.append(f"n_fdr_{fdr_threshold}")

        df_res = pd.DataFrame(index=group_list, columns=res_cols, dtype=np.float32)
        df_res.index.name = "group"

        df_fdr = pd.DataFrame(
            {"fdr": multipletests(df_reg["pval"].values, method="fdr_bh")[1]},
            index=df_reg.index,
        )

        for group in group_list:
            group_cell_list = list(df_reg.index[df_reg[group_col] == group])
            df_res.loc[group, ["n_cell", "n_ctrl"]] = [len(group_cell_list), n_ctrl]

            # Number of FDR < fdr_threshold cells in each group
            for fdr_threshold in fdr_thresholds:
                df_res.loc[group, f"n_fdr_{fdr_threshold}"] = (
                    df_fdr.loc[group_cell_list, "fdr"].values < fdr_threshold
                ).sum()
        
            # Association stats    
            score_q95 = np.quantile(df_reg.loc[group_cell_list, "norm_score"], 0.95)
            v_ctrl_score_q95 = np.quantile(
                df_reg.loc[group_cell_list, control_list], 0.95, axis=0
            )
            
            mc_p = ((v_ctrl_score_q95 >= score_q95).sum() + 1) / (
                v_ctrl_score_q95.shape[0] + 1
            )
            mc_z = (score_q95 - v_ctrl_score_q95.mean()) / v_ctrl_score_q95.std()
            df_res.loc[group, ["assoc_mcp", "assoc_mcz"]] = [mc_p, mc_z]
            
    
        dict_df_res[group_col] = df_res

    return modality, dict_df_res

