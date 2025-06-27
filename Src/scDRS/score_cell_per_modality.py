"""
Multi-modality scDRS scoring module.

This module defines functions for computing disease relevance scores at the single-cell level
across multiple data modalities using a modified version of scDRS score_cell function. It supports parallel computation, 
shared gene set intersection, and empirical p-value estimation using Monte Carlo control gene sets.

Main adaptaions:
- score_cell_per_modality: Main entry to compute cell-level disease scores across modalities.
- _intersect_genes: Align gene lists across modalities.
- _get_raw_stats: Compute raw and control scores for a given modality.
- _aggregate_modalities: Aggregate results across modalities by averaging.

Returns:
- Dictionary of DataFrames, each keyed by modality name (including combined scores),
  with raw scores, normalized scores, p-values, and z-scores for each cell.
"""

import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import anndata
import time
from typing import List, Dict, Tuple
from statsmodels.stats.multitest import multipletests
import scdrs
from multiprocessing import Pool


def score_cell_per_modality(
    data_dict, # adata_dict: Dict[str, Tuple[anndata.AnnData, List[str], List[float]]],
    ctrl_match_key="mean_var",
    n_cores=2,
    n_ctrl=1000,
    n_genebin=200,
    weight_opt="vs",
    copy=False,
    return_ctrl_raw_score=False,
    return_ctrl_norm_score=False,
    random_seed=0,
    verbose=False,
    save_intermediate=None,
):
    """
    Scores cells across multiple modalities using scDRS.

    Parameters:
        data_dict (dict): Dictionary where keys are modality names and values are tuples:
                          (AnnData, list of gene names, list of gene weights).
        ctrl_match_key (str): Gene-level statistic used to match control genes (default: "mean_var").
        n_cores (int): Number of CPU cores to use for parallel processing.
        n_ctrl (int): Number of Monte Carlo control sets.
        n_genebin (int): Number of gene bins for control gene sampling.
        weight_opt (str): Weighting method: one of {"uniform", "vs", "inv_std", "od"}.
        copy (bool): If True, copies AnnData before modification.
        return_ctrl_raw_score (bool): If True, includes all raw control scores in output.
        return_ctrl_norm_score (bool): If True, includes all normalized control scores in output.
        random_seed (int): Seed for reproducibility.
        verbose (bool): If True, print verbose output.
        save_intermediate (str | None): Optional path to save intermediate files.

    Returns:
        dict: Dictionary of pandas DataFrames indexed by cell ID. Each value contains:
              - raw_score
              - norm_score
              - mc_pval
              - pval
              - nlog10_pval
              - zscore
              - (optional) control scores
    """
    # Checking cell ids
    keys = list(data_dict.keys())
    ref_idx = data_dict[keys[0]][0].obs.index

    for key in keys[1:]:
        adata_idx = data_dict[key][0].obs.index
        assert ref_idx.equals(adata_idx), f"Cell IDs mismatch between '{keys[0]}' and '{key}'"
    
    print(f"Scoring cells for modalities: {', '.join(str(k) for k in data_dict.keys())}")

    # Find gene overlaps between modalities:
    _intersect_genes(data_dict)

    tasks = [
        (
        name,
        adata,
        gene_list,
        gene_weight,
        ctrl_match_key,
        n_ctrl,
        n_genebin,
        weight_opt,
        copy,
        random_seed,
        verbose,
        save_intermediate
        ) for name, (adata, gene_list, gene_weight) in data_dict.items()
    ]
    
    start = time.time()
    
    print(f"Processing {len(tasks)} modalities with {n_cores} cores...")
    with Pool(processes=n_cores) as pool:
        modality_results = pool.map(_get_raw_stats, tasks)

    end = time.time()
    elapsed = end - start
    minutes = int(elapsed // 60)
    seconds = int(elapsed % 60)
    print(f"Control cells generated for in {minutes} min {seconds} sec")
    
    print("Aggregating scores...")
    modality_results = _aggregate_modalities(modality_results)
    
    print("Gathering results...")
    all_df_res = {}
    for result in modality_results:
        
        name = result["name"]
        v_raw_score = result["v_raw_score"]
        mat_ctrl_raw_score = result["mat_ctrl_raw_score"]
        v_var_ratio_c2t = result["v_var_ratio_c2t"]
        
        v_norm_score, mat_ctrl_norm_score = _correct_background(
            v_raw_score,
            mat_ctrl_raw_score,
            v_var_ratio_c2t,
            save_intermediate=save_intermediate,
        )

        # Get p-values
        mc_p = (1 + (mat_ctrl_norm_score.T >= v_norm_score).sum(axis=0)) / (1 + n_ctrl)
        pooled_p = _get_p_from_empi_null(v_norm_score, mat_ctrl_norm_score.flatten())
        nlog10_pooled_p = -np.log10(pooled_p)
        pooled_z = -sp.stats.norm.ppf(pooled_p).clip(min=-10, max=10)

        # Gather results
        dic_res = {
            "raw_score": v_raw_score,
            "norm_score": v_norm_score,
            "mc_pval": mc_p,
            "pval": pooled_p,
            "nlog10_pval": nlog10_pooled_p,
            "zscore": pooled_z,
        }

        if return_ctrl_raw_score:
            for i in range(n_ctrl):
                dic_res["ctrl_raw_score_%d" % i] = mat_ctrl_raw_score[:, i]
        if return_ctrl_norm_score:
            for i in range(n_ctrl):
                dic_res["ctrl_norm_score_%d" % i] = mat_ctrl_norm_score[:, i]
        
        # Convert to dataframe
        df_res = pd.DataFrame(index=ref_idx, data=dic_res, dtype=np.float32)
    
        all_df_res[name] = df_res
    

    return all_df_res


def _intersect_genes(data_dict):
    """
    Intersects gene lists across all modalities and updates them to a shared gene set.

    Parameters:
        data_dict (dict): Dictionary of modalities, where each value is a tuple:
                          (AnnData, gene list, gene weights).

    Returns:
        dict: Updated data_dict where each modality's gene list and weights are restricted
              to the shared set.
    """

    gene_sets = []
    reference_gene_weights = None

    for _, gene_list, gene_weights in data_dict.values():
        gene_sets.append(set(gene_list))

        if reference_gene_weights is None:
            reference_gene_weights = {gene: weight for gene, weight in zip(gene_list, gene_weights)}

    # Compute shared genes across all sets
    shared_genes = set.intersection(*gene_sets)
    
    if len(shared_genes) == 0:
        raise ValueError("No shared genes found across the provided gene lists.")
    else:
        print(f"Found {len(shared_genes)} shared genes across modalities.")

    # Extract weights for shared genes from the reference set
    shared_gene_list = []
    shared_weights = []

    for gene in shared_genes:
        if gene in reference_gene_weights:
            shared_gene_list.append(gene)
            shared_weights.append(reference_gene_weights[gene])

    # Update the gene lists and weights for all modalities
    for modality in data_dict.keys():
        data_dict[modality][1] = shared_gene_list
        data_dict[modality][2] = shared_weights

    return data_dict
    

def _get_raw_stats(params):
    (   
        name,
        data,
        gene_list,
        gene_weight,
        ctrl_match_key,
        n_ctrl,
        n_genebin,
        weight_opt,
        copy,
        random_seed,
        verbose,
        save_intermediate
    ) = params
    """
    Compute raw and control scores for a single modality.

    Parameters:
        params (tuple): Contains modality name, AnnData, gene list, gene weights, control 
                        matching parameters, etc.

    Returns:
        dict: Dictionary with:
              - 'name': Modality name
              - 'v_raw_score': Raw score per cell
              - 'mat_ctrl_raw_score': Matrix of control raw scores [cells x n_ctrl]
              - 'v_var_ratio_c2t': Variance ratio for normalization
    """

    np.random.seed(random_seed)
    adata = data.copy() if copy else data
    n_cell, n_gene = adata.shape

    # Check preprocessing information
    assert (
        "SCDRS_PARAM" in adata.uns
    ), "adata.uns['SCDRS_PARAM'] not found, run `scdrs.pp.preprocess` first"

    # Check GENE_STATS from adata.uns["SCDRS_PARAM"]
    assert (
        "GENE_STATS" in adata.uns["SCDRS_PARAM"]
    ), "adata.uns['SCDRS_PARAM']['GENE_STATS'] not found, run `scdrs.pp.preprocess` first"

    gene_stats_set_expect = {"mean", "var", "var_tech"}
    gene_stats_set = set(adata.uns["SCDRS_PARAM"]["GENE_STATS"])
    assert (
        len(gene_stats_set_expect - gene_stats_set) == 0
    ), "One of 'mean', 'var', 'var_tech' not found in adata.uns['SCDRS_PARAM']['GENE_STATS'], run `scdrs.pp.preprocess` first"

    # Check if ctrl_match_key is in GENE_STATS
    assert ctrl_match_key in adata.uns["SCDRS_PARAM"]["GENE_STATS"], (
        "ctrl_match_key=%s not found in adata.uns['SCDRS_PARAM']['GENE_STATS']"
        % ctrl_match_key
    )

    # Check if weight_opt is legal
    assert weight_opt in {"uniform", "vs", "inv_std", "od"}, (
        "weight_opt=%s is not one of {'uniform', 'vs', 'inv_std', 'od}'" % weight_opt
    )

    if verbose:
        msg = "# score_modalities_per_cell summary:"
        msg += "\n    n_cell=%d, n_gene=%d," % (n_cell, n_gene)
        msg += "\n    n_disease_gene=%d," % len(gene_list)
        msg += "\n    n_ctrl=%d, n_genebin=%d," % (n_ctrl, n_genebin)
        msg += "\n    ctrl_match_key='%s'," % ctrl_match_key
        msg += "\n    weight_opt='%s'," % weight_opt
        msg += "\n    random_seed=%d, verbose=%s," % (random_seed, verbose)
        msg += "\n    save_intermediate=%s," % save_intermediate
        print(msg)

    # Load parameters
    flag_sparse = adata.uns["SCDRS_PARAM"]["FLAG_SPARSE"]
    flag_cov = adata.uns["SCDRS_PARAM"]["FLAG_COV"]

    df_gene = adata.uns["SCDRS_PARAM"]["GENE_STATS"].loc[adata.var_names].copy()
    df_gene["gene"] = df_gene.index
    df_gene.drop_duplicates(subset="gene", inplace=True)

    gene_list = list(gene_list)
    if gene_weight is not None:
        gene_weight = list(gene_weight)
    else:
        gene_weight = [1] * len(gene_list)

    # Overlap gene_list with df_gene["gene"]
    dic_gene_weight = {x: y for x, y in zip(gene_list, gene_weight)}
    gene_list = sorted(set(gene_list) & set(df_gene["gene"]))
    gene_weight = [dic_gene_weight[x] for x in gene_list]

    # Select control gene sets
    dic_ctrl_list, dic_ctrl_weight = _select_ctrl_geneset(
        df_gene, gene_list, gene_weight, ctrl_match_key, n_ctrl, n_genebin, random_seed
    )

    # Compute raw scores
    v_raw_score, v_score_weight = _compute_raw_score(
        adata, gene_list, gene_weight, weight_opt
    )

    # Compute control raw scores
    print(f"Generating {n_ctrl} control cells for {name}...")

    mat_ctrl_raw_score = np.zeros([n_cell, n_ctrl])
    mat_ctrl_weight = np.zeros([len(gene_list), n_ctrl])
    for i_ctrl in range(n_ctrl):
        v_ctrl_raw_score, v_ctrl_weight = _compute_raw_score(
            adata, dic_ctrl_list[i_ctrl], dic_ctrl_weight[i_ctrl], weight_opt
        )
        mat_ctrl_raw_score[:, i_ctrl] = v_ctrl_raw_score
        mat_ctrl_weight[:, i_ctrl] = v_ctrl_weight

    # Compute normalized scores
    v_var_ratio_c2t = np.ones(n_ctrl)
    if (ctrl_match_key == "mean_var") & (weight_opt in ["uniform", "vs", "inv_std"]):
        # For mean_var matched control genes and raw scores computed as weighted average,
        # estimate variance ratio assuming independence.
        for i_ctrl in range(n_ctrl):
            v_var_ratio_c2t[i_ctrl] = (
                df_gene.loc[dic_ctrl_list[i_ctrl], "var"]
                * mat_ctrl_weight[:, i_ctrl] ** 2
            ).sum()
        v_var_ratio_c2t /= (df_gene.loc[gene_list, "var"] * v_score_weight ** 2).sum()


    return {
        "name": name,
        "v_raw_score": v_raw_score,
        "mat_ctrl_raw_score": mat_ctrl_raw_score,
        "v_var_ratio_c2t": v_var_ratio_c2t
        }


def _aggregate_modalities(modality_results):
    """
    Averages raw and control scores across all modalities.

    Parameters:
        modality_results (list): List of dictionaries returned by _get_raw_stats,
                                 one for each modality.

    Returns:
        list: Updated modality_results with an additional entry "name" containing
              the combined average of all modalities.
    """
    
    # Assume all shapes are consistent
    n = len(modality_results)

    # Initialize variables
    mat_sum = np.zeros_like(modality_results[0]["mat_ctrl_raw_score"])
    v_raw_sum = np.zeros_like(modality_results[0]["v_raw_score"])
    v_var_sum = np.zeros_like(modality_results[0]["v_var_ratio_c2t"])
    name = ""

    for res in modality_results:
        mat_sum += res["mat_ctrl_raw_score"]
        v_raw_sum += res["v_raw_score"]
        v_var_sum += res["v_var_ratio_c2t"]
        name += f"+{res['name']}" if name else res["name"]
        
    # Compute the mean
    avg = {
        "name" : name,
        "v_raw_score" : v_raw_sum / n,
        "mat_ctrl_raw_score" : mat_sum / n,
        "v_var_ratio_c2t" : v_var_sum / n
    }
    
    modality_results.append(avg)
    return modality_results


def _select_ctrl_geneset(
    input_df_gene,
    gene_list,
    gene_weight,
    ctrl_match_key,
    n_ctrl,
    n_genebin,
    random_seed,
):

    """Subroutine for `scdrs.method.score_cell`. Select control gene sets that match
    the disease gene set by `ctrl_match_key`.

    It recognizes `ctrl_match_key` as categorical if the number of unique values is
    less than 10% of the total number of values, and otherwise continuous. For
    categorical `ctrl_match_key`, genes are matched within each category. For continuous
    `ctrl_match_key`, genes are divided into `n_genebin` bins and are matched within
    each bin. A matched control gene takes the same weight as the disease gene,

    Args
    ----
    input_df_gene : pd.DataFrame
        Gene-level statistics of shape (n_gene, n_stats).
    gene_list : list
        Disease gene list of length n_disease_gene.
    gene_weight : list
        Gene weights of length n_disease_gene for genes in the gene_list.
    ctrl_match_key : str
        Gene-level statistic used for matching control and disease genes;
        should be in `input_df_gene`.
    n_ctrl : int
        Number of control gene sets.
    n_genebin : int
        Number of bins for dividing genes by ctrl_match_key if
        `input_df_gene[ctrl_match_key]` is a continuous variable.
    random_seed : int
        Random seed.

    Returns
    -------
    dic_ctrl_list : dict of lists
        dic_ctrl_list[i]: the i-th control gene list
    dic_ctrl_weight : dict of lists
        dic_ctrl_weight[i]: weights for the i-th control gene list

    """

    np.random.seed(random_seed)
    df_gene = input_df_gene.copy()
    if "gene" not in df_gene:
        df_gene["gene"] = df_gene.index
    disease_gene_set = set(gene_list)
    dic_gene_weight = {x: y for x, y in zip(gene_list, gene_weight)}

    # Divide genes into equal-sized bins based on ctrl_match_key
    if df_gene[ctrl_match_key].unique().shape[0] < df_gene.shape[0] / 10:
        df_gene_bin = df_gene.groupby(ctrl_match_key).agg({"gene": set})
    else:
        df_gene["qbin"] = pd.qcut(
            df_gene[ctrl_match_key], q=n_genebin, labels=False, duplicates="drop"
        )
        df_gene_bin = df_gene.groupby("qbin").agg({"gene": set})

    # Find ctrl_match_key matched control genes
    dic_ctrl_list = {x: [] for x in range(n_ctrl)}
    dic_ctrl_weight = {x: [] for x in range(n_ctrl)}
    for bin_ in df_gene_bin.index:
        bin_gene = sorted(df_gene_bin.loc[bin_, "gene"])
        bin_disease_gene = sorted(df_gene_bin.loc[bin_, "gene"] & disease_gene_set)
        if len(bin_disease_gene) > 0:
            for i_list in np.arange(n_ctrl):
                dic_ctrl_list[i_list].extend(
                    np.random.choice(
                        bin_gene, size=len(bin_disease_gene), replace=False
                    )
                )
                dic_ctrl_weight[i_list].extend(
                    [dic_gene_weight[x] for x in bin_disease_gene]
                )

    return dic_ctrl_list, dic_ctrl_weight


def _compute_raw_score(adata, gene_list, gene_weight, weight_opt):
    """Compute raw score
        v_score_weight = gene_weight * {uniform/vs/inv_std}
        `SCDRS_PARAM` is assumed to have been computed using `sparse_reg_out`

    Parameters
    ----------
    adata : anndata.AnnData
        Single-cell data of shape (n_cell, n_gene). Assumed
        to be size-factor-normalized and log1p-transformed.
    gene_list : list
        Disease gene list of length n_disease_gene.
    gene_weight : list
        Gene weights of length n_disease_gene for genes in the gene_list.
    weight_opt : str
        Option for computing the raw score
        - 'uniform': average over the genes in the gene_list.
        - 'vs': weighted average with weights equal to 1/sqrt(technical_variance_of_logct).
        - 'inv_std': weighted average with weights equal to 1/std.
        - 'od': overdispersion score.

    Returns
    -------
    v_raw_score : np.ndarray
        Raw score of shape (n_cell,).
    v_score_weight : np.ndarray
        Gene weights of shape (n_disease_gene,).

    Notes
    -----

    """

    gene_list = list(gene_list)
    gene_weight = list(gene_weight)

    assert weight_opt in {"uniform", "vs", "inv_std", "od"}, (
        "weight_opt=%s is not one of {'uniform', 'vs', 'inv_std', 'od}'" % weight_opt
    )

    # Compute overdispersion score
    # (used only for benchmarking, do not support implicit covariate correction mode)
    if weight_opt == "od":
        return _compute_overdispersion_score(adata, gene_list, gene_weight)

    # Compute other weighted average scores
    assert (
        "SCDRS_PARAM" in adata.uns
    ), "adata.uns['SCDRS_PARAM'] not found, run `scdrs.pp.preprocess` first"

    df_gene = adata.uns["SCDRS_PARAM"]["GENE_STATS"]
    flag_sparse = adata.uns["SCDRS_PARAM"]["FLAG_SPARSE"]
    flag_cov = adata.uns["SCDRS_PARAM"]["FLAG_COV"]

    if weight_opt == "uniform":
        v_score_weight = np.ones(len(gene_list))
    if weight_opt == "vs":
        v_score_weight = 1 / np.sqrt(df_gene.loc[gene_list, "var_tech"].values + 1e-2)
    if weight_opt == "inv_std":
        v_score_weight = 1 / np.sqrt(df_gene.loc[gene_list, "var"].values + 1e-2)

    if gene_weight is not None:
        v_score_weight = v_score_weight * np.array(gene_weight)
    v_score_weight = v_score_weight / v_score_weight.sum()

    if flag_sparse and flag_cov:
        # Implicit covariate correction mode
        cell_list = list(adata.obs_names)
        cov_list = list(adata.uns["SCDRS_PARAM"]["COV_MAT"])
        cov_mat = adata.uns["SCDRS_PARAM"]["COV_MAT"].loc[cell_list, cov_list].values
        cov_beta = (
            adata.uns["SCDRS_PARAM"]["COV_BETA"].loc[gene_list, cov_list].values.T
        )
        gene_mean = adata.uns["SCDRS_PARAM"]["COV_GENE_MEAN"].loc[gene_list].values

        # Compute v_raw_score = transformed_X @ v_score_weight
        # where transformed_X = adata.X + cov_mat @ cov_beta + gene_mean
        v_raw_score = (
            adata[:, gene_list].X.dot(v_score_weight)
            + cov_mat @ (cov_beta @ v_score_weight)
            + gene_mean @ v_score_weight
        ).flatten()
    else:
        # Normal mode
        v_raw_score = adata[:, gene_list].X.dot(v_score_weight).reshape([-1])

    return v_raw_score, v_score_weight


def _compute_overdispersion_score(adata, gene_list, gene_weight):
    """Compute overdispersion score

        Raw weight: w_g_raw = gene_weight / \sigma_{tech,g}^2

        Normalized weight: w_g = w_g_raw / \sum_g w_g_raw

        Overdispersion score: s_c = \sum_g w_g * [(X_cg - \mu_g)^2 - \sigma_{tech,g}^2]

    Args
    ----
    adata : anndata.AnnData
        Single-cell data of shape (n_cell, n_gene). Assumed
        to be size-factor-normalized and log1p-transformed.
    gene_list : list
        Disease gene list of length n_disease_gene.
    gene_weight : list
        Gene weights of length n_disease_gene for genes in the gene_list.

    Returns
    -------
    v_raw_score : np.ndarray
        Raw score of shape (n_cell,).
    v_score_weight : np.ndarray
        Gene weights of shape (n_disease_gene,).
    """

    gene_list = list(gene_list)
    gene_weight = list(gene_weight)

    assert (
        "SCDRS_PARAM" in adata.uns
    ), "adata.uns['SCDRS_PARAM'] not found, run `scdrs.pp.preprocess` first"

    # Mode check
    flag_sparse = adata.uns["SCDRS_PARAM"]["FLAG_SPARSE"]
    flag_cov = adata.uns["SCDRS_PARAM"]["FLAG_COV"]
    if flag_sparse and flag_cov:
        cell_list = list(adata.obs_names)
        cov_list = list(adata.uns["SCDRS_PARAM"]["COV_MAT"])
        mat_X = (
            adata[:, gene_list].X.toarray()
            + adata.uns["SCDRS_PARAM"]["COV_MAT"]
            .loc[cell_list, cov_list]
            .values.dot(
                adata.uns["SCDRS_PARAM"]["COV_BETA"].loc[gene_list, cov_list].values.T
            )
            + adata.uns["SCDRS_PARAM"]["COV_GENE_MEAN"].loc[gene_list].values
        )
    else:
        mat_X = adata[:, gene_list].X

    v_mean = adata.uns["SCDRS_PARAM"]["GENE_STATS"].loc[gene_list, "mean"].values
    v_var_tech = (
        adata.uns["SCDRS_PARAM"]["GENE_STATS"].loc[gene_list, "var_tech"].values
    )

    v_w = 1 / (v_var_tech + 1e-2)
    if gene_weight is not None:
        v_w = v_w * np.array(gene_weight)
    v_w = v_w / v_w.sum()

    # Compute overdispersion score
    if sp.sparse.issparse(mat_X):
        v_raw_score = mat_X.power(2).dot(v_w).reshape([-1])  # Quadratic term
    else:
        v_raw_score = (mat_X ** 2).dot(v_w).reshape([-1])  # Quadratic term
    v_raw_score = v_raw_score - mat_X.dot(2 * v_w * v_mean).reshape([-1])  # Linear term
    v_raw_score = (
        v_raw_score + (v_w * (v_mean ** 2 - v_var_tech)).sum()
    )  # Constant term

    return v_raw_score, np.ones(len(gene_list))


def _correct_background(
    v_raw_score, mat_ctrl_raw_score, v_var_ratio_c2t, save_intermediate=None
):
    """Cell-wise and gene-wise background correction

    Args
    ----
    v_raw_score : np.ndarray
        Disease raw score of shape (n_cell,).
    mat_ctrl_raw_score : np.ndarray
        Disease raw control scores of shape (n_cell,n_ctrl).
    v_var_ratio_c2t : np.ndarray
        Ratio of independent variance between control scores and disease score,
        of shape (n_ctrl).
    save_intermediate : str
        File path prefix for saving intermediate results.

    Returns
    -------
    v_norm_score : np.ndarray
        Normalized disease score of shape (n_cell,)
    mat_ctrl_norm_score : np.ndarray
        Normalized control scores of shape (n_cell,n_ctrl).
    """

    if save_intermediate is not None:
        np.savetxt(
            save_intermediate + ".raw_score.tsv.gz",
            v_raw_score,
            fmt="%.9e",
            delimiter="\t",
        )
        np.savetxt(
            save_intermediate + ".ctrl_raw_score.tsv.gz",
            mat_ctrl_raw_score,
            fmt="%.9e",
            delimiter="\t",
        )

    # Zero-values are assigned the smallest values at the end
    ind_zero_score = v_raw_score == 0
    ind_zero_ctrl_score = mat_ctrl_raw_score == 0

    # First gene set alignment: mean 0 and same independent variance
    v_raw_score = v_raw_score - v_raw_score.mean()
    mat_ctrl_raw_score = mat_ctrl_raw_score - mat_ctrl_raw_score.mean(axis=0)
    mat_ctrl_raw_score = mat_ctrl_raw_score / np.sqrt(v_var_ratio_c2t)
    if save_intermediate is not None:
        np.savetxt(
            save_intermediate + ".raw_score.1st_gs_alignment.tsv.gz",
            v_raw_score,
            fmt="%.9e",
            delimiter="\t",
        )
        np.savetxt(
            save_intermediate + ".ctrl_raw_score.1st_gs_alignment.tsv.gz",
            mat_ctrl_raw_score,
            fmt="%.9e",
            delimiter="\t",
        )

    # Cell-wise standardization
    v_mean = mat_ctrl_raw_score.mean(axis=1)
    v_std = mat_ctrl_raw_score.std(axis=1)
    v_norm_score = v_raw_score.copy()
    v_norm_score = (v_norm_score - v_mean) / v_std
    mat_ctrl_norm_score = ((mat_ctrl_raw_score.T - v_mean) / v_std).T
    if save_intermediate is not None:
        np.savetxt(
            save_intermediate + ".raw_score.cellwise_standardization.tsv.gz",
            v_norm_score,
            fmt="%.9e",
            delimiter="\t",
        )
        np.savetxt(
            save_intermediate + ".ctrl_raw_score.cellwise_standardization.tsv.gz",
            mat_ctrl_norm_score,
            fmt="%.9e",
            delimiter="\t",
        )

    # Second gene set alignment: mean 0
    v_norm_score = v_norm_score - v_norm_score.mean()
    mat_ctrl_norm_score = mat_ctrl_norm_score - mat_ctrl_norm_score.mean(axis=0)
    if save_intermediate is not None:
        np.savetxt(
            save_intermediate + ".raw_score.2nd_gs_alignment.tsv.gz",
            v_norm_score,
            fmt="%.9e",
            delimiter="\t",
        )
        np.savetxt(
            save_intermediate + ".ctrl_raw_score.2nd_gs_alignment.tsv.gz",
            mat_ctrl_norm_score,
            fmt="%.9e",
            delimiter="\t",
        )

    # Set cells with raw_score=0 to the minimum norm_score value
    norm_score_min = min(v_norm_score.min(), mat_ctrl_norm_score.min())
    v_norm_score[ind_zero_score] = norm_score_min - 1e-3
    mat_ctrl_norm_score[ind_zero_ctrl_score] = norm_score_min
    if save_intermediate is not None:
        np.savetxt(
            save_intermediate + ".raw_score.final.tsv.gz",
            v_norm_score,
            fmt="%.9e",
            delimiter="\t",
        )
        np.savetxt(
            save_intermediate + ".ctrl_raw_score.final.tsv.gz",
            mat_ctrl_norm_score,
            fmt="%.9e",
            delimiter="\t",
        )

    return v_norm_score, mat_ctrl_norm_score


def _get_p_from_empi_null(v_t, v_t_null):
    """Compute p-value from empirical null
    For score T and a set of null score T_1,...T_N, the p-value is

        p= [1 + \Sigma_{i=1}^N 1_{ (T_i \geq T) }] / (1+N)

    If T, T_1, ..., T_N are i.i.d. variables following a null distritbuion,
    then p is super-uniform.

    The naive algorithm is N^2. Here we provide an O(N log N) algorithm to
    compute the p-value for each of the N elements in v_t

    Args
    ----
    v_t : np.ndarray
        Observed score of shape (M,).
    v_t_null : np.ndarray
        Null scores of shape (N,).

    Returns
    -------
    v_p: : np.ndarray
        P-value for each element in v_t of shape (M,).
    """

    v_t = np.array(v_t)
    v_t_null = np.array(v_t_null)

    v_t_null = np.sort(v_t_null)
    v_pos = np.searchsorted(v_t_null, v_t, side="left")
    v_p = (v_t_null.shape[0] - v_pos + 1) / (v_t_null.shape[0] + 1)
    return v_p