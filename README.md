# Improving scDRS with eQTL and Chromatin Accessibility Data

This repository contains notebooks and scripts related to my bioinformatics graduation project, which focused on enhancing the scDRS(single-cell Disease Relevance Score)<sup>[1]</sup> methodology by integrating expression Quantitative Trait Loci (eQTL) data and open chromatin information (scATAC-seq peaks and ArchR gene scores).

Below are brief usage instructions for the main scripts used in this project with the following functionality:
1. Gene scoring using PascalX<sup>[2]</sup>
2. Retrieving the 1000 most significant genes(lowest p-values) and creating a scDRS compatible gene set
3. Using scDRS on:
   - A single Anndata object (.h5ad file). Default scDRS usage, similar to the [tutorial](https://martinjzhang.github.io/scDRS/notebooks/quickstart.html) provided by the authors.
   - Multiple modalities (e.g. scATAC Gene Scores and scRNA-seq data) with a modified version of scDRS `score_cell` function. 

> [!NOTE]
> Some scripts are tailored for internal (in-house) datasets and may not generalize. However, they can typically be adapted with minimal code changes.

### References
1. Zhang*, Hou*, et al. (2022).  
   *Polygenic enrichment distinguishes disease associations of individual cells in single-cell RNA-seq data*.  
   *Nature Genetics*. https://www.nature.com/articles/s41588-022-01167-z

2. Krefl, D., Brandulas Cammarata, A., & Bergmann, S. (2023).  
   *PascalX: a python library for GWAS gene and pathway enrichment tests*.  
   *Bioinformatics*, btad296. https://doi.org/10.1093/bioinformatics/btad296

---

# Instructions

This guide outlines how to process GWAS summary statistics and run scDRS cell scoring using PascalX and scDRS. 


### Dependencies

Requires Linux OS for full functionality, but main scripts should work on MAC or Windows.

| Package      | Version |
|--------------|---------|
| **Python**   | 3.9.21  |
| scanpy       | 1.10.3  |
| scdrs        | 1.0.2   |
| pandas       | 2.2.3   |
| numpy        | 1.26.4  |
| matplotlib   | 3.9.4   |
| seaborn      | 0.13.2  |
| PascalX      | 0.0.5   |
---

## 1. Processing GWAS Summary Statistics with PascalX
`Src/PascalX/pascalx_gene_scoring.py`

This script scores genes using PascalX based on GWAS summary statistics.

> [!NOTE]
> External variant datasets (e.g., eQTLs) can optionally be included using the PascalX mapperfile feature.  
> Refer to the [PascalX documentation](https://bergmannlab.github.io/PascalX/index.html) for more details.

#### Output
- `<trait>_gscores.tsv` — Gene-level p-values and number of SNPs used per gene.
---
#### Script Arguments

| Argument            | Description                                                                                      | Required | Default        |
|---------------------|--------------------------------------------------------------------------------------------------|----------|----------------|
| `--summ_stats`      | Path to GWAS summary statistics file                                                             | ✅ Yes   | —              |
| `--ref_panel`       | Path to LD reference panel (chrromosome-separated PLINK format)                                                  | ✅ Yes   | —              |
| `--gene_ref`        | Path to gene annotation file                                                                     | ✅ Yes   | —              |
| `--out_dir`         | Output directory                                                                                 | ✅ Yes   | —              |
| `--rscol`           | Column index of variant ID (0-based) in summary statistics                                       | ✅ Yes   | —              |
| `--pcol`            | Column index of p-value (0-based) in summary statistics                                          | ✅ Yes   | —              |
| `--scorer`          | Scoring method to use (e.g., `chi2sum`)                                                          | ❌ No       | `chi2sum`      |
| `--window`          | Gene window size (base pairs)                                                                    | ❌ No       | `50000`        |
| `--cores`           | Number of CPU cores for parallel processing                                                      | ❌ No       | `1`            |
| `--keepfile`        | Optional file listing sample IDs to keep                                                         |❌ No       | `None`         |
| `--mapperfile`      | Optional SNP-to-gene mapping file (tab-delimited: SNP ID, ENSG ID)                               |❌ No       | `None`         |
| `--chr`             | Chromosomes to score; use `"all"` or comma-separated list (e.g., `"1,2,10"`)                     |❌ No       | `"all"`        |

> [!IMPORTANT]
> Gene reference file is currently hardcoded. Edit `Scorer.load_genome` column index parameters in the `main` function of `pascalx_gene_scoring.py` for your own gene location reference.


---

#### 💻 Example Command

```bash
python pascalx_gene_scoring.py \
    --summ_stats path/to/summary_stats.tsv \
    --ref_panel path/to/ld_reference \
    --gene_ref path/to/gene_reference.tsv \
    --rscol 4 \
    --pcol 2 \
    --out_dir results/pascalx_output/ \
    --cores 11 \
    --window 10000 \
    --mapperfile path/to/eqtl_mapper.tsv \
```
---

## 2. Create Gene Set File from PascalX Output
`Src/PascalX/pascalx_munge_gs.sh`

Generates a gene set file containing compatible with scDRS using the PascalX output file as input. Number of genes can optionally be adjusted within the script.

#### Output
- `output_geneset.gs` - Top `n`(default 1000) genes - gene set compatable with scDRS.  

#### Usage

```bash
bash Src/PascalX/pascalx_munge_gs.sh /path/to/{file_name}_gscores.tsv /path/to/output_geneset.gs TRAIT_NAME
```
---
## 3. Run scDRS cell scoring
This section covers how to compute single-cell disease relevance scores using scDRS, based on GWAS-derived gene sets. The scripts support both single-modality (e.g., scRNA-seq) and multi-modality (e.g., scRNA + scATAC) input.  

> [!IMPORTANT]
> All input `.h5ad` files in this section are expected to contain `AnnData.X` matrices that are **log1p-transformed** and **total count normalized**

### 1. Scoring on a single adata object only
`Src/scDRS/scdrs_score_cells.py`

This script performs single-modality scDRS scoring using an input `.h5ad` file and a GWAS trait-associated gene set. This script computes
individual-cell level scDRS scoress, performs downstream group analysis, and generates UMAP plots of the normalized scDRS scores.

#### Main functionality:
- Loads an AnnData object containing single-cell data
- Preprocesses data (optionally correcting for covariates)
- Loads a single gene set (from .gs format)
- Computes scDRS scores per cell
- Performs group-level analysis across user-specified cell annotations
- Saves all output files to a structured output directory
- Generates UMAP plots of scDRS scores

#### Outputs
- `<trait>_scores.csv`  
  Normalized scDRS scores per cell (indexed by cell ID)

- `<trait>_<group>_stats.csv`  
  Group-level trait association and FDR results for each group

- `umap_<trait>_<size>.png`  
  UMAP plots showing per-cell scDRS scores, generated at multiple point sizes (`s=2`, `2.5`, `3`, `3.5`) using the `RdBu_r` colormap. 


#### Recommended Resources
- 40 GB RAM  
- 1 CPU core
- About 2 hours and 30 minutes runtime

---

#### Script Arguments

| Argument         | Description                                                        | Required | Default                   |
|------------------|--------------------------------------------------------------------|----------|---------------------------|
| `-d`, `--adata`  | Path to input AnnData `.h5ad` file with gene expression data       | ✅ Yes      | —                         |
| `-gs`, `--geneset` | Path to gene set file                                             | ✅ Yes      | —                         |
| `-o`, `--outdir` | Output directory                                                   | ❌ No       | `scdrs_results`            |
| `-t`, `--title`  | Title for UMAP plot (uses trait name in gene set if not provided) | ❌ No       | None                      |
| `-cov`, `--covariates` | List of covariates to include (space-separated)              | ❌ No       | `sex_code age n_genes batch_code` |

---

#### 💻 Example Command

```bash
python3 scdrs_score_cells.py \
    -d path/to/adata.h5ad \
    -gs path/to/geneset.gs \
    -o output_directory \
    -t "Trait Name" \
    -cov sex_code age n_genes batch_code
```
---
	
### 2. Scoring cells across multiple modalities (e.g. scRNA-seq and corresponding scATAC Gene scores)
`Src/scDRS/scdrs_score_modalities.py`

Performs multi-modal scDRS analysis using gene expression and ATAC-derived gene scores from single-cell data. This script computes
individual-cell level scDRS scores for scRNA-seq and scATAC Gene Score data, and the combination of both.
Performs downstream group-level enrichment analysis, and generates UMAPs of normalized scDRS scores.

#### Main functionality:
- Loads RNA and ATAC `.h5ad` files
- Loads a single gene set file
- Computes scDRS scores for RNA and ATAC and for a combination of both using a modified version of scDRS `score_cell` function
- Performs group analysis for RNA, ATAC and the combination
- Saves normalized scDRS scores and group analysis results
- Generates UMAP plots of normalized scDRS scores for RNA and ATAC. Also generates a UMAP for the combination if concat-adata is provided.

#### Outputs
For scRNA-seq, scATAC Gene Scores and the combination, generates:

- `<trait>_scores.csv`  
  Normalized scDRS scores per cell (indexed by cell ID)

- `<trait>_<group>_stats.csv`  
  Group-level trait association and FDR results for each group

- `umap_<trait>_<size>.png`  
  UMAP plots showing per-cell scDRS scores, generated at multiple point sizes (`s=2`, `2.5`, `3`, `3.5`) using the `RdBu_r` colormap.
  Only generates UMAPs for the combined modalities if concat-adata argument is provided.

#### Recommended Resources
- 60 GB RAM  
- 3 CPU cores  
- About 3 hours runtime

---
#### Script Arguments

| **Argument**         | **Description**                                                                                       | **Required** | **Example Path / Value**                                                 |
|----------------------|-------------------------------------------------------------------------------------------------------|--------------|---------------------------------------------------------------------------|
| `-rna`, `--rna_file`         | Preprocessed scRNA-seq data (.h5ad) with gene expression                                              | ✅ Yes       | `Data/SCDRS/multiome/sc_gene_expression_pp.h5ad`                        |
| `-atac`,   `--atac_file`        | Preprocessed scATAC-seq data (.h5ad) with gene scores                                                 | ✅ Yes       | `Data/SCDRS/multiome/sc_gene_scores_pp.h5ad`                            |
| `-gs`, `--geneset`          | Gene set file in `.gs` format (e.g., from PascalX output)                                             | ✅ Yes       | `Data/SCDRS/genesets/px_BMI_w10kb.gs`                                     |
| `-o`, `--outdir`           | Output directory where results will be saved                                                          | ❌ No        | `Results/scdrs_multi/` (default: `scdrs_results/`)                        |
| `-cov`, `--covariates`       | (Optional) Covariate file with dataset keys and covariates (no header)                               | ❌ No        | `Data/SCDRS/multiome/covariates.csv`                                     |
| `-concat`,   `--concat-adata`     | (Optional) Concatonated RNA and ATAC `.h5ad` file for plotting scDRS score UMAPs                                     | ❌ No        | `Data/SCDRS/multiome/atac_rna_concat.h5ad`                                |
| `-t`, `--title`            | (Optional) Title used for UMAP plot annotation (overrides trait name in geneset file)                 | ❌ No        | `"Body Mass Index"`                                                      |
| `-c`, `--cores`            | Number of cores used for parallel processing                                                          | ❌ No        | `3` (default)                                                             |

---

#### 💻 Example Command

```bash
python3 scdrs_score_modalities.py \
    --rna_file Data/SCDRS/multiome/ATAC_gene_expression_pp.h5ad \
    --atac_file Data/SCDRS/multiome/ATAC_gene_scores_pp.h5ad \
    --geneset Data/SCDRS/genesets/px_BMI_w10kb.gs \
    --covariates Data/SCDRS/multiome/covariates.csv \
    --concat-adata Data/SCDRS/multiome/atac_rna_concat.h5ad \
    --outdir Results/multiomics/ \
    --title "Body Mass Index" \
    --cores 3
```
