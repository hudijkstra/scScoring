"""
This module provides utilities for running and managing gene-level association analyses using the 
MAGMA tool. It includes functions for:

- Running gene annotation with sliding windows of varying sizes.
- Annotating and analyzing genes with reference and GWAS summary statistics.
- Submitting SLURM batch jobs for MAGMA gene-level analyses.

Usage:
# Annotate using multiple window sizes
result_df = run_window_analysis(
    snploc_file=Path("gwas_snps.loc"),
    geneloc_file=Path("genes.loc"),
    window_sizes=[0, 10_000, 20_000]
)
print(result_df)

# Run a local MAGMA gene analysis
run_gene_analysis(
    reference="1000G_EUR",
    gene_anno_file="my_annot.genes.annot",
    pval_file="gwas.tsv",
    out_file="output/gene_results",
    trait="height",
    N=500000
)

# Submit batch jobs to SLURM
batch_gene_analysis(
    trait="height",
    samples="N=500000",
    n_batches=10,
    gene_anno_file=Path("my_annot.genes.annot"),
    pval_file=Path("gwas.tsv"),
    out_dir=Path("magma_out"),
    reference=Path("1000G_EUR"),
)


Author: Hessel Dijkstra
"""

import subprocess as sp
import pandas as pd
import os
import re
from pathlib import Path

def run_window_analysis(
    snploc_file: Path,
    geneloc_file: Path,
    window_sizes: list,
):
    """ """

    annotation_results = pd.DataFrame(
        columns=["Window size", "Mapped SNPs", "Mapped genes (total genes)"]
    )

    print(f"Running analysis for {snploc_file}")
    for window_size in window_sizes:
        command = f"magma --annotate window={window_size} --snp-loc {snploc_file} --gene-loc {geneloc_file} --out window_test"
        # Run the process and capture the output
        result = sp.run(command, shell=True, stdout=sp.PIPE, stderr=sp.PIPE, text=True)
        
        # Initialize variables to store results
        mapped_snps, mapped_genes, snp_percentage, total_genes = None, None, None, None

        # Process the output after the command has completed
        for line in result.stdout.splitlines():
            # Capture SNPs mapped so far
            snp_match = re.search(r"SNPs mapped so far: (\d+)", line)
            if snp_match:
                mapped_snps = snp_match.group(
                    1
                )  # Capture the number of SNPs mapped so far

            # Capture SNP percentage mapped to at least one gene
            snp_percentage_match = re.search(
                r"(\d+)\s?\((-?\d+\.\d+%)\)\s?mapped to at least one gene", line
            )
            if snp_percentage_match:
                mapped_snps = snp_percentage_match.group(
                    1
                )  # Capture the number of SNPs
                snp_percentage = snp_percentage_match.group(2)  # Capture the percentage

            # Capture the number of genes mapped
            gene_match = re.search(
                r"at least one SNP mapped to each of a total of (\d+)\s+genes \(out of (\d+)\)",
                line,
            )
            if gene_match:
                mapped_genes = gene_match.group(1)  # Number of genes with mapped SNPs
                total_genes = gene_match.group(2)  # Total number of genes

        # Check if both SNPs and genes were mapped, then create a DataFrame
        if mapped_snps and mapped_genes:
            # Format the numbers with commas as thousand separators
            mapped_snps = f"{int(mapped_snps):,}"
            mapped_genes = f"{int(mapped_genes):,}"
            total_genes = f"{int(total_genes):,}"

            snp_combined = (
                f"{mapped_snps} ({snp_percentage})" if snp_percentage else mapped_snps
            )
            mapped_genes_combined = f"{mapped_genes} ({total_genes})"

            window_data = pd.DataFrame(
                {
                    "Window size": [window_size],
                    "Mapped SNPs": [snp_combined],
                    "Mapped genes (total genes)": [mapped_genes_combined],
                }
            )
            annotation_results = pd.concat(
                [annotation_results, window_data], ignore_index=True
            )

    for file in ["window_test.genes.annot.txt", "window_test.log"]:
        path = Path(file)
        if path.exists() and path.is_file():
            path.unlink()

    return annotation_results


def run_gene_annotation(
    window: int, snploc_file: str, geneloc_file: str, out_file: str
):
    """
    Runs the MAGMA annotation command with the specified parameters.

    :param window: Window size for annotation.
    :param snploc_file: Path to SNP location file.
    :param geneloc_file: Path to gene location file.
    :param out_file: Output file name.
    :return: Command output as a string.
    """
    command = f"magma --annotate window={window} --snp-loc {snploc_file} --gene-loc {geneloc_file} --out {out_file}"
    result = sp.run(command, shell=True, text=True, capture_output=True)
    return result.stdout


def run_gene_analysis(
    reference: str,
    gene_anno_file: str,
    pval_file: str,
    out_file: str,
    trait: str,
    ncol: str = None,
    N: int = None,
    batch: list = None,
) -> None:
    """
    Runs the MAGMA analysis command with the specified parameters and prints the output line by line.

    :param reference: Path to reference file.
    :param gene_anno_file: Path to gene annotation file.
    :param pval_file: Path to p-value file.
    :param out_file: Output file name.
    :param trait: Trait name to prefix output lines.
    """

    # Check inputs
    if ncol == None and N == None:
        raise ValueError("Either 'ncol' or 'N' must be specified.")
    else:
        samples = "ncol=" + ncol if ncol else "N=" + str(N)
    if batch:
        if not isinstance(batch, list):
            raise ValueError("batch must be a list")
        if len(batch) != 2:
            raise ValueError("batch must contain [batch_nr, n_batches]")
        else:
            batch = f"--batch {int(batch[0])} {int(batch[1])}"
    else:
        batch = ""

    # Run gene analysis
    print(f"Processing trait: {trait}...")
    command = f"magma --bfile {reference} --gene-annot {gene_anno_file} --pval {pval_file} use=1,2 {samples} {batch} --out {out_file}"
    process = sp.Popen(command, shell=True, stdout=sp.PIPE, stderr=sp.PIPE, text=True)

    for line in iter(process.stdout.readline, ""):
        print(f"{trait}  {line}", end="")

    process.stdout.close()
    process.wait()

#TODO: Possibly change to chr batch method
#TODO: Add reference file check
def batch_gene_analysis(trait: str, 
                        samples: str,
                        n_batches: int,
                        gene_anno_file: Path,
                        pval_file: Path,
                        out_dir: Path,
                        reference: Path,
                        time: str = "2:00:00",
                        mem: str = "1500MB"
                       ):
    """
    Submits MAGMA analysis jobs to SLURM in a loop without writing script files to disk.

    Parameters:
        trait (str): The trait name used for input and output file naming.
        samples (str): Sample size argument for MAGMA. must be either "N=" for total samples, 
                       or "ncol=" for sample size column name (e.g., "N=500000", "ncol=sample_size").
        n_batches (int): Number of batches to process.
        data_dir (Path): Base directory where all data is stored.
        out_dir (Path): Output directory
        reference (Path): Path to the reference file for MAGMA.
        time (str): Time limit for the SLURM job (default: "2:00:00").
        mem (str): Memory allocation for the SLURM job (default: "1500MB").
    """
    

    if not gene_anno_file.exists():
        raise FileNotFoundError(f"Gene annotation file does not exist: {gene_anno_file}")

    if not pval_file.exists():
        raise FileNotFoundError(f"P-value file does not exist: {pval_file}")

    if not samples.startswith(("ncol=", "N=")):
        raise ValueError('Samples must be in the format "ncol={column_name}" or "N={total_samples}".') 
    
    # Add dir for trait to ouput directory 
    out_dir = out_dir / trait
    # Dir for job logs
    log_dir = out_dir / "job_logs"
    # Add trait as prefix for outfile
    out_file = out_dir / trait
    
    # Ensure necessary directories exist
    out_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    # Loop over the batch numbers and submit jobs
    for batch_num in range(1, n_batches + 1):
        batch_script = f"""#!/bin/bash
#SBATCH --job-name=magma_analysis
#SBATCH --output={log_dir}/job_output_{batch_num}.log
#SBATCH --error={log_dir}/job_error_{batch_num}.log
#SBATCH --time={time}
#SBATCH --mem={mem}

magma --bfile {reference} --gene-annot {gene_anno_file} --pval {pval_file} use=1,2 {samples} --batch {batch_num} {n_batches} --out {out_file}
"""
        # Submit job
        process = sp.run(["sbatch"], input=batch_script, text=True, capture_output=True)

        # Print job submission output
        print(f"Submitted batch {batch_num}: {process.stdout.strip()}")

