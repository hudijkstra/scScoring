"""
This module provides functionality for submitting batch-based gene-level association analyses 
using the MAGMA software through SLURM job scheduling. The core function, `batch_gene_analysis`, 
dynamically constructs SLURM batch scripts in memory and submits them without writing any 
temporary script files to disk.

Example:
--------
batch_gene_analysis(
    trait="schizophrenia",
    samples="N=500000",
    n_batches=10,
    gene_anno_file=Path("/path/to/gene_annotation.genes.annot"),
    pval_file=Path("/path/to/gwas_results.tsv"),
    out_dir=Path("/output/magma_results"),
    reference=Path("/path/to/reference_data")
)

Author: Hessel Dijkstra

"""

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