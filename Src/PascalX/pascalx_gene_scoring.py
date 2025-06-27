"""
This script provides a command-line interface for scoring genes using the PascalX framework.
It supports flexible configuration of the gene scoring procedure with various inputs:
reference LD panels, GWAS summary statistics, optional SNP-to-gene mappings, and chromosome filtering.

For more information, please visit the official PascalX documentation: https://bergmannlab.github.io/PascalX/

### IMPORTANT ###
The gene location file format is hardcoded. Modify the column indices in `Scorer.load_genome` in the main function if using a different format.


Key Features:

- Supports gene scoring using the `chi2sum` method (currently the only available method).
- Handles reference LD panels and gene location references.
- Accepts custom GWAS summary statistics and optional SNP-gene mapper files.
- Allows scoring across all chromosomes or a user-defined subset.
- Supports parallel execution via the `--cores` argument.
- Saves results as TSV files with gene scores and optional error reports.

Command-Line Arguments:
-----------------------
- `--scorer`        : Scoring method (default: "chi2sum").
- `--window`        : Genomic window size around each gene (default: 50000).
- `--cores`         : Number of CPU cores to use for parallel scoring.
- `--ref_panel`     : Path to the LD reference panel.
- `--gene_ref`      : Gene location reference file.
- `--summ_stats`    : Path to the GWAS summary statistics file.
- `--rscol`         : Column index of rsID in summary statistics (0-based).
- `--pcol`          : Column index of p-values in summary statistics (0-based).
- `--keepfile`      : (Optional) File with individual IDs to include from the LD panel.
- `--mapperfile`    : (Optional) SNP-to-gene mapping file (TSV format: SNP ID, gene ID).
- `--out_dir`       : Directory where output files are saved.
- `--chr`           : Chromosomes to process ("all" or comma-separated list, e.g. "1,2,5,20").

Output:
-------
- `{trait}_gscores.tsv`: Gene scores with p-values and number of SNPs used.
- `{trait}_gscores.err`: Optional error file listing genes that failed scoring.

Example Usage:
--------------
```bash
python3 pascalx_gene_scoring.py \
  --ref_panel data/1000G_EUR.ld \
  --gene_ref data/gene_locations.tsv \
  --summ_stats data/height_gwas.tsv \
  --rscol 0 --pcol 2 \
  --out_dir results/height \
  --cores 4 --window 10000 \


Author: Hessel Dijkstra

"""

import argparse
import os
import csv
from PascalX import genescorer

#NOTE: only one scorer supported for now
def get_scorer(scorer: str = "chi2sum", window: int = 50000):
    genescorers = {
        "chi2sum": genescorer.chi2sum(window = window)
    }
    if scorer not in genescorers:
        raise ValueError(f"Invalid scorer '{scorer}'. Available options: {list(genescorers.keys())}")
    return genescorers[scorer]


def chromosome_list(value):
    """
    Argument parser helper for validating and parsing chromosome input.

    Parameters:
    -----------
    value : str
        A string specifying chromosomes to process.
        Acceptable formats:
        - "all": Score all chromosomes (chr1 to chr22)
        - Comma-separated list of integers, e.g., "1,2,5,20"

    Returns:
    --------
    str or list of int
        Returns "all" if specified, otherwise a list of chromosome numbers as integers.

    """
    try:
        if value == "all":
            return value
        else:
            chroms = [int(ch) for ch in value.split(",")]
            if not all(1 <= ch <= 22 for ch in chroms):
                raise ValueError
            return chroms
    except ValueError:
        raise argparse.ArgumentTypeError("Chromosomes must be a comma-separated list of numbers between 1 and 22.")


def score_genes(chromosomes: str, Scorer: genescorer, cores = int):
    """
    Performs gene scoring using the provided PascalX scorer object for selected chromosomes.

    Parameters:
    -----------
    chromosomes : str or list of int
        Either the string "all" (to score all chromosomes), or a list of chromosome numbers (e.g., [1, 2, 22]).
    
    Scorer : genescorer
        An instance of a PascalX gene scoring object, such as one returned by `genescorer.chi2sum(...)`.
        This object must have `.score_all()` and `.score_chr()` methods implemented.

    cores : int
        Number of CPU cores to use for parallel processing.

    Returns:
    --------
    tuple
        A tuple `R = (scored_genes, warnings, errors)` returned by PascalX, where:
        - `scored_genes`: list of [gene, pvalue, nsnps]
        - `warnings`: implementation-dependent
        - `errors`: list of [gene, reason] for genes that failed to score
    """
    if chromosomes == "all":
        return Scorer.score_all(parallel = cores, nobar=True)
    else:
        return Scorer.score_chr(chrs = chromosomes, parallel = cores, nobar=True)


def save_score_data(file_name, R, out_dir):
    """
    Saves PascalX gene scoring results to disk.

    Parameters:
    -----------
    file_name : str
        The name of the summary statistics file used to derive a base output name.
        Only the base (without extension) is used for naming the output files.
    
    R : tuple
        A tuple returned by PascalX scoring, containing:
        - R[0]: List of [gene, pvalue, nsnps] rows for successfully scored genes.
        - R[1]: (Unused, typically warnings or internal logs)
        - R[2]: List of [gene, reason] rows for genes that failed to score.
    
    out_dir : str
        The output directory where the result files will be written.
        The function will create this directory if it doesn't exist.

    Output:
    -------
    - {base_name}_gscores.tsv : Tab-separated file with columns [gene, pvalue, nsnps].
    - {base_name}_gscores.err : (Optional) Tab-separated file with columns [gene, reason],
                                only created if there are failed genes.
    """
    os.makedirs(out_dir, exist_ok=True)
    
    base_name = os.path.basename(file_name).split(".")[0]
    output_file = os.path.join(out_dir,base_name)
    gene_score_file = output_file + "_gscores.tsv"
    
    # Write gene scores
    with open(gene_score_file, 'w') as out:
        csv_out = csv.writer(out, delimiter='\t')
        # Write the header
        csv_out.writerow(['gene', 'pvalue', 'nsnps'])
        for row in R[0]:
            csv_out.writerow(row)

    print(f"Scores saved to {gene_score_file}\n")

    # Write error file if any genes failed to score
    if R[2]:
        error_file = output_file + "_gscores.err"
        with open(error_file, 'w') as out:
            csv_out = csv.writer(out, delimiter='\t')
            # Write the header
            csv_out.writerow(["gene", "reason"])
            for row in R[2]:
                csv_out.writerow(row)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Score genes using PascalX.")
    
    parser.add_argument("--scorer", type=str, default="chi2sum", help="Scoring method (default: chi2sum)")
    parser.add_argument("--window", type=int, default=50000, help="Gene window (default: 50000)")
    parser.add_argument("--cores", type=int, default=1, help="Number of cores to use (default: 1)")
    parser.add_argument("--ref_panel", type=str, required=True, help="Path to LD reference file")
    parser.add_argument("--gene_ref", type=str, required=True, help="Path to gene reference file")
    parser.add_argument("--summ_stats", type=str, required=True, help="Summary statistics file path")
    parser.add_argument("--rscol", type=int, required=True, help="Summary statistics variant column index (0-based)")
    parser.add_argument("--pcol", type=int, required=True, help="Summary statistics p-value column index(0-based)")
    parser.add_argument("--keepfile", type=str, default=None, help="Optional file with sample IDs to keep")
    parser.add_argument("--mapperfile", type=str, default=None, help="Optional SNP-gene mapping tsv file. First column should be SNP ID, second column ensembl gene id.")
    parser.add_argument("--out_dir", type=str, required=True, help="Output directory")
    parser.add_argument(
        "--chr", type=chromosome_list, default="all",
        help='Chromosomes to score, default is "all". '
            'Specify as a comma-separated list (e.g., "1,2,5,10").'
            )

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
        print(f"{arg.ljust(max_key_length)} : {value}")  # Align keys dynamically
    print("###########################\n")


def main():
    args = parse_arguments()
    print_args(args)
    
    # Get scoring method
    Scorer = get_scorer(args.scorer, args.window)
    
    # Load LD reference panel
    Scorer.load_refpanel(
        filename = args.ref_panel,
        parallel=args.cores,
        keepfile=args.keepfile
        )
    
    # Load gene location reference
    print("Loading gene locations...")
    Scorer.load_genome(
        file = args.gene_ref,
        ccol = 1,       # Column containing chromosome number
        cid = 4,        # Column containing gene id
        csymb = 0,      # Column containing gene symbol
        cstx = 2,       # Column containing transcription start
        cetx = 3,       # Column containing transcription end
        cs = 6,         # Column containing strand
        chrStart = 0,   # Number of leading characters to skip in ccol
        splitchr = "\t",
        header = True
        )
    
    # Process GWAS summary statistics
    print(f"Processing GWAS: {args.summ_stats}...")
    Scorer.load_GWAS(
        args.summ_stats, 
        rscol = args.rscol,   # rsID 
        pcol = args.pcol,     # p-values 
        header=True, 
        delimiter="\t"
    )
    
    # Load SNP-gene mapping file if provided
    if args.mapperfile:
        Scorer.load_mapping(args.mapperfile, rcol=0, gcol=1, delimiter='\t', header=True, joint=True)

    print("Scoring genes...")
    R = score_genes(chromosomes = args.chr, Scorer = Scorer, cores = args.cores)
    
    # Save scores to file
    save_score_data(args.summ_stats, R, args.out_dir)


if __name__ == "__main__":
    main()