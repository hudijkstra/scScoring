"""
bed_to_mapperfile.py

This script extracts variant IDs (e.g., rsIDs) and gene IDs (e.g., Ensembl IDs) from a BED-like 
(tab-delimited) file and creates a mapper file suitable for use in PascalX.

The script expects column indices (0-based) for both the variant and gene ID fields, 
and outputs a TSV file with two columns: `rsID` and `ENSG`.

Example usage:
    python bed_to_mapperfile.py input.bed 3 7 output_mapper.tsv

Arguments:
    bed_file     : Path to the input BED file (tab-delimited with no header).
    variant_col  : 0-based index of the column containing variant IDs (e.g., rsIDs).
    gene_col     : 0-based index of the column containing gene IDs (e.g., Ensembl gene IDs).
    output_file  : Path to save the resulting TSV mapper file with header.

Output:
    A TSV file with columns:
        rsID    ENSG

This format can be used as the --mapper argument in PascalX to map variants to genes.

Author: Hessel Dijkstra
"""
import pandas as pd
import argparse
import os

def process_bed_file(bed_file, variant_col, gene_col, output_file):
    # Read the BED file
    df = pd.read_csv(bed_file, sep="\t", header=None)

    df_variant_gene = df[[variant_col, gene_col]]
    df_variant_gene.columns = ['rsID', 'ENSG']

    # Save to file with header
    df_variant_gene.to_csv(output_file, sep="\t", index=False, header=True)

    print(f"Processed data saved to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract variant and gene IDs from a BED file.")
    parser.add_argument("bed_file", type=str, help="Path to the BED file (tab-delimited).")
    parser.add_argument("variant_col", type=int, help="0-based column index for variant ID.")
    parser.add_argument("gene_col", type=int, help="0-based column index for gene ID.")
    parser.add_argument("output_file", type=str, help="Path to the output TSV file.")

    args = parser.parse_args()

    process_bed_file(args.bed_file, args.variant_col, args.gene_col, args.output_file)
