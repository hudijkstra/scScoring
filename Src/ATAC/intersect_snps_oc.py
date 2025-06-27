"""
intersect_snps_with_peaks.py

This script intersects a set of SNPs with open chromatin regions (e.g., ATAC-seq peaks),
optionally filtering peaks by signal intensity before the intersection. The output is a BED 
file containing only the SNPs that fall within filtered peak regions.

Functionality:
- Loads an ATAC peak BED file and filters peaks based on a signal threshold (column 8).
- Loads a SNP BED file.
- Intersects SNPs with the filtered open chromatin peaks using pybedtools.
- Saves the intersected SNPs to an output BED file.

Usage:
    python intersect_snps_with_peaks.py --snps snps.bed --peaks atac_peaks.bed --out snps_in_open_chromatin.bed

Arguments:
    --snps   : Path to a BED file containing SNP positions.
    --peaks  : Path to a BED file containing ATAC-seq peaks, where column 8 (0-based index 7) contains signal values.
    --out    : Path to the output BED file where intersected SNPs will be written.

Notes:
- The script strips the 'chr' prefix and converts chromosome 'X' and 'Y' to '23' and '24' respectively.
- Output BED file has no header and contains only the SNPs that fall within open chromatin regions.


Author: Hessel Dijkstra
"""

import argparse
import pandas as pd
import pybedtools
import os


def filter_peaks(peaks_df, threshold=0.001):
    "Filters region where atleast 1/1000 cells 
    # Filter peaks with signal > threshold in column 8 (V8)
    print(f"Total peaks before filter: {len(peaks_df)}")
    filtered_peaks = peaks_df[peaks_df[7] > threshold].copy()
    print(f"Peaks after filter: {len(filtered_peaks)}")
    
    return filtered_peaks

def get_open_regions(peaks_df):
    # CHR, start, end
    region_data = peaks_df[[0, 1, 2]].copy()
    
    # make sure data types are correct
    region_data[0] = region_data[0].astype(int)
    region_data[1] = region_data[1].astype(int)
    region_data[2] = region_data[2].astype(int)

    # convert to BedTool object
    open_regions = pybedtools.BedTool.from_dataframe(region_data)

    return open_regions


def main(snp_bed_path, peak_bed_path, output_file):
    # Read and filter open chromatin peaks
    peaks_df = pd.read_csv(peak_bed_path, sep="\t", compression="infer", comment="#", header=None)
    
    # Make sure format is correct 
    peaks_df[0] = peaks_df[0].str.replace("^chr", "", regex=True)
    peaks_df[0] = peaks_df[0].replace({"X": "23", "Y": "24"})
   
    # Filter peaks based on signal value
    filtered_peaks = filter_peaks(peaks_df)

    # Get open chromatine regions
    open_regions = get_open_regions(filtered_peaks)

    # Load SNPs as BedTool
    snps = pybedtools.BedTool(snp_bed_path)
    initial_snp_count = len(snps)

    # Intersect SNPs with open chromatin
    intersected = snps.intersect(open_regions, wa=True)
    mapped_snp_count = len(intersected)

    print(f"Initial SNPs: {initial_snp_count}, SNPs in open chromatin: {mapped_snp_count}")

    # Write results
    df_mapped = intersected.to_dataframe()
    df_mapped.to_csv(output_file, sep="\t", index=False, header=False)

    print(f"Saved intersected SNPs to: {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Intersect SNPs with open chromatin regions.")
    parser.add_argument("--snps", required=True, help="Path to SNP BED file")
    parser.add_argument("--peaks", required=True, help="Path to ATAC peak BED file")
    parser.add_argument("--out", required=True, help="Output file")

    args = parser.parse_args()
    main(args.snps, args.peaks, args.out)
