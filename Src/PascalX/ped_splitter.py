"""
This script extracts sample IDs from a PED-like file where population ancestry is annotated.
Specifically, it selects individuals labeled as "EUR" (European ancestry) in the 7th column
and writes their Family ID and Individual ID to a new file in tab-separated format.

Usage:
------
    python3 ped_splitter.py <input_file> <output_file>

Arguments:
----------
    <input_file>   : Path to the input PED or metadata file (whitespace-delimited).
    <output_file>  : Path where the filtered sample IDs will be saved.


Author : Hessel Dijkstra
"""

import sys

def extract_samples(file_path, output_file):
    with open(file_path, 'r') as file, open(output_file, 'w') as out_file:
        for line in file:
            columns = line.strip().split()
            # Ensure there are at least 7 columns and the 7th column is "EUR"
            if len(columns) >= 7 and columns[6] == "EUR":
                # Write both Family ID (columns[0]) and Individual ID (columns[1])
                out_file.write(columns[0] + "\t" + columns[1] + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 ped_splitter.py <input_file> <output_file")
        sys.exit(1)
    extract_samples(sys.argv[1], sys.argv[2])
