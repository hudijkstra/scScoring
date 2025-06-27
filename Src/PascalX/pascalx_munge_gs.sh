#!/bin/bash

# Usage: create_geneset.sh input_file.tsv output_file.tsv trait_name
# Description: Creates a gene set file for scDRS from a p-value file using scdrs munge-gs.
# Expects input file to have two columns: gene_id and p-value.

source /groups/umcg-franke-scrna/tmp04/users/umcg-hdijkstra/scscore/bin/activate

INPUT_FILE="$1"
OUTPUT_FILE="$2"
TRAIT="$3"

if [[ -z "$INPUT_FILE" || -z "$OUTPUT_FILE" || -z "$TRAIT" ]]; then
    echo "Usage: $0 input_file.tsv output_file.tsv trait_name"
    exit 1
fi

echo "Processing $INPUT_FILE..."

# Create a temporary file for the p-value subset
PTMP=$(mktemp)

# Clean up temp file on exit
cleanup() {
    rm -f "$PTMP"
}
trap cleanup EXIT

# Read and modify header
HEADER=$(head -n 1 "$INPUT_FILE")
COL1=$(echo "$HEADER" | cut -f1)

# Write new header to temp file
echo -e "${COL1}\t${TRAIT}" > "$PTMP"

# Add first two columns from data (excluding header)
tail -n +2 "$INPUT_FILE" | cut -f1,2 >> "$PTMP"

scdrs munge-gs \
    --out-file "$OUTPUT_FILE" \
    --pval_file "$PTMP" \
    --weight zscore \
    --n-max 1000

echo "Gene set file created: $OUTPUT_FILE"
