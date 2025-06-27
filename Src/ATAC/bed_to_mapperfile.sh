#!/bin/bash

# Script: run_extract_ids.sh
# Usage: ./run_extract_ids.sh <BED_FILE> <VARIANT_COL> <GENE_COL> <OUTPUT_FILE>

source /groups/umcg-franke-scrna/tmp02/users/umcg-hdijkstra/scscore/bin/activate

BED_FILE="/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/sc-scoring-student-project/scScoring/Data/ATAC/mapped_snps/oc_eqtl_intersect_mono.tsv"
VARIANT_COL="3"
GENE_COL="4"
OUTPUT_FILE="/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/sc-scoring-student-project/scScoring/Data/PascalX/mapperfiles/gwas_eqtl_oc_intersect.tsv"

# Run the Python script with the provided arguments
python3 bed_to_mapperfile.py "$BED_FILE" "$VARIANT_COL" "$GENE_COL" "$OUTPUT_FILE"
