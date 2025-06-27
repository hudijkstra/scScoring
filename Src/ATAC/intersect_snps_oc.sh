#!/bin/bash

ml BEDTools/2.30.0-GCCcore-11.3.0
source /groups/umcg-franke-scrna/tmp02/users/umcg-hdijkstra/scscore/bin/activate


SNP_FILE="/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/sc-scoring-student-project/scScoring/Data/ATAC/unmapped_snps/gwas_eqtl_intersect.bed"
PEAKS_FILE="/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/sc-scoring-student-project/scScoring/Data/ATAC/C_peaks/mo_peaks_lane1to80_wrna_UT_monocyte.bed.gz"
OUT_PREF="/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/sc-scoring-student-project/scScoring/Data/ATAC/mapped_snps/oc_eqtl_intersect_mono.tsv"

echo "Intersecting SNPs with peaks..."

python3 intersect_snps_oc.py \
    --snps $SNP_FILE \
    --peaks $PEAKS_FILE \
    --out $OUT_PREF