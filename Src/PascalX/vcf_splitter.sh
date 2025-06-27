#!/bin/bash

# Check if a VCF file path argument is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <path_to_vcf_file>"
    exit 1
fi

# Get the VCF file path from the argument
VCF_FILE=$1

# Check if the provided VCF file exists
if [ ! -f "$VCF_FILE" ]; then
    echo "Error: File '$VCF_FILE' not found!"
    exit 1
fi

# Load bcftools module
module load BCFtools/1.19-GCCcore-11.3.0

# Split the VCF file into individual chromosome files (1-22)
for chr in {1..22}
do
    # Define the output file name for each chromosome
    OUT_FILE="${VCF_FILE%.vcf.gz}.chr$chr.vcf.gz"

    # Check if the chromosome file already exists
    if [ ! -f "$OUT_FILE" ]; then
        echo "Creating chromosome $chr file: $OUT_FILE"
        # If the file doesn't exist, create it
        bcftools view -r $chr "$VCF_FILE" -Oz -o "$OUT_FILE"
    else
        echo "Chromosome $chr file already exists: $OUT_FILE"
    fi
done

echo "Chromosome splitting complete."
