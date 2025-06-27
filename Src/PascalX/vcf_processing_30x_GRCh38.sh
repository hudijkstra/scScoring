#!/bin/bash
# This script automates preprocessing of chromosome-wise VCF data into PLINK transposed format (.tped/.tfam).
# It filters individuals using a custom sample list (`custom_keep.tsv`) and supports both EUR-only and full-population output.
# It also handles optional download of PLINK if needed

# Make DIR and copy over
mkdir -p $1/
cp $5 $1/custom_keep.tsv
cd $1/

# Download plink
if [ $4 = "tped" ];
then
    wget -nc http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20201019.zip
    unzip -o plink_linux_x86_64_20201019.zip
fi

function pwait() {
    while [ $(jobs -p | wc -l) -ge $3 ]; do
        sleep 1
    done
}

if [ $4 = "tped" ];
then
    if [ $2 = "EUR" ];
    then
        for i in {1..22}
        do
            ./plink --recode 12 transpose --vcf-half-call missing --vcf 30x-GRCh38.chr$i.vcf.gz -keep custom_keep.tsv --out EUR.30x-GRCh38.chr$i &
            pwait $3
        done
    else
        for i in {1..22}
        do

            ./plink --recode 12 transpose --vcf-half-call missing --vcf 30x-GRCh38.chr$i.vcf.gz -keep custom_keep.tsv --out ALL.30x-GRCh38.chr$i &
            pwait $3
        done
    fi

    wait
    gzip *.tped
else
    for i in {1..22}
    do
        mv 30x-GRCh38.chr$i.vcf.gz $2.30x-1KG.GRCh38.chr$i.vcf.gz
    done
fi
