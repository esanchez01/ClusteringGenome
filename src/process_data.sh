#!/bin/bash


create_list () {
   ls $1 | \
   grep .vcf.gz | \
   grep -v ".tbi" | \
   sed -e $2 > $3
}


header () {
    zgrep '^\#' "$(head -1 data/temp/input.list)" | \
    gzip > data/temp/full_chroms.vcf.gz;
}


concatenate() {
    zgrep -v "^\#" $(cat data/temp/input.list) | \
    gzip >> data/temp/full_chroms.vcf.gz
}


filter() {
    plink2 \
    --vcf data/temp/full_chroms.vcf.gz \
    --snps-only \
    --maf $1 \
    --geno $2 \
    --mind $3 \
    --recode \
    --allow-extra-chr \
    --make-bed \
    --out data/temp/chromosomes
}


pca() {
    plink2 \
    --bfile data/temp/chromosomes \
    --pca $1 \
    $2 $3 \
    --allow-extra-chr \
    --out $4
}


"$@"
