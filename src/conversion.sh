#!/bin/bash


mapper() {
   bwa mem \
   -R $1 \
   -p $2 \
   $3 > $4
}


fasta_index() {
    samtools faidx $1
}


sam_bam() {
    gatk SortSam \
    -I $1 \
    -O $2 \
    -SO coordinate
}


index_bam() {
    samtools index $1 $2.bai
}


seq_dict() {
    gatk CreateSequenceDictionary \
    -R $1\
    -O $2/$3.dict
}


haplotype() {
    gatk HaplotypeCaller \
    -R $1 \
    -I $2 \
    -O $3
}


copier() {
    cp $1 $2
}


"$@"
