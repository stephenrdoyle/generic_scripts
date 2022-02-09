#!/bin/bash

# GFF/FASTA chunker

# to run: ./gff_fasta_chunker.sh <FASTA> <GFF> <WINDOW_SIZE>
# eg. ./gff_fasta_chunker.sh REF.fa GFF.gff 1000000


# input files
FASTA=$1
GFF=$2
WINDOW=$3

# index reference and get chromosome sizes
samtools faidx ${FASTA}
cut -f1,2 ${FASTA}.fai > REF.genome

# make bed file of windows
bedtools makewindows -g REF.genome -w ${WINDOW} | awk '{print $1,$2+1,$3}' OFS="\t" > REF.${WINDOW}.bed


# chunk fasta into individual files
while read CHR START END; do
     samtools faidx ${FASTA} ${CHR}:${START}-${END} > ${CHR}_${START}_${END}.fa ;
     done < REF.${WINDOW}.bed


# chunk gff into individual files
while read CHR START END; do
     echo "##gff-version 3" > ${CHR}_${START}_${END}.gff
     grep "^${CHR}" ${GFF} | awk -v start=${START} -v end=${END} '{if($4>=start && $4<end) print $1,$2,$3,$4-start+1,$5-start+1,$6,$7,$8,$9}' OFS="\t" >> ${CHR}_${START}_${END}.gff;
     done < REF.${WINDOW}.bed
