#!/bin/bash

# GFF/FASTA chunker

FASTA=$1
GFF=$2
WINDOW=$3



# make bed file of windows

samtools faidx ${FASTA}

cut -f1,2 REF.fa.fai > REF.genome

bedtools makewindows -g REF.genome -w ${WINDOW} > REF.${WINDOW}.bed



# chunk fasta into individual files

while read CHR START END; do
     samtools faidx ${FASTA} ${CHR}:${START}-${END} > ${CHR}_${START}_${END}.fa ;
     done < REF.${WINDOW}.bed


# chunk gff into individual files

while read CHR START END; do
     echo "##gff-version 3" > ${CHR}_${START}_${END}.gff
     grep "^${CHR}" ${GFF} | awk -v start=${START} -v end=${END} '{if($4>=start && $4<end) print $1,$2,$3,$4-start,$5-start,$6,$7,$8,$9}' OFS="\t" >> ${CHR}_${START}_${END}.gff;
     done < REF.${WINDOW}.bed
