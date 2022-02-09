#!/usr/bin/env bash

# ========================================================================================
# run_spliceleader_finder.sh
#
# Workflow to find known splice leader sequences in RNAseq data
#
# Reference:
# Software:
#
# Usage: ~sd21/bash_scripts/run_spliceleader_finder.sh <PREFIX> <REFERENCE> <GFF> <splice leader sequence> <SL_MIN_LENGTH> <R1.fastq> <R2.fastq>
#
# @authors
# Stephen Doyle <sd21@sanger.ac.uk>
#
# ----------------------------------------------------------------------------------------

# Requirements in path
#--- samtools-1.6
#--- bedtools v2.17.0
#--- hisat2


# Notes
#--- use bsub to run
#--- relies on a Augustus-like GFF to detect overlap of SL with transcription start sites, and may fail to to the last step if not compatible. Can hack if needed.


# Modifications
#- 180712: added a minimum length of SL parameter for cutadapt to work with - originally set to 10 bp, but probably should be set to ~14 based on Alans script to find SL ßrelated seq by chance in genome
#- 180712: changed length of upstream window from start codon from 50 to 100 bp when looking for overlap between SL sequence and transcripts


########################################################################

PREFIX=$1
REFERENCE=$2
GFF=$3
SL_SEQ=$4
SL_MIN_LENGTH=$5
R1=$6
R2=$7


if [ "$#" -eq 0 ] || [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
	echo ""
    echo "Usage: ~sd21/bash_scripts/run_spliceleader_finder.sh <PREFIX> <REFERENCE> <GFF> <splice leader sequence> <SL_MIN_LENGTH> <R1.fastq> <R2.fastq>"
    echo ""
    exit 0
fi

if [ -d "${PREFIX}_SL_ANALYSIS_out" ]; then
	echo -e "\nThere is already a directory started with this sample name. Rename and start again\n"
    exit 0
fi


mkdir ${PREFIX}_SL_ANALYSIS_out
cd ${PREFIX}_SL_ANALYSIS_out
ln -s ${R1}
ln -s ${R2}
ln -s ${REFERENCE}
ln -s ${GFF}


# cutadapt - http://cutadapt.readthedocs.io/en/stable/

# R1 reads with mates
cutadapt --trimmed-only --prefix=${PREFIX}_SL_trimmed_ --overlap=${SL_MIN_LENGTH} --error-rate=0.1 --info-file=info1.tmp -g ${SL_SEQ} -o ${PREFIX}.trimmed.1_1.tmp.fastq
-p ${PREFIX}.trimmed.1_2.tmp.fastq ${R1} ${R2} > ${PREFIX}_cutadapt_trim.summary


# get the reverse reads by switching R1 and R2 around
cutadapt --trimmed-only --prefix=${PREFIX}_SL_trimmed_ --overlap=${SL_MIN_LENGTH} --error-rate=0.1 --info-file=info2.tmp -g ${SL_SEQ} -o ${PREFIX}.trimmed.2_2.tmp.fastq
-p ${PREFIX}.trimmed.2_1.tmp.fastq ${R2} ${R1} >> ${PREFIX}_cutadapt_trim.summary


# extract adapter hit sequences, sort and count, make a sequence logo
#--- extract adapters
cat info1.tmp info2.tmp | awk '{if($2!="-1") print $0}' | cut -f 6 | sort | uniq -c | sort -k1,1n > ${PREFIX}_adapter_hits.counts
#--- make a fasta
while read count sequence; do echo -e ">adapter_${count}\n$sequence" >> ${PREFIX}_adapters.fa ; done < ${PREFIX}_adapter_hits.counts
#--- extend variable sequence to a fixed length (n=25bp) using "n"s to pad the sequence, and generate a sequence logo
awk '$1~">"{print $0}$1!~">"{tmp="";for(i=1;i<25-length($0)+1;i++){tmp=tmp"n"};print tmp""$0}' ${PREFIX}_adapters.fa >> ${PREFIX}_adapters.fixedlength.fa
~sd21/lustre118_link/software/bin/weblogo --format pdf --datatype fasta --ignore-lower-case < ${PREFIX}_adapters.fixedlength.fa > ${PREFIX}_adapters.fixedlength.pdf


# merge the two datasets together
cat ${PREFIX}.trimmed.1_1.tmp.fastq ${PREFIX}.trimmed.2_1.tmp.fastq > ${PREFIX}.trimmed.merged_R1.tmp.fastq
cat ${PREFIX}.trimmed.1_2.tmp.fastq ${PREFIX}.trimmed.2_2.tmp.fastq > ${PREFIX}.trimmed.merged_R2.tmp.fastq


# map reads using histat
/nfs/users/nfs_a/ar11/hisat2-2.0.0-beta/hisat2-build -p 1 ${REFERENCE} ${PREFIX}.REFidx.tmp
/nfs/users/nfs_a/ar11/hisat2-2.0.0-beta/hisat2 ${PREFIX}.REFidx.tmp -1 ${PREFIX}.trimmed.merged_R1.tmp.fastq -2 ${PREFIX}.trimmed.merged_R2.tmp.fastq -q -S ${PREFIX}.tmp
.sam


# extract trimmed and mapped SL reads
samtools-1.6 view ${PREFIX}.tmp.sam -H > ${PREFIX}.tmp.sam.header
grep "${PREFIX}_SL_trimmed_" ${PREFIX}.tmp.sam > ${PREFIX}.tmp.sam.body
cat ${PREFIX}.tmp.sam.header ${PREFIX}.tmp.sam.body > ${PREFIX}.SL-only.tmp.sam
samtools-1.6 view -bS ${PREFIX}.SL-only.tmp.sam > ${PREFIX}.SL-only.tmp.bam
samtools-1.6 sort ${PREFIX}.SL-only.tmp.bam -o ${PREFIX}.SL-only.sorted.bam
samtools index -b ${PREFIX}.SL-only.sorted.bam


# extract coordinates, correcting for read direction, and output in bed format
bedtools bamtobed -i ${PREFIX}.SL-only.sorted.bam | grep "${PREFIX}_SL_trimmed_" | awk '{if($6=="+" && $5> 0) print $1,$2-1,$2+1,$4,$5,$6; else if($6=="-" && $5> 0) prin
t $1,$3-1,$3+1,$4,$5,$6}' OFS="\t" > ${PREFIX}.SL-only.tmp.bed


# extract dinuclotide at splice site
bedtools getfasta -s -fi ${REFERENCE} -bed ${PREFIX}.SL-only.tmp.bed -fo ${PREFIX}.dinucleotide.tmp.fasta
cat ${PREFIX}.dinucleotide.tmp.fasta | grep -v ">" > ${PREFIX}.dinucleotide.tmp.seq


# extract twenty-mer at splice site
bedtools bamtobed -i ${PREFIX}.SL-only.sorted.bam | grep "${PREFIX}_SL_trimmed_" | awk '{if($6=="+" && $5> 0) print $1,$2-10,$2+10,$4,$5,$6; else if($6=="-" && $5> 0) pr
int $1,$3-10,$3+10,$4,$5,$6}' OFS="\t" > ${PREFIX}.twentymer.tmp.bed
bedtools getfasta -s -fi  ${REFERENCE} -bed ${PREFIX}.twentymer.tmp.bed -fo ${PREFIX}.twentymer.fasta
cat ${PREFIX}.twentymer.fasta | grep -v ">" >  ${PREFIX}.twentymer.tmp.seq
~sd21/lustre118_link/software/bin/weblogo --format pdf < ${PREFIX}.twentymer.fasta > ${PREFIX}.twentymer.weblogo.pdf


# combine data
paste ${PREFIX}.SL-only.tmp.bed ${PREFIX}.dinucleotide.tmp.seq ${PREFIX}.twentymer.tmp.seq > ${PREFIX}.SL-only.bed


# calculate overlap with transcripts

awk '{if($3=="mRNA" || $3=="transcript") print $9}' ${GFF} | sed -e 's/ID=//g' -e 's/Parent=//g' -e 's/Name=//g' -e 's/;/\t/g' > ${PREFIX}.gene_transcript.list

#--- extract transcripts from GFF, and generate a flaking window around the start site 50 bp upstream and 20 downstream

while read -r transcript gene name other; do \
	grep "Parent=${transcript}" ${GFF} | \
	grep "CDS" | sort -k1,1 -k4,4n | \
	awk -v transcript_id="${transcript}" -v gene_id="${gene}" -v name_id="${name}" ' $7=="+" {print $1,$4-200,$4+30,"ID="transcript_id";Parent="gene_id";Name="name_i
d,".",$7; p=1 ; exit } END{if (p != 1) {print $1,$5-30,$5+200,"ID="transcript_id";Parent="gene_id";Name="name_id,".",$7}}' OFS="\t" >> ${PREFIX}_transcript_start_windows
.tmp.bed ; \
	done < ${PREFIX}.gene_transcript.list

# clean up a little
#--- make sure all bed coords have stranded information - those that dont are from features without CDSs and therefore not needed
#--- fix and start coordinates in bed file that are less than 0. Replace with start coordinate of 1.

awk '{if($6=="+" || $6=="-") print $0}' ${PREFIX}_transcript_start_windows.tmp.bed | sed 's/-[0-9][0-9]*/1/' > ${PREFIX}_transcript_start_windows.bed

#--- calculate coverage of SL sequences in transcript start site windows
bedtools coverage -s -b ${PREFIX}_transcript_start_windows.bed -a ${PREFIX}.SL-only.bed | sort -k1,1 -k2,2n > ${PREFIX}_transcript_start_windows.SL.coverage








#--- extract internal CDS sequences from GFF, and generate a flaking window around the start site 10 bp upstream and 10 downstream

while read -r transcript gene name other; do \
	grep "Parent=${transcript}$" ${GFF} | \
	grep "CDS" | sort -k1,1 -k4,4n | \
	sed '1d;$d' | \
	awk -v transcript_id="${transcript}" -v gene_id="${gene}" -v name_id="${name}" '{print $1,$4-10,$5+10,"ID="transcript_id";Parent="gene_id";Name="name_id,".",$7}'
 OFS="\t" >> ${PREFIX}_internal_cds_windows.tmp.bed ; \
	done < ${PREFIX}.gene_transcript.list

awk '{if($6=="+" || $6=="-") print $0}' ${PREFIX}_internal_cds_windows.tmp.bed | sed 's/-[0-9][0-9]*/1/' > ${PREFIX}_internal_cds_windows.bed

bedtools coverage -s -b ${PREFIX}_internal_cds_windows.bed -a ${PREFIX}.SL-only.bed | sort -k1,1 -k2,2n > ${PREFIX}_internal_cds_windows.SL.coverage




#--- calculate coverage


awk '{if($7>1) print $4}' ${PREFIX}_transcript_start_windows.SL.coverage | sed -e 's/;/\t/g' -e 's/ID=//g' | cut -f1 > ${PREFIX}_transcript_start_windows.SL.genelist

awk '{if($7>1) print $4}' ${PREFIX}_internal_cds_windows.SL.coverage | sed -e 's/;/\t/g' -e 's/ID=//g' | cut -f1 | sort | uniq -c > ${PREFIX}_internal_cds_windows.SL.gen
elist



# final cleanup
#rm *.tmp.*














# while read -r transcript gene name other; do \
# 	grep "Parent=${transcript}$" ${GFF} | grep "CDS" | \
# 	# if 1st exon < 20 bp on the positive strand, skip and calculate of the second exon
# 	if [[ $(head -n1 | awk '{print $5-$4}') -le "20" ]] && [[ $(head -n1 | awk '{print $7}') == "+" ]]
# 	then
# 		sed -n '2p' | awk -v transcript_id="${transcript}" -v gene_id="${gene}" -v name_id="${name}" '{print $1,$4-100,$4+20,"ID="transcript_id";Parent="gene_id"
;Name="name_id,".",$7}' OFS="\t" >> ${PREFIX}_transcript_start_windows.tmp.bed
# 	# if 1st exon > 20 bp on the positive strand, ok
# 	elif [[ $(head -n1 | awk '{print $5-$4}') -gt "20" ]] && [[ $(head -n1 | awk '{print $7}') == "+" ]]
# 	then
# 			sed -n '1p' | awk -v transcript_id="${transcript}" -v gene_id="${gene}" -v name_id="${name}" '{print $1,$4-100,$4+20,"ID="transcript_id";Parent="
gene_id";Name="name_id,".",$7}' OFS="\t" >> ${PREFIX}_transcript_start_windows.tmp.bed
#
# 	elif [[ $(tail -n1 | awk '{print $5-$4}') -le "20" ]] && [[ $(tail -n1 | awk '{print $7}') == "-" ]]
# 	then
# 			sed 'x;$!d' | awk -v transcript_id="${transcript}" -v gene_id="${gene}" -v name_id="${name}" '{print $1,$5-20,$5+100,"ID="transcript_id";Parent="
gene_id";Name="name_id,".",$7}' OFS="\t" >> ${PREFIX}_transcript_start_windows.tmp.bed
#
# 	else [[ $(tail -n1 | awk '{print $5-$4}') -gt "20" ]] && [[ $(tail -n1 | awk '{print $7}') == "-" ]]
# 	 		sed '$!d' | awk -v transcript_id="${transcript}" -v gene_id="${gene}" -v name_id="${name}" '{print $1,$5-20,$5+100,"ID="transcript_id";Parent="ge
ne_id";Name="name_id,".",$7}' OFS="\t" >> ${PREFIX}_transcript_start_windows.tmp.bed
# 	fi; done < <(head -n10 ${PREFIX}.gene_transcript.list)
