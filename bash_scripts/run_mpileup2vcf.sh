#!/usr/bin/env bash

# ========================================================================================
# run_mpileup2vcf.sh
#
# Tool will run an mpileup on all bams in a specified directory, and split the job up by sequence in the reference to parallelise the job. Finally, it shoudl merge all the split jobs into a single file.
# mpileup command based on https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1024-y
#
# Reference:
# Software:
#
# Usage: ~sd21/bash_scripts/run_mpileup2vcf.sh <PREFIX> <REFERENCE_FASTA> <BAMLIST>
#
# @authors
# Stephen Doyle <sd21@sanger.ac.uk>
#
# ----------------------------------------------------------------------------------------


# load modules
module load \
samtools/1.6--h244ad75_4 \
fastaq/3.17.0-docker3 \
bcftools/1.9--h68d8f2e_9



(($# == 3)) || { echo -e "\nUsage: $0 <prefix> <reference sequence> <bamlist_file>\n\n"; exit; }

prefix=$1
ref=$2
bamlist=$3

#--- Step 1: prepare input files
cp $ref ref.tmp
samtools faidx ref.tmp
fastaq to_fasta -l 0 ref.tmp ref.tmp2
samtools faidx ref.tmp2
grep ">" ref.tmp | cut -f1  -d" " | sed -e 's/>//g' | cat -n > ref_seq.tmp.list

while read number sequences; do grep -A1 "$sequences" ref.tmp2 > $sequences.tmp.fasta; done < ref_seq.tmp.list


#while read name; do
#	if [ ! -f ${name}.bai ]; then
#	 samtools index -b ${name}
#	fi; done < ${bamlist}

while read number sequences; do
echo -e "samtools mpileup --ignore-RG -Ou -t DP,SP,AD,ADF,INFO/AD -b $bamlist -r $sequences -f ref.tmp -F0.25 -d500 -E | bcftools call -vm -Oz -o $number.$sequences.tmp.vcf.gz" > run_mpileup2vcf.tmp.$number;
done < ref_seq.tmp.list
chmod a+x run_mpileup2vcf.tmp*


jobs=$( wc -l ref_seq.tmp.list | cut -f1 -d" " )
bsub -q long -R'span[hosts=1] select[mem>10000] rusage[mem=10000]' -M10000 -J mpileup2vcf_array[1-$jobs] -e mpileup2vcf_array[1-$jobs].e -o mpileup2vcf_array[1-$jobs].o ./run_mpileup2vcf.tmp.\$LSB_JOBINDEX



#--- Step 2: bring it together


echo -e "ls -1v *.tmp.vcf.gz > vcffiles.fofn && vcf-concat -f vcffiles.fofn | gzip -c >  ${prefix}.raw.vcf.gz" > run_mp2vcf_combine
chmod a+x run_mp2vcf_combine

bsub -q normal -w "done(mpileup2vcf_array)"  -R'span[hosts=1] select[mem>1000] rusage[mem=1000]' -n1 -M1000 -J mp2vcf_combine -e mp2vcf_combine.e -o mp2vcf_combine.o ./run_mp2vcf_combine ${prefix}
