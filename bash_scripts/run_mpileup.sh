#!/usr/bin/env bash

# ========================================================================================
# run_mpileup.sh
#
# Tool will run an mpileup on all bams in a specified directory, and split the job up by sequence in the reference to parallelise the job. Finally, it shoudl merge all the split jobs into a single file.
# mpileup command based on https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1024-y
#
# Reference:
# Software:
#
# Usage: ~sd21/bash_scripts/run_mpileup2popoolation2.sh <PREFIX> <REFERENCE_FASTA> <BAMLIST>
#
# @authors
# Stephen Doyle <sd21@sanger.ac.uk>
#
# ----------------------------------------------------------------------------------------


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

#ls -1 $bam_dir/*.bam > bam.list

while read name; do
	if [ ! -f ${name%.bam}.bai ]; then
	 samtools index -b ${name}
	fi; done < ${bamlist}

while read number sequences; do \
echo -e "samtools mpileup -b $bamlist -r $sequences -f ref.tmp -t DP,SP,AD,ADF,INFO/AD -F0.25 -d500 -E -o $number.$sequences.tmp.mpileup" > run_mpileup.tmp.$number;
done < ref_seq.tmp.list
chmod a+x run_mpileup.tmp*

#while read number sequences; do \
#echo -e "\
#samtools-1.3 mpileup -b $bamlist -r $sequences -f ref.tmp -t DP,SP,AD,ADF,INFO/AD -F0.25 -d500 -E -o $number.$sequences.tmp.mpileup" > run_mp_$sequences.$number; done < ref_seq.tmp.list
jobs=$( wc -l ref_seq.tmp.list | cut -f1 -d" " )
bsub -q long -R'span[hosts=1] select[mem>10000] rusage[mem=10000]' -M10000 -J mpileup_array[1-$jobs] -e mpileup_array[1-$jobs].e -o mpileup_array[1-$jobs].o ./run_mpileup.tmp.\$LSB_JOBINDEX



#--- Step 2: run mpileup on split chromosomes/sequences
#chmod a+x run_mpileup.tmp
#./run_mpileup.tmp


#--- Step 3: check to see all mpileup jobs are complete. If so, them merge into a single mpileup file
#a=$( ls -1 *.FINISHED | wc -l )
#b=$( cat ref_seq.tmp.list | wc -l )
#while [[ $a != $b ]]
#	do
#	sleep 100
#	a=$( ls -1 *.FINISHED | wc -l )
#	b=$( cat ref_seq.tmp.list | wc -l )
#	done
#cat $(find ./ -name "*.mpileup" | sort -V) > ${prefix}.mpileup; rm *tmp*



echo -e "cat \$(find ./ -name "*.mpileup" | sort -V) > \${prefix}.mpileup; touch mpileup_FINISHED" > run_mp_combine
chmod a+x run_mp_combine

bsub -q normal -w "done(mpileup_array)"  -R'span[hosts=1] select[mem>1000] rusage[mem=1000]' -n1 -M1000 -J mp_combine -e mp_combine.e -o mp_combine.o ./run_mp_combine

# if final mpileup throws strange error in popoolation2 for example, can be fixed using the following
#sed 's/\t\t/\t!\t!/g' final.mpileup > final.mpileup2
