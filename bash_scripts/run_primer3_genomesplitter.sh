#!/usr/bin/env bash

# ========================================================================================
# run_primer3_genomesplitter.sh
#
# Primer design tool to automate genome-wide primer design, either whole genome, or in specific coordinates. Will design primers evently spaced throguhout the genome or region using a defined window size between primers
#
# Reference:
# Software:
#
# Usage: ~sd21/bash_scripts/run_primer3_genomesplitter.sh <prefix> <reference.fa> <window_size> <window_step> [optional: specific coordinates]
#
# @authors
# Stephen Doyle <sd21@sanger.ac.uk>
#
# ----------------------------------------------------------------------------------------

if [ "$#" -eq 0 ] || [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
	echo ""
    echo "Usage: ./run_primer3_genomesplitter.sh <prefix> <reference.fa> <window_size> <window_step> [optional: specific coordinates]"
    echo ""
    echo "Usage example: Hc_XQTL_IVM_chr5_45-47Mb HAEM_V4_final.chr.fa 10000 200000"
    echo "#>> will result in a set of primers every ~200kb throughout the genome"
    echo "Usage example2: Hc_XQTL_IVM_chr5_45-47Mb HAEM_V4_final.chr.fa 10000 200000 hcontortus_chr5_Celeg_TT_arrow_pilon:45000000-47000000"
    echo "#>> will result in a set of primers every ~200kb throughout chromosome 5 region between 45-47 Mb"
    echo ""
    exit 0
fi


PREFIX=$1
REFERENCE=$2
WINDOW_SIZE=$3
WINDOW_STEP=$4
COORDS=$5


# step1 - get reference

#cp ${REFERENCE} REF.tmp.fa
#samtools-1.3 faidx REF.tmp.fa
#cut -f1,2 REF.tmp.fa.fai > REF.tmp.genome

if [ "$#" -eq  "5" ]
then
    cp ${REFERENCE} REF.tmp.${PREFIX}.tmp.fa
    samtools-1.3 faidx REF.tmp.${PREFIX}.tmp.fa ${COORDS} > REF.tmp.${PREFIX}.fa
    samtools-1.3 faidx REF.tmp.${PREFIX}.fa
    cut -f1,2 REF.tmp.${PREFIX}.fa.fai > REF.tmp.${PREFIX}.genome
else
   	cp ${REFERENCE} REF.tmp.${PREFIX}.fa
	samtools-1.3 faidx REF.tmp.${PREFIX}.fa
	cut -f1,2 REF.tmp.${PREFIX}.fa.fai > REF.tmp.${PREFIX}.genome
fi




# step2 - design windows
bedtools-2 makewindows -g REF.tmp.${PREFIX}.genome -w ${WINDOW_SIZE} -s ${WINDOW_STEP} > REF.tmp.${PREFIX}.${WINDOW_SIZE}bp_windows_${WINDOW_STEP}bp_apart.bed


# step3 - get fastas of windows
bedtools-2 getfasta -fi REF.tmp.${PREFIX}.fa -bed REF.tmp.${PREFIX}.${WINDOW_SIZE}bp_windows_${WINDOW_STEP}bp_apart.bed > REF.tmp.${PREFIX}.${WINDOW_SIZE}bp_windows_${WINDOW_STEP}bp_apart.fasta
grep ">" REF.tmp.${PREFIX}.${WINDOW_SIZE}bp_windows_${WINDOW_STEP}bp_apart.fasta | sed 's/>//g' > REF.tmp.${PREFIX}.${WINDOW_SIZE}bp_windows_${WINDOW_STEP}bp_apart.sequences.list

# step4 - run primer3


while read NAME SEQUENCE; do
grep -A 1 -w ${NAME} REF.tmp.${PREFIX}.${WINDOW_SIZE}bp_windows_${WINDOW_STEP}bp_apart.fasta | awk '{ORS=(NR%2?FS:RS)}1' OFS="\t" | sed 's/>//g' > ${NAME}.tmp.sequence; \
done < REF.tmp.${PREFIX}.${WINDOW_SIZE}bp_windows_${WINDOW_STEP}bp_apart.sequences.list

n=1
for i in *.sequence; do \
NAME=`cut -f1 -d " " ${i}`
SEQUENCE=`cut -f2 -d " " ${i}`
echo -e "\
SEQUENCE_ID=${NAME}
SEQUENCE_TEMPLATE=${SEQUENCE}
PRIMER_TASK=pick_pcr_primers
PRIMER_OPT_SIZE=19
PRIMER_MIN_SIZE=17
PRIMER_MAX_SIZE=21
PRIMER_MAX_NS_ACCEPTED=1
PRIMER_MIN_TM=58
PRIMER_MAX_TM=62
PRIMER_PAIR_MAX_DIFF_TM=5
PRIMER_PRODUCT_SIZE_RANGE=250-300
PRIMER_GC_CLAMP=1
PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/nfs/users/nfs_s/sd21/lustre118_link/software/primer3/src/primer3_config/
PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1
PRIMER_MAX_SELF_ANY_TH=40
PRIMER_PAIR_MAX_COMPL_ANY_TH=40
PRIMER_MAX_SELF_END_TH=40
PRIMER_MAX_HAIRPIN_TH=40
PRIMER_MAX_POLY_X=4
P3_FILE_FLAG=0
PRIMER_EXPLAIN_FLAG=1
=" > ${NAME}_p3.config.tmp.${n};
let "n+=1"; done



# Step 4 - extract primer data
# to do
#--- extract primer coordinates in window, and adjust based on genome to actual genome coordinates
#make files
echo -e "#experiment_prefix\tpcr_id\tprimer1_name\tprimer1_sequence\tprimer1_tm\tprimer2_name\tprimer2_sequence\tprimer2_tm\tproduct_size" > primers.${PREFIX}.${WINDOW_SIZE}bp_windows_${WINDOW_STEP}bp_apart.database
> primers.${PREFIX}.${WINDOW_SIZE}bp_windows_${WINDOW_STEP}bp_apart.fasta
> primers.${PREFIX}.${WINDOW_SIZE}bp_windows_${WINDOW_STEP}bp_apart.bed
> primers.${PREFIX}_MISSING_primers.list

# extract data

n=1
for i in ` ls -1 *.config* | sort -V `; do

~sd21/lustre118_link/software/bin/primer3_core < ${i} > ${i}.runlog ;
	if [ -z "$(grep PRIMER_LEFT_0 ${i}.runlog)" ]
	then
		echo "No suitable primers found in ${PREFIX} ${i}" >> primers.${PREFIX}_MISSING_primers.list
	else
		if [ "$#" -eq  "5" ]
		then
			CHR=`echo $COORDS | cut -f 1 -d ":" `
			WINDOW_START1=`grep "SEQUENCE_ID=" ${i}.runlog | cut -f 2 -d ":" | cut -f1 -d "-" `
			WINDOW_START2=`grep "SEQUENCE_ID=" ${i}.runlog | cut -f 3 -d ":" | cut -f1 -d "-" `
			WINDOW_START=$((WINDOW_START1+WINDOW_START2))

		else
			CHR=`grep "SEQUENCE_ID=" ${i}.runlog | sed 's/SEQUENCE_ID=//g' | cut -f 1 -d ":" `
			WINDOW_START=`grep "SEQUENCE_ID=" ${i}.runlog | sed 's/SEQUENCE_ID=//g' | awk -F":" '{print $NF}' | cut -f1 -d "-" `
		fi

		WINDOW_END=`grep "SEQUENCE_ID=" ${i}.runlog | sed 's/SEQUENCE_ID=//g' | awk -F":" '{print $NF}' | cut -f2 -d "-" `;

		P1_START=`grep "PRIMER_LEFT_0=" ${i}.runlog | sed 's/PRIMER_LEFT_0=//g' | cut -f1 -d "," `;
		P1_END=`grep "PRIMER_LEFT_0=" ${i}.runlog | sed 's/PRIMER_LEFT_0=//g' | awk -F"," '{print $1+$2}' `;
		P2_START=`grep "PRIMER_RIGHT_0=" ${i}.runlog | sed 's/PRIMER_RIGHT_0=//g' | cut -f1 -d "," `;
		P2_END=`grep "PRIMER_RIGHT_0=" ${i}.runlog | sed 's/PRIMER_RIGHT_0=//g' | awk -F"," '{print $1+$2}' `;

		P1_SEQ=`grep "PRIMER_LEFT_0_SEQUENCE=" ${i}.runlog | sed 's/PRIMER_LEFT_0_SEQUENCE=//g'`;
		P2_SEQ=`grep "PRIMER_RIGHT_0_SEQUENCE=" ${i}.runlog | sed 's/PRIMER_RIGHT_0_SEQUENCE=//g'`;
		P1_TM=`grep "PRIMER_LEFT_0_TM=" ${i}.runlog | sed 's/PRIMER_LEFT_0_TM=//g'`;
		P2_TM=`grep "PRIMER_RIGHT_0_TM=" ${i}.runlog | sed 's/PRIMER_RIGHT_0_TM=//g'`;
		PRODUCT_SIZE=`grep "PRIMER_PAIR_0_PRODUCT_SIZE=" ${i}.runlog | sed 's/PRIMER_PAIR_0_PRODUCT_SIZE=//g'`;

		# collate output files
		echo -e ">${CHR}:$((WINDOW_START+P1_START))-$((WINDOW_START+P1_END))_PCR${n}_F\n${P1_SEQ}\n>${CHR}:$((WINDOW_START+P2_START))-$((WINDOW_START+P2_END))_PCR${n}_R\n${P2_SEQ}" >> primers.${PREFIX}.${WINDOW_SIZE}bp_windows_${WINDOW_STEP}bp_apart.fasta;
		echo -e "${PREFIX}\tPCR${n}\t${CHR}:$((WINDOW_START+P1_START))-$((WINDOW_START+P1_END))_PCR${n}_F\t${P1_SEQ}\t${P1_TM}\t${CHR}:$((WINDOW_START+P2_START))-$((WINDOW_START+P2_END))_PCR${n}_R\t${P2_SEQ}\t${P2_TM}\t${PRODUCT_SIZE}" >> primers.${PREFIX}.${WINDOW_SIZE}bp_windows_${WINDOW_STEP}bp_apart.database;
		echo -e "${CHR}\t$((WINDOW_START+P1_START))\t$((WINDOW_START+P1_END))\t${PREFIX}:${CHR}:$((WINDOW_START+P1_START))-$((WINDOW_START+P1_END))_PCR${n}_F\n${CHR}\t$((WINDOW_START+P2_START))\t$((WINDOW_START+P2_END))\t${PREFIX}:${CHR}:$((WINDOW_START+P2_START))-$((WINDOW_START+P2_END))_PCR${n}_R" >> primers.${PREFIX}.${WINDOW_SIZE}bp_windows_${WINDOW_STEP}bp_apart.bed;
		let "n+=1";
	fi;
done

mkdir LOGFILES.${PREFIX}.${WINDOW_SIZE}bp_windows_${WINDOW_STEP}bp_apart
mv *tmp* LOGFILES.${PREFIX}.${WINDOW_SIZE}bp_windows_${WINDOW_STEP}bp_apart


# check primer dimers: http://www.primer-dimer.com/
