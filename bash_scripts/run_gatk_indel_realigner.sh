#!/usr/bin/env bash

# ========================================================================================
# run_gatk_indel_realigner.sh
#
# Realign indels using GATK for bams that will either be used in popoolation2 or unified genotyper. This does not have to be done if using GATK haplotyp caller.
#
# Reference:
# Software:
#
# Usage: ~sd21/bash_scripts/run_gatk_indel_realigner.sh <REFERENCE_FASTA> <BAMLIST>
#
# @authors
# Stephen Doyle <sd21@sanger.ac.uk>
#
# ----------------------------------------------------------------------------------------

REFERENCE=$1
BAM_LIST=$2

module load gatk/3.7.0

samtools faidx ${REFERENCE}
samtools dict  ${REFERENCE} >  ${REFERENCE%.fa}.dict

echo -e "/software/pathogen/etc/gatk/3.7.0/wrappers/gatk -T RealignerTargetCreator --num_threads 7 -R ${REFERENCE} \\" > run_gatk_RealignerTargetCreator
while read SAMPLE; do
echo -e "--input_file ${SAMPLE} \\" >> run_gatk_RealignerTargetCreator;
done < ${BAM_LIST}
echo -e "-o all_samples.intervals" >> run_gatk_RealignerTargetCreator

chmod a+x run_gatk_RealignerTargetCreator
bsub -q long -n 15 -R'span[hosts=1] select[mem>10000] rusage[mem=10000]' -M10000 -J step_01_gatk_indel_realigner -e step_01_gatk_indel_realigner.e -o step_01_gatk_indel_realigner.o ./run_gatk_RealignerTargetCreator

while read SAMPLE; do
echo -e "bsub -q long -w "step_01_gatk_indel_realigner" -R'span[hosts=1] select[mem>10000] rusage[mem=10000]' -M10000 -J step_02_gatk_indel_realigner -e step_02_gatk_indel_realigner.e -o  step_02_gatk_indel_realigner.o /software/pathogen/etc/gatk/3.7.0/wrappers/gatk -T IndelRealigner \
-R ${REFERENCE} \\
-I ${SAMPLE} \\
-targetIntervals all_samples.intervals \\
-o ${SAMPLE%.bam}.realigned.bam" >> run_gatk_IndelRealigner; done < ${BAM_LIST}


chmod a+x run_gatk_IndelRealigner
./run_gatk_IndelRealigner
