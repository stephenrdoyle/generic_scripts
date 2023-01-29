#!/bin/bash
#-------------------------------------------------------------------------------
# run_mapping.sh
#-------------------------------------------------------------------------------

### author: Stephen Doyle, stephen.doyle[at]sanger.ac.uk

# Input required
PREFIX=STOP_Trichuris
REF=/nfs/users/nfs_s/sd21/lustre118_link/STH/STOP_TRICHURIS/REF/STH_multispecies.renamed.fa
WORKING_DIR=/nfs/users/nfs_s/sd21/lustre118_link/STH/STOP_TRICHURIS
THREADS=20


# should already exist 
DATA_DIR=${WORKING_DIR}/DATA

mkdir ${PREFIX}_MAPPING
MAPPING_DIR=${WORKING_DIR}/${PREFIX}_MAPPING


# save current script in run folder to reproduce the exact output and make a progress log
cp ${PWD}/run_mapping.sh ${MAPPING_DIR}/commands.$(date -Iminutes).txt
touch ${MAPPING_DIR}/progress.log

# load GATK 4.1.4.1
module load gatk/4.1.4.1
module load samtools/1.14--hb421002_0
module load minimap2/2.16=h84994c4_1-c1
#conda install -c bioconda sambamba


# reference 
#---  setup directory a reference and get reference 
[ -d ${MAPPING_DIR}/REF ] || mkdir ${MAPPING_DIR}/REF
REF_DIR=${MAPPING_DIR}/REF

cp ${REF} ${REF_DIR}/REF.fa

#--- CreateSequenceDictionary
[ -f ${REF_DIR}/REF.fa.fai ] || samtools faidx ${REF_DIR}/REF.fa

[ -f ${REF_DIR}/REF.dict ] || gatk --java-options "-Xmx8g -Xms4g -Djava.io.tmpdir=${TMP_DIR}" CreateSequenceDictionary \
     --REFERENCE ${REF_DIR}/REF.fa \
     --OUTPUT ${REF_DIR}/REF.dict \
     --spark-runner LOCAL


# Mapping 
#--- interate through samples in data directory
n=0
for FASTQ in `ls -1 ${DATA_DIR}/*_1.fastq.gz`; do

n=$((n + 1))
count=$(ls -1 ${DATA_DIR}/*_1.fastq.gz | wc -l)

DATE=$(date -Iminutes)

SAMPLE_ID=$(basename ${FASTQ} _1.fastq.gz)

echo "Started mapping" ${SAMPLE_ID} "(sample: "${n}" of "${count}")" ${DATE}  >> ${MAPPING_DIR}/progress.log

if [ -d "${MAPPING_DIR}/${PREFIX}_mapping" ]; then
     echo -e "\nThere is already a run started with this sample name. Rename and start again\n"
     exit 0
fi

mkdir -p ${MAPPING_DIR}/${SAMPLE_ID}_mapping ${MAPPING_DIR}/${SAMPLE_ID}_mapping/.tmp

PLATFORM_UNIT=$( zcat ${FASTQ} | head -n1 | awk '{for(i=1;i<=NF;i++){if($i~/\:/){a=$i}} print a}' | cut -f1 -d ":" | sed 's/@//g')
DATE=$(date -Iminutes)

# FastqToSam: convert fastq to uBAM
gatk --java-options "-Xmx8g -Xms4g -Djava.io.tmpdir=${MAPPING_DIR}/${SAMPLE_ID}_mapping/.tmp -Dsamjdk.compression_level=1" FastqToSam \
     --FASTQ ${DATA_DIR}/${SAMPLE_ID}_1.fastq.gz \
     --FASTQ2 ${DATA_DIR}/${SAMPLE_ID}_2.fastq.gz \
     --OUTPUT ${MAPPING_DIR}/${SAMPLE_ID}_mapping/${SAMPLE_ID}.unaligned.bam \
     --SAMPLE_NAME ${SAMPLE_ID} \
     --LIBRARY_NAME ${SAMPLE_ID}.lib \
     --READ_GROUP_NAME ${SAMPLE_ID} \
     --SORT_ORDER queryname \
     --TMP_DIR ${MAPPING_DIR}/${SAMPLE_ID}_mapping/.tmp \
     --spark-runner LOCAL

# MarkIlluminaAdapters: mark adapters in the raw reads
gatk --java-options "-Xmx8g -Xms4g -Djava.io.tmpdir=${SAMPLE_ID}_mapping/.tmp -Dsamjdk.compression_level=1" MarkIlluminaAdapters \
     --INPUT ${MAPPING_DIR}/${SAMPLE_ID}_mapping/${SAMPLE_ID}.unaligned.bam \
     --OUTPUT ${MAPPING_DIR}/${SAMPLE_ID}_mapping/${SAMPLE_ID}.adaptMarked.bam \
     --METRICS ${MAPPING_DIR}/${SAMPLE_ID}_mapping/${SAMPLE_ID}.adaptMarked.metrics.txt \
     --TMP_DIR ${MAPPING_DIR}/${SAMPLE_ID}_mapping/.tmp \
     --spark-runner LOCAL

# SamToFastq: make interleaved fastq for mapping
gatk --java-options "-Xmx8g -Xms4g -Djava.io.tmpdir=${MAPPING_DIR}/${SAMPLE_ID}_mapping/.tmp -Dsamjdk.compression_level=1" SamToFastq \
     --INPUT ${MAPPING_DIR}/${SAMPLE_ID}_mapping/${SAMPLE_ID}.adaptMarked.bam \
     --FASTQ ${MAPPING_DIR}/${SAMPLE_ID}_mapping/${SAMPLE_ID}.interleaved.fastq.gz \
     --CLIPPING_ATTRIBUTE XT --CLIPPING_ACTION 2 --INTERLEAVE true --INCLUDE_NON_PF_READS true \
     --TMP_DIR ${MAPPING_DIR}/${SAMPLE_ID}_mapping/.tmp \
     --spark-runner LOCAL

# mapping
# bwa mem -K 100000000 -v 3 -t ${THREADS} -Y -p ${REF_DIR}/REF.fa ${SAMPLE_ID}_mapping/${SAMPLE_ID}.interleaved.fastq.gz | samtools view -h -b- > ${SAMPLE_ID}_mapping/${SAMPLE_ID}.aligned.bam

minimap2 -t ${THREADS} -a -Y -x sr ${REF_DIR}/REF.fa ${MAPPING_DIR}/${SAMPLE_ID}_mapping/${SAMPLE_ID}.interleaved.fastq.gz | samtools view -h -b - > ${MAPPING_DIR}/${SAMPLE_ID}_mapping/${SAMPLE_ID}.aligned.bam

# MergeBamAlignment: sort alignment
gatk --java-options "-Xmx8g -Xms4g -Djava.io.tmpdir=${MAPPING_DIR}/${SAMPLE_ID}_mapping/.tmp -Dsamjdk.compression_level=1" MergeBamAlignment \
     --REFERENCE_SEQUENCE ${REF_DIR}/REF.fa \
     --UNMAPPED_BAM ${MAPPING_DIR}/${SAMPLE_ID}_mapping/${SAMPLE_ID}.unaligned.bam \
     --ALIGNED_BAM ${MAPPING_DIR}/${SAMPLE_ID}_mapping/${SAMPLE_ID}.aligned.bam \
     --OUTPUT ${MAPPING_DIR}/${SAMPLE_ID}_mapping/${SAMPLE_ID}.alnMerged.bam \
     --CREATE_INDEX false --ADD_MATE_CIGAR true --CLIP_ADAPTERS true --CLIP_OVERLAPPING_READS true --INCLUDE_SECONDARY_ALIGNMENTS true --MAX_INSERTIONS_OR_DELETIONS -1 --PRIMARY_ALIGNMENT_STRATEGY BestMapq --ATTRIBUTES_TO_RETAIN XS \
     --spark-runner LOCAL

# MarkDuplicatesSpark: mark duplicates
#gatk --java-options "-Xmx8g -Xms4g -Djava.io.tmpdir=${MAPPING_DIR}/${SAMPLE_ID}_mapping/.tmp -Dsamjdk.compression_level=5" MarkDuplicatesSpark \
#     --input ${MAPPING_DIR}/${SAMPLE_ID}_mapping/${SAMPLE_ID}.alnMerged.bam \
#     --output ${MAPPING_DIR}/${SAMPLE_ID}_mapping/${SAMPLE_ID}.deduped.bam \
#     --metrics-file ${MAPPING_DIR}/${SAMPLE_ID}_mapping/${SAMPLE_ID}.deduped.metrics.txt \
#     --spark-runner LOCAL

sambamba markdup --nthreads ${THREADS} ${MAPPING_DIR}/${SAMPLE_ID}_mapping/${SAMPLE_ID}.alnMerged.bam ${MAPPING_DIR}/${SAMPLE_ID}_mapping/${SAMPLE_ID}.deduped.bam

samtools flagstat ${MAPPING_DIR}/${SAMPLE_ID}_mapping/${SAMPLE_ID}.deduped.bam > ${MAPPING_DIR}/${SAMPLE_ID}_mapping/${SAMPLE_ID}.deduped.bam.flagstat
samtools stats ${MAPPING_DIR}/${SAMPLE_ID}_mapping/${SAMPLE_ID}.deduped.bam > ${MAPPING_DIR}/${SAMPLE_ID}_mapping/${SAMPLE_ID}.deduped.bam.stats;

# cleanup
rm ${MAPPING_DIR}/${SAMPLE_ID}_mapping/${SAMPLE_ID}.unaligned.bam
rm ${MAPPING_DIR}/${SAMPLE_ID}_mapping/${SAMPLE_ID}.adaptMarked.bam
rm ${MAPPING_DIR}/${SAMPLE_ID}_mapping/${SAMPLE_ID}.aligned.bam
rm ${MAPPING_DIR}/${SAMPLE_ID}_mapping/${SAMPLE_ID}.alnMerged.bam

echo "Finished mapping" ${SAMPLE_ID} "(sample: "${n}" of "${count}")" ${DATE}  >> ${MAPPING_DIR}/progress.log

done

