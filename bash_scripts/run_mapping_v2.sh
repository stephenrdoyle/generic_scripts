#!/bin/bash
#-------------------------------------------------------------------------------
# run_mapping.sh
#-------------------------------------------------------------------------------

### author: Stephen Doyle, stephen.doyle[at]sanger.ac.uk

# Input required
PREFIX=mapping_hcontortus_qtl
REF=/lustre/scratch125/pam/teams/team333/sd21/haemonchus_contortus/QTL/01_REFERENCE/HAEM_V4_final.chr.fa
DATA_DIR=/lustre/scratch125/pam/teams/team333/sd21/haemonchus_contortus/QTL/02_RAW
LANES_IDS=/lustre/scratch125/pam/teams/team333/sd21/haemonchus_contortus/QTL/lanes_samples.list

THREADS=20

[ -d MAPPING_${PREFIX} ] || mkdir MAPPING_${PREFIX}
WORKING_DIR=${PWD}/MAPPING_${PREFIX}

# save current script in run folder to reproduce the exact output and make a progress log
cp ${PWD}/run_mapping.sh ${WORKING_DIR}/commands.$(date -Iminutes).txt
cp ${LANES_IDS} ${WORKING_DIR}/samples.$(date -Iminutes).txt

touch ${WORKING_DIR}/progress.log


# load GATK 4.1.4.1
module load gatk/4.1.4.1
module load samtools/1.14--hb421002_0
module load minimap2/2.16=h84994c4_1-c1
#conda install -c bioconda sambamba


# reference 
#---  setup directory a reference and get reference 
[ -d ${WORKING_DIR}/REF ] || mkdir ${WORKING_DIR}/REF
REF_DIR=${WORKING_DIR}/REF

cp ${REF} ${REF_DIR}/REF.fa

#--- CreateSequenceDictionary
[ -f ${REF_DIR}/REF.fa.fai ] || samtools faidx ${REF_DIR}/REF.fa

[ -f ${REF_DIR}/REF.dict ] || gatk --java-options "-Xmx8g -Xms4g -Djava.io.tmpdir=${TMP_DIR}" CreateSequenceDictionary \
     --REFERENCE ${REF_DIR}/REF.fa \
     --OUTPUT ${REF_DIR}/REF.dict \
     --spark-runner LOCAL

#---  setup directory for storing mapping job scrips
[ -d ${WORKING_DIR}/MAPPING_JOBS ] || mkdir ${WORKING_DIR}/MAPPING_JOBS
MAPPING_JOBS=${WORKING_DIR}/MAPPING_JOBS



# Mapping 
#--- interate through samples in data directory
n=0
#for FASTQ in `ls -1 ${DATA_DIR}/*_1.fastq.gz`; do

while read -r LANE_ID SAMPLE_ID; do

n=$((n + 1))
count=$( cat ${LANES_IDS} | wc -l )

#SAMPLE_ID=$(basename ${FASTQ} _1.fastq.gz)


if [ -d "${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping" ]; then
     echo -e "\nThere is already a run started with the sample name ${SAMPLE_ID}_${LANE_ID}. Rename and start again\n"
     else

mkdir -p ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/.tmp

FASTQ=`ls -1 ${DATA_DIR}/${LANE_ID}_1.fastq.gz`

PLATFORM_UNIT=$( zcat ${FASTQ} | head -n1 | awk '{for(i=1;i<=NF;i++){if($i~/\:/){a=$i}} print a}' | cut -f1 -d ":" | sed 's/@//g')
DATE=$(date -Iminutes)

echo -e "
# mapping script for ${SAMPLE_ID}

module load gatk/4.1.4.1
module load samtools/1.14--hb421002_0
module load minimap2/2.16=h84994c4_1-c1

DATE=\$(date -Iminutes)

echo -e \""Started mapping" ${SAMPLE_ID} "\(sample: "${n}" of "${count}"\)" \${DATE} \" >> ${WORKING_DIR}/progress.log


# FastqToSam: convert fastq to uBAM
gatk --java-options \"-Xmx8g -Xms4g -Djava.io.tmpdir=${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/.tmp -Dsamjdk.compression_level=1\" FastqToSam \
     --FASTQ ${DATA_DIR}/${LANE_ID}_1.fastq.gz \
     --FASTQ2 ${DATA_DIR}/${LANE_ID}_2.fastq.gz \
     --OUTPUT ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.unaligned.bam \
     --SAMPLE_NAME ${SAMPLE_ID} \
     --LIBRARY_NAME ${SAMPLE_ID}.lib \
     --READ_GROUP_NAME ${SAMPLE_ID}_${LANE_ID} \
     --SORT_ORDER queryname \
     --TMP_DIR ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/.tmp \
     --spark-runner LOCAL

# MarkIlluminaAdapters: mark adapters in the raw reads
gatk --java-options \"-Xmx8g -Xms4g -Djava.io.tmpdir=${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/.tmp -Dsamjdk.compression_level=1\" MarkIlluminaAdapters \
     --INPUT ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.unaligned.bam \
     --OUTPUT ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.adaptMarked.bam \
     --METRICS ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.adaptMarked.metrics.txt \
     --TMP_DIR ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/.tmp \
     --spark-runner LOCAL

# SamToFastq: make interleaved fastq for mapping
gatk --java-options \"-Xmx8g -Xms4g -Djava.io.tmpdir=${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/.tmp -Dsamjdk.compression_level=1\" SamToFastq \
     --INPUT ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.adaptMarked.bam \
     --FASTQ ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.interleaved.fastq.gz \
     --CLIPPING_ATTRIBUTE XT --CLIPPING_ACTION 2 --INTERLEAVE true --INCLUDE_NON_PF_READS true \
     --TMP_DIR ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/.tmp \
     --spark-runner LOCAL

# mapping
# bwa mem -K 100000000 -v 3 -t ${THREADS} -Y -p ${REF_DIR}/REF.fa ${SAMPLE_ID}_mapping/${SAMPLE_ID}.interleaved.fastq.gz | samtools view -h -b- > ${SAMPLE_ID}_mapping/${SAMPLE_ID}.aligned.bam

minimap2 -t ${THREADS} -a -Y -x sr ${REF_DIR}/REF.fa ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.interleaved.fastq.gz | samtools view -h -b - > ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.aligned.bam

# MergeBamAlignment: sort alignment
gatk --java-options \"-Xmx8g -Xms4g -Djava.io.tmpdir=${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/.tmp -Dsamjdk.compression_level=1\" MergeBamAlignment \
     --REFERENCE_SEQUENCE ${REF_DIR}/REF.fa \
     --UNMAPPED_BAM ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.unaligned.bam \
     --ALIGNED_BAM ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.aligned.bam \
     --OUTPUT ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.alnMerged.bam \
     --CREATE_INDEX false --ADD_MATE_CIGAR true --CLIP_ADAPTERS true --CLIP_OVERLAPPING_READS true --INCLUDE_SECONDARY_ALIGNMENTS true --MAX_INSERTIONS_OR_DELETIONS -1 --PRIMARY_ALIGNMENT_STRATEGY BestMapq --ATTRIBUTES_TO_RETAIN XS \
     --spark-runner LOCAL

# MarkDuplicatesSpark: mark duplicates
#gatk --java-options "-Xmx8g -Xms4g -Djava.io.tmpdir=${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/.tmp -Dsamjdk.compression_level=5" MarkDuplicatesSpark \
#     --input ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.alnMerged.bam \
#     --output ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.deduped.bam \
#     --metrics-file ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.deduped.metrics.txt \
#     --spark-runner LOCAL

sambamba markdup --nthreads ${THREADS} ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.alnMerged.bam ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.deduped.bam

samtools flagstat ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.deduped.bam > ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.prefilter.flagstat
samtools stats ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.deduped.bam > ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.prefilter.stats;

samtools view --threads 4 -F 12 -b ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.deduped.bam -o ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.bam
samtools view --threads 4 -f 12 ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.deduped.bam -o ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.unmapped.bam

samtools index -b ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.bam

# cleanup
rm ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.unaligned.bam
rm ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.adaptMarked.bam
rm ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.aligned.bam
rm ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.alnMerged.bam
rm ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.interleaved.fastq.gz
rm ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.deduped.bam.bai
rm ${WORKING_DIR}/${SAMPLE_ID}_${LANE_ID}_mapping/${SAMPLE_ID}.deduped.bam

echo -e \" "Finished mapping" ${SAMPLE_ID} "\(sample: "${n}" of "${count}"\)" \${DATE} \"  >> ${WORKING_DIR}/progress.log

" > ${MAPPING_JOBS}/map_${SAMPLE_ID}.job_${n}

fi;

done < ${LANES_IDS}

chmod a+x ${MAPPING_JOBS}/map*

# setup job conditions
JOBS=$( ls -1 ${MAPPING_JOBS}/map* | wc -l )

#submit job array to call variants put scaffold / contig
bsub -q normal -R'span[hosts=1] select[mem>10000] rusage[mem=10000]' -n ${THREADS} -M10000 -J "mapping_${PREFIX}_[1-$JOBS]%20" -e "${WORKING_DIR}/mapping_${PREFIX}_[1-$JOBS].e" -o "${WORKING_DIR}/mapping_${PREFIX}_[1-$JOBS].o" "${MAPPING_JOBS}/map_*.job_\$LSB_JOBINDEX"


