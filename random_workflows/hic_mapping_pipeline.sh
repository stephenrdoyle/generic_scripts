#!/bin/bash

#-------------------------------------------------------------------------------
# run_hic_pipeline.sh
#-------------------------------------------------------------------------------

# usage: ./run_hic_pipeline.sh

# author: Stephen R. Doyle
# contact: stephen.doyle[at]sanger.ac.uk
# date: 2022/11/16

# input required
export PREFIX=tcircumcinca_dnazoo
export REFERENCE=/nfs/users/nfs_s/sd21/lustre118_link/teladorsagia_circumcincta/RAPID_CURATION/CUSTOM_HIC/T_circumcincta.14.0.ec.cg.pg_purgehaplotigs_HiC.fasta
export HIC_READS_R1=/nfs/users/nfs_s/sd21/lustre118_link/teladorsagia_circumcincta/RAPID_CURATION/CUSTOM_HIC/42782_8_1/42782_8#1_1.fastq.gz
export HIC_READS_R2=/nfs/users/nfs_s/sd21/lustre118_link/teladorsagia_circumcincta/RAPID_CURATION/CUSTOM_HIC/42782_8_1/42782_8#1_2.fastq.gz
export PACBIO_READS=/nfs/users/nfs_s/sd21/lustre118_link/teladorsagia_circumcincta/raw/pacbio/2017_test_noWGA/all_subreads.fasta
export TELOMERE_SEQ=gcctaa


# requirements
# adapter from https://github.com/ArimaGenomics/mapping_pipeline
# arima tools:
# cooler:
# pretext:
# picard:
# seqtk: module load seqtk/1.3--ha92aebf_0
# samtools:
# fastaq:
# juicer: https://github.com/aidenlab/juicer/releases/tag/1.6



## file locations
export TMP_DIR="$PWD/hicpipe_${PREFIX}_out/tmp"
export LOG_FILES="$PWD/hicpipe_${PREFIX}_out/00_LOG_FILES"
export REF_DIR="$PWD/hicpipe_${PREFIX}_out/01_REFERENCE"
export MAPPING_DIR="$PWD/hicpipe_${PREFIX}_out/02_HIC_MAPPING"
export PRETEXT_DIR="$PWD/hicpipe_${PREFIX}_out/03_PRETEXT"
export HIGLASS_DIR="$PWD/hicpipe_${PREFIX}_out/04_HIGLASS"
export COVERAGE_DIR="$PWD/hicpipe_${PREFIX}_out/05_PACBIO_COVERAGE"
export GAPS_DIR="$PWD/hicpipe_${PREFIX}_out/06_GAPS"
export REPEATS_DIR="$PWD/hicpipe_${PREFIX}_out/07_REPEATS"
export TELOMERES_DIR="$PWD/hicpipe_${PREFIX}_out/08_TELOMERES"
export JUICER_DIR="$PWD/hicpipe_${PREFIX}_out/09_JUICER"


## tool locations
export FILTER='/nfs/users/nfs_s/sd21/lustre118_link/software/HIC/mapping_pipeline/filter_five_end.pl'
export COMBINER='/nfs/users/nfs_s/sd21/lustre118_link/software/HIC/mapping_pipeline/two_read_bam_combiner.pl'
export STATS='/nfs/users/nfs_s/sd21/lustre118_link/software/HIC/mapping_pipeline/get_stats.pl'
export PICARD='/nfs/users/nfs_s/sd21/lustre118_link/software/picard-tools-2.5.0/picard.jar'
export JUICER='/nfs/users/nfs_s/sd21/lustre118_link/software/HIC/juicer-1.6/scripts/common/juicer_tools.jar'


## Other
export CPU=30
export MAPQ_FILTER=10
export FASTQ_BN=$(basename ${HIC_READS_R1} _1.fastq.gz)



### Check output directories exist & create them as needed
[ -d ${REF_DIR} ] || mkdir -p ${REF_DIR}
[ -d $MAPPING_DIR ] || mkdir -p $MAPPING_DIR
[ -d ${PRETEXT_DIR} ] || mkdir -p ${PRETEXT_DIR}
[ -d ${HIGLASS_DIR} ] || mkdir -p ${HIGLASS_DIR}
[ -d ${TMP_DIR} ] || mkdir -p ${TMP_DIR}
[ -d ${COVERAGE_DIR} ] || mkdir -p ${COVERAGE_DIR}
[ -d ${LOG_FILES} ] || mkdir -p ${LOG_FILES}
[ -d ${GAPS_DIR} ] || mkdir -p ${GAPS_DIR}
[ -d ${REPEATS_DIR} ] || mkdir -p ${REPEATS_DIR}
[ -d ${TELOMERES_DIR} ] || mkdir -p ${TELOMERES_DIR}
[ -d ${JUICER_DIR} ] || mkdir -p ${JUICER_DIR}


# save current script in run folder to reproduce the exact output
cp ${PWD}/run_hic_pipeline.sh ${PWD}/hicpipe_${PREFIX}_out/commands.txt

# make a progress file
> ${PWD}/hicpipe_${PREFIX}.progress.log


#-------------------------------------------------------------------------------
### Prepare reference files
#-------------------------------------------------------------------------------
func_build_reference() {

echo "Started: building reference files:" $(date -u) >> hicpipe_${PREFIX}.progress.log

# copy reference to reference directory
cp ${REFERENCE} ${REF_DIR}/REF.fa

# samtools - index
samtools faidx ${REF_DIR}/REF.fa

# bwa - index for mapping
#bwa index -a bwtsw ${REF_DIR}/REF.fa

# minimap - index for mapping
minimap2 -t ${CPU} -d ${REF_DIR}/REF.mmi ${REF_DIR}/REF.fa

# make a genome file for bedtools
cut -f1,2 ${REF_DIR}/REF.fa.fai | sed 's/-/_/g'| sort -k2,2 -nr > ${REF_DIR}/REF.genome


# check if finished successfully
if [[ -s ${REF_DIR}/REF.fa ]] && [[ -s ${REF_DIR}/REF.fa.fai ]] && [[ -s ${REF_DIR}/REF.mmi ]] && [[ -s ${REF_DIR}/REF.genome ]]
then
     echo "Complete: building reference files:" $(date -u) >> hicpipe_${PREFIX}.progress.log
     touch ${REF_DIR}/COMPLETE_reference_files
else
     echo "PROBLEM: building reference files:" $(date -u) >> hicpipe_${PREFIX}.progress.log
     touch ${REF_DIR}/PROBLEM_reference_files
fi

}

export -f func_build_reference


#if [ -s ${REF_DIR}/REF.fa  ];
#then
#    echo "REF.fa is ready"
#else
#     # copy reference to reference directory
#    cp ${REFERENCE} ${REF_DIR}/REF.fa
#fi


#-------------------------------------------------------------------------------
### Mapping HiC illumina reads to reference
#-------------------------------------------------------------------------------
func_map_hic_reads(){

echo "Started: Mapping HiC illumina reads to reference:" $(date -u) >> hicpipe_${PREFIX}.progress.log

#- R1 FASTQ to BAM
#bwa mem -t ${CPU} ${REF_DIR}/REF.fa ${HIC_READS_R1} |\
#     samtools view --threads ${CPU} -Sb - \
#     > ${MAPPING_DIR}/${PREFIX}_1.bam

#- R2 FASTQ to BAM
#bwa mem -t ${CPU} ${REF_DIR}/REF.fa ${HIC_READS_R2} |\
#     samtools view --threads ${CPU} -Sb - \
#     > ${MAPPING_DIR}/${PREFIX}_2.bam

#- Filter 5' end (R1)
#samtools view -h ${MAPPING_DIR}/${PREFIX}_1.bam |\
#     perl $FILTER |\
#     samtools view -Sb - -o ${MAPPING_DIR}/${PREFIX}_1.filtered.bam

#- Filter 5' end (R2)
#samtools view -h ${MAPPING_DIR}/${PREFIX}_2.bam |\
#     perl $FILTER |\
#     samtools view -Sb - -o ${MAPPING_DIR}/${PREFIX}_2.filtered.bam


#- map R1 using minimap and filter 5' ends
minimap2 -x sr -a -t 20 ${REF_DIR}/REF.fa ${HIC_READS_R1} |\
     samtools view -h |\
     sed 's/\/1//1' |\
     perl $FILTER |\
     samtools view -Sb - -o ${MAPPING_DIR}/${PREFIX}_1.filtered.bam

#- map R2 using minimap and filter 5' ends
minimap2 -x sr -a -t 20 ${REF_DIR}/REF.fa ${HIC_READS_R2} |\
     samtools view -h |\
     sed 's/\/2//1' |\
     perl $FILTER |\
     samtools view -Sb - -o ${MAPPING_DIR}/${PREFIX}_2.filtered.bam

#- combine read pairs and filter on mapping quality
perl $COMBINER ${MAPPING_DIR}/${PREFIX}_1.filtered.bam ${MAPPING_DIR}/${PREFIX}_2.filtered.bam samtools $MAPQ_FILTER |\
     samtools view -bS -t ${REF_DIR}/REF.fa.fai - |\
     samtools sort --threads ${CPU} -o ${MAPPING_DIR}/${PREFIX}.combined_reads.bam -

#- add read group to BAM
java -Xmx4G -Djava.io.tmpdir=${TMP_DIR} -jar $PICARD AddOrReplaceReadGroups INPUT=${MAPPING_DIR}/${PREFIX}.combined_reads.bam OUTPUT=${MAPPING_DIR}/${PREFIX}.combined_reads.rg.bam ID=${FASTQ_BN} LB=${FASTQ_BN} SM=${PREFIX} PL=ILLUMINA PU=none

#- mark duplicates
java -Xmx20G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar $PICARD MarkDuplicates INPUT=${MAPPING_DIR}/${PREFIX}.combined_reads.rg.bam OUTPUT=${MAPPING_DIR}/${PREFIX}.combined_reads.md.bam METRICS_FILE=${MAPPING_DIR}/markdups.metrics.txt TMP_DIR=${TMP_DIR} ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

# index BAM
samtools index ${MAPPING_DIR}/${PREFIX}.combined_reads.md.bam

perl $STATS ${MAPPING_DIR}/${PREFIX}.combined_reads.md.bam > ${MAPPING_DIR}/${PREFIX}.combined_reads.md.bam.stats


# check if finished successfully
if [[ -s ${MAPPING_DIR}/${PREFIX}.combined_reads.bam ]] && [[ -s ${MAPPING_DIR}/${PREFIX}.combined_reads.md.bam ]] && [[ -s ${MAPPING_DIR}/${PREFIX}.combined_reads.md.bam.stats ]]
then
     echo "Complete: mapping HiC reads to reference:" $(date -u) >> hicpipe_${PREFIX}.progress.log
     touch ${MAPPING_DIR}/COMPLETE_hic_mapping

     #rm ${MAPPING_DIR}/${PREFIX}_1.bam ${MAPPING_DIR}/${PREFIX}_2.bam ${MAPPING_DIR}/${PREFIX}_1.filtered.bam ${MAPPING_DIR}/${PREFIX}_2.filtered.bam ${MAPPING_DIR}/${PREFIX}.combined_reads.bam ${MAPPING_DIR}/${PREFIX}.combined_reads.rg.bam
else
     echo "PROBLEM: mapping HiC reads to reference:" $(date -u) >> hicpipe_${PREFIX}.progress.log
     touch ${MAPPING_DIR}/PROBLEM_hic_mapping
fi

# clean up


}

export -f func_map_hic_reads


#-------------------------------------------------------------------------------
### Making pretext files for visualisation
#-------------------------------------------------------------------------------
func_pretext(){

echo "Started: Making pretext files for visualisation:" $(date -u) >> hicpipe_${PREFIX}.progress.log

# run PretextMap
samtools view -h ${MAPPING_DIR}/${PREFIX}.combined_reads.md.bam | PretextMap -o ${PRETEXT_DIR}/${PREFIX}.pretext.map --sortby length --sortorder descend --mapq ${MAPQ_FILTER} --highRes

samtools view --threads 4 -u -F0x400 ${MAPPING_DIR}/${PREFIX}.combined_reads.md.bam | bedtools bamtobed | sort -k4 --parallel=8 -S50G > ${PRETEXT_DIR}/${PREFIX}.pretext.part

# check if finished successfully
if [[ -s ${PRETEXT_DIR}/${PREFIX}.pretext.map ]] && [[ -s ${PRETEXT_DIR}/${PREFIX}.pretext.part ]]
then
     echo "Complete: Making pretext files for visualisation:" $(date -u) >> hicpipe_${PREFIX}.progress.log
     touch ${PRETEXT_DIR}/COMPLETE_pretext
else
     echo "PROBLEM: Making pretext files for visualisation:" $(date -u) >> hicpipe_${PREFIX}.progress.log
     touch ${PRETEXT_DIR}/PROBLEM_pretext
fi

}

export -f func_pretext



#-------------------------------------------------------------------------------
### Making HiGlass files for visualisation
#-------------------------------------------------------------------------------
func_higlass() {

echo "Started: Making HiGlass files for visualisation:" $(date -u) >> hicpipe_${PREFIX}.progress.log

samtools view --threads 4 -u -F0x400 ${MAPPING_DIR}/${PREFIX}.combined_reads.md.bam | bedtools bamtobed | sort -k4 --parallel=8 -S50G > ${HIGLASS_DIR}/merge.mkdup.bed

paste -d '\t' - - < ${HIGLASS_DIR}/merge.mkdup.bed | awk 'BEGIN {FS="\t"; OFS="\t"} {if ($1 > $7) {print substr($4,1,length($4)-2),$12,$7,$8,"16",$6,$1,$2,"8",$11,$5} else { print substr($4,1,length($4)-2),$6,$1,$2,"8",$12,$7,$8,"16",$5,$11} }' | tr '\-+' '10'  | sort --parallel=8 -S10G -k3,3d -k7,7d > ${HIGLASS_DIR}/pre.bed

# run cooler to make .cool and .mcool files for HiGlass
cooler cload pairs -0 -c1 3 -p1 4 -c2 7 -p2 8 ${REF_DIR}/REF.genome:1000 ${HIGLASS_DIR}/pre.bed ${HIGLASS_DIR}/${PREFIX}.higlass.cool

cooler zoomify -o ${HIGLASS_DIR}/${PREFIX}.higlass.mcool ${HIGLASS_DIR}/${PREFIX}.higlass.cool



# check if finished successfully
if [[ -s ${HIGLASS_DIR}/${PREFIX}.higlass.cool ]] && [[ -s ${HIGLASS_DIR}/${PREFIX}.higlass.mcool ]]
then
     echo "Complete: Making HiGlass files for visualisation:" $(date -u) >> hicpipe_${PREFIX}.progress.log
     touch ${HIGLASS_DIR}/COMPLETE_higlass
else
     echo "PROBLEM: Making HiGlass files for visualisation:" $(date -u) >> hicpipe_${PREFIX}.progress.log
     touch ${HIGLASS_DIR}/PROBLEM_higlass
fi

}

export -f func_higlass







#-------------------------------------------------------------------------------
### Making coverage files from pacbio reads to help with visualisation
#-------------------------------------------------------------------------------
func_coverage() {

echo "Started: Making coverage files from pacbio reads to help with visualisation:" $(date -u) >> hicpipe_${PREFIX}.progress.log

# minimap - mapping
minimap2 -x map-pb --MD -t ${CPU} -a ${REF_DIR}/REF.fa ${PACBIO_READS} | samtools view -Sb - > ${COVERAGE_DIR}/minimap.bam

# minimap - sort bam
samtools sort -T ${COVERAGE_DIR}/sort_tmp -o ${COVERAGE_DIR}/minimap.sorted.bam ${COVERAGE_DIR}/minimap.bam

# minimap - pri_bam
samtools view -b -hF 256 ${COVERAGE_DIR}/minimap.sorted.bam > ${COVERAGE_DIR}/${PREFIX}.minimap.sorted.pri.bam

# index bam
samtools index ${COVERAGE_DIR}/${PREFIX}.minimap.sorted.pri.bam

# bam to bed
bedtools bamtobed -i ${COVERAGE_DIR}/${PREFIX}.minimap.sorted.pri.bam | sort -k1,1 -T ${TMP_DIR} > ${COVERAGE_DIR}/${PREFIX}.geval.bed

# bed coverage
bedtools genomecov -bga -split -i ${COVERAGE_DIR}/${PREFIX}.geval.bed -g ${REF_DIR}/REF.genome | bedtools sort > ${COVERAGE_DIR}/${PREFIX}.coverage.bed

# make_bigwig
bedGraphToBigWig ${COVERAGE_DIR}/${PREFIX}.coverage.bed ${REF_DIR}/REF.genome ${COVERAGE_DIR}/${PREFIX}.coverage.bw


# check if finished successfully
if [[ -s ${COVERAGE_DIR}/${PREFIX}.coverage.bed ]] && [[ -s ${COVERAGE_DIR}/${PREFIX}.coverage.bw ]]
then
     echo "Complete: Making coverage files from pacbio reads to help with visualisation:" $(date -u) >> hicpipe_${PREFIX}.progress.log
     touch ${COVERAGE_DIR}/COMPLETE_coverage

     # clean up
     #rm ${COVERAGE_DIR}/minimap.bam ${COVERAGE_DIR}/minimap.sorted.bam
else
     echo "PROBLEM: Making coverage files from pacbio reads to help with visualisation:" $(date -u) >> hicpipe_${PREFIX}.progress.log
     touch ${COVERAGE_DIR}/PROBLEM_coverage
fi



}

export -f func_coverage




#-------------------------------------------------------------------------------
### Making gaps files from genome to help with visualisation
#-------------------------------------------------------------------------------
func_gaps() {

echo "Started: Making gaps files from genome to help with visualisation:" $(date -u) >> hicpipe_${PREFIX}.progress.log

seqtk cutN -n 1 -g ${REF_DIR}/REF.fa > ${GAPS_DIR}/${PREFIX}.gaps.bed

bedtools genomecov -bg -i ${GAPS_DIR}/${PREFIX}.gaps.bed -g ${REF_DIR}/REF.genome > ${GAPS_DIR}/${PREFIX}.gaps.bedgraph

# check if finished successfully
if [[ -s ${GAPS_DIR}/${PREFIX}.gaps.bed ]] && [[ -s ${GAPS_DIR}/${PREFIX}.gaps.bedgraph ]]
then
     echo "Complete: Making gaps files from genome to help with visualisation:" $(date -u) >> hicpipe_${PREFIX}.progress.log
     touch ${GAPS_DIR}/COMPLETE_gaps
else
     echo "PROBLEM: Making gaps files from genome to help with visualisation:" $(date -u) >> hicpipe_${PREFIX}.progress.log
     touch ${GAPS_DIR}/PROBLEM_gaps
fi

}

export -f func_gaps




#-------------------------------------------------------------------------------
### Make repeat files from genome to help with visualisation
#-------------------------------------------------------------------------------

func_repeats() {

echo "Started: Make repeat files from genome to help with visualisation:" $(date -u) >> hicpipe_${PREFIX}.progress.log

mkdir ${REPEATS_DIR}/red_out

/nfs/users/nfs_s/sd21/lustre118_link/software/RED_REPEAT_DETECTOR/redUnix64/Red -gnm ${REF_DIR} -rpt ${REPEATS_DIR}/red_out -frm 2

cat ${REPEATS_DIR}/red_out/*.rpt | sed 's/>//g' > ${REPEATS_DIR}/${PREFIX}.repeats.bed

bedtools makewindows -g ${REF_DIR}/REF.genome -w 10000 > ${REPEATS_DIR}/REF.genome.10000bp.bed

bedtools coverage -a ${REPEATS_DIR}/REF.genome.10000bp.bed -b ${REPEATS_DIR}/${PREFIX}.repeats.bed | awk '{print $1, $2, $3, $5}' OFS="\t" | bedtools sort > ${REPEATS_DIR}/${PREFIX}.repeats.bedgraph

bedGraphToBigWig ${REPEATS_DIR}/${PREFIX}.repeats.bedgraph ${REF_DIR}/REF.genome ${REPEATS_DIR}/${PREFIX}.repeats.bw


# check if finished successfully
if [[ -s ${REPEATS_DIR}/${PREFIX}.repeats.bedgraph ]] && [[ -s ${REPEATS_DIR}/${PREFIX}.repeats.bw ]]
then
     echo "Complete: Make repeat files from genome to help with visualisation:" $(date -u) >> hicpipe_${PREFIX}.progress.log
     touch ${REPEATS_DIR}/COMPLETE_repeats
else
     echo "PROBLEM: Make repeat files from genome to help with visualisation:" $(date -u) >> hicpipe_${PREFIX}.progress.log
     touch ${REPEATS_DIR}/PROBLEM_repeats
fi

}

export -f func_repeats


#-------------------------------------------------------------------------------
### Find telomeres for visualisation
#-------------------------------------------------------------------------------

func_telomeres() {

echo "Started: Find telomeres for visualisation:" $(date -u) >> hicpipe_${PREFIX}.progress.log

fastaq search_for_seq ${REF_DIR}/REF.fa ${TELOMERES_DIR}/telomere_dimer.pos ${TELOMERE_SEQ}${TELOMERE_SEQ}

cat ${TELOMERES_DIR}/telomere_dimer.pos | awk '{print $1,$2,$2+1}' OFS="\t" > ${TELOMERES_DIR}/telomere_dimer.bed

bedtools makewindows -g ${REF_DIR}/REF.genome -w 10000 > ${TELOMERES_DIR}/REF.genome.10000bp.bed

bedtools coverage -b ${TELOMERES_DIR}/telomere_dimer.bed -a ${TELOMERES_DIR}/REF.genome.10000bp.bed  | awk '{print $1, $2, $3, $5}' OFS="\t" | bedtools sort > ${TELOMERES_DIR}/${PREFIX}.telomere_coverage.bedgraph

bedGraphToBigWig ${TELOMERES_DIR}/${PREFIX}.telomere_coverage.bedgraph ${REF_DIR}/REF.genome ${TELOMERES_DIR}/${PREFIX}.telomere_coverage.bw

# check if finished successfully
if [[ -s ${TELOMERES_DIR}/${PREFIX}.telomere_coverage.bedgraph ]] && [[ -s ${TELOMERES_DIR}/${PREFIX}.telomere_coverage.bw ]]
then
     echo "Complete: Find telomeres for visualisation:" $(date -u) >> hicpipe_${PREFIX}.progress.log
     touch ${TELOMERES_DIR}/COMPLETE_telomeres
else
     echo "PROBLEM: Find telomeres for visualisation:" $(date -u) >> hicpipe_${PREFIX}.progress.log
     touch ${TELOMERES_DIR}/PROBLEM_telomeres
fi

}

export -f func_telomeres



#-------------------------------------------------------------------------------
### Making .hic files for juicer
#-------------------------------------------------------------------------------

func_juicer() {

echo "Started: Making .hic files for juicer:" $(date -u) >> hicpipe_${PREFIX}.progress.log

samtools view --threads 4 -u -F0x400 ${MAPPING_DIR}/${PREFIX}.combined_reads.md.bam | bedtools bamtobed | sort -k4 --parallel=8 -S50G > ${JUICER_DIR}/merge.mkdup.bed

paste -d '\t' - - < ${JUICER_DIR}/merge.mkdup.bed | awk 'BEGIN {FS="\t"; OFS="\t"} {if ($1 > $7) {print substr($4,1,length($4)-2),$12,$7,$8,"16",$6,$1,$2,"8",$11,$5} else { print substr($4,1,length($4)-2),$6,$1,$2,"8",$12,$7,$8,"16",$5,$11} }' | tr '\-+' '10'  | sort --parallel=8 -S10G -k3,3d -k7,7d > ${JUICER_DIR}/pre.bed

java -jar ${JUICER} pre -q ${MAPQ_FILTER} ${JUICER_DIR}/pre.bed ${JUICER_DIR}/${PREFIX}.juicer.hic ${REF_DIR}/REF.genome

awk -f /nfs/users/nfs_s/sd21/lustre118_link/software/HIC/3d-dna/utils/generate-assembly-file-from-fasta.awk ${REF_DIR}/REF.fa > ${JUICER_DIR}/${PREFIX}.juicer.assembly

# check if finished successfully
if [[ -s ${JUICER_DIR}/pre.bed ]] && [[ -s ${JUICER_DIR}/${PREFIX}.juicer.hic ]]
then
     echo "Complete: Find telomeres for visualisation:" $(date -u) >> hicpipe_${PREFIX}.progress.log
     touch ${JUICER_DIR}/COMPLETE_juicer
else
     echo "PROBLEM: Find telomeres for visualisation:" $(date -u) >> hicpipe_${PREFIX}.progress.log
     touch ${JUICER_DIR}/PROBLEM_juicer
fi

}

export -f func_juicer





#-------------------------------------------------------------------------------
### Finishing up
#-------------------------------------------------------------------------------

func_finish() {

     # check if finished successfully
     if [[ -f ${REF_DIR}/COMPLETE_reference_files ]] && [[ -f ${MAPPING_DIR}/COMPLETE_hic_mapping ]] && [[ -f ${PRETEXT_DIR}/COMPLETE_pretext ]] && [[ -f ${HIGLASS_DIR}/COMPLETE_higlass ]] && [[ -f ${COVERAGE_DIR}/COMPLETE_coverage ]] && [[ -f ${GAPS_DIR}/COMPLETE_gaps ]] && [[ -f ${REPEATS_DIR}/COMPLETE_repeats ]] && [[ -f ${JUICER_DIR}/COMPLETE_juicer ]] && [[ -f ${TELOMERES_DIR}/COMPLETE_telomeres ]]
     then
          echo "HiC mapping pipeline is COMPLETE!:" $(date -u) >> hicpipe_${PREFIX}.progress.log
          touch ${PWD}/hicpipe_${PREFIX}_COMPLETE
     else
          echo "PROBLEM: HiC mapping pipeline failed along the way. Revise and resubmit:" $(date -u) >> hicpipe_${PREFIX}.progress.log
          touch ${PWD}/hicpipe_${PREFIX}_ERROR
     fi

}

export -f func_finish



#-------------------------------------------------------------------------------
# running the pipeline
#-------------------------------------------------------------------------------

bsub -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>7000] rusage[mem=7000]" -n ${CPU} -M7000 -o ${LOG_FILES}/01_build_reference.o -e ${LOG_FILES}/01_build_reference.e -J 01_build_reference_${PREFIX} func_build_reference

bsub -w "01_build_reference_${PREFIX}" -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>30000] rusage[mem=30000]" -q long -n 20 -M30000 -o ${LOG_FILES}/02_map_hic.o -e ${LOG_FILES}/02_map_hic.e -J 02_map_hic_${PREFIX} func_map_hic_reads

bsub -w "02_map_hic_${PREFIX}" -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>60000] rusage[mem=60000]" -n 4 -M60000 -o ${LOG_FILES}/03_pretext.o -e ${LOG_FILES}/03_pretext.e -J 03_pretext_${PREFIX} func_pretext

bsub -w "02_map_hic_${PREFIX}" -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>60000] rusage[mem=60000]" -n 4 -M60000 -o ${LOG_FILES}/04_higlass.o -e ${LOG_FILES}/04_higlass.e -J 04_higlass_${PREFIX} func_higlass

bsub -w "01_build_reference_${PREFIX}" -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>15000] rusage[mem=15000]" -q long -n 20 -M15000 -o ${LOG_FILES}/05_coverage.o -e ${LOG_FILES}/05_coverage.e -J 05_coverage_${PREFIX} func_coverage

bsub -w "01_build_reference_${PREFIX}" -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>1000] rusage[mem=1000]" -M1000 -o ${LOG_FILES}/06_gaps.o -e ${LOG_FILES}/06_gaps.e -J 06_gaps_${PREFIX} func_gaps

bsub -w "01_build_reference_${PREFIX}" -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>5000] rusage[mem=5000]" -M5000 -o ${LOG_FILES}/07_repeats.o -e ${LOG_FILES}/07_repeats.e -J 07_repeats_${PREFIX} func_repeats

bsub -w "01_build_reference_${PREFIX}" -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>1000] rusage[mem=1000]" -M1000 -o ${LOG_FILES}/08_telomeres.o -e ${LOG_FILES}/08_telomeres.e -J 08_telomeres_${PREFIX} func_telomeres

bsub -w "02_map_hic_${PREFIX}" -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>60000] rusage[mem=60000]" -n 4 -M60000 -o ${LOG_FILES}/09_juicer.o -e ${LOG_FILES}/09_juicer.e -J 09_juicer_${PREFIX} func_juicer

# wrap up
bsub -w "01_build_reference_${PREFIX} && 02_map_hic_${PREFIX} && 03_pretext_${PREFIX} && 04_higlass_${PREFIX} && 05_coverage_${PREFIX} && 06_gaps_${PREFIX} && 07_repeats_${PREFIX} && 08_telomeres_${PREFIX} && 09_juicer_${PREFIX}" -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>1000] rusage[mem=1000]" -q small -M1000 -o ${LOG_FILES}/finish.o -e ${LOG_FILES}/finish.e -J finish_${PREFIX} func_finish
