#!/bin/bash

#-------------------------------------------------------------------------------
# run_arrow_pacbio_genome_polish.sh
#-------------------------------------------------------------------------------

# stephen doyle
# Jan 2023


export PREFIX=tc_ASSEMBLY_230308
export REFERENCE=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ASSEMBLY/POST_CANU_IMPROVEMENT/PRE_POLISH_COMPLETE/ASSEMBLY_230308.fa
export PB_READ_DATA_FOFN=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/POLISH/subread_bams.fofn

# load required modules
module load pacbio-smrttools/7.0.1.66768
module load samtools/1.14--hb421002_0

## file locations
export LOG_FILES="$PWD/pb_arrow_polish_${PREFIX}_out/LOG_FILES"
export REFERENCE_FILES="$PWD/pb_arrow_polish_${PREFIX}_out/REFERENCE_FILES"
export MAPPING="$PWD/pb_arrow_polish_${PREFIX}_out/MAPPING"
export SPLIT_FASTAS="$PWD/pb_arrow_polish_${PREFIX}_out/SPLIT_FASTAS"
export SPLIT_BAMS="$PWD/pb_arrow_polish_${PREFIX}_out/SPLIT_BAMS"
export SPLIT_POLISHED_FASTAS="$PWD/pb_arrow_polish_${PREFIX}_out/SPLIT_POLISHED_FASTAS"





### Check output directories exist & create them as needed
[ -d ${LOG_FILES} ] || mkdir -p ${LOG_FILES}
[ -d ${REFERENCE_FILES} ] || mkdir -p ${REFERENCE_FILES}
[ -d ${MAPPING} ] || mkdir -p ${MAPPING}
[ -d ${SPLIT_FASTAS} ] || mkdir -p ${SPLIT_FASTAS}
[ -d ${SPLIT_BAMS} ] || mkdir -p ${SPLIT_BAMS}
[ -d ${SPLIT_POLISHED_FASTAS} ] || mkdir -p ${SPLIT_POLISHED_FASTAS}



# save current script in run folder to reproduce the exact output
cp ${PWD}/run_arrow_pacbio_genome_polish.sh ${PWD}/pb_arrow_polish_${PREFIX}_out/commands.$(date -Iminutes).txt

# make a progress file
#> ${PWD}/pb_arrow_polish_${PREFIX}.progress.log





#-------------------------------------------------------------------------------
### 01. Prepare reference files
#-------------------------------------------------------------------------------
# pbmm2: https://github.com/PacificBiosciences/pbmm2/ 

func_build_reference() {

if [ -f "${REFERENCE_FILES}/REF.fa" ]; then
        echo "Reference is already setup. Moving on."
        exit 0
    else

    cp ${REFERENCE} ${REFERENCE_FILES}/REF.fa

    pbmm2 index ${REFERENCE_FILES}/REF.fa ${REFERENCE_FILES}/REF.mmi

    samtools faidx ${REFERENCE_FILES}/REF.fa > ${REFERENCE_FILES}/REF.fa.fai 

    cat ${REFERENCE_FILES}/REF.fa.fai | cut -f1 > ${REFERENCE_FILES}/sequences.list

    while read SEQUENCE; do 
       echo "samtools faidx ${REFERENCE_FILES}/REF.fa ${SEQUENCE} > ${SPLIT_FASTAS}/${SEQUENCE}.fa";
    done < ${REFERENCE_FILES}/sequences.list | parallel -j20
fi
}

export -f func_build_reference



#-------------------------------------------------------------------------------
### 02. Map pacbio reads to reference
#-------------------------------------------------------------------------------
func_map_reads () {
    # check if bam file exits - if yes, then exit. Else, run mapping
    if [ -s "${MAPPING}/REF.bam" ]; then
            echo "Bam file is already setup. Moving on."
            exit 0
        else

        pbmm2 align ${REFERENCE_FILES}/REF.fa ${PB_READ_DATA_FOFN} ${MAPPING}/REF.bam --preset SUBREAD --sort -j 20 -J 20

        pbindex ${MAPPING}/REF.bam

        samtools index ${MAPPING}/REF.bam
    fi
}

export -f func_map_reads



#-------------------------------------------------------------------------------
### 03. Split bams per reference
#-------------------------------------------------------------------------------
func_split_bams () {

while read SEQUENCE; do 
    echo "samtools view -b --reference ${SPLIT_FASTAS}/${SEQUENCE}.fa -o ${SPLIT_BAMS}/${SEQUENCE}.bam ${MAPPING}/REF.bam ${SEQUENCE}"; 
    done < ${REFERENCE_FILES}/sequences.list  | parallel -j20

while read SEQUENCE; do 
    echo "pbindex ${SPLIT_BAMS}/${SEQUENCE}.bam";
    done < ${REFERENCE_FILES}/sequences.list | parallel -j20

while read SEQUENCE; do 
    echo "samtools index ${SPLIT_BAMS}/${SEQUENCE}.bam";
    done < ${REFERENCE_FILES}/sequences.list | parallel -j20

}

export -f func_split_bams


#-------------------------------------------------------------------------------
### 04. Polish split bams and fasta sequences using arrow
#-------------------------------------------------------------------------------
func_arrow_polish () {

while read SEQUENCE; do 
    echo "arrow --algorithm arrow --threaded -j 20 --referenceFilename ${SPLIT_FASTAS}/${SEQUENCE}.fa -o ${SPLIT_POLISHED_FASTAS}/${SEQUENCE}.arrow.fa ${SPLIT_BAMS}/${SEQUENCE}.bam";
    done < ${REFERENCE_FILES}/sequences.list | parallel -j20

ls -v ${SPLIT_POLISHED_FASTAS}/*.arrow.fa | xargs cat > $PWD/pb_arrow_polish_${PREFIX}_out/${PREFIX}.arrow.fa

} 

export -f func_arrow_polish



#-------------------------------------------------------------------------------
# running the pipeline
#-------------------------------------------------------------------------------

# func_build_reference
bsub -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>7000] rusage[mem=7000]" -M7000 -o ${LOG_FILES}/01_build_reference.o -e ${LOG_FILES}/01_build_reference.e -J 01_build_reference_${PREFIX} func_build_reference

# func_map_reads
bsub -w "01_build_reference_${PREFIX}" -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>50000] rusage[mem=50000]" -q long -M50000 -n20 -o ${LOG_FILES}/02_map_reads.o -e ${LOG_FILES}/02_map_reads.e -J 02_map_reads_${PREFIX} func_map_reads

# func_split_bams
bsub -w "02_map_reads_${PREFIX}" -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>10000] rusage[mem=10000]" -q long -M10000 -o ${LOG_FILES}/03_split_bams.o -e ${LOG_FILES}/03_split_bams.e -J 03_split_bams_${PREFIX} func_split_bams

# func_arrow_polish
bsub -w "03_split_bams_${PREFIX}" -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>200000] rusage[mem=200000]" -q hugemem -n 20 -M200000 -o ${LOG_FILES}/04_arrow_polish.o -e ${LOG_FILES}/04_arrow_polish.e -J 04_arrow_polish_${PREFIX} func_arrow_polish


