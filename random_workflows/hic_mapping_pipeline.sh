#! /bin/bash



# adapter from https://github.com/ArimaGenomics/mapping_pipeline 

## input required
FQ_BASENAME='basename_of_fastq_files'
EXP_NAME='overall_exp_name'
PREFIX='bwa_index_name'
IN_DIR='/path/to/gzipped/fastq/files'
REF='/path/to/reference_sequences/reference_sequeneces.fa'
CPU=12

## file locations
REF_DIR='$PWD/ref'
RAW_DIR='$PWD/bams'
FILT_DIR='$PWD/filtered_bams'
PAIR_DIR='$PWD/paired_bams'
REP_DIR='$PWD/deduplicated_files'
TMP_DIR='$PWD/tmp'


## tool locations
FILTER='/nfs/users/nfs_s/sd21/lustre118_link/software/HIC/mapping_pipeline/filter_five_end.pl'
COMBINER='/nfs/users/nfs_s/sd21/lustre118_link/software/HIC/mapping_pipeline/two_read_bam_combiner.pl'
STATS='/nfs/users/nfs_s/sd21/lustre118_link/software/HIC/mapping_pipeline/get_stats.pl'
PICARD='/nfs/users/nfs_s/sd21/lustre118_link/software/picard-tools-2.5.0/picard.jar'


## Other
REP_LABEL=$EXP_NAME\_rep1
MAPQ_FILTER=10

---

echo "### Step 0: Check output directories exist & create them as needed"
[ -d $REF_DIR ] || mkdir -p $REF_DIR
[ -d $RAW_DIR ] || mkdir -p $RAW_DIR
[ -d $FILT_DIR ] || mkdir -p $FILT_DIR
[ -d $TMP_DIR ] || mkdir -p $TMP_DIR
[ -d $PAIR_DIR ] || mkdir -p $PAIR_DIR
[ -d $REP_DIR ] || mkdir -p $REP_DIR
[ -d $MERGE_DIR ] || mkdir -p $MERGE_DIR



echo "### Step 0.A: Get reference"
cp $REF $REF_DIR/REF.fasta
samtools faidx $REF_DIR/REF.fasta

echo "### Step 0.B: Index reference"
bwa index -a bwtsw -p $PREFIX $REF_DIR/REF.fasta


echo "### Step 1.A: FASTQ to BAM (1st)"
bwa mem -t $CPU $REF_DIR/REF.fasta $IN_DIR/$FQ_BASENAME\_1.fastq.gz |\
     $SAMTOOLS view -@ $CPU -Sb - \
     > $RAW_DIR/$FQ_BASENAME\_1.bam


echo "### Step 1.B: FASTQ to BAM (2nd)"
bwa mem -t $CPU $REF_DIR/REF.fasta $IN_DIR/$FQ_BASENAME\_2.fastq.gz |\
     $SAMTOOLS view -@ $CPU -Sb - \
     > $RAW_DIR/$FQ_BASENAME\_2.bam



echo "### Step 2.A: Filter 5' end (1st)"
samtools view -h $RAW_DIR/$FQ_BASENAME\_1.bam |\
     perl $FILTER |\
     samtools view -Sb - \
     > $FILT_DIR/$FQ_BASENAME\_1.bam

echo "### Step 2.B: Filter 5' end (2nd)"
samtools view -h $RAW_DIR/$FQ_BASENAME\_2.bam |\
     perl $FILTER |\
     samtools view -Sb - \
     > $FILT_DIR/$FQ_BASENAME\_2.bam




echo "### Step 3A: Pair reads & mapping quality filter"
perl $COMBINER $FILT_DIR/$FQ_BASENAME\_1.bam $FILT_DIR/$FQ_BASENAME\_2.bam samtools $MAPQ_FILTER |\
     samtools view -bS -t $REF_DIR/REF.fasta.fai - |\
     samtools sort -@ $CPU -o $TMP_DIR/$FQ_BASENAME.bam -

echo "### Step 3.B: Add read group"
java -Xmx4G -Djava.io.tmpdir=temp/ -jar $PICARD AddOrReplaceReadGroups INPUT=$TMP_DIR/$FQ_BASENAME.bam OUTPUT=$PAIR_DIR/$FQ_BASENAME.bam ID=$FQ_BASENAME LB=$FQ_BASENAME SM=$EXP_NAME PL=ILLUMINA PU=none


echo "### Step 4: Mark duplicates"
java -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar $PICARD MarkDuplicates INPUT=$PAIR_DIR/$FQ_BASENAME.bam OUTPUT=$REP_DIR/$REP_LABEL.bam METRICS_FILE=$REP_DIR/metrics.$REP_LABEL.txt TMP_DIR=$TMP_DIR ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

samtools index $REP_DIR/$REP_LABEL.bam

perl $STATS $REP_DIR/$REP_LABEL.bam > $REP_DIR/$REP_LABEL.bam.stats

echo "Finished Mapping Pipeline through Duplicate Removal"
