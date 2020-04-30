# S. mansoni 10X single cell RNASeq analysis
# -- initial prep of data for Carmen - reference prep, mapping



# working directory
cd /nfs/users/nfs_s/sd21/lustre118_link/sm/10X/CELL_RANGER



# Step 1 -  make reference
bsub.py --queue yesterday --threads 7 20 01_10X_CR_mkref "~sd21/lustre118_link/software/10X_GENOMICS/cellranger-2.1.0/cellranger mkref --genome=Sm_v7.1 --fasta=Smansoni_v7.fa --genes=Sm_v7.1.gtf --nthreads=7"



# Step 2 - get Carmens data
#Lane	dataset
#24056_3	FUGI_R_D7119553
#24056_4	FUGI_R_D7119554
#24346_5	FUGI_R_D7159524
#24346_6	FUGI_R_D7159525


mkdir RAW_FASTQ
cd RAW_FASTQ

kinit
icd /seq/24056/cellranger/fastq
ils | grep "24056_[34]" | while read -r data; do iget /seq/24056/cellranger/fastq/$data ./; done &
icd /seq/24346/cellranger/fastq
ils | grep "24346_[56]" | while read -r data; do iget /seq/24346/cellranger/fastq/$data ./; done &


# fastqs were not availabel for the latest batch of data, so had to download crams, and convert to fastqs
icd /seq/25786
ils  | grep cram$ | grep -v phi | grep 25786_[1234] | grep -v "888" | grep -v 0 | while read -r data; do iget /seq/25786/${data} ./ ; done &
ils  | grep cram$ | grep -v phi | grep 25786_[4] | grep -v "888" | grep -v 0 | while read -r data; do iget /seq/25786/${data} ./ ; done &

#--- made a "samples.list" file containing sample metadata, which could be used to fill in while loop
# 25786_1#1.cram	25786_1	FUGI_R_D7462031	257	TCTTAAAG	Adult_F
# 25786_1#2.cram	25786_1	FUGI_R_D7462031	258	CGAGGCTC	Adult_F
# 25786_1#3.cram	25786_1	FUGI_R_D7462031	259	GTCCTTCT	Adult_F
# 25786_1#4.cram	25786_1	FUGI_R_D7462031	260	AAGACGGA	Adult_F
# 25786_2#1.cram	25786_2	FUGI_R_D7462032	261	CTGTAACT	Adult_F
# 25786_2#2.cram	25786_2	FUGI_R_D7462032	262	TCTAGCGA	Adult_F
# 25786_2#3.cram	25786_2	FUGI_R_D7462032	263	AGAGTGTG	Adult_F
# 25786_2#4.cram	25786_2	FUGI_R_D7462032	264	GACCCTAC	Adult_F
# 25786_3#1.cram	25786_3	FUGI_R_D7462033	265	GCGCAGAA	Adult_F
# 25786_3#2.cram	25786_3	FUGI_R_D7462033	266	ATCTTACC	Adult_F
# 25786_3#3.cram	25786_3	FUGI_R_D7462033	267	TATGGTGT	Adult_F
# 25786_3#4.cram	25786_3	FUGI_R_D7462033	268	CGAACCTG	Adult_F
# 25786_4#1.cram	25786_4	FUGI_R_D7462034	269	AGGAGATG	Adult_F
# 25786_4#2.cram	25786_4	FUGI_R_D7462034	270	GATGTGGT	Adult_F
# 25786_4#3.cram	25786_4	FUGI_R_D7462034	271	CTACATCC	Adult_F
# 25786_4#4.cram	25786_4	FUGI_R_D7462034	272	TCCTCCAA	Adult_F
# 25872_2#1.cram	25872_2	FUGI_R_D7465035	341	GCGAGAGT 	Sporocysts
# 25872_2#2.cram	25872_2 FUGI_R_D7465035	342	TACGTTCA	Sporocysts
# 25872_2#3.cram	25872_2	FUGI_R_D7465035	343	AGTCCCAC	Sporocysts
# 25872_2#4.cram	25872_2	FUGI_R_D7465035	344	CTATAGTG	Sporocysts
# 25872_3#1.cram	25872_3 FUGI_R_D7465036	345 TTATCGTT	Sporocysts
# 25872_3#2.cram	25872_3 FUGI_R_D7465036	346 AGCAGAGC	Sporocysts
# 25872_3#3.cram	25872_3 FUGI_R_D7465036	347 CATCTCCA	Sporocysts
# 25872_3#4.cram	25872_3 FUGI_R_D7465036	348 GCGGATAG	Sporocysts
# 25872_4#1.cram	25872_4	FUGI_R_D7465037	349 GGCGAGTA	Sporocysts
# 25872_4#2.cram	25872_4	FUGI_R_D7465037	350 ACTTCTAT	Sporocysts
# 25872_4#3.cram	25872_4	FUGI_R_D7465037	351 CAAATACG	Sporocysts
# 25872_4#4.cram	25872_4	FUGI_R_D7465037	352 TTGCGCGC	Sporocysts
# 25892_2#1.cram	25892_2 FUGI_R_D7465034	337 AAGCGCTG	Sporocysts
# 25892_2#2.cram	25892_2 FUGI_R_D7465034	338 CGTTTGAT	Sporocysts
# 25892_2#3.cram	25892_2 FUGI_R_D7465034	339 GTAGCACA	Sporocysts
# 25892_2#4.cram	25892_2 FUGI_R_D7465034	340 TCCAATGC	Sporocysts











while read cram lane id tag_id tag; do
/software/solexa/pkg/samtools/1.5/bin/samtools fastq \
-1 ${lane}_${tag}_S${tag_id}_L00${lane#*_}_R1_001.fastq.gz \
-2 ${lane}_${tag}_S${tag_id}_L00${lane#*_}_R2_001.fastq.gz \
--i1 ${lane}_${tag}_S${tag_id}_L00${lane#*_}_I1_001.fastq.gz \
--index-format "i8" -n -i $cram; done < samples.list &



cd ../



# Step 3 - mapping data to reference

# run count
bsub.py --queue long --threads 16 30 02_10X_CR_count_Lane3 "/nfs/users/nfs_s/sd21/lustre118_link/software/10X_GENOMICS/cellranger-2.1.0/cellranger count \
--id=FUGI_R_D7119553 \
--transcriptome=/nfs/users/nfs_s/sd21/lustre118_link/sm/10X/CELL_RANGER/Sm_v7.1 \
--fastqs=/nfs/users/nfs_s/sd21/lustre118_link/sm/10X/CELL_RANGER/RAW_FASTQ/24056_3 \
--sample=24056_3_AAGGGTGA,24056_3_CGTTAATC,24056_3_GCCACGCT,24056_3_TTACTCAG \
--localcores=16 \
--localmem=30"



bsub.py --queue long --threads 16 30 02_10X_CR_count_Lane4 "/nfs/users/nfs_s/sd21/lustre118_link/software/10X_GENOMICS/cellranger-2.1.0/cellranger count \
--id=FUGI_R_D7119554 \
--transcriptome=/nfs/users/nfs_s/sd21/lustre118_link/sm/10X/CELL_RANGER/Sm_v7.1 \
--fastqs=/nfs/users/nfs_s/sd21/lustre118_link/sm/10X/CELL_RANGER/RAW_FASTQ/24056_4 \
--sample=24056_4_ATTACTTC,24056_4_CAGCGGAA,24056_4_GCATTCGG,24056_4_TGCGAACT \
--localcores=16 \
--localmem=30"



bsub.py --queue long --threads 16 30 02_10X_CR_count_Lane5 "/nfs/users/nfs_s/sd21/lustre118_link/software/10X_GENOMICS/cellranger-2.1.0/cellranger count \
--id=FUGI_R_D7159524 \
--transcriptome=/nfs/users/nfs_s/sd21/lustre118_link/sm/10X/CELL_RANGER/Sm_v7.1 \
--fastqs=/nfs/users/nfs_s/sd21/lustre118_link/sm/10X/CELL_RANGER/RAW_FASTQ/24346_5 \
--sample=24346_5_ACATTACT,24346_5_CAGCCCAC,24346_5_GGCAATGG,24346_5_TTTGGGTA \
--localcores=16 \
--localmem=30"



bsub.py --queue long --threads 16 30 02_10X_CR_count_Lane6 "/nfs/users/nfs_s/sd21/lustre118_link/software/10X_GENOMICS/cellranger-2.1.0/cellranger count \
--id=FUGI_R_D7159525 \
--transcriptome=/nfs/users/nfs_s/sd21/lustre118_link/sm/10X/CELL_RANGER/Sm_v7.1 \
--fastqs=/nfs/users/nfs_s/sd21/lustre118_link/sm/10X/CELL_RANGER/RAW_FASTQ/24346_6 \
--sample=24346_6_AGGTATTG,24346_6_CTCCTAGT,24346_6_GATGCCAA,24346_6_TCAAGGCC \
--localcores=16 \
--localmem=30"


bsub.py --queue long --threads 16 30 02_10X_CR_count_25786_1_1 "/nfs/users/nfs_s/sd21/lustre118_link/software/10X_GENOMICS/cellranger-2.1.0/cellranger count \
--id=FUGI_R_D7462031 \
--transcriptome=/nfs/users/nfs_s/sd21/lustre118_link/sm/10X/CELL_RANGER/Sm_v7.1 \
--fastqs=/nfs/users/nfs_s/sd21/lustre118_link/sm/10X/CELL_RANGER/RAW_FASTQ/25786_1 \
--sample=25786_1_TCTTAAAG,25786_1_CGAGGCTC,25786_1_GTCCTTCT,25786_1_AAGACGGA \
--localcores=16 \
--localmem=30"

bsub.py --queue long --threads 16 30 02_10X_CR_count_25786_2 "/nfs/users/nfs_s/sd21/lustre118_link/software/10X_GENOMICS/cellranger-2.1.0/cellranger count \
--id=FUGI_R_D7462032 \
--transcriptome=/nfs/users/nfs_s/sd21/lustre118_link/sm/10X/CELL_RANGER/Sm_v7.1 \
--fastqs=/nfs/users/nfs_s/sd21/lustre118_link/sm/10X/CELL_RANGER/RAW_FASTQ/25786_2 \
--sample=25786_2_CTGTAACT,25786_2_TCTAGCGA,25786_2_AGAGTGTG,25786_2_GACCCTAC \
--localcores=16 \
--localmem=30"

bsub.py --queue long --threads 16 30 02_10X_CR_count_25786_3 "/nfs/users/nfs_s/sd21/lustre118_link/software/10X_GENOMICS/cellranger-2.1.0/cellranger count \
--id=FUGI_R_D7462033 \
--transcriptome=/nfs/users/nfs_s/sd21/lustre118_link/sm/10X/CELL_RANGER/Sm_v7.1 \
--fastqs=/nfs/users/nfs_s/sd21/lustre118_link/sm/10X/CELL_RANGER/RAW_FASTQ/25786_3 \
--sample=25786_3_GCGCAGAA,25786_3_ATCTTACC,25786_3_TATGGTGT,25786_3_CGAACCTG \
--localcores=16 \
--localmem=30"


bsub.py --queue long --threads 16 30 02_10X_CR_count_25786_4 "/nfs/users/nfs_s/sd21/lustre118_link/software/10X_GENOMICS/cellranger-2.1.0/cellranger count \
--id=FUGI_R_D7462034 \
--transcriptome=/nfs/users/nfs_s/sd21/lustre118_link/sm/10X/CELL_RANGER/Sm_v7.1 \
--fastqs=/nfs/users/nfs_s/sd21/lustre118_link/sm/10X/CELL_RANGER/RAW_FASTQ/25786_4 \
--sample=25786_4_AGGAGATG,25786_4_GATGTGGT,25786_4_CTACATCC,25786_4_TCCTCCAA \
--localcores=16 \
--localmem=30"




bsub.py --queue long --threads 16 30 02_10X_CR_count_25872_2 "/nfs/users/nfs_s/sd21/lustre118_link/software/10X_GENOMICS/cellranger-2.1.0/cellranger count \
--id=FUGI_R_D7465035 \
--transcriptome=/nfs/users/nfs_s/sd21/lustre118_link/sm/10X/CELL_RANGER/Sm_v7.1 \
--fastqs=/nfs/users/nfs_s/sd21/lustre118_link/sm/10X/CELL_RANGER/RAW_FASTQ/25872_2 \
--sample=25872_2_GCGAGAGT,25872_2_TACGTTCA,25872_2_AGTCCCAC,25872_2_CTATAGTG \
--localcores=16 \
--localmem=30"



bsub.py --queue long --threads 16 30 02_10X_CR_count_25872_3 "/nfs/users/nfs_s/sd21/lustre118_link/software/10X_GENOMICS/cellranger-2.1.0/cellranger count \
--id=FUGI_R_D7465036 \
--transcriptome=/nfs/users/nfs_s/sd21/lustre118_link/sm/10X/CELL_RANGER/Sm_v7.1 \
--fastqs=/nfs/users/nfs_s/sd21/lustre118_link/sm/10X/CELL_RANGER/RAW_FASTQ/25872_3 \
--sample=25872_3_TTATCGTT,25872_3_AGCAGAGC,25872_3_CATCTCCA,25872_3_GCGGATAG \
--localcores=16 \
--localmem=30"

bsub.py --queue long --threads 16 30 02_10X_CR_count_25872_4 "/nfs/users/nfs_s/sd21/lustre118_link/software/10X_GENOMICS/cellranger-2.1.0/cellranger count \
--id=FUGI_R_D7465037 \
--transcriptome=/nfs/users/nfs_s/sd21/lustre118_link/sm/10X/CELL_RANGER/Sm_v7.1 \
--fastqs=/nfs/users/nfs_s/sd21/lustre118_link/sm/10X/CELL_RANGER/RAW_FASTQ/25872_4 \
--sample=25872_4_GGCGAGTA,25872_4_ACTTCTAT,25872_4_CAAATACG,25872_4_TTGCGCGC \
--localcores=16 \
--localmem=30"



bsub.py --queue long --threads 16 30 02_10X_CR_count_25872_2 "/nfs/users/nfs_s/sd21/lustre118_link/software/10X_GENOMICS/cellranger-2.1.0/cellranger count \
--id=FUGI_R_D7465034 \
--transcriptome=/nfs/users/nfs_s/sd21/lustre118_link/sm/10X/CELL_RANGER/Sm_v7.1 \
--fastqs=/nfs/users/nfs_s/sd21/lustre118_link/sm/10X/CELL_RANGER/RAW_FASTQ/25872_2 \
--sample=25872_2_CTATAGTG,25872_2_AGTCCCAC,25872_2_TACGTTCA,25872_2_GCGAGAGT \
--localcores=16 \
--localmem=30"
