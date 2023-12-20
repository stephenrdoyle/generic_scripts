

# get reference files for mouse

wget  wget https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

wget https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz

for i in *gz; do gunzip $i; done

# Step 1 -  make reference
 module load cellranger-7/7.0.1

bsub.py --queue yesterday --threads 7 20 01_10X_CR_mkref "cellranger mkref --genome=GRCm39 --fasta=Mus_musculus.GRCm39.dna.primary_assembly.fa --genes=Mus_musculus.GRCm39.110.gtf --nthreads=7"



# data 
cat data.txt
4672STDY6814755
4672STDY6814756
4672STDY7140017
4672STDY7140018
4672STDY7140019
4672STDY7140020
4672STDY7631188
4672STDY7631189
4672STDY7631190
4672STDY7631191
4672STDY7631192
4672STDY7631193
4672STDY7623907
4672STDY7623908
4672STDY7623909
4672STDY7623910
4672STDY7623911
4672STDY7623912
4672STDY7623913
4672STDY7623914
4672STDY8112878
4672STDY8112879
4672STDY8112880
4672STDY8112881
4672STDY8112882
4672STDY8112883
4672STDY8112884
4672STDY8112885
4672STDY8112974
4672STDY8112975
4672STDY8112976
4672STDY8112977
4672STDY8113070
4672STDY8112979
4672STDY8112980
4672STDY8112981



pf supplementary --type file --file-id-type sample --id data.txt > pf.supplementary.txt


ls -1 *_1.fastq.gz |\
sed -e 's/_1.fastq.gz//g' -e 's/_/\t/g' |\
while read RUN LANE SAMPLE; do 
    cp ${RUN}_${LANE}_${SAMPLE}_1.fastq.gz ${RUN}_${LANE}_S${SAMPLE}_L00${LANE}_R1_001.fastq.gz; 
    cp ${RUN}_${LANE}_${SAMPLE}_2.fastq.gz ${RUN}_${LANE}_S${SAMPLE}_L00${LANE}_R2_001.fastq.gz; 
    done



# Step 2 - run cellranger
bsub.py --queue normal --threads 16 30 02_10X_CR_count "cellranger count \
--id=test \
--transcriptome=/nfs/users/nfs_s/sd21/lustre_link/trichuris_muris/scRNAseq/CELLRANGER/GRCm39 \
--fastqs=/nfs/users/nfs_s/sd21/lustre_link/trichuris_muris/scRNAseq/pathfind_data.txt \
--sample=31438_2 \
--localcores=16 \
--localmem=30"



SampleName_S1_L001_R1_001.fastq.gz