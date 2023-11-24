# load kraken 2
module load kraken2/2.1.2

cd /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/GENOME/TRANSCRIPTOME/RAW

# run kraken on the modern PE trimmed reads
for i in $(ls -1 *_1.fastq.gz | sed 's/_1.fastq.gz//g'); do 
     bsub.py 10 kraken2 "kraken2 --db /nfs/users/nfs_s/sd21/lustre_link/databases/kraken/kraken2-microbial-fatfree \
     --report ${i}.kraken2report \
     --paired ${i}_1.fastq.gz ${i}_2.fastq.gz";
done



# once the kraken runs have completed, run multiqc .
multiqc *kraken2report --title kraken