# clust analysis of gene expression profiles

cd ~/lustre_link/schistosoma_mansoni/RNASEQ_SINGLE_WORM/data

while read samplename lane1 lane2; do 
    pf data --type lane --id $lane1 --filetype fastq --symlink ./ ; 
    done < samplename_lane1_lane2.txt


while read samplename lane1 lane2; do 
    pf data --type lane --id $lane2 --filetype fastq --symlink ./ ; 
    done < samplename_lane1_lane2.txt






### Run Kallisto
```shell
cd ~/lustre_link/schistosoma_mansoni/RNASEQ_SINGLE_WORM/mapping


wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS18/species/schistosoma_mansoni/PRJEA36577/schistosoma_mansoni.PRJEA36577.WBPS18.genomic.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS18/species/schistosoma_mansoni/PRJEA36577/schistosoma_mansoni.PRJEA36577.WBPS18.annotations.gff3.gz
gunzip *gz

# make a transcripts fasta
gffread -x TRANSCRIPTS.fa -g schistosoma_mansoni.PRJEA36577.WBPS18.genomic.fa schistosoma_mansoni.PRJEA36577.WBPS18.annotations.gff3


# index the transcripts
module load kallisto/0.46.2--h4f7b962_1

kallisto index --index TRANSCRIPTS.ixd TRANSCRIPTS.fa

# run kallisto
ln -s ../data/samplename_lane1_lane2.txt

###------- placed in "run_kalliso.sh"
while read samplename lane1 lane2; do
kallisto quant \
     --bias \
     --index TRANSCRIPTS.ixd \
     --output-dir ${samplename}_out \
     --bootstrap-samples 100 \
     --threads 7 \
     --fusion \
     ../data/${lane1}_1.fastq.gz ../data/${lane1}_2.fastq.gz ../data/${lane2}_1.fastq.gz ../data/${lane2}_2.fastq.gz;
done < samplename_lane1_lane2.txt
###-------

bsub.py --queue long --threads 7 20 kallisto_quant "bash ./run_kalliso.sh"



mkdir KALLISTO_MAPPED_SAMPLES
mv kallisto_* KALLISTO_MAPPED_SAMPLES/





### generate the tpm data
# extract TPMs per sample
for i in ` ls -1d *out `; do 
    echo $i > ${i}.tpm ; cat ${i}/abundance.tsv | cut -f5 | sed '1d' >> ${i}.tpm; 
    done

# generate a "transcripts" list, taken from the TRANSCRIPTS.fa file

echo "ID" > transcripts.list; grep ">" TRANSCRIPTS.fa | cut -f1 -d" " | sed -e 's/>//g' >> transcripts.list



# make a data frame containing all TMP values from all samples
paste transcripts.list \
X_Eggs_R1_out.tpm \
X_Eggs_R2_out.tpm \
X_Eggs_R3_out.tpm \
X_Eggs_R4_out.tpm \
X_Eggs_R5_out.tpm \
X_Miracidia_R1_out.tpm \
X_Miracidia_R2_out.tpm \
X_Miracidia_R3_out.tpm \
X_Miracidia_R4_out.tpm \
X_Miracidia_R5_out.tpm \
X_1d_Sporocysts_R1_out.tpm \
X_1d_Sporocysts_R2_out.tpm \
X_1d_Sporocysts_R3_out.tpm \
X_1d_Sporocysts_R4_out.tpm \
X_1d_Sporocysts_R5_out.tpm \
X_5d_Sporocysts_R1_out.tpm \
X_5d_Sporocysts_R2_out.tpm \
X_5d_Sporocysts_R3_out.tpm \
X_5d_Sporocysts_R4_out.tpm \
X_5d_Sporocysts_R5_out.tpm \
M_32d_Sporocysts_R1_out.tpm \
M_32d_Sporocysts_R2_out.tpm \
M_32d_Sporocysts_R3_out.tpm \
M_32d_Sporocysts_R4_out.tpm \
M_32d_Sporocysts_R5_out.tpm \
M_Cercariae_R1_out.tpm \
M_Cercariae_R2_out.tpm \
M_Cercariae_R3_out.tpm \
M_Cercariae_R4_out.tpm \
M_Cercariae_R5_out.tpm \
M_2d_Somules_R1_out.tpm \
M_2d_Somules_R2_out.tpm \
M_2d_Somules_R3_out.tpm \
M_2d_Somules_R4_out.tpm \
M_2d_Somules_R5_out.tpm \
M_26d_Juveniles_R1_out.tpm \
M_26d_Juveniles_R2_out.tpm \
M_26d_Juveniles_R3_out.tpm \
M_26d_Juveniles_R4_out.tpm \
M_26d_Juveniles_R5_out.tpm \
> kallisto_male.tpm.table


paste transcripts.list \
X_Eggs_R1_out.tpm \
X_Eggs_R2_out.tpm \
X_Eggs_R3_out.tpm \
X_Eggs_R4_out.tpm \
X_Eggs_R5_out.tpm \
X_Miracidia_R1_out.tpm \
X_Miracidia_R2_out.tpm \
X_Miracidia_R3_out.tpm \
X_Miracidia_R4_out.tpm \
X_Miracidia_R5_out.tpm \
X_1d_Sporocysts_R1_out.tpm \
X_1d_Sporocysts_R2_out.tpm \
X_1d_Sporocysts_R3_out.tpm \
X_1d_Sporocysts_R4_out.tpm \
X_1d_Sporocysts_R5_out.tpm \
X_5d_Sporocysts_R1_out.tpm \
X_5d_Sporocysts_R2_out.tpm \
X_5d_Sporocysts_R3_out.tpm \
X_5d_Sporocysts_R4_out.tpm \
X_5d_Sporocysts_R5_out.tpm \
F_32d_Sporocysts_R1_out.tpm \
F_32d_Sporocysts_R2_out.tpm \
F_32d_Sporocysts_R3_out.tpm \
F_32d_Sporocysts_R4_out.tpm \
F_32d_Sporocysts_R5_out.tpm \
F_Cercariae_R1_out.tpm \
F_Cercariae_R2_out.tpm \
F_Cercariae_R3_out.tpm \
F_Cercariae_R4_out.tpm \
F_Cercariae_R5_out.tpm \
F_2d_Somules_R1_out.tpm \
F_2d_Somules_R2_out.tpm \
F_2d_Somules_R3_out.tpm \
F_2d_Somules_R4_out.tpm \
F_2d_Somules_R5_out.tpm \
F_26d_Juveniles_R1_out.tpm \
F_26d_Juveniles_R2_out.tpm \
F_26d_Juveniles_R3_out.tpm \
F_26d_Juveniles_R4_out.tpm \
F_26d_Juveniles_R5_out.tpm \
> kallisto_female.tpm.table



# female replicates table - placed in "females.replicates.txt"
kallisto_female.tpm.table	eggs	X_Eggs_R1_out	X_Eggs_R2_out	X_Eggs_R3_out	X_Eggs_R4_out	X_Eggs_R5_out
kallisto_female.tpm.table	miracidia	X_Miracidia_R1_out	X_Miracidia_R2_out	X_Miracidia_R3_out	X_Miracidia_R4_out	X_Miracidia_R5_out
kallisto_female.tpm.table	1d-sporocysts	X_1d_Sporocysts_R1_out	X_1d_Sporocysts_R2_out	X_1d_Sporocysts_R3_out	X_1d_Sporocysts_R4_out	X_1d_Sporocysts_R5_out
kallisto_female.tpm.table	5d-sporocysts	X_5d_Sporocysts_R1_out	X_5d_Sporocysts_R2_out	X_5d_Sporocysts_R3_out	X_5d_Sporocysts_R4_out	X_5d_Sporocysts_R5_out
kallisto_female.tpm.table	32d-sporocysts	F_32d_Sporocysts_R1_out	F_32d_Sporocysts_R2_out	F_32d_Sporocysts_R3_out	F_32d_Sporocysts_R4_out	F_32d_Sporocysts_R5_out
kallisto_female.tpm.table	cercariae	F_Cercariae_R1_out	F_Cercariae_R2_out	F_Cercariae_R3_out	F_Cercariae_R4_out	F_Cercariae_R5_out
kallisto_female.tpm.table	2d-somules	F_2d_Somules_R1_out	F_2d_Somules_R2_out	F_2d_Somules_R3_out	F_2d_Somules_R4_out	F_2d_Somules_R5_out
kallisto_female.tpm.table	26d-juveniles	F_26d_Juveniles_R1_out	F_26d_Juveniles_R2_out	F_26d_Juveniles_R3_out	F_26d_Juveniles_R4_out	F_26d_Juveniles_R5_out


# run clust 
module load clust/1.17.0--pyhfa5458b_0

clust kallisto_female.tpm.table -r females.replicates.txt