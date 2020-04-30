



~sd21/lustre118_link/software/GENOME_IMPROVEMENT/Metassembler/bin/metassemble


#!/bin/bash
#### Metassembler - merging multiple genome assemblies

#### Job submission information
queue=hugemem
memory=30000
threads=16

#### Parameters to be defined
# experiment name
exp=phusion2-canu-metassemblr

# assemblies to be merged
assembly_1_name=canu
assembly_1_fasta=/nfs/users/nfs_s/sd21/sd21_lustre_link/tc/hybrid_assembly/masucra-canu-merge/canu.fasta
assembly_2_name=phusion
assembly_2_fasta=/nfs/users/nfs_s/sd21/sd21_lustre_link/tc/scaffolding/zemin/worm-phusion.fasta

# short read dataset
shortread_R1=/lustre/scratch108/parasites/sd21/tc/raw/illumina/Tc_1M_427732_WGA2/15063_2#54_nodup_PE1.fastq.gz
shortread_R2=/lustre/scratch108/parasites/sd21/tc/raw/illumina/Tc_1M_427732_WGA2/15063_2#54_nodup_PE2.fastq.gz
shortread_maxins=800
shortread_minins=150
shortread_mean_ins_size=350
shortread_SD=200

###############################################################################################
#### Generate config file
echo -e "\
#######################################################\n\
###   Metassemble configuration file suggested template\n\
#######################################################\n\
\n\
[global]\n\
\n\
#Mate-pair mapping parameters:\n\
bowtie2_threads=${threads}\n\
bowtie2_read1=${shortread_R1}\n\
bowtie2_read2=${shortread_R2}\n\
bowtie2_maxins=${shortread_maxins}\n\
bowtie2_minins=${shortread_minins}\n\
\n\
#CE-stat computation parameters:\n\
#mateAn_A=<int>\n\
#mateAn_B=<int>\n\
#Or:\n\
mateAn_s=${shortread_SD}\n\
mateAn_m=${shortread_mean_ins_size}\n\
\n\
#Whole Genome Alignment parametesr:\n\
nucmer_l=50\n\
nucmer_c=300\n\
\n\
[1]\n\
\n\
fasta=${assembly_1_fasta}\n\
ID=${assembly_1_name}\n\
#mateAn_file=<path>\n\
\n\
[2]\n\
\n\
fasta=${assembly_2_fasta}\n\
ID=${assembly_2_name}\n\
#mateAn_file=<path>\n\
" > ${exp}.config

bsub -q ${queue} -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>${memory}] rusage[mem=${memory}] span[hosts=1]" -M${memory} -n ${threads} -o out.${exp}.%
J.o -e error.${exp}.%J.e -J ${exp} \
/nfs/users/nfs_s/sd21/software/Metassembler/bin/metassemble \
--conf ${exp}.config \
--outd ${exp}_out
