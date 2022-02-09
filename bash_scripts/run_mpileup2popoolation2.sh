#!/usr/bin/env bash

# ========================================================================================
# run_mpileup2popoolation2.sh
#
# Wrapper to convert a mpileup file to popoolation2 FST and FET datasets,
#
# Reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3232374/pdf/btr589.pdf
# Software: https://sourceforge.net/p/popoolation2/wiki/Main/
#
# Usage: ~sd21/bash_scripts/run_mpileup2popoolation2.sh <PREFIX> <REFERENCE_FASTA> <MPILEUP> <POOLSIZE> <WINDOW_SIZE>
#
# @authors
# Stephen Doyle <sd21@sanger.ac.uk>
#
# ----------------------------------------------------------------------------------------

PREFIX=$1
FASTA=$2
MPILEUP=$3
POOL_SIZE=$4
WINDOW=$5

ID="U$(date +%s)"

grep ">" ${FASTA} | sed -e 's/>//g' | cat -n > ref.fa.sequence-names.tmp


#--- make sync
echo -e "\
java -jar ~sd21/lustre118_link/software/POOLSEQ/popoolation2_1201/mpileup2sync.jar \
--input ${MPILEUP} \
--output ${PREFIX}.raw.sync \
--min-qual 20 \
--threads 7 \
--fastq-type sanger

while read NUMBER NAME; do
grep "\${NAME}" ${PREFIX}.raw.sync > \${NUMBER}.\${NAME}.raw.sync.tmp
done < ref.fa.sequence-names.tmp" > run_make_syncronised_file.tmp.${ID}

chmod a+x run_make_syncronised_file.tmp.${ID}


#---- make fst and fet arrays
while read NUMBER NAME; do
echo -e "\
perl ~sd21/lustre118_link/software/POOLSEQ/popoolation2_1201/fisher-test.pl --input ${NUMBER}.${NAME}.raw.sync.tmp --output ${NUMBER}.${NAME}.fet.tmp --min-count 4 --min-coverage 30 --max-coverage 2% --suppress-noninformative

perl ~sd21/lustre118_link/software/POOLSEQ/popoolation2_1201/fst-sliding.pl --pool-size ${POOL_SIZE} --window-size ${WINDOW} --step-size ${WINDOW} --min-count 4 --min-coverage 30 --max-coverage 2% --input ${NUMBER}.${NAME}.raw.sync.tmp --output ${NUMBER}.${NAME}.raw.fst.tmp" > run_pp2_split.tmp.${ID}.${NUMBER};
done < ref.fa.sequence-names.tmp
chmod a+x *run_pp2_split*

echo -e "\
PREFIX=\$1
cat \$(find ./ -name \"*.fet.tmp\" | sort -V) | sed -e \"s/=/\\\t/g\" > \${PREFIX}.merged.fet
cat \$(find ./ -name \"*.fst.tmp\" | sort -V) | sed -e \"s/=/\\\t/g\" > \${PREFIX}.merged.fst
#rm *.tmp*" > run_pp2_combine.tmp.${ID}
chmod a+x run_pp2_combine.tmp.${ID}




jobs=$( wc -l ref.fa.sequence-names.tmp | cut -f1 -d" " )

bsub -q yesterday -R'span[hosts=1] select[mem>10000] rusage[mem=10000]' -M10000 -n 7 -e 01_mp2pp2_makesync.e -o 01_mp2pp2_makesync.o -J 01_mp2pp2_makesync ./run_make_syncronised_file.tmp.${ID}

bsub -q long -w "done(01_mp2pp2_makesync)" -R'span[hosts=1] select[mem>1000] rusage[mem=1000]' -M1000 -J 02_mp2pp2_array[1-$jobs] -e 02_pp2_array[1-$jobs].e -o 02_pp2_array[1-$jobs].o ./run_pp2_split.tmp.${ID}.\$LSB_JOBINDEX

bsub -q normal -w "done(02_mp2pp2_array[1-$jobs])" -R'span[hosts=1] select[mem>1000] rusage[mem=1000]' -n1 -M1000 -J 03_mp2pp2_combine -e 03_pp2_combine.e -o 03_pp2_combine.o ./run_pp2_combine.tmp.${ID} ${PREFIX}
