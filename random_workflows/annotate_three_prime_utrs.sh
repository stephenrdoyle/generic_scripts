#!/bin/bash

# ========================================================================================
# annotate_three_prime_utrs.sh
#
# Generate a new gff of 3' utrs for a list of transcripts from a gff file. If a
# utr is missing or is below a minimum length, a fake utr can be added at a defined length
#
#
# Usage: annotate_three_prime_utrs.sh <PREFIX> <GFF> <TRANSCRIPT_LIST> <MINIMUM_UTR_LENGTH> <FAKE_UTR_LENGTH>
#
# Example: ./annotate_three_prime_utrs.sh HCON_V4_WBP15 HCON_V4_curated_20200422_WBPS15.gff3 transcript_ids.list 5 1000
#
# @authors
# Stephen Doyle <sd21@sanger.ac.uk>
#
# ----------------------------------------------------------------------------------------

# read inputs from command line and set variables
PREFIX=$1
GFF=$2
TRANSCRIPT_LIST=$3
min_utr_length=$4
fake_utr_length=$5

# example of how a "TRANSCRIPT_LIST" file can be generated
#awk -F'[\t;]' '$3=="mRNA"  {print $9}' HCON_V4_curated_20200422_WBPS15.gff3 | sed 's/Name=//g' > transcript_ids.list


# make the final output file to write the results to.
> ${PREFIX}_three_prime_utr_database.gff

# loop over the transcript list
while read mRNA_ID; do

     grep ${mRNA_ID} ${GFF} > ${mRNA_ID}.utrtmp

     mRNA_direction=$(awk '{print $7}' ${mRNA_ID}.utrtmp | head -n1)

     # determine stand of transcript, so that utr can be placed correctly.
     if [ $mRNA_direction == "+" ];
     then

          chr=$(cut -f1 ${mRNA_ID}.utrtmp | head -n1 )
          cds_last=$(grep "CDS" ${mRNA_ID}.utrtmp | sort -k5 | tail -n1 | cut -f5 )
          exon_last=$(grep "exon" ${mRNA_ID}.utrtmp | sort -k5 | tail -n1 | cut -f5 )
          let  "utr_length=${exon_last}-${cds_last}"

     # check length of utr relative to minimum length threshold. Keep if above threshold, or use fake utr length if too short
     if [ $utr_length -le $min_utr_length ];
     then
          utr_length=${fake_utr_length} ;
     elif [ $utr_length -gt 1000 ];
     then
          utr_length=${fake_utr_length} ;
     else
          utr_length=${utr_length};
     fi

          let "utr_start=${cds_last}+1"
          let "utr_end=${utr_start}+${utr_length}-1"

          # print data into new file
          printf ${chr}"\t""HC_WBP15_3UTRs""\t""three_prime_utr""\t"${utr_start}"\t"${utr_end}"\t"".""\t"${mRNA_direction}"\t"${utr_length}"\t""Name="${mRNA_ID}"_three_prime_utr;ID="${mRNA_ID}"_three_prime_utr;Parent="${mRNA_ID}"\n" >> ${PREFIX}_three_prime_utr_database.gff;

     else

     # if not on the positive (+) strand, transcript must be on the antisense strand (-).
     chr=$(cut -f1 ${mRNA_ID}.utrtmp | head -n1 )
     cds_last=$(grep "CDS" ${mRNA_ID}.utrtmp | sort -k4 | head -n1 | cut -f4 )
     exon_last=$(grep "exon" ${mRNA_ID}.utrtmp | sort -k4 | head -n1 | cut -f4 )
     let  "utr_length=${cds_last}-${exon_last}"

     # check length of utr relative to minimum length threshold. Keep if above threshold, or use fake utr length if too short
     if [ $utr_length -le $min_utr_length ];
     then
          utr_length=${fake_utr_length} ;
     elif [ $utr_length -gt 1000 ];
     then
          utr_length=${fake_utr_length} ;
     else
          utr_length=${utr_length};
     fi

     let "utr_start=${cds_last}-1"
     let "utr_end=${utr_start}-${utr_length}+1"

     # print data into new file
     printf ${chr}"\t""HC_WBP15_3UTRs""\t""three_prime_utr""\t"${utr_end}"\t"${utr_start}"\t"".""\t"${mRNA_direction}"\t"${utr_length}"\t""Name="${mRNA_ID}"_three_prime_utr;ID="${mRNA_ID}"_three_prime_utr;Parent="${mRNA_ID}"\n" >> ${PREFIX}_three_prime_utr_database.gff;
     fi;

# remove the temporary transcript file
rm *.utrtmp

done < ${TRANSCRIPT_LIST}




# to get fasta sequences after generating the gff, could do something like the following:
# cat HCON_V4_WBP15_three_prime_utr_database.gff | sed 's/Name=//g' | awk -F'[\t;]' '{print $1,$4,$5,$9,"1",$7}' OFS="\t" > HCON_V4_WBP15_three_prime_utr_database.bed
# bedtools getfasta -fi ../../../REF/HAEM_V4_final.chr.fa -fo HCON_V4_WBP15_three_prime_utr_database.fasta -bed HCON_V4_WBP15_three_prime_utr_database.bed -name -s
