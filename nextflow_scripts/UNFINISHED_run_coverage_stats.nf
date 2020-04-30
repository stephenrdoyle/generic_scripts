#!/usr/bin/env nextflow
/*
* ========================================================================================
* run_coverage_stats.nf
*
* Generate coverage statistics per contig and per window for one or more bam files. Bam files are
* automatically detected, and a window size of 10kb is used, however, can be set manually on the
* command line (see example).
*
*    Usage: nextflow ~sd21/bash_scripts/run_coverage_stats.nf [optional (default = 10000): --window_length=XXX]
*
* @authors
* Stephen Doyle <sd21@sanger.ac.uk>
* ----------------------------------------------------------------------------------------
*/


params.window_length = 10000
params.bams = "$PWD/*.bam"
params.bams_index = "$PWD/*.bai"
params.outdir = 'bam_coverage_results'

bam_files = Channel.fromPath(params.bams)
bai_files = Channel.fromPath(params.bams_index)

bam_files.into {
  bam_files1
  bam_files2
}

bai_files.into {
  bai_files1
  bai_files2
}

process get_coverage_window {
     publishDir params.outdir, mode:'copy'

     clusterOptions '-M 5000 -R "select[mem>5000] rusage[mem=5000]"'

     input:
     path bam from bam_files2
     path bai from bai_files2
     val window_length from params.window_length

     output:
     file "${bam.baseName}*.cov"

     module 'samtools/1.6--h244ad75_4:bamtools/2.5.1--he860b03_5:bedtools/2.29.0--hc088bd4_3'

     shell:
     '''
     bamtools header -in !{bam} | grep "^@SQ" | awk -F'[:\\t]' '{print $3,$5}' OFS="\\t" > !{bam}.chr.genome
     bedtools makewindows -g !{bam}.chr.genome -w !{window_length} > !{bam}.!{window_length}_window.bed
     samtools bedcov -Q 20 !{bam}.!{window_length}_window.bed !{bam} | awk -F'\\t' '{print $1,$2,$3,$4,$4/($3-$2)}' OFS="\\t" > !{bam}.!{window_length}_window.cov
     '''
}

process get_coverage_chromosome {
     publishDir params.outdir, mode:'copy'

     clusterOptions '-M 5000 -R "select[mem>5000] rusage[mem=5000]"'

     input:
     path bam from bam_files1
     path bai from bai_files1

     output:
     file "${bam.baseName}*.cov"

     module 'samtools/1.6--h244ad75_4:bamtools/2.5.1--he860b03_5:bedtools/2.29.0--hc088bd4_3'

     shell:
     '''
     bamtools header -in !{bam} | grep "^@SQ" | awk -F'[:\\t]' '{print $3,"1",$5}' OFS="\\t" > !{bam}.chromosome.bed
     samtools bedcov -Q 20 !{bam}.chromosome.bed !{bam} | awk -F'\\t' '{print $1,$2,$3,$4,$4/($3-$2)}' OFS="\\t" > !{bam}.chromosome.cov
     '''
}





process cleanup {
     input:


     output:


     script:
     """

     """
}

















##########################################################################################
# run_cov_stats
##########################################################################################

# Usage: ~sd21/bash_scripts/run_cov_stats < window size >


WINDOW=$1

for i in *.bam
do

bamtools header -in ${i} | grep "^@SQ" | awk -F'[:\t]' '{printf $3"\t"1"\t"$5"\n"}' OFS="\t" > ${i%.bam}.chr.bed
bamtools header -in ${i} | grep "^@SQ" | awk -F'[:\t]' '{printf $3"\t"$5"\n"}' OFS="\t" > ${i%.bam}.chr.genome

bedtools makewindows -g ${i%.bam}.chr.genome -w ${WINDOW} > ${i%.bam}.${WINDOW}_window.bed

samtools-1.6 bedcov -Q 20 ${i%.bam}.chr.bed ${i} | awk -F'\t' '{printf $1"\t"$2"\t"$3"\t"$4"\t"$4/($3-$2)"\n"}' OFS="\t" > ${i%.bam}.chr.cov
samtools-1.6 bedcov -Q 20 ${i%.bam}.${WINDOW}_window.bed ${i} | awk -F'\t' '{printf $1"\t"$2"\t"$3"\t"$4"\t"$4/($3-$2)"\n"}' OFS="\t" > ${i%.bam}.${WINDOW}_window.cov

rm ${i%.bam}.chr.bed ${i%.bam}.${WINDOW}_window.bed ${i%.bam}.chr.genome

done

for i in *.chr.cov; do printf "${i}\n" > ${i}.tmp | awk '{print $5}' OFS="\t" ${i} >> ${i}.tmp; done
paste *.tmp > coverage_stats.summary
rm *.tmp#!/usr/bin/env nextflow
