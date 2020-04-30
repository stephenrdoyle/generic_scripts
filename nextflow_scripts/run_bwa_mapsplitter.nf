#!/usr/bin/env nextflow
/*
* ========================================================================================
* run_bwa_mapsplitter.nf
* @authors
* Stephen Doyle <sd21@sanger.ac.uk>
* ----------------------------------------------------------------------------------------
*/


reference_ch = Channel.fromPath(params.ref)
raw_dir_ch = Channel.fromPath(params.rawreads)
params.reads = "$PWD/*_{1,2}.fastq.gz"
params.outdir = 'mapsplitter_results'


Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs_ch }

reference_ch.into {
      reference1
      reference2
}

log.info """\
         MAPSPLITTER   P I P E L I N E
         =============================
         reference : ${params.ref}
         raw_reads : ${raw_dir_ch}
         outdir: ${params.outdir}
         =============================
         """
         .stripIndent()



process index_ref {

     input:
     file reference from reference1

     output:
     file "*" into bwa_index_out

     script:
     """
     bwa index -b 100000000 ${reference}
     """

}


process split_reads {

     tag "$pair_id"

     input:
     tuple val(pair_id), path(reads) from read_pairs_ch

     output:
     tuple val(pair_id), path('*_R1_tmp_*'), path('*_R2_tmp_*') into split_reads_out

     script:
     """
     zcat ${reads[0]} | split -d -a 3 -l 4000000 - ${pair_id}_R1_tmp_
     zcat ${reads[1]} | split -d -a 3 -l 4000000 - ${pair_id}_R2_tmp_
     """
}



process_map_split_reads {
     tag "$pair_id"

     input:
     tuple val(pair_id), path(R1), path(R2) from split_reads_out
     file reference from reference2

     output:
     file "*" into mapped_bams_out

     script:
     """
     bwa mem -t 4 -R '@RG\tID:${pair_id}\tSM:${pair_id}' -Y -M ${reference} ${reads[0]} ${reads[1]}  | samtools view --threads 4 -b - | samtools sort --threads 4 -o ${pair_id}.tmp.sort.bam
     """
}



process_combine_bams {

     input:
     file all_mapped_bams_out from mapped_bams_out.collect()

     output:
     file "*" into combined_bams_out

     script:
     """
     ls -1 *.tmp.sort.bam > bam.fofn
     samtools merge --threads 4 -cpf -b bam.fofn tmp.merged.sorted.bam
     """
}



process_mark_duplicates {

     input:
     file combined_bam from combined_bams_out

     output:
     file deduped_bam into combined_bams_deduped_out

     script:
     """
     java -Xmx20g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/picard-tools-2.5.0/picard.jar MarkDuplicates INPUT=tmp.merged.sorted.bam OUTPUT=merged.sorted.marked.bam METRICS_FILE=tmp.merged.sorted.marked.metrics TMP_DIR=$PWD/tmp
     """
}


process_finish_bam {
     publishDir params.outdir, mode:'copy'

     input:
     file deduped_bam from combined_bams_deduped_out
     output:
     path "*"

     script:
     """
     samtools flagstat merged.sorted.marked.bam > $"{sample_name}".merged.sorted.marked.flagstat
     bamtools stats -in merged.sorted.marked.bam > $"{sample_name}".merged.sorted.marked.bamstats
     samtools view --threads 4 -F 12 -b merged.sorted.marked.bam -o $"{sample_name}".merged.sorted.marked.bam
     samtools view --threads 4 -f 12 merged.sorted.marked.bam -o $"{sample_name}".unmapped.bam
     samtools index -b $"{sample_name}".merged.sorted.marked.bam
     """
}
