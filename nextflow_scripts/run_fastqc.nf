#!/usr/bin/env nextflow
/*
* ========================================================================================
* run_fastqc.nf
* @authors
* Stephen Doyle <sd21@sanger.ac.uk>
* ----------------------------------------------------------------------------------------
*/

params.reads = "$PWD/*_{1,2}.fastq.gz"
params.outdir = 'multiqc_results'

log.info """\
         FASTQC/MULTIQC   P I P E L I N E
         =============================
         reads : ${params.reads}
         outdir: ${params.outdir}
         =============================
         """
         .stripIndent()


Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs_ch }


process fastqc {
     tag "$pair_id"

     input:
     tuple val(pair_id), path(reads) from read_pairs_ch

     output:
     file "*.{zip,html}" into fastqc_results

     module 'fastqc/0.11.8-c2:ISG/jre/1.8.0_131'

     """
     fastqc ${reads[0]} ${reads[1]}
     """
}

process multiqc {
    publishDir params.outdir, mode:'copy'

    input:
    file all_fastqc_results from fastqc_results.collect()

    output:
    path 'multiqc_report.html'
    path 'multiqc_data'

    module 'multiqc/1.8:ISG/jre/1.8.0_131'

    """
    multiqc .
    """
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}

workflow.onComplete {

    def msg = """\
        Pipeline execution summary: run_fastqc_
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()

    sendMail(to: 'sd21@sanger.ac.uk', subject: 'My pipeline execution', body: msg)
}
