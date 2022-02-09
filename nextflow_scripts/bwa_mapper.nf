#!/usr/bin/env nextflow
/*
 * Authors:
 * - Stephen Doyle <stephen.doyle@sanger.ac.uk>
 *
 * inspired by Anderson labs https://github.com/AndersenLab/nil-ril-nf/blob/master/main.nf
 */

 // update to nextflow 20+ DSL2
  if( !nextflow.version.matches('21.0+') ) {
     println "This workflow requires Nextflow version 21.0 or greater -- You are running version $nextflow.version"
     //exit 1
 }


nextflow.enable.dsl=2


// defaults
date = new Date().format( 'yyyyMMdd' )
params.debug = false
params.cores = 4
params.tmpdir = "/tmp"
params.relative = true
params.out = "./MAPPING.${date}"


// debug
if (params.debug == true) {
println """
    ***Using debug mode***
"""
params.fastq_fofn = "/nfs/users/nfs_s/sd21/lustre118_link/play/NEXTFLOW/test_data/sample_sheet.txt"
params.reference = "/nfs/users/nfs_s/sd21/lustre118_link/play/NEXTFLOW/test_data/hcontortus_chr_mtDNA_arrow_pilon.fa"
} else {
    params.fastq_fofn = "(required)"
    params.reference = "(required)"
}



// checks
if (params.reference == "(required)" || params.fastq_fofn == "(required)") {
    println """
    Error: Reference, and FQ sheet are required for analysis. Please check parameters.
    Reference: ${params.reference}
    fastqs: ${params.fastq_fofn}
    """
    System.exit(1)
}

if (!reference.exists()) {
    println """
    Error: Reference ${params.reference} does not exist.
    """
    System.exit(1)
}

if (!fq_file.exists()) {
    println """
    Error: fastq sheet ${params.fastq_fofn} does not exist.
    """
    System.exit(1)
}



// reference
File reference = new File("${params.reference}")
if (params.reference != "(required)") {
   reference_handle = reference.getAbsolutePath();
} else {
    reference_handle = "(required)"
}

// Define FASTQ files
File fq_file = new File("${params.fastq_fofn}")
if (params.relative) {
    fq_file_prefix = fq_file.getParentFile().getAbsolutePath();
    fastq_input = Channel.from(fq_file.collect { it.tokenize( ',' ) })
                 .map { SM, ID, LB, read1, read2 -> [SM, ID, LB, file("${fq_file_prefix}/${read1}"), file("${fq_file_prefix}/${read2}")] }
} else {
    fastq_input = Channel.from(fq_file.collect { it.tokenize( ',' ) })
             .map { SM, ID, LB, read1, read2 -> [SM, ID, LB, file("${read1}"), file("${read1}")] }
}



param_summary = '''
worm_mapper

''' + """
    parameters           description                    Set/Default
    ==========           ===========                    =======

    --debug              Set to 'true' to test          ${params.debug}
    --cores              Number of cores                ${params.cores}
    --out                Directory to output results    ${params.out}
    --fastq_input        fastq file (see help)          ${params.fastq_fofn}
    --reference          Reference Genome               ${reference_handle}
    --tmpdir             A temporary directory          ${params.tmpdir}
    --email              Email to be sent results       ${params.email}
    HELP: TO DO
  ---------------------------------------------------------------------------
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
"""

println param_summary



// Generate workflow
workflow {

    // alignment
    fastq_input | perform_alignment
    perform_alignment.out.sample_aligned_bams.groupTuple() | merge_bam
    merge_bam.out.merged_SM | SM_coverage
    SM_coverage.out.toSortedList() | SM_coverage_merge
    perform_alignment.out.aligned_bams
    perform_alignment.out.aligned_bams | fq_coverage
    fq_coverage.out
        .toSortedList() | fq_coverage_merge



    // stats
    merge_bam.out.merged_SM | SM_bam_stats
    SM_bam_stats.out
        .toSortedList() | combine_SM_bam_stats
    merge_bam.out.merged_SM | idx_stats_SM
    idx_stats_SM.out
        .toSortedList() | combine_idx_stats
    perform_alignment.out.aligned_bams | fq_idx_stats
    fq_idx_stats.out
        .toSortedList() | fq_combine_idx_stats
    perform_alignment.out.aligned_bams | fq_bam_stats
    fq_bam_stats.out
        .toSortedList() | combine_bam_stats


}



/*
    ===============
    Fastq alignment
    ===============
*/

process perform_alignment {

    cpus params.cores

    tag { ID }

    input:
        tuple val(SM), val(ID), val(LB), val(read1), val(read2)

    output:
        tuple val(SM), val(ID), val(LB), file("${ID}.bam"), file("${ID}.bam.bai"), emit: aligned_bams
        tuple val(SM), file("${ID}.bam"), emit: sample_aligned_bams

    """
        bwa mem -t ${task.cpus} -R '@RG\\tID:${ID}\\tLB:${LB}\\tSM:${SM}' ${reference_handle} ${read1} ${read2} | \\
        sambamba view --nthreads=${task.cpus} --show-progress --sam-input --format=bam --with-header /dev/stdin | \\
        sambamba sort --nthreads=${task.cpus} --show-progress --tmpdir=${params.tmpdir} --out=${ID}.bam /dev/stdin 2>&1
        sambamba index --nthreads=${task.cpus} ${ID}.bam
        if [[ ! \$(samtools view ${ID}.bam | head -n 10) ]]; then
            exit 1;
        fi
    """
}

/*
    fq idx stats
*/

process fq_idx_stats {

    tag { ID }

    input:
        tuple val(SM), val(ID), val(LB), file("${ID}.bam"), file("${ID}.bam.bai")

    output:
        file 'fq_idxstats'

    """
        samtools idxstats ${ID}.bam | awk '{ print "${ID}\\t" \$0 }' > fq_idxstats
    """
}

process fq_combine_idx_stats {

    publishDir params.out + "/fq", mode: 'copy'

    input:
        val bam_idxstats

    output:
        file("fq_bam_idxstats.tsv")

    """
        echo -e "SM\\treference\\treference_length\\tmapped_reads\\tunmapped_reads" > fq_bam_idxstats.tsv
        cat ${bam_idxstats.join(" ")} >> fq_bam_idxstats.tsv
    """

}

/*
    fq bam stats
*/

process fq_bam_stats {

    tag { ID }

    input:
        tuple val(SM), val(ID), val(LB), file("${ID}.bam"), file("${ID}.bam.bai")

    output:
        file 'bam_stat'

    """
        cat <(samtools stats ${ID}.bam | grep ^SN | cut -f 2- | awk '{ print "${ID}\t" \$0 }' | sed 's/://g') > bam_stat
    """
}

process combine_bam_stats {

    publishDir params.out + "/fq", mode: 'copy'

    input:
        val stat_files

    output:
        file("fq_bam_stats.tsv")

    """
    echo -e "fq_pair_id\\tvariable\\tvalue\\tcomment" > fq_bam_stats.tsv
    cat ${stat_files.join(" ")} >> fq_bam_stats.tsv
    """
}

/*
    Fastq coverage
*/
process fq_coverage {

    tag { ID }

    input:
        tuple val(SM), val(ID), val(LB), file("${ID}.bam"), file("${ID}.bam.bai")
    output:
        file("${ID}.coverage.tsv")


    """
    samtools faidx ${reference_handle}
    awk '{print \$1,1,\$2}' OFS="\t" ${reference_handle}.fai > ${reference_handle}.bed
    samtools bedcov ${reference_handle}.bed ${ID}.bam > ${ID}.coverage.tsv
    rm ${reference_handle}.bed
    """
}

process fq_coverage_merge {

    publishDir params.out + "/fq", mode: 'copy'

    input:
        val fq_set

    output:
        file("fq_coverage.full.tsv")
        path "fq_coverage.tsv", emit: fq_coverage_plot

    """
        echo -e 'bam\\tcontig\\tstart\\tend\\tproperty\\tvalue' > fq_coverage.full.tsv
        cat ${fq_set.join(" ")} >> fq_coverage.full.tsv
        cat <(echo -e 'fq\\tcoverage') <( cat fq_coverage.full.tsv | grep 'genome' | grep 'depth_of_coverage' | cut -f 1,6) > fq_coverage.tsv
    """
}


process merge_bam {

    echo true

    storeDir params.out + "/bam"

    tag { SM }

    cpus params.cores

    memory { 16.GB * task.attempt }

    errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }

    input:
        tuple val(SM), val(bam)

    output:
        tuple val(SM), file("${SM}.final.bam"), file("${SM}.final.bam.bai"), emit: merged_SM

    """
    count=`echo ${bam.join(" ")} | tr ' ' '\\n' | wc -l`
    if [ "\${count}" -eq "1" ]; then
        ln -s ${bam[0]} ${SM}.merged.bam
        ln -s ${bam[0]}.bai ${SM}.merged.bam.bai
    else
        sambamba merge --nthreads=${task.cpus} --show-progress ${SM}.merged.bam ${bam.join(" ")}
        sambamba index --nthreads=${task.cpus} ${SM}.merged.bam
    fi
	samtools collate --threads=${task.cpus} -o ${SM}.namecollate.bam ${SM}.merged.bam
	samtools fixmate --threads=${task.cpus} -m ${SM}.namecollate.bam ${SM}.fixmate.bam
	samtools sort --threads=${task.cpus} -o ${SM}.positionsort.bam ${SM}.fixmate.bam
	samtools markdup --threads=${task.cpus} ${SM}.positionsort.bam ${SM}.final.bam
	sambamba index --nthreads=${task.cpus} ${SM}.final.bam

    """
}


/*
 SM idx stats
*/

process idx_stats_SM {

    tag { SM }

    input:
        tuple val(SM), file("${SM}.bam"), file("${SM}.bam.bai")
    output:
        file 'SM_bam_idxstats'

    """
        samtools idxstats ${SM}.bam | awk '{ print "${SM}\\t" \$0 }' > SM_bam_idxstats
    """
}

process combine_idx_stats {

    publishDir params.out +"/SM", mode: 'copy'

    input:
        val bam_idxstats

    output:
        file("SM_bam_idxstats.tsv")

    """
        echo -e "SM\\treference\\treference_length\\tmapped_reads\\tunmapped_reads" > SM_bam_idxstats.tsv
        cat ${bam_idxstats.join(" ")} >> SM_bam_idxstats.tsv
    """

}


/*
    SM bam stats
*/

process SM_bam_stats {

    tag { SM }

    input:
        tuple val(SM), file("${SM}.bam"), file("${SM}.bam.bai")

    output:
        file 'bam_stat'

    """
        cat <(samtools stats ${SM}.bam | grep ^SN | cut -f 2- | awk '{ print "${SM}\t" \$0 }' | sed 's/://g') > bam_stat
    """
}

process combine_SM_bam_stats {

    publishDir params.out + "/SM", mode: 'copy'

    input:
        val stat_files

    output:
        file("SM_bam_stats.tsv")

    """
        echo -e "fq_pair_id\\tvariable\\tvalue\\tcomment" > SM_bam_stats.tsv
        cat ${stat_files.join(" ")} >> SM_bam_stats.tsv
    """
}





/*
    Coverage Bam
*/
process SM_coverage {

    tag { SM }

    input:
        tuple val(SM), file("${SM}.bam"), file("${SM}.bam.bai")

    output:
        file("${SM}.coverage.tsv")

    """
	samtools faidx ${reference_handle}
        awk '{print \$1,1,\$2}' OFS="\t" ${reference_handle}.fai > ${reference_handle}.bed
        samtools bedcov ${reference_handle}.bed ${SM}.bam > ${SM}.coverage.tsv
    """
}


process SM_coverage_merge {

    publishDir params.out + "/SM", mode: 'copy'


    input:
        val sm_set

    output:
        file("SM_coverage.full.tsv")
        path "SM_coverage.tsv", emit: SM_coverage_plot

    """
        echo -e 'SM\\tcontig\\tstart\\tend\\tproperty\\tvalue' > SM_coverage.full.tsv
        cat ${sm_set.join(" ")} >> SM_coverage.full.tsv
        # Generate condensed version
        cat <(echo -e 'strain\\tcoverage') <(cat SM_coverage.full.tsv | grep 'genome' | grep 'depth_of_coverage' | cut -f 1,6) > SM_coverage.tsv
    """

}



workflow.onComplete {

    summary = """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
    Parameters
    ----------
    debug              Set to 'true' to test          ${params.debug}
    cores              Number of cores                ${params.cores}
    out                Directory to output results    ${params.out}
    fastq_input        fastq file (see help)          ${params.fastq_fofn}
    reference          Reference Genome               ${reference_handle}
    """

    println summary


}
