#!/usr/bin/env nextflow
/*
* ========================================================================================
* run_busco_cegma.nf
* @authors
* Stephen Doyle <sd21@sanger.ac.uk>
*
*
*    Usage: nextflow ~sd21/bash_scripts/run_busco_metazoa.nf --ref '<fasta>' OR --ref '*.fasta'
*
* ----------------------------------------------------------------------------------------
*/


reference_ch = Channel.fromPath(params.ref)
params.outdir = 'busco_cegma_results'
params.busco_lineage = '/nfs/users/nfs_s/sd21/databases/busco/metazoa_odb9'

reference_ch.into {
  reference1
  reference2
}

log.info """\
         BUSCO / CEGMA   P I P E L I N E
         =============================
         reference : ${params.ref}
         outdir: ${params.outdir}
         =============================
         """
         .stripIndent()


process busco {
     publishDir params.outdir, mode:'copy'

     cpus 15
     queue 'long'
     clusterOptions '-M 8000 -R "select[mem>8000] rusage[mem=8000]"'

     input:
     file reference from reference1

     output:
     path 'run*_busco3.02.metazoa'


     """
     export AUGUSTUS_CONFIG_PATH=/nfs/users/nfs_s/sd21/software/augustus-3.2.1/config
     export PATH=$PATH:/nfs/users/nfs_s/sd21/software/augustus-3.2.1/bin/

     /nfs/users/nfs_s/sd21/lustre118_link/software/ASSEMBLY_QC/busco_v3/scripts/run_BUSCO.py \
     --in ${reference} \
     --out ${reference}_busco3.02.metazoa \
     --mode genome \
     --lineage_path /nfs/users/nfs_s/sd21/databases/busco/metazoa_odb9 \
     --species caenorhabditis \
     --cpu 15 --tarzip --force --restart --long --blast_single_core \
     --tmp_path ${reference}.tmp
     """
     }

process cegma {
     publishDir params.outdir, mode:'copy'

     input:
     file reference from reference2

     output:
     path 'cegma2.5_*_out'

     """
     perl /nfs/users/nfs_s/sd21/lustre118_link/software/ASSEMBLY_QC/CEGMA_v2.5/bin/cegma \
     --genome ${reference} \
     --output cegma2.5_${reference}_out
     --quiet

     """
     }
