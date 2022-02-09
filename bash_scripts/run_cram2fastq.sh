#!/usr/bin/env bash

# ========================================================================================
# run_cram2fastq.sh
#
# Convert CRAM file(s) to paired end FASTAQ files(s). Run in directory with CRAM files
#
# Reference:
# Software:
#
# Usage: ~sd21/bash_scripts/run_cram2fastq.sh
#
# @authors
# Stephen Doyle <sd21@sanger.ac.uk>
#
# ----------------------------------------------------------------------------------------



for i in *cram; do samtools view -ub --threads 4 ${i} | samtools sort -n - | samtools fastq -1 ${i%.cram}_1.fastq.gz -2 ${i%.cram}_2.fastq.gz - ; done
