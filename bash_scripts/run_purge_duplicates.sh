#!/usr/bin/env bash

# ========================================================================================
# run_purge_duplicates.sh
#
# Remove duplicate haplotypes and heterozygous overlaps from a PacBio assembly.
#
# Reference: https://www.biorxiv.org/content/10.1101/729962v1.full.pdf
# Github: https://github.com/dfguan/purge_dups
#
# Usage: ~sd21/bash_scripts/run_purge_duplicates.sh <PREFIX> <REFERENCE_FASTA> <Pacbio File of file names>
#
# @authors
# Stephen Doyle <sd21@sanger.ac.uk>
#
# ----------------------------------------------------------------------------------------


PREFIX=$1
REFERENCE=$2
PB_FOFNs=$3

pd_config.py ${REFERENCE} ${PB_FOFNs}

run_purge_dups.py config.json ~sd21/lustre118_link/software/GENOME_IMPROVEMENT/purge_dups/bin ${PREFIX}
