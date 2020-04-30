#!/usr/bin/env bash

# ========================================================================================
# run_gmato_microsat_pipeline.sh
#
# Wrapper script to use the tool GMATo to mine microsatellite sequences from a genome reference file
# Simple takes a fasta file containing th reference sequence as input, and returns the type and frequency of microsatellites
# as well as their position in the genome. This information could be used to design new microsatellite genetic markers
#
# Reference:
# Software: https://sourceforge.net/projects/gmato/files/GMAToV1.2/
# Manual: https://sourceforge.net/projects/gmato/files/GMAToV1.2/help_20130106.pdf
#
# Usage: ~sd21/bash_scripts/run_SSR_microsat_pipeline.sh <PREFIX> <REFERENCE>
#
# @authors
# Stephen Doyle <sd21@sanger.ac.uk>
#
# ----------------------------------------------------------------------------------------

prefix=$1
reference=$2

export PATH=$PATH:/nfs/users/nfs_s/sd21/sd21_lustre_link/software/MICROSATELLITES/gmatoV1.2Build20130106
### note: this is setup to run from my (sd21) path. This is fine providing you have access to the path.


mkdir ${prefix}_gmato_out

gmat.pl -r 5 -m 2 -x 10 -s 0 -i $reference 2&>${prefix}.gmato.log

mv *fms* ${prefix}_gmato_out/
mv *ssr* ${prefix}_gmato_out/
