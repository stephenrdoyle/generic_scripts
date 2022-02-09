#!/usr/bin/env bash

# ========================================================================================
# run_minimap2syntenyplot.sh
#
# Quick synteny plots comparing two genomes using minimap
#
# Reference:
# Software: https://github.com/zeeev/minimap
#
# Usage: ~sd21/bash_scripts/run_minimap2syntenyplot.sh <TARGET> <QUERY>
#
# @authors
# Stephen Doyle <sd21@sanger.ac.uk>
#
# ----------------------------------------------------------------------------------------


target=$1
query=$2

# run minimap
/nfs/users/nfs_s/sd21/lustre118_link/software/GENOME_ASSEMBLY/minimap/minimap $target $query > ${target%.fa*}_vs_${query%.fa*}.mini

# parse minimap output and plot comparison
cat ${target%.fa*}_vs_${query%.fa*}.mini | /nfs/users/nfs_s/sd21/lustre118_link/software/GENOME_ASSEMBLY/minimap/utils/bin/layout > ${target%.fa*}_vs_${query%.fa*}.layout.txt

/software/R-3.5.0/bin/Rscript --vanilla /nfs/users/nfs_s/sd21/lustre118_link/software/GENOME_ASSEMBLY/minimap/utils/plot/plotLayout.R -f ${target%.fa*}_vs_${query%.fa*}.layout.txt -p ${target%.fa*}_vs_${query%.fa*}.miniplot.pdf
