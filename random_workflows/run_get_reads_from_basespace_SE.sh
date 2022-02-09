#!/usr/bin/env bash

# ========================================================================================
# run_get_reads_from_basespace_SE
#
# Get single end fastq reads from basespace
#
# Reference:
# Software:
#
# Usage: ~sd21/bash_scripts/run_get_reads_from_basespace_SE <INPUT>
#
# @authors
# Stephen Doyle <sd21@sanger.ac.uk>
#
# ----------------------------------------------------------------------------------------




rm *tmp
touch get_reads.tmp

while read name R1_id; do \
echo -e "\
wget -nc -O ${name}_R1.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/"${R1_id}"/content?access_token=7a496770d4294732b7011a3de5ee5e7d'\n\
" >> get_reads.tmp; done < ${1}

chmod a+x get_reads.tmp
bsub.py 1 get_reads ./get_reads.tmp


# example "INPUT" file - a tab sepearet file containing a sample name and run ID
# Tc_RAD_IDX01_run1	15988991
# Tc_RAD_IDX02_run1	15988992
# Tc_RAD_IDX03_run1	15988993
# Tc_RAD_IDX07_run1	15988994
# Tc_RAD_IDX09_run1	15988995
# Tc_RAD_IDX11_run1	15988996
# Tc_RAD_IDX01_run2	16002030
# Tc_RAD_IDX02_run2	16002031
# Tc_RAD_IDX03_run2	16002032
# Tc_RAD_IDX07_run2	16002033
# Tc_RAD_IDX09_run2	16002034
# Tc_RAD_IDX11_run2	16002035
