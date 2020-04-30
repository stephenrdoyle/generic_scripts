#!/usr/bin/env bash

# ========================================================================================
# run_SSR_microsat_pipeline.sh
#
# Tool to identify microsatellites in short read sequencing data
#
# Reference: https://academic.oup.com/jhered/article/104/6/881/798273
# Software: https://pubs.usgs.gov/ds/778/
#
# Usage: ~sd21/bash_scripts/run_SSR_microsat_pipeline.sh <PREFIX> <READ_1> <READ_2>
#
# @authors
# Stephen Doyle <sd21@sanger.ac.uk>
#
# ----------------------------------------------------------------------------------------



prefix=$1
read1=$2
read2=$3

SSR_dir=/nfs/users/nfs_s/sd21/sd21_lustre_link/software/MICROSATELLITES/SSR_pipeline_src_and_docs


export PATH=$PATH:/nfs/users/nfs_s/sd21/sd21_lustre_link/software/MICROSATELLITES/SSR_pipeline_src_and_docs


echo -e "\
nSets=1\n\
#set1:\n\
#comments can be placed on a line after a '#' sign\n\
\n\
read_set_1 = ['$read1']\n\
read_set_2 = ['$read2']\n\
#put an 'r' in front of the path name on Windows machines\n\
###UPDATE 'path_to_files' FOR YOUR COMPUTER\n\
path_to_files = r'$PWD'\n\
generic_output_name = '${prefix}_SSR_out'\n\
keepZipped = 'T'\n\
alignment_Parameters = [10, 0.25, 70]\n\
microsat_search_params = [(2, 7, 40),\n\
                          (3, 6, 40),\n\
                          (4, 5, 40),\n\
                          (5, 4, 40),\n\
                          (6, 4, 40)]\n\
#end of information for a set of files is indicated by a semicolon (';')\n\
;" > $prefix.config_file

python $SSR_dir/SSR_pipeline.py $prefix.config_file
