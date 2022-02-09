#!/usr/bin/env bash

# ========================================================================================
# run_chromosomer_scaffolder.sh
#
# Scaffold a draft genome assembly against a complete reference assembly.
#
# Reference: https://doi.org/10.1186/s13742-016-0141-6
# Software: https://github.com/gtamazian/chromosomer
#
# Usage: ~sd21/bash_scripts/run_chromosomer_scaffolder.sh <PREFIX> <REFERENCE_ASSEMBLY> <DRAFT_ASSEMBLY>
#
# @authors
# Stephen Doyle <sd21@sanger.ac.uk>
#
# ----------------------------------------------------------------------------------------


PREFIX=$1
REF=$2
FRAG_ASSEMBLY=$3


# chromosomer is installing in the anaconda bin. Need to load it into my path
export PATH="/nfs/users/nfs_s/sd21/lustre118_link/software/anaconda2/bin:$PATH"



# step 1 - calculate the sequence lengths of the fragmented assembly
chromosomer fastalength ${FRAG_ASSEMBLY} ${PREFIX}_fragments.length

# step 2 - make a blast database of the reference assembly
makeblastdb -in ${REF} -dbtype nucl -out ${PREFIX}


# step 3 - blast the fragmented assembly against reference assembly
blastn -query ${FRAG_ASSEMBLY} -db ${PREFIX} -outfmt 6 -out ${PREFIX}_fragment_alignments.txt

# step 4 - make a fragment map, introducing 100bp gaps where needed. This can be changed
chromosomer fragmentmap ${PREFIX}_fragment_alignments.txt 100 ${PREFIX}_fragments.length ${PREFIX}_fragment_map.txt


# step 5 - use the fragment map to output a new genome
chromosomer assemble ${PREFIX}_fragment_map.txt ${FRAG_ASSEMBLY} ${PREFIX}_assembled_chromosomes.fa
