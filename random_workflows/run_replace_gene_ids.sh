# Replacing gene names in GFF with Species ID and incremental gene specific ID

gff=$1
species_prefix=$2

# TESTING
# wd: /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/BRAKER_CHR/POST_BRAKER
#gff=augustus.filtered.gff3
#species_prefix=HCONnew

# copy gff to a tmp file - just to avoid doing something silly like overwriting original
cat ${gff} | grep -v "#" > ${gff}.tmp


# STEP 1 - generate using IDs for each gene
# - extract gene name from lines in gff matching gene|GENE
# - remove surrounding characters
# - sort by gene ID and remove duplicates
# - print new gene id (species prefix with 8 digit ID that increases by 10 for each gene), and old gene id

awk -F'[\t;]' '{if($3=="gene" || $3=="GENE") print $9}' ${gff}.tmp | sed -e 's/ID=//g' -e 's/\;//g' | sort -V | uniq | awk -v species_prefix="$species_prefix" '{fmt=species_prefix"%08d\t%s\n"; printf fmt,NR*10,$0}' > genes_renames.list

# pasa dependent fix
cat genes_renames.list | sed 's/gene/mRNA/g' >  mRNA_renames.list


# STEP 2 - replace gene IDs
# - read new / old gene IDs from above
# - replace IDs, recognising differences in placement of gene, mRNA, and progeny
# NOTE: match - equal sign at beginning, followed by one of [.;$] where $ is the end of the line, which is what AUGUSTSUS / BRAKER produce - this might have to be tweaked for other outputs
# NOTE: eg.  =XXX[. ;]

#while read new_gene old_gene; do 
#grep "$old_gene[.;$]" ${gff}.tmp | sed -i -e "s/=${old_gene}[.]/=${new_gene}\./g" -e "s/=${old_gene}[;]/=${new_gene}\;/g" -e "s/=${old_gene}$/=${new_gene}/g" -e "s/Name=.*$/Name=${new_gene}/g"; 
#done < genes_renames.list > ${gff%.gff*}.renamed.tmp


# while read new_mRNA old_mRNA; do 
# grep "$old_mRNA[.;$]" ${gff%.gff*}.genesrenamed.tmp | grep -v "gene" | sed -e "s/=${old_mRNA}[.]/=${new_mRNA}\./g" -e "s/=${old_mRNA}[;]/=${new_mRNA}\;/g" -e "s/=${old_mRNA}$/=${new_mRNA}/g" ; 
# done < mRNA_renames.list > ${gff%.gff*}.allrenamed.tmp


while read new_gene old_gene; do 
sed -i -e "s/=${old_gene}[.]/=${new_gene}\./g" -e "s/=${old_gene}[;]/=${new_gene}\;/g" -e "s/=${old_gene}$/=${new_gene}/g" -e "s/Name=.*$/Name=${new_gene}/g" ${gff}.tmp; 
done < genes_renames.list



while read new_mRNA old_mRNA; do 
sed -i -e "s/=${old_mRNA}[.]/=${new_mRNA}\./g" -e "s/=${old_mRNA}[;]/=${new_mRNA}\;/g" -e "s/=${old_mRNA}$/=${new_mRNA}/g" ${gff}.tmp; 
done < mRNA_renames.list

# STEP 3 - clean up

#echo "##gff-version 3" > ${species_prefix}.renamed.gff3; cat ${gff%.gff*}.allrenamed.tmp | sort -k1,1 -k4,4n >> ${species_prefix}.renamed.gff3

#rm *tmp*