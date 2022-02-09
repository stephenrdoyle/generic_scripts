
# Transfer mRNA NAMES to IDs for Apollo generated GFF files

- Apollo give as unique number/letter sting as an ID to mRNA trascnripts, which is the parent ID for all subsequent features
- It gives the GENE name as mRNA name with unique transcript ID
    eg
      Gene: HCON_00000001
      mRNA: HCON_00000001-00001

- Need to do two things
    - give mRNA transcripts a unquique ID if they dont already have done
    - get mRNA ID, and replcate it with mRNA names






```shell
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/TRANSCRIPTOME_CURATION
mkdir TESTING
cd TESTING

ln -s ../HCON_V4_WBP11plus_190114.gff3
head -n 500 HCON_V4_WBP11plus_190114.gff3 > testing.gff
cd TESTING

# get mRNAs IDs and NAMES
awk '$3=="mRNA" {print $0}' OFS="\t"  testing.gff | sed -e 's/Note=Manually dissociate transcript from gene;//g' | cut -f3,5 -d ";" | sed -e 's/ID=//g' -e 's/;Name=/\t/g' > mRNA_IDs_NAMEs.txt

# remove transcript extensions
cat mRNA_IDs_NAMEs.txt | awk  '{ $2 = substr($2, 1,13); print }' | head


cat mRNA_IDs_NAMEs.txt | awk  '{ $2 = substr($2, 1,13); print }' OFS="\t" | sort -k2 > mRNA_IDs_NAMEs_trimmed-sorted.txt

# make a unique list
cut -f2  mRNA_IDs_NAMEs_trimmed-sorted.txt | uniq > mRNA_IDs_NAMEs_unique.txt

# add unique ID to each transcript
>mRNA_IDs_NAMEs_transcriptIDs.txt
while read NAME; do grep -w ${NAME} mRNA_IDs_NAMEs_trimmed-sorted.txt | cat -n | awk '{print $2,$3,$3"-0000"$1}' OFS="\t" >> mRNA_IDs_NAMEs_transcriptIDs.txt; done < mRNA_IDs_NAMEs_unique.txt &


while read OLD GENE NEW; do sed -i "/$OLD/ s//$NEW/g" testing.gff; done < mRNA_IDs_NAMEs_transcriptIDs.txt &

# run the real substitution
while read OLD GENE NEW; do sed -i "/$OLD/ s//$NEW/g" HCON_V4_WBP11plus_190114.gff3; done < mRNA_IDs_NAMEs_transcriptIDs.txt &

```
