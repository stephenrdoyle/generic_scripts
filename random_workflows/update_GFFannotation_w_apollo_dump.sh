# update GFF with Apollo dump


# GFF fixes
- Apollo dumps
- for some reason, Apollo has a bug that splits a interpro records over multiple lines.
eg.
##gff-version 3
##sequence-region hcontortus_chr3_Celeg_TT_arrow_pilon 1 43560352
hcontortus_chr3_Celeg_TT_arrow_pilon	iris	gene	6120	8670	.	+	.	biotype=protein_coding;owner=irisadmin@local.host;description=UniProtKB/TrEMBL;ID=057ebf0f-52bc-4e95-9850-e23df56a0bff;date_last_modified=2019-02-27;Name=HCON_00064935;date_creation=2019-02-27
hcontortus_chr3_Celeg_TT_arrow_pilon	iris	mRNA	6120	8670	.	+	.	owner=irisadmin@local.host;Parent=057ebf0f-52bc-4e95-9850-e23df56a0bff;ID=aa8f17ce-add4-4491-b5e7-96660beeac84;date_last_modified=2019-02-27;Name=HCON_00064935;info=method:InterPro accession:IPR000306 description:FYVE zinc finger
method:InterPro accession:IPR011011 description:Zinc finger%2C FYVE/PHD-type
method:InterPro accession:IPR013083 description:Zinc finger%2C RING/FYVE/PHD-type
method:InterPro accession:IPR017455 description:Zinc finger%2C FYVE-related;date_creation=2019-02-27
hcontortus_chr3_Celeg_TT_arrow_pilon	iris	exon	6819	6998	.	+	.	Parent=aa8f17ce-add4-4491-b5e7-96660beeac84;ID=ec68fe21-08df-4b88-ae32-dcb3295257be;Name=ec68fe21-08df-4b88-ae32-dcb3295257be
hcontortus_chr3_Celeg_TT_arrow_pilon	iris	exon	7190	7357	.	+	.	Parent=aa8f17ce-add4-4491-b5e7-96660beeac84;ID=36373354-cb22-40da-b7e5-f2c87e01b3b2;Name=36373354-cb22-40da-b7e5-f2c87e01b3b2
hcontortus_chr3_Celeg_TT_arrow_pilon	iris	exon	6394	6483	.	+	.	Parent=aa8f17ce-add4-4491-b5e7-96660beeac84;ID=f65191f9-ec3a-4463-8b48-341ec70ca011;Name=f65191f9-ec3a-4463-8b48-341ec70ca011


- to fix this, need to run the following
```shell
gunzip Annotations.gff3.gz

cat Annotations.gff3 | while read LINE; do if [[ $LINE =~ (^hcontortus*|\#) ]]; then echo -ne "\n${LINE} "; else echo -n "${LINE} ";fi; done > Annotations.v2.gff3


```



# get mRNAs IDs and NAMES
awk '$3=="mRNA" {print $0}' OFS="\t" Annotations.v2.gff3 | sed -e 's/Note=Manually dissociate transcript from gene;//g' | cut -f3,5 -d ";" | sed -e 's/ID=//g' -e 's/;Name=/\t/g' > mRNA_IDs_NAMEs.txt

# remove Names from mRNA lines
cat Annotations.v2.gff3 | sed -e '/mRNA/ s/Name=.*;//' -e 's/owner=irisadmin@local.host;//' > Annotations.v2.2.gff3

# replace name with ID
cat Annotations.v2.2.gff3 | sed 's/ *$//' | awk -F '[\t; ]' '{if($3=="mRNA") print $0";"$10; else print $0}' | sed -e 's/ID=/Name=/2' > Annotations.v3.gff3



APOLLO_GFF=Annotations.v3.gff3

# CHECK - are the names the right lenght?
#--- count length
cat mRNA_IDs_NAMEs.txt | cut -f2 | awk '{ print length }' | sort | uniq -c
# find IDs with strang lengths, eg, should be 19, but print those with 20
cat mRNA_IDs_NAMEs.txt | cut -f2 | awk '{if(length!=19) print }'
# check IDs are unique
cat mRNA_IDs_NAMEs.txt | cut -f2 | sort | uniq -c | sort





# remove transcript extensions
cat mRNA_IDs_NAMEs.txt | awk  '{ $2 = substr($2, 1,13); print }' OFS="\t" | sort -k2 > mRNA_IDs_NAMEs_trimmed-sorted.txt

# make a unique list
cut -f2  mRNA_IDs_NAMEs_trimmed-sorted.txt | uniq > mRNA_IDs_NAMEs_unique.txt

# add unique ID to each transcript
>mRNA_IDs_NAMEs_transcriptIDs.txt
while read NAME; do grep -w ${NAME} mRNA_IDs_NAMEs_trimmed-sorted.txt | cat -n | awk '{print $2,$3,$3"-0000"$1}' OFS="\t" >> mRNA_IDs_NAMEs_transcriptIDs.txt; done < mRNA_IDs_NAMEs_unique.txt

# run the real substitution
export PATH="/nfs/users/nfs_s/sd21/lustre118_link/software/anaconda2/bin:$PATH"
# fsed - https://github.com/wroberts/fsed
awk '{print $1,$3}' OFS="\t" mRNA_IDs_NAMEs_transcriptIDs.txt > mRNA_IDs_NAMEs_transcriptIDs.2.txt

fsed --pattern-format=tsv --output ${APOLLO_GFF}.renamed  mRNA_IDs_NAMEs_transcriptIDs.2.txt ${APOLLO_GFF} &



# remove Names from mRNA lines

cat old_annotation.gff | sed -e 's/ *$//' -e '/mRNA/ s/Name=.*;//' -e 's/owner=irisadmin@local.host;//' -e 's/owner=sd21@sanger.ac.uk;//g'| awk -F '[\t;]' '{if($3=="mRNA") print $0";"$10; else print $0}' | sed -e 's/ID=/Name=/2' > old_annotation.v2.gff

APOLLO_GFF=Annotations.v3.gff3
OLD_GFF=old_annotation.v2.gff

#dump updated genes from old gff


# remove overlapping models which have been updated
bedtools intersect -s -v  -f 0.1 -b ${APOLLO_GFF}.renamed -a ${OLD_GFF} > ${OLD_GFF}.filtered

# make a list of updated genes
sed -e 's/Note=Manually.*.gene;//g' -e 's/Name=//g' ${APOLLO_GFF}.renamed | awk '{if($3~"gene") print $9}' | cut -d ";" -f3 | grep -v "PB" > updated_genes.list

sed -e 's/Note=Manually.*.gene;//g' -e 's/Name=//g' ${APOLLO_GFF}.renamed | awk '{if($3~"repeat") print $9}' | cut -d ";" -f3 | grep -v "PB" >> updated_genes.list



# make a list of original genes
awk '{if($3=="mRNA") print $9}' ${OLD_GFF}.filtered | sed -e 's/Note=Manually.*.gene;//g' -e 's/ID=//g' -e 's/Name=//g' | cut -d ";" -f2 | cut -c-13 | sort | uniq > original_genes.list
# --- from GAG output
#awk '{if($3=="mRNA") print $9}' ${OLD_GFF}.filtered | sed -e 's/Note=Manually dissociate transcript from gene;//g' -e 's/ID=//g' | cut -d ";" -f1 | cut -c-13 | sort | uniq > original_genes.list

# identify genes in original file that need to be kept as they have not been updated
diff -u <(sort updated_genes.list) <(sort original_genes.list) | grep "^+" | sed 's/+//g' > retained_original_genes.list


screen
APOLLO_GFF=Annotations.v3.gff3
OLD_GFF=old_annotation.v2.gff
mkdir NEW_MODELS

# make new gene models of retained genes
while read GENE; do
CHR=$( grep "${GENE}" ${OLD_GFF}.filtered | cut -f1 | head -n1);
>NEW_MODELS/${CHR}_${GENE}.model;
grep "${GENE}" ${OLD_GFF} >> NEW_MODELS/${CHR}_${GENE}.model;
echo "###" >> NEW_MODELS/${CHR}_${GENE}.model;
done < retained_original_genes.list &

# make new gene models of new curated genes from apollo
while read GENE; do
CHR=$( grep "${GENE}" ${APOLLO_GFF}.renamed | cut -f1 | head -n1);
>NEW_MODELS/${CHR}_${GENE}.model;
grep "${GENE}" ${APOLLO_GFF}.renamed >> NEW_MODELS/${CHR}_${GENE}.model;
echo "###" >> NEW_MODELS/${CHR}_${GENE}.model;
done < updated_genes.list &


# make sequence region file and collate genes per region

grep "sequence-region" ${APOLLO_GFF}  > sequence_regions.list

while read NAME CHR START END; do
>${CHR}.genes
echo -e $NAME"\t"$CHR"\t"$START"\t"$END >${CHR}.genes;
cat NEW_MODELS/${CHR}* | sort | uniq | sort -k1,1 -k4,4n >> ${CHR}.genes;
done < sequence_regions.list


# bring it all together into a multisequence merged GFF
echo -e "##gff-version 3" > UPDATED_annotation.gff3
cat *.genes >> UPDATED_annotation.gff3

sed -i -e 's/iris\t/WSI_sd21_apollo\t/g' -e 's/owner=irisadmin@local.host;//g' UPDATED_annotation.gff3

now=$(date +'%Y%m%d')
ln -s UPDATED_annotation.gff3 HCON_V4_curated_${now}.gff3
gag.py -f ../../../REF/HAEM_V4_final.chr.fa -g UPDATED_annotation.gff3 --fix_start_stop -o UPDATED_annotation_GAGout
# note - this will remove non-standard

# once completed,
```
