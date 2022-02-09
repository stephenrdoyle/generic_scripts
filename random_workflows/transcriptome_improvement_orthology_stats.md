# Generating stats to show genome improvement - orthology




## Orthology <a name="orthology"></a>
```bash
# working dir:
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/TRANSCRIPTOME_CURATION/CURRENT/ORTHOLOGY

# get data
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS11/species/haemonchus_placei/PRJEB509/haemonchus_placei.PRJEB509.WBPS11.protein.fa.gz
gunzip haemonchus_placei.PRJEB509.WBPS11.protein.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS11/species/caenorhabditis_elegans/PRJNA13758/caenorhabditis_elegans.PRJNA13758.WBPS11.protein.fa.gz
gunzip caenorhabditis_elegans.PRJNA13758.WBPS11.protein.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS11/species/haemonchus_contortus/PRJNA205202/haemonchus_contortus.PRJNA205202.WBPS11.protein.fa.gz
gunzip haemonchus_contortus.PRJNA205202.WBPS11.protein.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS10/species/haemonchus_contortus/PRJEB506/haemonchus_contortus.PRJEB506.WBPS10.protein.fa.gz
gunzip haemonchus_contortus.PRJEB506.WBPS10.protein.fa.gz



ln -s ../UPDATED_annotation.gff3 new_annotation.gff
gffread new_annotation.gff -g /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/REF/HAEM_V4_final.chr.fa -y new_annotation.proteins.fa


# curate data for input into orthofinder
# get one coding sequence per gene - certainly Ce and HcV4 has multiple  isoforms, and therefore multiple coding sequneces per gene
fastaq to_fasta -l0 caenorhabditis_elegans.PRJNA13758.WBPS11.protein.fa caenorhabditis_elegans.PRJNA13758.WBPS11.protein.fa2
fastaq to_fasta -l0 haemonchus_placei.PRJEB509.WBPS11.protein.fa haemonchus_placei.PRJEB509.WBPS11.protein.fa2
fastaq to_fasta -l0 haemonchus_contortus.PRJEB506.WBPS10.protein.fa haemonchus_contortus.PRJEB506.WBPS10.protein.fa2
fastaq to_fasta -l0 haemonchus_contortus.PRJNA205202.WBPS11.protein.fa haemonchus_contortus.PRJNA205202.WBPS11.protein.fa2
fastaq to_fasta -l0 new_annotation.proteins.fa  new_annotation.proteins.fa2

# fix fasta header - likely only c elegans that had excessive information in the header, but fixed to make them all consistent
awk '{if($1 ~ /^>/) print ">"$3,$2; else print $0}' caenorhabditis_elegans.PRJNA13758.WBPS11.protein.fa2 > ce.proteins.fa
awk '{if($1 ~ /^>/) print ">"$3,$2; else print $0}' haemonchus_placei.PRJEB509.WBPS11.protein.fa2 > hp.proteins.fa
awk '{if($1 ~ /^>/) print ">"$3,$2; else print $0}' haemonchus_contortus.PRJEB506.WBPS10.protein.fa2 > hc_V1.proteins.fa
awk '{if($1 ~ /^>/) print ">"$3,$2; else print $0}' haemonchus_contortus.PRJNA205202.WBPS11.protein.fa2 > hc_McM.proteins.fa
cat new_annotation.proteins.fa2 > hc_new_annotation.proteins.fa

#remove stop codons from hc_V4
sed -i 's/[A-Z]*\.\.*[A-Z]//g' hc_new_annotation.proteins.fa
sed -i 's/\.$//g' hc_new_annotation.proteins.fa

# make unique gene lists
grep ">" ce.proteins.fa | cut -f 1 -d " " | sort | uniq > unique.Ce_genes.list
#> 20208 unique.Ce_genes.list

grep ">" hp.proteins.fa | cut -f 1 -d " " | sort | uniq > unique.Hp_genes.list
#> 21928 unique.Hp_genes.list

grep ">" hc_McM.proteins.fa | cut -f 1 -d " " | sort | uniq > unique.Hc_McM_genes.list
#> 23610 unique.Hc_McM_genes.list

grep ">" hc_V1.proteins.fa | cut -f 1 -d " " | sort | uniq > unique.Hc_V1_genes.list
#> 21869 unique.Hc_V1_genes.list


grep ">" hc_new_annotation.proteins.fa | sort | uniq | sed 's/>//g' | sed 's/transcript://g' | sort | uniq | cut -c-13 | sort | uniq > unique.hc_new_annotation_genes.list
#> 19610 unique.hc_new_annotation_genes.list


while read -r gene; do grep -m1 -A1 ${gene} ce.proteins.fa; done < unique.Ce_genes.list > ce.proteins.unique.fa &
while read -r gene; do grep -m1 -A1 ${gene} hp.proteins.fa; done < unique.Hp_genes.list > hp.proteins.unique.fa &
while read -r gene; do grep -m1 -A1 ${gene} hc_McM.proteins.fa; done < unique.Hc_McM_genes.list > hc_McM.proteins.unique.fa &
while read -r gene; do grep -m1 -A1 ${gene} hc_V1.proteins.fa; done < unique.Hc_V1_genes.list > hc_V1.proteins.unique.fa &
while read -r gene; do grep -m1 -A1 ${gene} hc_new_annotation.proteins.fa; done < unique.hc_new_annotation_genes.list > hc_new_annotation.proteins.unique.fa &


# run OrthoFinder
#--- setup data
mkdir PROTEIN_FASTAs
mv *.unique.fa PROTEIN_FASTAs
cd PROTEIN_FASTAs
for i in *.fa; do cut -f1 -d " " $i | sed '/^$/d' > tmp; mv tmp $i; done
cd ../

# run OF
bsub.py --queue small --threads 20 20 orthofinder_2.5.2 "/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/OrthoFinder/./orthofinder -t 20 -a 20 -S diamond -M msa -A mafft -f PROTEIN_FASTAs/"




# extract gene start/end and cds start end data for calculating UTR lengths
> coding.length
for GENE in ` awk -F '[\t]' '$3=="gene" {print $9}' UPDATED_annotation.gff3 | grep -Po "Name=HCON_........" | sed 's/Name=//g'`; do
     chr=$( grep ${GENE} UPDATED_annotation.gff3 | grep -v "#" |  head -n1 | cut -f1 );
     gene_start=$(grep ${GENE} UPDATED_annotation.gff3 | grep -v "#" | awk '$3=="gene" {print}' | head -n1 | awk '{print $4}') ;
     gene_end=$(grep ${GENE} UPDATED_annotation.gff3 | grep -v "#" | awk '$3=="gene" {print}' | head -n1 | awk '{print $5}');
     cds1=$(grep ${GENE} UPDATED_annotation.gff3 | grep -v "#" | awk '$3=="CDS" {print}' | head -n1 | awk '{print $4}');
     cds2=$(grep ${GENE} UPDATED_annotation.gff3 | grep -v "#" | awk '$3=="CDS" {print}' | tail -n1 | awk '{print $5}'); strand=$(grep ${GENE} UPDATED_annotation.gff3 | grep -v "#" | head -n1 | cut -f7 );
     echo -e "$chr\t$gene_start\t$gene_end\t$cds1\t$cds2\t$GENE\t.\t$strand" >> coding.length; done &

# calculated difference between gene and cds starts and ends - note this doesnt take into account direction, so not specifically 5 and 3 in order, just a difference
cat coding.length | awk '{print $4-$2, $3-$5}' OFS="\t" > utr.lengths

awk '{if($1>0) print }' utr.lengths | wc -l
#> 12807

awk '{if($2>0) print }' utr.lengths | wc -l
#> 12938

19612 utr.lengths
39224 total UTRs

25745 UTRs  annotated =  65.6%


bedtools-2 subtract -s -a genes.bed -b coding.length > utrs.bed

awk '{print $3-$2}' utrs.bed > utr.length

#
ls -lrt







R

library(tidyverse)

data <- read.delim("WBPS15.assemblystats.data")


ggplot(data) +
     geom_point(aes(log10(N50),log10(N50n),size=log2(total_length),col=log2(Gaps))) +
     theme_bw() +
     labs(x="Genome contiguity (log10[N50])", y="Genome fragmentation (log10[N50n])", size = "Genome size (log2[total length])", Colour="Number of gaps (log2[gaps])")
