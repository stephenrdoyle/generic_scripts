


cd /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/ORTHOLOGY/SYNTENY

ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/SELECTION/PROTEIN_FASTAs/Results_Jan25/Orthologues_Jan25/Orthologues/Orthologues_hc_V4.proteins.unique/hc_V4.proteins.unique__v__hp.proteins.unique.csv
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/SELECTION/PROTEIN_FASTAs/Results_Jan25/Orthologues_Jan25/Orthologues/Orthologues_hc_V4.proteins.unique/hc_V4.proteins.unique__v__ce.proteins.unique.csv

ln -fs /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/TRANSCRIPTOME_CURATION/HCON_V4_WBP11plus_190125.ips.gff3 hc.gff3
ln -fs /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/TRANSCRIPTOME_SUMMARY_STATS/CELEGANS/wormbase.20240.complete.gff3 ce.gff3
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS12/species/haemonchus_placei/PRJEB509/haemonchus_placei.PRJEB509.WBPS12.annotations.gff3.gz 
gunzip haemonchus_placei.PRJEB509.WBPS12.annotations.gff3.gz
ln -sf haemonchus_placei.PRJEB509.WBPS12.annotations.gff3 hp.gff3




# make lists of 1 to 1 orthologs
awk 'NF==3 {print $1,$2,$3}' OFS="\t" hc_V4.proteins.unique__v__hp.proteins.unique.csv > hc_hp_1to1.list
#> 9970 hc_hp_1to1.list

awk 'NF==3 {print $2,$3}' OFS="\t" hc_V4.proteins.unique__v__ce.proteins.unique.csv > hc_ce_1to1.list
#> 7361 hc_ce_1to1.list


awk -F '[\t;]' '$3=="gene" {print $1,$4,$5,$7,$10}' OFS="\t" ce.gff3 | sed 's/Name=//g' > ce.gff.coords
awk -F '[\t;]' '$3=="gene" {print $1,$4,$5,$7,$10}' OFS="\t" hp.gff3 | sed 's/Name=//g' > hp.gff.coords
awk -F '[\t;]' '$3=="mRNA" {print $1,$4,$5,$7,$9}' OFS="\t" hc.gff3 | sed 's/ID=//g' > hc.gff.coords

while read HC_ID CE_ID; do COORDS1=$( grep ${HC_ID} hc.gff.coords); COORDS2=$( grep ${CE_ID} ce.gff.coords); echo -e "${COORDS1}\t${COORDS2}"; done <hc_ce_1to1.list > hc_ce_1to1.hccoords
sort -k1,1 -k2,2n hc_ce_1to1.hccoords | awk '{if($2+0==$2 && $7+0==$7) print $0}' > hc_ce_1to1.coords.chrsorted



# count the number of concordant vs non concordant chromosomal orthologs

cut -f1,6 hc_ce_1to1.coords.chrsorted | sort | uniq -c | awk '{print $2,$3,$1}' OFS="\t" > hc_ce_orthologcount_perchr.txt





# load
R-3.5.0
library(ggplot2)
library(patchwork)

# plot number of genes per chromosome
data<-read.table("hc_ce_orthologcount_perchr.txt",header=F)
plot1<- ggplot(data,aes(data$V2,data$V3,group=data$V1))+
	geom_bar(aes(fill=data$V2), stat = "identity", position = "dodge")+
	theme_bw()+facet_grid(.~data$V1)

# plot scatter between hc and ce orthologs by position, by hc chromosome, coloured by ce chromosome
data2<-read.table("hc_ce_1to1.coords.chrsorted",header=F)
plot2<-ggplot(data2,aes(data2$V2,data2$V7,group=data2$V1))+
	geom_point(aes(col=data2$V6))+
	facet_grid(.~data2$V1)+
	theme_bw()

plot1 + plot2 + plot_layout(ncol=1)
