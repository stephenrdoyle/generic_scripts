# example of converting orthofinder output to upsetR to compare shared 1-to-1 and >1 orthogroups 


# you need
#--- "Orthogroups.GeneCount.tsv", which looks something like:

sd21@farm5-head2:~/lustre118_link/hc/GENOME/TRANSCRIPTOME/TRANSCRIPTOME_CURATION/20200421/ORTHOFINDER/PROTEIN_FASTAs/OrthoFinder/Results_Apr30_1/Orthogroups$ head Orthogroups.GeneCount.tsv
Orthogroup	HCON_WBPS15.proteins.unique	ce.proteins.unique	hc_V1.proteins.unique	hc_V4.proteins.unique	Total
OG0000000	0	0	116	0	116
OG0000001	0	87	0	0	87
OG0000002	0	87	0	0	87
OG0000003	0	78	1	0	79
OG0000004	2	0	71	2	75
OG0000005	0	73	0	0	73
OG0000006	0	65	0	0	65
OG0000007	0	59	0	0	59
OG0000008	0	59	0	0	59



#--- UpSetR (https://www.rdocumentation.org/packages/UpSetR/versions/1.4.0)


#---------------------------------------------------------------------------------------------
# Step 1: get one-to-ones

head -n 1 Orthogroups.GeneCount.tsv | awk '{print $1,$2,$3,$4,$5}' OFS="\t" > orthosingles.data
awk '{if($2<=1 && $3 <=1 && $4<=1 && $5<=1) print $1,$2,$3,$4,$5}' OFS="\t" Orthogroups.GeneCount.tsv >> orthosingles.data

# note1: only want the first 5 columns of the table above, hence the $1,$2,$3,$4,$5 - you will need to look and adjust based on your data
# note2: this bit - $2<=1 && $3 <=1 && $4<=1 && $5<=1 - selects for data where the number of genes is 0 or 1 - hence missing or in single copy - again, you will need to adjust this based on the columns you have




# Step 2: get rest, convert to right format
awk '{print $1,$2,$3,$4,$5}' OFS="\t" Orthogroups.GeneCount.tsv | awk 'NR>1 {for(i=2;i<=NF;i++)if($i>0)$i=1}1' OFS="\t" > orthogroups.data

# note: 'NR>1 {for(i=2;i<=NF;i++)if($i>0)$i=1}1'  simply converts any number, ie gene number, greater than 1 to a value of 1. Needed for upsetR



# load R
R
library(UpSetR)

singles <- read.table("orthosingles.data",header=T,comment.char="")
groups <- read.table("orthogroups.data",header=T,comment.char="")


pdf("one-to-one_orthogroups_plot.upsetr.pdf",height=5,width=10,useDingbats=FALSE)
upset(singles)
dev.off()

pdf("orthogroups_plot.upsetr.pdf",height=5,width=10,useDingbats=FALSE)
upset(groups)
dev.off()
