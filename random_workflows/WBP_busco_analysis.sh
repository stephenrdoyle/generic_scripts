#WBP BUSCO / CEGMA ANALYSIS

wget -r "ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS11/species/" -P ./ -A "*genomic.fa.gz"
gunzip *.gz

wget -r "ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/" -P ./ -A "*protein.fa.gz"
gunzip *.gz




# BUSCO ARRAY
export PATH="/nfs/users/nfs_s/sd21/lustre118_link/software/hmmer-3.1b2-linux-intel-x86_64/binaries:$PATH"
export PATH="/nfs/users/nfs_s/sd21/software/augustus-3.2.1/bin:$PATH"
export PATH="/nfs/users/nfs_s/sd21/software/augustus-3.2.1/scripts:$PATH"
export AUGUSTUS_CONFIG_PATH="/nfs/users/nfs_s/sd21/software/augustus-3.2.1/config"
export LD_LIBRARY_PATH=/software/gcc-4.9.2/lib64/:$LD_LIBRARY_PATH

num=1
while read nematode; do
echo -e "~sd21/lustre118_link/software/ASSEMBLY_QC/busco_v3/scripts/run_BUSCO.py --in ${nematode}.WBPS11.genomic.fa --out ${nematode}_BUSCO_out --mode genome --lineage_path /nfs/users/nfs_s/sd21/databases/busco/metazoa_odb9 --species caenorhabditis --cpu 12 --force --restart --tmp_path ${nematode}.tmp --long --blast_single_core" > run_BUSCO.nematode.tmp.${num} ;
let num++ ; done < nematode.list

chmod a+x run_BUSCO.nematode.tmp.*
JOBS=$( ls -1 run_BUSCO.nematode.tmp.* | wc -l )
ID="U$(date +%s)"

bsub -q long -R'span[hosts=1] select[mem>20000] rusage[mem=20000]' -n 12 -M20000 -J BUSCO_nematode_${ID}_[1-$JOBS] -e BUSCO_nematode_${ID}_[1-$JOBS].e -o BUSCO_nematode_${ID}_[1-$JOBS].o "./run_BUSCO.nematode.tmp.\$LSB_JOBINDEX"


num=1
while read platy; do
echo -e "~sd21/lustre118_link/software/ASSEMBLY_QC/busco_v3/scripts/run_BUSCO.py --in ${platy}.WBPS11.genomic.fa --out ${platy}_BUSCO_out --mode genome --lineage_path /nfs/users/nfs_s/sd21/databases/busco/metazoa_odb9 --species schistosoma2 --cpu 12 --force --restart --tmp_path ${platy}.tmp --long --blast_single_core" > run_BUSCO.platy.tmp.${num} ;
let num++ ; done < platy.list

chmod a+x run_BUSCO.platy.tmp.*
JOBS=$( ls -1 run_BUSCO.platy.tmp.* | wc -l )
ID="U$(date +%s)"

bsub -q long -R'span[hosts=1] select[mem>20000] rusage[mem=20000]' -n 12 -M20000 -J BUSCO_platy_${ID}_[1-$JOBS] -e BUSCO_platy_${ID}_[1-$JOBS].e -o BUSCO_platy_${ID}_[1-$JOBS].o "./run_BUSCO.platy.tmp.\$LSB_JOBINDEX"



# make plots

sed -e 's/Complete/1/g' -e 's/Duplicated/2/g' -e 's/Fragmented/0.5/g' -e 's/Missing/0/g' clade1.txt > clade1_coded.txt
sed -e 's/Complete/1/g' -e 's/Duplicated/2/g' -e 's/Fragmented/0.5/g' -e 's/Missing/0/g' clade3.txt > clade3_coded.txt
sed -e 's/Complete/1/g' -e 's/Duplicated/2/g' -e 's/Fragmented/0.5/g' -e 's/Missing/0/g' clade4.txt > clade4_coded.txt
sed -e 's/Complete/1/g' -e 's/Duplicated/2/g' -e 's/Fragmented/0.5/g' -e 's/Missing/0/g' clade5.txt > clade5_coded.txt
sed -e 's/Complete/1/g' -e 's/Duplicated/2/g' -e 's/Fragmented/0.5/g' -e 's/Missing/0/g' platy.txt > platy_coded.txt


R-3.4.0
library(gplots)

clade1<-as.matrix(read.table("clade1_coded.txt",header=T,row.names=1))
clade3<-as.matrix(read.table("clade3_coded.txt",header=T,row.names=1))
clade4<-as.matrix(read.table("clade4_coded.txt",header=T,row.names=1))
clade5<-as.matrix(read.table("clade5_coded.txt",header=T,row.names=1))
platy<-as.matrix(read.table("platy_coded.txt",header=T,row.names=1))

heatmap.2(clade1,trace="none",margins=c(16,3),labRow="",cexCol=0.75)
heatmap.2(clade3,trace="none",margins=c(16,3),labRow="",cexCol=0.75)
heatmap.2(clade4,trace="none",margins=c(16,3),labRow="",cexCol=0.75)
heatmap.2(clade5,trace="none",margins=c(16,3),labRow="",cexCol=0.75)
heatmap.2(platy,trace="none",margins=c(16,3),labRow="",cexCol=0.75)
