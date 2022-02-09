# Repeat analysis of S. mansoni
- comparison of repeats from Grevelding paper to V9 assembly
- trying to reconcile SBs repeats with Greveldings repeats


```bash
# working directory
cd /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/REPEATS

# copied W-specific repeats from Supplementary figure 2 in
# https://doi.org/10.1093/gbe/evab204 into a file called "w-repeats.fa"

# had to fix the file to remove spaces
sed -i 's/ //g' w-repeats.fa


# get both V7 and V9 genomes

cp ../../papers/genome_improvement/smansoni/WBP_V7/schistosoma_mansoni.PRJEA36577.WBPS15.genomic.fa sm_v7.fa
cp ../../papers/genome_improvement/smansoni/buddenborg_2021/SM_V9_ENA.genomic.fa sm_v9.fa


# make a list of w-repeats
grep ">" w-repeats.fa | sed 's/>//g' > w-repeats.list.txt



# run nucmer to identify repeats in both V7 and V9 genomes
nucmer --prefix w-repeats-v7 --maxmatch w-repeats.fa sm_v7.fa
show-coords -lTHc -I 80 w-repeats-v7.delta > w-repeats-v7.coords

nucmer --prefix w-repeats-v9 --maxmatch w-repeats.fa sm_v9.fa
show-coords -lTHc -I 80 w-repeats-v9.delta > w-repeats-v9.coords

# extract some basic count information of repeat hits in both V7 and v9 genomes
printf "REPEAT\tV7_count\tV9_count\n" > w-repeats.counts.txt
while read REPEAT; do
  # v7 count
  v7_count=$(grep ${REPEAT} w-repeats-v7.coords | sort -k3 | wc -l)
  # v9 count
  v9_count=$(grep ${REPEAT} w-repeats-v9.coords | sort -k3 | wc -l)
  #echo -e "${REPEAT}\t${v7_count}\t${v9_count}";
  printf "${REPEAT}\t${v7_count}\t${v9_count}\n" >> w-repeats.counts.txt;
  done < w-repeats.list.txt



# counts of repeat per chromosome
#--- V7
while read REPEAT; do
  grep "$REPEAT" w-repeats-v9.coords | awk '{print $12,$13}' OFS="\t" | sort | uniq -c ;
  done < w-repeats.list.txt

#--- V9
while read REPEAT; do
  grep "$REPEAT" w-repeats-v9.coords | awk '{print $12,$13}' OFS="\t" | sort | uniq -c ;
  done < w-repeats.list.txt


```



```R
library(tidyverse)


# load the V9 data
a<- read.table("w-repeats-v9.coords")

# get rid of the highly abundant repeat for on all chromosomes
b<- filter(a,V12!="W36.2_2_SM_V7_W001_335_bp")
ggplot(b,aes(V3,fill=V12))+geom_histogram(bins=1000)+facet_grid(V13~.)

# pull of repeats only on WSR
c <- filter(b, V13=="SM_V9_WSR")
ggplot(c,aes(V3,fill=V12))+geom_histogram(bins=1000)+facet_grid(V13~.)

```
