###


```bash
# download all genomes from WBP
wget -r "ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS16/species/" -P ./ -A "*genomic.fa.gz"
find . -name "*.gz" -exec mv {} ./ \;
gunzip *.gz

# generate assembly stats for all genomes
module load assembly-stats/1.0.1

assembly-stats -t *fa > WBP16_assembly_stats.txt


# clean it up to make it compatible with the quality stat data
cat WBP16_assembly_stats.txt | sed -e 's/.WBPS16.genomic.fa//g' -e 's/\.P/_P/g' | tr '[:upper:]' '[:lower:]' > tmp; mv tmp WBP16_assembly_stats.txt

```

```bash
# download BUSCO and CEGMA data from WBP
ls -1 *.genomic.fa | sed -e 's/.WBPS16.genomic.fa//g' -e 's/\./_/g' | tr '[:upper:]' '[:lower:]' |\
     while read -r NAME; do
     DATA=$( wget -q --header='Content-type:application/json' "https://parasite.wormbase.org/rest-16/info/quality/${NAME}" -O - | jq --sort-keys --tab -c ".busco" | sed -e 's/"//g' -e 's/,/\t/g' -e 's/:/\t/g' -e 's/{//g' -e 's/}//g' ) ;
     echo -e ${NAME}"\t"${DATA} ;
     done > WBP16_quality_stats.txt

awk '{print $1,$3,$5,$7,$9}' OFS="\t" WBP16_quality_stats.txt > tmp; mv tmp WBP16_quality_stats.txt
```


```R
library(tidyverse)
library(ggrepel)
library(viridis)

data_assembly <-read.table("WBP16_assembly_stats.txt", header=T)
data_quality <- read.table("WBP16_quality_stats.txt", header=F, sep="\t")
colnames(data_quality) <- c("filename", "complete",  "duplicated", "fragmented", "missing")

data <- full_join(data_assembly,data_quality, by="filename")

ggplot(data, aes(log10(n50), log10(number), size=total_length, col=as.numeric(complete), label=paste0(filename, " (",number,")" ))) +
     geom_point() +
     scale_colour_viridis(option="plasma", direction=-1, limits=c(0,100)) +
     theme_bw() +
     labs(x="N50 (log10[bp])", y="Scaffold count (log10[n])", size="Genome size (log10[bp])", col="BUSCO: Complete (%)") +
     #scale_y_reverse() +
     geom_text_repel(data=subset(data, number <200 ),size=2, box.padding = 0.5, max.overlaps = Inf) +
     scale_x_continuous(expand = expansion(mult = 0.15)) +
     scale_y_reverse() +
     scale_size()


```

### WBP17 data from Dio
```R
library(tidyverse)
library(ggrepel)
library(viridis)
library(patchwork)

data <- read.table("WBP17_stats.txt", header=T)

plot1 <- ggplot(data, aes(log10(n50), log10(scaffold_count), col=as.numeric(assembly_busco3_complete), label=paste0(species, " (",scaffold_count,")" ))) +
     geom_point() +
     scale_colour_viridis(option="plasma", direction=-1, limits=c(0,100)) +
     theme_bw() +
     labs(x="N50 (log10[bp])", y="Scaffold count (log10[n])", size="Genome size (log10[bp])", col="BUSCO: Complete (%)") +
     geom_text_repel(data=subset(data, scaffold_count <200 ),size=2, box.padding = 0.5, max.overlaps = Inf) +
     scale_x_continuous(expand = expansion(mult = 0.15)) +
     scale_y_reverse() +
     scale_size()


plot2 <- ggplot(data, aes(log10(n50), log10(scaffold_count), col=as.numeric(annotation_busco3_complete), label=paste0(species, " (",scaffold_count,")" ))) +
     geom_point() +
     scale_colour_viridis(option="plasma", direction=-1, limits=c(0,100)) +
     theme_bw() +
     labs(x="N50 (log10[bp])", y="Scaffold count (log10[n])", size="Genome size (log10[bp])", col="BUSCO: Complete (%)") +
     geom_text_repel(data=subset(data, scaffold_count <200 ),size=2, box.padding = 0.5, max.overlaps = Inf) +
     scale_x_continuous(expand = expansion(mult = 0.15)) +
     scale_y_reverse() +
     scale_size()

plot1 + plot2 + plot_layout(ncol=1)

```









- showing genome improvement over time
```R
library(tidyverse)
library(ggrepel)
library(viridis)

data <- read.table("species_year_genome_improvement.txt", header=T, sep="\t")

ggplot(data, aes(Year, log10(N50.scaffolds), col=Species, size=log10(as.numeric(Genome_size)))) +
     geom_line(size=1) +
     geom_point() +
     scale_colour_viridis(discrete=TRUE, option="plasma") +
     theme_bw() + theme(legend.position="none") +
     geom_text_repel(aes(label = Species),
               data = data %>% group_by(Species) %>% filter(Year == max(Year)),
               nudge_x = 0.5,
               box.padding = 0.5, max.overlaps = Inf,
               size = 4) +
     scale_x_continuous(expand = expansion(mult = 0.2)) +
     scale_y_continuous(expand = expansion(mult = 0.2)) +
     labs(y="N50 (log10[bp])", x="Year assembly published")
```


```bash
#hcontortus BUSCO
conda activate busco_5.3.0

cd ~/lustre118_link/papers/genome_improvement/hcontortus

ln -s laing_2013/haemonchus_contortus.PRJEB506.WBPS1.genomic.fa hcontortus_laing2013_genome.fa
ln -s doyle_2020/haemonchus_contortus.PRJEB506.WBPS14.genomic.fa hcontortus_doyle2020_genome.fa
ln -s haemonchus_contortus.PRJEB506.WBPS14.genomic.fa hcontortus_doyle2020_wbp17_genome.fa

for i in *genome.fa; do
     bsub.py 10 --threads 7 busco_${i%.fa} \
     "busco --lineage /nfs/users/nfs_s/sd21/lustre118_link/databases/busco/eukaryota_odb10 --in ${i} --out ${i%.fa}.busco5genome --mode genome -f --augustus --long --cpu 7";
     done


ln -s laing_2013/proteins.fa hcontortus_laing2013_proteins.fa
ln -s doyle_2020/proteins.fa hcontortus_doyle2020_proteins.fa
ln -s WBP17/proteins.fa hcontortus_doyle2020_wbp17_proteins.fa

for i in *proteins.fa; do
     bsub.py 10 --threads 7 busco_${i%.fa} \
     "busco --lineage /nfs/users/nfs_s/sd21/lustre118_link/databases/busco/eukaryota_odb10 --in ${i} --out ${i%.fa}.busco5protein --mode protein -f";
     done
```





# genome assembly - cumulative size curves per species

module load seqtk/1.3--ha92aebf_0

find ../../ -name "*genome.fa" | while read -r NAME; do cp ${NAME} .; done


# extract scaffold lengths and total sequecne length for plotting
>data.seqtk_out.txt
for i in *.fa; do
	SPECIES=`echo ${i} | cut -f1 -d '_'`;
	seqtk comp ${i} | sort -k2,2nr | awk -v SPECIES=${SPECIES} -v NAME=${i%_genome.fa} '{print $1,$2,NAME,SPECIES}' OFS="\t" > ${i}.seqtk_out.txt;
	done

for i in *.seqtk_out.txt; do
	SUM=$(cat ${i} | datamash sum 2 );
	awk -v SUM=${SUM} '{print $1,$2,$3,$4,SUM,$2/SUM}' OFS="\t" ${i} > ${i}.seqtk_out2.txt;
	done

cat *.seqtk_out2.txt > data.seqtk_out.txt
rm *.seqtk_out2.txt

# plot in R
```R
library(tidyverse)

data <- read.table("data.seqtk_out.txt", sep="\t")
data <- data %>% group_by(V3) %>% mutate(csum = cumsum(V6))



ggplot(data, aes(csum*100, log10(V2), group=V3, colour = V3)) +
	geom_point() +
	geom_line() +
	labs(x="Cumulative percentage of contig/scaffold lengths relative to the total genome assembly length", y="Contig/scaffold length (log10[bp])", colour="Genome version\nby species") +
	facet_grid(.~V4) +
	theme_bw()
```


### Extract and count UTR sequences
```bash

# V4 genome - WBP11+
GFF=HCON_V4_WBP11plus_190125.ips.gff3
> coding.length

for MRNA in ` awk -F '[\t]' '$3=="mRNA" {print $9}' ${GFF} | grep -Po "ID=HCON_.............." | sed -e 's/ID=//g'`; do
     chr=$( grep ${MRNA} ${GFF} | grep -v "#" |  head -n1 | cut -f1 );
     gene_start=$(grep ${MRNA} ${GFF} | grep -v "#" | awk '$3=="mRNA" {print}' | head -n1 | awk '{print $4}') ;
     gene_end=$(grep ${MRNA} ${GFF} | grep -v "#" | awk '$3=="mRNA" {print}' | head -n1 | awk '{print $5}');
     strand=$(grep ${MRNA} ${GFF} | grep -v "#" | head -n1 | cut -f7 );
     if [[ ${strand} == "+" ]]; then
     cds1=$(grep ${MRNA} ${GFF} | grep -v "#" | awk '$3=="CDS" {print}' | head -n1 | awk '{print $4}');
     cds2=$(grep ${MRNA} ${GFF} | grep -v "#" | awk '$3=="CDS" {print}' | tail -n1 | awk '{print $5}');
     else
     cds1=$(grep ${MRNA} ${GFF} | grep -v "#" | awk '$3=="CDS" {print}' | tail -n1 | awk '{print $4}');
     cds2=$(grep ${MRNA} ${GFF} | grep -v "#" | awk '$3=="CDS" {print}' | head -n1 | awk '{print $5}');
     fi;
     five_utr=$(echo "$cds1 - $gene_start" | bc);
     three_utr=$(echo "$gene_end - $cds2" | bc);
     echo -e "$chr\t$gene_start\t$gene_end\t$cds1\t$cds2\t$five_utr\t$three_utr\t$MRNA\t.\t$strand" >> coding.length; done &

number=$(wc -l coding.length); count1=$(awk '{if($6>0) print}' coding.length | wc -l); count2=$(awk '{if($7>0) print}' coding.length | wc -l); echo -e "$number\t$count1\t$count2"

# 21320 coding.length	14412	14488
# = (14412+14488)/(21320*2)
# = 0.677767355

# V4+ genome
GFF=haemonchus_contortus.PRJEB506.WBPS16.annotations.gff3
> coding.length

for MRNA in ` awk -F '[\t]' '$3=="mRNA" {print $9}' ${GFF} | grep -Po "ID=transcript:HCON_.............." | sed -e 's/ID=transcript://g'`; do
     chr=$( grep ${MRNA} ${GFF} | grep -v "#" |  head -n1 | cut -f1 );
     gene_start=$(grep ${MRNA} ${GFF} | grep -v "#" | awk '$3=="mRNA" {print}' | head -n1 | awk '{print $4}') ;
     gene_end=$(grep ${MRNA} ${GFF} | grep -v "#" | awk '$3=="mRNA" {print}' | head -n1 | awk '{print $5}');
     strand=$(grep ${MRNA} ${GFF} | grep -v "#" | head -n1 | cut -f7 );
     if [[ ${strand} == "+" ]]; then
     cds1=$(grep ${MRNA} ${GFF} | grep -v "#" | awk '$3=="CDS" {print}' | head -n1 | awk '{print $4}');
     cds2=$(grep ${MRNA} ${GFF} | grep -v "#" | awk '$3=="CDS" {print}' | tail -n1 | awk '{print $5}');
     else
     cds1=$(grep ${MRNA} ${GFF} | grep -v "#" | awk '$3=="CDS" {print}' | tail -n1 | awk '{print $4}');
     cds2=$(grep ${MRNA} ${GFF} | grep -v "#" | awk '$3=="CDS" {print}' | head -n1 | awk '{print $5}');
     fi;
     five_utr=$(echo "$cds1 - $gene_start" | bc);
     three_utr=$(echo "$gene_end - $cds2" | bc);
     echo -e "$chr\t$gene_start\t$gene_end\t$cds1\t$cds2\t$five_utr\t$three_utr\t$MRNA\t.\t$strand" >> coding.length; done &

number=$(wc -l coding.length); count1=$(awk '{if($6>0) print}' coding.length | wc -l); count2=$(awk '{if($7>0) print}' coding.length | wc -l); echo -e "$number\t$count1\t$count2"

# 20987 coding.length	9089	9134
# = (9089+9134)/(20987*2)
# = 0.434149712

# V1 genome
GFF=haemonchus_contortus.PRJEB506.WBPS8.annotations.gff3
> coding.length

for MRNA in ` awk -F '[\t]' '$3=="mRNA" {print $9}' ${GFF} | grep -Po "ID=transcript:HCOI..........." | sed -e 's/ID=transcript://g'`; do
     chr=$( grep ${MRNA} ${GFF} | grep -v "#" |  head -n1 | cut -f1 );
     gene_start=$(grep ${MRNA} ${GFF} | grep -v "#" | awk '$3=="mRNA" {print}' | head -n1 | awk '{print $4}') ;
     gene_end=$(grep ${MRNA} ${GFF} | grep -v "#" | awk '$3=="mRNA" {print}' | head -n1 | awk '{print $5}');
     strand=$(grep ${MRNA} ${GFF} | grep -v "#" | head -n1 | cut -f7 );
     if [[ ${strand} == "+" ]]; then
     cds1=$(grep ${MRNA} ${GFF} | grep -v "#" | awk '$3=="CDS" {print}' | head -n1 | awk '{print $4}');
     cds2=$(grep ${MRNA} ${GFF} | grep -v "#" | awk '$3=="CDS" {print}' | tail -n1 | awk '{print $5}');
     else
     cds1=$(grep ${MRNA} ${GFF} | grep -v "#" | awk '$3=="CDS" {print}' | tail -n1 | awk '{print $4}');
     cds2=$(grep ${MRNA} ${GFF} | grep -v "#" | awk '$3=="CDS" {print}' | head -n1 | awk '{print $5}');
     fi;
     five_utr=$(echo "$cds1 - $gene_start" | bc);
     three_utr=$(echo "$gene_end - $cds2" | bc);
     echo -e "$chr\t$gene_start\t$gene_end\t$cds1\t$cds2\t$five_utr\t$three_utr\t$MRNA\t.\t$strand" >> coding.length; done &

number=$(wc -l coding.length); count1=$(awk '{if($6>0) print}' coding.length | wc -l); count2=$(awk '{if($7>0) print}' coding.length | wc -l); echo -e "$number\t$count1\t$count2"









```
