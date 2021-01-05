
# playtime



# for mintie

transcript	chrom	exonStarts	exonEnds	gene
CHS.1.1	chr1	11873,12612,13220	12227,12721,14409	DDX11L1
CHS.2.1	chr1	14361,14969,15795,16606,16857,17232,17605,17914,18267,24737,29320	14829,15038,15947,16765,17055,17368,17742,18061,18366,24891,29370	WASH7P



awk -F'[\t;]' '{if($3=="mRNA") print $9}' HCON_V4_curated_20200422_WBPS15.gff3 | sed 's/Name=//g' >  transcript.list

#test
#echo "HCON_00000020-00001" > transcript.list

printf transcript"\t"chrom"\t"exonStarts"\t"exonEnds"\t"gene"\n" > HAEM_V4_final.info

while read TRANSCRIPT; do
grep "${TRANSCRIPT}" HCON_V4_curated_20200422_WBPS15.gff3 > ${TRANSCRIPT}.tmp;
chrom=$(cat ${TRANSCRIPT}.tmp | cut -f1 | head -n1)
exonStarts=$(cat ${TRANSCRIPT}.tmp | awk '{if($3=="exon") print $4}' | sort | xargs | sed -e 's/ /,/g')
exonEnds=$(cat ${TRANSCRIPT}.tmp | awk '{if($3=="exon") print $5}' | sort | xargs | sed -e 's/ /,/g')
gene=$( echo ${TRANSCRIPT} | cut -c-13 )

printf ${TRANSCRIPT}"\t"${chrom}"\t"${exonStarts}"\t"${exonEnds}"\t"${gene}"\n" >> HAEM_V4_final.info
rm *.tmp
done < transcript.list

sed -i 's/^[ \t]*//' HAEM_V4_final.info

cat transcript.list | awk '{print $1,substr($1,1,13)}' > tx2gene.txt








---------------------------------------------------------------------------------------------------------------------
#popooolation2 TE

working dir:
/nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/PP2_TE


# genome
ln -s ../../01_REFERENCE/HAEM_V4_final.chr.fa

# LTR digest repeats
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/REPEATS/CHR_LTR/HAEM_V4_final.chr.ltrdigest.gff.2.fa

#repeat modeller output
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/REPEATS/CHR/RM_43484.TueMar271017292018/consensi.fa.classified

cat HAEM_V4_final.chr.ltrdigest.gff.2.fa consensi.fa.classified > TE.fa



module load tetools/1.1-c3

# run repeat masker
bsub.py 10 --threads 4 RM "RepeatMasker -gccalc -s -cutoff 200 -no_is -nolow -norna -gff -u -pa 4 -lib  TE.fa HAEM_V4_final.chr.fa"

# merge masked reference and TE reference
cat HAEM_V4_final.chr.fa.masked TE.fa > HAEM_V4_final.masked.TE.merge.fa


# make a "te-heirarchy" file, containing repeat name, repeat family, repeat type
printf "id\tfamily\torder\n" > te-hierarchy.txt
cat  HAEM_V4_final.chr.fa.out | awk '{print $10,$11,$11}' OFS="\t"  | sed -e 's/\/[^\/]*//2g' -e '1,3d' | sort | uniq  >> te-hierarchy.txt

module load bwa/0.7.17=pl5.22.0_2
module load samtools/1.6--h244ad75_4

~sd21/bash_scripts/run_bwamem_splitter MHCO3_P0_L3_n200_01 $PWD/HAEM_V4_final.masked.TE.merge.fa /lustre/scratch118/infgen/pathogen/pathpipe/helminths/seq-pipelines/Haemonchus/contortus/TRACKING/2679/2679STDY7527671/SLX/22599742/26582_8#1/26582_8#1_1.fastq.gz /lustre/scratch118/infgen/pathogen/pathpipe/helminths/seq-pipelines/Haemonchus/contortus/TRACKING/2679/2679STDY7527671/SLX/22599742/26582_8#1/26582_8#1_2.fastq.gz &
~sd21/bash_scripts/run_bwamem_splitter MHCO18_P0_L3_n200_IVM_01 $PWD/HAEM_V4_final.masked.TE.merge.fa /lustre/scratch118/infgen/pathogen/pathpipe/helminths/seq-pipelines/Haemonchus/contortus/TRACKING/2679/2679STDY7527672/SLX/22599754/26582_8#2/26582_8#2_1.fastq.gz /lustre/scratch118/infgen/pathogen/pathpipe/helminths/seq-pipelines/Haemonchus/contortus/TRACKING/2679/2679STDY7527672/SLX/22599754/26582_8#2/26582_8#2_2.fastq.gz &

~sd21/bash_scripts/run_bwamem_splitter XQTL_F3_L3_n200_IVM_post_01 $PWD/HAEM_V4_final.masked.TE.merge.fa /lustre/scratch118/infgen/pathogen/pathpipe/helminths/seq-pipelines/Haemonchus/contortus/TRACKING/2679/2679STDY6583135/SLX/18074181/21395_2#2/21395_2#2_1.fastq.gz /lustre/scratch118/infgen/pathogen/pathpipe/helminths/seq-pipelines/Haemonchus/contortus/TRACKING/2679/2679STDY6583135/SLX/18074181/21395_2#2/21395_2#2_2.fastq.gz &




bsub.py 10 pp2_te \
java -Xmx10g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/POOLSEQ/popte2-v1.10.04.jar ppileup \
--bam /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/PP2_TE/MHCO3_P0_L3_n200_01_bwasplitter_out/MHCO3_P0_L3_n200_01.merged.sorted.marked.bam \
--bam /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/PP2_TE/MHCO18_P0_L3_n200_IVM_01_bwasplitter_out/MHCO18_P0_L3_n200_IVM_01.merged.sorted.marked.bam \
--bam /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/PP2_TE/XQTL_F3_L3_n200_IVM_post_01_bwasplitter_out/XQTL_F3_L3_n200_IVM_post_01.merged.sorted.marked.bam \
--map-qual 15 --hier te-hierarchy.txt --output output.ppileup \
--detailed-log

mv output.ppileup output.ppileup.gz

bsub.py 10 pp2_te_identifySignatures "java -Xmx10g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/POOLSEQ/popte2-v1.10.04.jar identifySignatures --ppileup output.ppileup.gz --mode joint --output ISEvUGA.signatures --min-count 10 --signature-window minimumSampleMedian --min-valley minimumSampleMedian"
bsub.py 10 pp2_te_freq "java -Xmx10g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/POOLSEQ/popte2-v1.10.04.jar frequency --ppileup output.ppileup.gz --signature ISEvUGA.signatures --output ISEvUGA.freqsig"
bsub.py 10 pp2_te_PUsigs "java -Xmx10g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/POOLSEQ/popte2-v1.10.04.jar pairupSignatures --signature ISEvUGA.freqsig --ref-genome HAEM_V4_final.masked.TE.merge.fa --hier te-hierarchy.txt --min-distance -200 --max-distance 300 --output ISEvUGA.teinsertions"












# map RNAseq reads to reference using STAR
for i in $(ls *R1_001.fastq.gz | rev | cut -c 17- | rev) ; do \
/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/STAR/bin/Linux_x86_64/STAR \
--twopassMode Basic \
--runThreadN 8 \
--genomeDir /nfs/users/nfs_s/sd21/lustre118_link/REFERENCE_SEQUENCES/haemonchus_contortus/hc_v4chr_rnaseq_star_index \
--readFilesIn ${i}_R1_001.fastq.gz ${i}_R2_001.fastq.gz \
--readFilesCommand zcat \
--outSAMstrandField intronMotif \
--outSAMtype BAM Unsorted \
--outFileNamePrefix ${i}; \
done







conda activate py37
export LD_LIBRARY_PATH=/nfs/users/nfs_s/sd21/lustre118_link/software/anaconda2/envs/py37/lib/


#!/bin/sh
cases=`ls cases/*fastq.gz`
controls=`ls controls/*fastq.gz`

/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/MINTIE/tools/bin/bpipe \
run @/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/MINTIE/params.txt \
/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/MINTIE/MINTIE.groovy $cases $controls
#kisssplice
/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/STAR/bin/Linux_x86_64/STARlong --genomeDir /nfs/users/nfs_s/sd21/lustre118_link/REFERENCE_SEQUENCES/haemonchus_contortus/hc_v4chr_rnaseq_star_index --readFilesIn *type_1.fa --outSAMunmapped Within --outFileNamePrefix type1_mapped
/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/STAR/bin/Linux_x86_64/STARlong --genomeDir /nfs/users/nfs_s/sd21/lustre118_link/REFERENCE_SEQUENCES/haemonchus_contortus/hc_v4chr_rnaseq_star_index --readFilesIn *type_2.fa --outSAMunmapped Within --outFileNamePrefix type2_mapped
/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/STAR/bin/Linux_x86_64/STARlong --genomeDir /nfs/users/nfs_s/sd21/lustre118_link/REFERENCE_SEQUENCES/haemonchus_contortus/hc_v4chr_rnaseq_star_index --readFilesIn *type_3.fa --outSAMunmapped Within --outFileNamePrefix type3_mapped
/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/STAR/bin/Linux_x86_64/STARlong --genomeDir /nfs/users/nfs_s/sd21/lustre118_link/REFERENCE_SEQUENCES/haemonchus_contortus/hc_v4chr_rnaseq_star_index --readFilesIn *type_4.fa --outSAMunmapped Within --outFileNamePrefix type4_mapped






-rw-r--r-- 1 sd21  225 Oct 17  2018 ZAI.population.list
-rw-r--r-- 1 sd21   66 Oct 17  2018 US.population.list
-rw-r--r-- 1 sd21  135 Oct 17  2018 STO.population.list
-rw-r--r-- 1 sd21  660 Oct 17  2018 STA.population.list
-rw-r--r-- 1 sd21   75 Oct 17  2018 POR.population.list
-rw-r--r-- 1 sd21  165 Oct 17  2018 PK.population.list
-rw-r--r-- 1 sd21  225 Oct 17  2018 NAM.population.list
-rw-r--r-- 1 sd21  195 Oct 17  2018 MOR.population.list
-rw-r--r-- 1 sd21  135 Oct 17  2018 ITA.population.list
-rw-r--r-- 1 sd21  105 Oct 17  2018 IND.population.list
-rw-r--r-- 1 sd21  528 Oct 17  2018 GB.population.list
-rw-r--r-- 1 sd21  345 Oct 17  2018 FRG.population.list
-rw-r--r-- 1 sd21 1005 Oct 17  2018 FRA.population.list
-rw-r--r-- 1 sd21   77 Oct 17  2018 CH.population.list
-rw-r--r-- 1 sd21  210 Oct 17  2018 CAP.population.list
-rw-r--r-- 1 sd21  225 Oct 17  2018 BRA.population.list
-rw-r--r-- 1 sd21  195 Oct 17  2018 BEN.population.list
-rw-r--r-- 1 sd21  225 Oct 17  2018 AUS.population.list


#P167
hcontortus_chr1_Celeg_TT_arrow_pilon 7029569
#P198
hcontortus_chr1_Celeg_TT_arrow_pilon 7029788
# P200
hcontortus_chr1_Celeg_TT_arrow_pilon 7029790
# btub iso2
hcontortus_chr2_Celeg_TT_arrow_pilon 13435823



# for i in $(ls *.population.list); do \
#      vcftools \
#      --gzvcf ../../1.hcontortus_chr1_Celeg_TT_arrow_pilon.cohort.vcf.gz \
#      --gzvcf ../../2.hcontortus_chr2_Celeg_TT_arrow_pilon.cohort.vcf.gz \
#      --keep ${i} \
#      --positions btub.positions \
#      --freq --out global_${i}; \
# done


for i in $(ls *.population.list); do \
zcat ../../1.hcontortus_chr1_Celeg_TT_arrow_pilon.cohort.vcf.gz ../../2.hcontortus_chr2_Celeg_TT_arrow_pilon.cohort.vcf.gz | \
     vcftools \
     --vcf - \
     --keep ${i} \
     --positions btub.positions \
     --freq --out global_${i}; \
done

#iso2 p200
hcontortus_chr2_Celeg_TT_arrow_pilon 13435829

for i in $(ls *.population.list); do \
     vcftools \
     --gzvcf ../../2.hcontortus_chr2_Celeg_TT_arrow_pilon.cohort.vcf.gz \
     --keep ${i} \
     --positions btub.positions2 \
     --freq --out global_iso2P200.${i}; \
done
#------










# plot of backcross data for Erik Andersen
cavr <- ggplot(chr5,aes(V2,V21))+
     geom_point(size=0.5)+
     xlim(30e6,45e6)+
     geom_vline(xintercept=c(36.8e6,41.3e6))+
     geom_hline(yintercept=mean(chr5$V21)+3*sd(chr5$V21))+
     geom_vline(xintercept=c(37.2e6,37.5e6),col="red")


wrs <- ggplot(chr5,aes(V2,V25))+
     geom_point(size=0.5)+
     xlim(30e6,45e6)+
     geom_vline(xintercept=c(36.8e6,42.4e6))+
     geom_hline(yintercept=mean(chr5$V25)+3*sd(chr5$V25))+
     geom_vline(xintercept=c(37.2e6,37.5e6),col="red")


cavr + wrs + plot_layout(ncol=1)








# forqs
LANG=/usr/lib/locale/en_US
export LC_ALL=C; unset LANGUAGE




# XQTL supplementary data - looking a fst distribution, and cutoffs

cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_CONTROL


R
# load libraries
library(ggplot2)
library(ggridges)
library(dplyr)
library(patchwork)

# get data
control <- read.table("XQTL_CONTROL/XQTL_CONTROL.merged.fst",header=F)
bz <- read.table("XQTL_BZ/XQTL_BZ.merged.fst",header=F)
lev <- read.table("XQTL_LEV/XQTL_LEV.merged.fst",header=F)
ivm <- read.table("XQTL_IVM/XQTL_IVM.merged.fst",header=F)

# add labels before merging
control$id <- "1. Control"
bz$id <- "2. Benzimidazole treated"
lev$id <- "3. Levamisole treated"
ivm$id <- "4. Ivermectin treated"

control <- select(control,c(V1,V2,V13,id))
bz <- select(bz,c(V1,,V2,V13,id))
lev <- select(lev,c(V1,V2,V13,id))
ivm <- select(ivm,c(V1,V2,V13,id))

data <- bind_rows(control, bz, lev, ivm)
data <- data[data$V1!="hcontortus_chr_mtDNA_arrow_pilon",]


vline.data <- data %>%
              group_by(id) %>%
              summarize(mean_fst_3sd = mean(V13)+3*sd(V13))

# plot
fst_distribution_plot <- ggplot(data,aes(V13,id)) +
     geom_density_ridges(aes(fill = id),scale = 1) + xlim(0,0.05) +
     geom_vline(aes(xintercept = mean_fst_3sd, col = id), vline.data, size = 1, linetype = "dashed") +
     scale_y_discrete(limits = rev(data$id)) + theme_bw() +
     labs(title="A", x="FST")


# proportion of values above the control threshold
nrow(control[control$V13>0.0235,])/nrow(control)*100
= 1.033201

nrow(bz[bz$V13>0.0235,])/nrow(bz)*100
= 4.210414

nrow(lev[lev$V13>0.0235,])/nrow(lev)*100
= 11.21552

nrow(ivm[ivm$V13>0.0235,])/nrow(ivm)*100
= 2.399702




# positions of variants gt fst+3sd
bz_high <- bz %>% filter(V13 > mean(V13)+3*sd(V13))
control_bz_high <- dplyr::inner_join(control, bz_high, by = c("V1","V2"))
#> control_bz_high = 712

# false positives
control_bz_high_control_high <- control_bz_high %>% filter(V13.x > mean(control$V13)+3*sd(control$V13))
#> 21 gt control

nrow(control_bz_high_control_high)/nrow(control_bz_high)
#> 0.02949438

# positions of variants gt fst+3sd
lev_high <- lev %>% filter(V13 > mean(V13)+3*sd(V13))
control_lev_high <- dplyr::inner_join(control, lev_high, by = c("V1","V2"))
nrow(control_lev_high)
#= 715

# false positives
control_lev_high_control_high <- control_lev_high %>% filter(V13.x > mean(control$V13)+3*sd(control$V13))
nrow(control_lev_high_control_high)
#> 33 gt control

nrow(control_lev_high_control_high)/nrow(control_lev_high)
#> 0.04615385

# positions of variants gt fst+3sd
ivm_high <- ivm %>% filter(V13 > mean(V13)+3*sd(V13))
control_ivm_high <- dplyr::inner_join(control, ivm_high, by = c("V1","V2"))
nrow(control_ivm_high)
#= 698

# false positives
control_ivm_high_control_high <- control_ivm_high %>% filter(V13.x > mean(control$V13)+3*sd(control$V13))
nrow(control_ivm_high_control_high)
#> 39 gt control

nrow(control_ivm_high_control_high)/nrow(control_ivm_high)
#> 0.05587393


data2 <- bind_rows(control_bz_high, control_lev_high, control_ivm_high)



fp_plot <- ggplot(data2,aes(V13.x,1)) +
     geom_jitter(aes(color = V13.x > mean(control$V13)+3*sd(control$V13)), alpha=0.4)+ xlim(0,0.05) +
     scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"))+
     theme_bw() +
     labs(title = "B", x="FST", colour="False positive")+
     facet_grid(id.y~.)

fst_distribution_plot + fp_plot + plot_layout(ncol = 1)
