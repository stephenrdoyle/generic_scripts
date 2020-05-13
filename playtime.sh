# playtime

extracting allele frequencies of beta tubulin P167,P198, P200 to make figure


working dir:
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_BZ


grep "pre" bam.list  > bz_pretreatment_samples.list
grep "post" bam.list  > bz_posttreatment_samples.list

grep "pre" ../XQTL_CONTROL/bam.list  > control_pretreatment_samples.list
grep "post" ../XQTL_CONTROL/bam.list  > control_posttreatment_samples.list

vcftools --vcf XQTL_BZ.raw.snpeff.vcf --keep bz_pretreatment_samples.list --positions btub.positions --extract-FORMAT-info AD --out bz_pretreatment
vcftools --vcf XQTL_BZ.raw.snpeff.vcf --keep bz_posttreatment_samples.list --positions btub.positions --extract-FORMAT-info AD --out bz_posttreatment

vcftools --gzvcf ../XQTL_CONTROL/1.hcontortus_chr1_Celeg_TT_arrow_pilon.tmp.vcf.gz --keep control_pretreatment_samples.list --positions btub.positions --extract-FORMAT-info AD --out control_pretreatment
vcftools --gzvcf ../XQTL_CONTROL/1.hcontortus_chr1_Celeg_TT_arrow_pilon.tmp.vcf.gz --keep control_posttreatment_samples.list --positions btub.positions --extract-FORMAT-info AD --out control_posttreatment

#where "btub.positions" contains:
#P167
#hcontortus_chr1_Celeg_TT_arrow_pilon 7029569
#P198
#hcontortus_chr1_Celeg_TT_arrow_pilon 7029788
# P200
#hcontortus_chr1_Celeg_TT_arrow_pilon 7029790



for i in `ls bz*AD.FORMAT`; do
      grep "^hcon" ${i} | awk -F '[\t,]' '{print $1,$2,$4/($3+$4),$6/($5+$6),$8/($7+$8),$10/($9+$10)}' OFS="\t" > ${i%.AD.FORMAT}.ADfreq;
done

for i in `ls control*AD.FORMAT`; do
      grep "^hcon" ${i} | awk -F '[\t,]' '{print $1,$2,$4/($3+$4),$6/($5+$6),$8/($7+$8)}' OFS="\t" > ${i%.AD.FORMAT}.ADfreq;
done

# hcontortus_chr1_Celeg_TT_arrow_pilon	7029569	0.125	0.197368	0.223301	0.264463	0.15	0.128571	0.131148	0.155172
# hcontortus_chr1_Celeg_TT_arrow_pilon	7029790	0.721311	0.805556	0.866071	0.926316	0.534884	0.477612	0.446281	0.347368
#
#
#
# awk '{print $1,$2,$7,$8,$9,$10,"bz_pre"}' OFS="\t" bz.freq > bz.freq2
# awk '{print $1,$2,$3,$4,$5,$6,"bz_post"}' OFS="\t" bz.freq >> bz.freq2
#
#
# awk '{print $1,$2,$7,$3,"BZ_treated","R1"}' OFS="\t" bz.freq > bz.freq2
# awk '{print $1,$2,$8,$4,"BZ_treated","R2"}' OFS="\t" bz.freq >> bz.freq2
# awk '{print $1,$2,$9,$5,"BZ_treated","R3"}' OFS="\t" bz.freq >> bz.freq2
# awk '{print $1,$2,$10,$6,"BZ_treated","R4"}' OFS="\t" bz.freq >> bz.freq2
#
#
#
# ggplot(a) + geom_segment(aes(x="1", xend="2", y=V3, yend=V4,col=factor(V2),group=V2), size=.75)+facet_grid(V5~V2)
#
#
#
#
# grep "^hcon" out.AD.FORMAT | awk -F '[\t,]' '{for(i=3;i<=NF;i+=2) print $1,$2,(i+=1)/($i+(i+=1))}'
#
#
# for ((i=3,j=4;i<=j;i+=2,j+=2)); do \
#      grep "^hcon" out.AD.FORMAT | \
#      awk -v i=$i -v j=$j -F '[\t,]' '{if(i<NF) print $1,$2,$j/($i+$j)}'
# done


```R
R
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(rstatix)

bz_pre <- read.table("bz_pretreatment.ADfreq")
colnames(bz_pre) <- c("CHR","POS","R1","R1.2","R2","R3")
bz_pre <- melt(bz_pre, id = c("CHR", "POS"), variable.name = "SAMPLE_ID")

bz_post <- read.table("bz_posttreatment.ADfreq")
colnames(bz_post) <- c("CHR","POS","R1","R1.2","R2","R3")
bz_post <- melt(bz_post, id = c("CHR", "POS"), variable.name = "SAMPLE_ID")

bz_data <- dplyr::full_join(bz_pre, bz_post, by = c("CHR","POS","SAMPLE_ID"))
bz_data$TREATMENT <- "BZ_treated"
colnames(bz_data) <- c("CHR","POS","SAMPLE_ID","PRE_TREATMENT","POST_TREATMENT","TREATMENT")


control_pre <- read.table("control_pretreatment.ADfreq")
colnames(control_pre) <- c("CHR","POS","R1","R2","R3")
control_pre <- melt(control_pre, id = c("CHR", "POS"), variable.name = "SAMPLE_ID")

control_post <- read.table("control_posttreatment.ADfreq")
colnames(control_post) <- c("CHR","POS","R1","R2","R3")
control_post <- melt(control_post, id = c("CHR", "POS"), variable.name = "SAMPLE_ID")

control_data <- dplyr::full_join(control_pre, control_post, by = c("CHR","POS","SAMPLE_ID"))
control_data$TREATMENT <- "Untreated_control"
colnames(control_data) <- c("CHR","POS","SAMPLE_ID","PRE_TREATMENT","POST_TREATMENT","TREATMENT")


# bring datasets together
data <- dplyr::bind_rows(control_data, bz_data)

# change the labels
data <- data %>%
  mutate(POS = str_replace(POS, c("7029569","7029790"), c("Phe167Tyr (SNP:7029569)","Phe200Tyr (SNP:7029790)")))

# make the plot
plot <- ggplot(data) +
     geom_segment(aes(x="1.PRE", xend="2.POST", y=PRE_TREATMENT, yend=POST_TREATMENT,col=factor(SAMPLE_ID),group=POS), size=1) +
     labs(title="Beta tubulin isotype 1 (HCON_00005260): \nResistant allele frequency pre- and post-treatment",x="Sampling time-point",y="Resistant allele frequency",col="Replicate ID") +
     ylim(0,1)+
     facet_grid(TREATMENT~POS)+
     theme_bw()
```


# perform pairwise t tests between pre/post for each SNP on BZ treated samples
bz_data_stats <- bz_data %>%
  gather(key = "TREATMENT", value = "FREQ", PRE_TREATMENT, POST_TREATMENT)

bz_stat.test <- bz_data_stats %>%
    group_by(POS) %>%
    pairwise_t_test(
      FREQ ~ TREATMENT, paired = TRUE,
      p.adjust.method = "bonferroni"
      ) %>%
    select(-df, -statistic, -p) # Remove details
bz_stat.test

bz_stat.test$TREATMENT <- "BZ_treated"
bz_stat.test <- bz_stat.test %>%
  mutate(POS = str_replace(POS, c("7029569","7029790"), c("Phe167Tyr (SNP:7029569)","Phe200Tyr (SNP:7029790)")))

# perform pairwise t tests between pre/post for each SNP on control samples
control_data_stats <- control_data %>%
  gather(key = "TREATMENT", value = "FREQ", PRE_TREATMENT, POST_TREATMENT)

control_stat.test <- control_data_stats %>%
    group_by(POS) %>%
    pairwise_t_test(
      FREQ ~ TREATMENT, paired = TRUE,
      p.adjust.method = "bonferroni"
      ) %>%
    select(-df, -statistic, -p) # Remove details
control_stat.test

control_stat.test$TREATMENT <- "Untreated_control"
control_stat.test <- control_stat.test %>%
  mutate(POS = str_replace(POS, c("7029569","7029790"), c("Phe167Tyr (SNP:7029569)","Phe200Tyr (SNP:7029790)")))

p.data <- dplyr::bind_rows(control_stat.test, bz_stat.test)


# make new plot with p values annotated on it
plot_pvalues <- plot +
     geom_text(data=p.data, aes(x="1.PRE", y=0.95, group=POS, label = paste('P = ',p.adj)))


#-----------------------------------------------------------------------------------------

# isotype 2 in US farm data

working dir:
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/US_FIELD/VCF

echo "hcontortus_chr2_Celeg_TT_arrow_pilon 13435823" > btub.positions

cat bam.list > samples.list

vcftools --vcf 2.hcontortus_chr2_Celeg_TT_arrow_pilon.snpeff.vcf --keep samples.list --positions btub.positions --extract-FORMAT-info AD --out us_farms_btub2

for i in `ls *AD.FORMAT`; do
      grep "^hcon" ${i} | awk -F '[\t,]' '{print $1,$2,$4/($3+$4),$6/($5+$6),$8/($7+$8),$10/($9+$10),$12/($11+$12),$14/($13+$14),$16/($15+$16),$18/($17+$18),$20/($19+$20),$22/($21+$22)}' OFS="\t" > ${i%.AD.FORMAT}.ADfreq;
done

R
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(rstatix)

us_btub2 <- read.table("us_farms_btub2.ADfreq")
colnames(us_btub2) <- c("CHR","POS","UGA Susceptible","Ober","Histon","Krushing","Alday","Mulligan","Rau-Shelby Acres","Pena","679-238-USDA","Strickland")
us_btub2 <- melt(us_btub2, id = c("CHR", "POS"), variable.name = "SAMPLE_ID")
