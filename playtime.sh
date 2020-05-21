
# playtime




#--------------------------------------------------------------------------------
# Figure 2 - Fst plot for parent strains

working dir:
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/PARENTS

parents <- read.table("XQTL_PARENTS.merged.fst",header=F)
parents <- parents[parents$V1!="hcontortus_chr_mtDNA_arrow_pilon",]
parents <- dplyr::select(parents, V1, V2, V7)
colnames(parents) <- c("CHR","POS","FST")
parents$LABEL <- "P0: MHco3(ISE) vs MHco18(UGA)"
parents$ROW_ID <- 1:nrow(parents)
colnames(parents) <- c("CHR","POS","FST","LABEL","ROW_ID")

data <- parents

# make the parental plot

# genome wide signficance per sample
data_gws <- data %>%
    group_by(LABEL) %>%
    summarise(GWS = mean(FST)+3*sd(FST))


plot_a <- ggplot(data)+
     geom_hline(data = data_gws, aes(yintercept=GWS), linetype="dashed",col="black")+
     geom_point(aes(ROW_ID*5000, FST, colour = CHR, group=LABEL),size=0.1)+
     ylim(0,1)+xlim(0,300e6)+
     labs(title="A",x="Chromosome position (5 kbp window)", y="Genetic differentiation (Fst)")+
     scale_color_manual(values=chr_colours)+
     scale_x_continuous(breaks=seq(0,3e8,0.5e8),limits = c(0,300e6))+
     theme_bw()+theme(legend.position="none",text = element_text(size=10))+
     facet_grid(LABEL~.)



control <- read.table("XQTL_CONTROL.merged.fst",header=F)
control <- control[control$V1!="hcontortus_chr_mtDNA_arrow_pilon",]
control <- dplyr::select(control, V1, V2, V11)
colnames(control) <- c("CHR","POS","FST")
control$LABEL <- "1. F2: Control"
control$ROW_ID <- 1:nrow(control)
colnames(control) <- c("CHR","POS","FST","LABEL","ROW_ID")

bz <- read.table("XQTL_BZ.merged.fst",header=F)
bz <- bz[bz$V1!="hcontortus_chr_mtDNA_arrow_pilon",]
bz <- dplyr::select(bz, V1, V2, V13)
colnames(bz) <- c("CHR","POS","FST")
bz$LABEL <- "2. F2: Benzimidazole"
bz$ROW_ID <- 1:nrow(bz)
colnames(bz) <- c("CHR","POS","FST","LABEL","ROW_ID")

data <- dplyr::bind_rows(control,bz)

lev <- read.table("XQTL_LEV.merged.fst",header=F)
lev <- lev[lev$V1!="hcontortus_chr_mtDNA_arrow_pilon",]
lev <- dplyr::select(lev, V1, V2, V13)
colnames(lev) <- c("CHR","POS","FST")
lev$LABEL <- "3. F2: Levamisole"
lev$ROW_ID <- 1:nrow(lev)
colnames(lev) <- c("CHR","POS","FST","LABEL","ROW_ID")

data <- dplyr::bind_rows(data,lev)

ivm <- read.table("XQTL_IVM.merged.fst",header=F)
ivm <- ivm[ivm$V1!="hcontortus_chr_mtDNA_arrow_pilon",]
ivm <- dplyr::select(ivm, V1, V2, V13)
colnames(ivm) <- c("CHR","POS","FST")
ivm$LABEL <- "4. F2: Ivermectin"
ivm$ROW_ID <- 1:nrow(ivm)
colnames(ivm) <- c("CHR","POS","FST","LABEL","ROW_ID")

data <- dplyr::bind_rows(data,ivm)

# make the XQTL plot
chr_colours<-c("blue","cornflowerblue","blue","cornflowerblue","blue","cornflowerblue")

# genome wide signficance per sample
data_gws <- data %>%
    group_by(LABEL) %>%
    summarise(GWS = mean(FST)+3*sd(FST))

plot_b <- ggplot(data)+
     geom_hline(data = data_gws, aes(yintercept=GWS), linetype="dashed",col="black")+
     geom_point(aes(ROW_ID*5000, FST, colour = CHR, group=LABEL),size=0.1)+
     ylim(0,0.1)+
     labs(title="B",x="Chromosome position (5 kbp window)", y="Genetic differentiation (Fst)")+
     scale_color_manual(values=chr_colours)+
     scale_x_continuous(breaks=seq(0,3e8,0.5e8),limits = c(0,300e6))+
     theme_bw()+theme(legend.position="none",text = element_text(size=10))+
     facet_grid(LABEL~.)


# combine the parental and XQTL plots into a single figure
plot_a + plot_b + plot_layout(ncol=1, heights = c(1,4))
ggsave("genomewide_fst_plots.pdf",useDingbats=FALSE,width=170,height=200,units="mm")

<!-- ```R
R
library(ggplot2)

# get data
data <- read.table("XQTL_PARENTS.merged.fst",header=F)
data <- data[data$V1!="hcontortus_chr_mtDNA_arrow_pilon",]

data$LABEL <- "Parents: MHco3(ISE) vs MHco18(UGA)"
# genome wide levels of significance


chr_colours<-c("blue","cornflowerblue","blue","cornflowerblue","blue","cornflowerblue")

fst_column <- data$V7

genomewide_sig <- mean(fst_column)+(3*sd(fst_column))

plot_a <- ggplot(data)+
     geom_hline(yintercept=genomewide_sig, linetype="dashed",col="black")+
     geom_point(aes(1:nrow(data)*5000, fst_column, colour = V1),size=0.1)+
     ylim(0,1)+
     labs(title="A",x="Chromosome position (5 kbp window)", y="Genetic differentiation (Fst)")+
     scale_color_manual(values=chr_colours)+
     scale_x_continuous(breaks=seq(0,3e8,0.5e8))+
     theme_bw()+theme(legend.position="none",text = element_text(size=10))+
     facet_grid(LABEL~.)


#ggsave("XQTL_parents_fst.pdf",useDingbats=FALSE,width=170,height=50,units="mm")

```

```R
R
library(ggplot2)

# get data
data <- read.table("XQTL_CONTROL.merged.fst",header=F)
data <- data[data$V1!="hcontortus_chr_mtDNA_arrow_pilon",]
data$LABEL <- "Control"
# genome wide levels of significance


chr_colours<-c("blue","cornflowerblue","blue","cornflowerblue","blue","cornflowerblue")

# Replicates
# 1:4 = V11
# 2:5 = V21
# 3:6 = V29

fst_column <- data$V11

genomewide_sig <- mean(fst_column)+(3*sd(fst_column))

plot_b <- ggplot(data)+
     geom_hline(yintercept=genomewide_sig, linetype="dashed",col="black")+
     geom_point(aes(1:nrow(data)*5000, fst_column, colour = V1),size=0.1)+
     ylim(0,0.1)+
     labs(title="B",x="Chromosome position (5 kbp window)", y="Genetic differentiation (Fst)")+
     scale_color_manual(values=chr_colours)+
     scale_x_continuous(breaks=seq(0,3e8,0.5e8))+
     theme_bw()+theme(legend.position="none",text = element_text(size=10))+
     facet_grid(LABEL~.)

```
```R
R
library(ggplot2)

# get data
data <- read.table("XQTL_BZ.merged.fst",header=F)
data <- data[data$V1!="hcontortus_chr_mtDNA_arrow_pilon",]
data$LABEL <- "Benzimidazole"
# genome wide levels of significance


chr_colours<-c("blue","cornflowerblue","blue","cornflowerblue","blue","cornflowerblue")

# Replicates
# 1:5 = V13
# 2:6 = V27
# 3:7 = V39
# 4:8 = V49


fst_column <- data$V13

genomewide_sig <- mean(fst_column)+(3*sd(fst_column))

plot_c <- ggplot(data)+
     geom_hline(yintercept=genomewide_sig, linetype="dashed",col="black")+
     geom_point(aes(1:nrow(data)*5000, fst_column, colour = V1),size=0.1)+
     ylim(0,0.1)+
     labs(title="C",x="Chromosome position (5 kbp window)", y="Genetic differentiation (Fst)")+
     scale_color_manual(values=chr_colours)+
     scale_x_continuous(breaks=seq(0,3e8,0.5e8))+
     theme_bw()+theme(legend.position="none",text = element_text(size=10))+
     facet_grid(LABEL~.)
```

```R
R
library(ggplot2)

# get data
data <- read.table("XQTL_LEV.merged.fst",header=F)
data <- data[data$V1!="hcontortus_chr_mtDNA_arrow_pilon",]
data$LABEL <- "Levamisole"
# genome wide levels of significance


chr_colours<-c("blue","cornflowerblue","blue","cornflowerblue","blue","cornflowerblue")

# Replicates
# 1:5 = V13
# 2:6 = V27
# 3:7 = V39
# 4:8 = V49


fst_column <- data$V13

genomewide_sig <- mean(fst_column)+(3*sd(fst_column))

plot_d <- ggplot(data)+
     geom_hline(yintercept=genomewide_sig, linetype="dashed",col="black")+
     geom_point(aes(1:nrow(data)*5000, fst_column, colour = V1),size=0.1)+
     ylim(0,0.1)+
     labs(title="D",x="Chromosome position (5 kbp window)", y="Genetic differentiation (Fst)")+
     scale_color_manual(values=chr_colours)+
     scale_x_continuous(breaks=seq(0,3e8,0.5e8))+
     theme_bw()+theme(legend.position="none",text = element_text(size=10))+
     facet_grid(LABEL~.)
```

## Ivermectin
```R
R
library(ggplot2)

# get data
data <- read.table("XQTL_IVM.merged.fst",header=F)
data <- data[data$V1!="hcontortus_chr_mtDNA_arrow_pilon",]
data$LABEL <- "Ivermectin"
# genome wide levels of significance


chr_colours<-c("blue","cornflowerblue","blue","cornflowerblue","blue","cornflowerblue")

# Replicates
# 1:5 = V13
# 2:6 = V27
# 3:7 = V39
# 4:8 = V49


fst_column <- data$V13

genomewide_sig <- mean(fst_column)+(3*sd(fst_column))

plot_e <- ggplot(data)+
     geom_hline(yintercept=genomewide_sig, linetype="dashed",col="black")+
     geom_point(aes(1:nrow(data)*5000, fst_column, colour = V1),size=0.1)+
     ylim(0,0.1)+
     labs(title="E",x="Chromosome position (5 kbp window)", y="Genetic differentiation (Fst)")+
     scale_color_manual(values=chr_colours)+
     scale_x_continuous(breaks=seq(0,3e8,0.5e8))+
     theme_bw()+theme(legend.position="none",text = element_text(size=10))+
     facet_grid(LABEL~.)
```

# bring it together
```R
library(patchwork)
plot_a + plot_b + plot_c + plot_d + plot_e + plot_layout(ncol=1)


``` -->




#-------------------------------------------------------------------------------
# benzimidaole figure - panel A - chromosome 1
```R
R
library(ggplot2)

# import fst data
xqtl_bz_fst <- read.table("XQTL_BZ.merged.fst",header=F)

xqtl_bz_fst_chr1 <- xqtl_bz_fst[xqtl_bz_fst$V1=="hcontortus_chr1_Celeg_TT_arrow_pilon",]

# calculate a genome wide significance cutoff
genomewide_sig <- mean(xqtl_bz_fst_chr1$V13)+(3*sd(xqtl_bz_fst_chr1$V13))

xqtl_bz_fst_chr1_peaks <- xqtl_bz_fst_chr1[(xqtl_bz_fst_chr1$V2 >= peaks$PEAK_START_COORD) & (xqtl_bz_fst_chr1$V2 <= peaks$PEAK_END_COORD),]

# get predicted peak windows

peaks <- read.table("peak.windows.bed",header=T)

peak_subset <- xqtl_bz_fst_chr1[(xqtl_bz_fst_chr1$V2 >= peaks$PEAK_START_COORD) & (xqtl_bz_fst_chr1$V2 <= peaks$PEAK_END_COORD),]

# make the plot
plot_a <- ggplot(xqtl_bz_fst_chr1)+
     geom_hline(yintercept=genomewide_sig, linetype="dashed",col="black")+
     geom_vline(xintercept=7029790,linetype="dashed",col="grey")+
     geom_point(aes(V2,V13,group=V1), col="cornflowerblue",size=0.5)+
     geom_point(data = subset(xqtl_bz_fst_chr1,(xqtl_bz_fst_chr1$V2 >= peaks$PEAK_START_COORD[1]) & (xqtl_bz_fst_chr1$V2 <= peaks$PEAK_END_COORD[1]) & (V13 > genomewide_sig)),aes(V2,V13),col="red",size=1)+
     geom_point(data = subset(xqtl_bz_fst_chr1,(xqtl_bz_fst_chr1$V2 >= peaks$PEAK_START_COORD[2]) & (xqtl_bz_fst_chr1$V2 <= peaks$PEAK_END_COORD[2]) & (V13 > genomewide_sig)),aes(V2,V13),col="red",size=1)+
     geom_point(data = subset(xqtl_bz_fst_chr1,(xqtl_bz_fst_chr1$V2 >= peaks$PEAK_START_COORD[3]) & (xqtl_bz_fst_chr1$V2 <= peaks$PEAK_END_COORD[3]) & (V13 > genomewide_sig)),aes(V2,V13),col="red",size=1)+
     ylim(0,0.1)+xlim(0,50e6)+
     theme_bw()+theme(legend.position="none",text = element_text(size=10))+
     labs(title="A",x="Chromosome position (5 kbp window)", y="Genetic differentiation (Fst)")+
     facet_grid(.~V1)
```



#-------------------------------------------------------------------------------

# benzimidaole figure - panel B - beta tubulin isotype 1 data


#extracting allele frequencies of beta tubulin P167,P198, P200 to make figure

```shell
#working dir:
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

```


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
bz_data$TREATMENT <- "Benzimidazole"
colnames(bz_data) <- c("CHR","POS","SAMPLE_ID","PRE_TREATMENT","POST_TREATMENT","TREATMENT")


control_pre <- read.table("control_pretreatment.ADfreq")
colnames(control_pre) <- c("CHR","POS","R1","R2","R3")
control_pre <- melt(control_pre, id = c("CHR", "POS"), variable.name = "SAMPLE_ID")

control_post <- read.table("control_posttreatment.ADfreq")
colnames(control_post) <- c("CHR","POS","R1","R2","R3")
control_post <- melt(control_post, id = c("CHR", "POS"), variable.name = "SAMPLE_ID")

control_data <- dplyr::full_join(control_pre, control_post, by = c("CHR","POS","SAMPLE_ID"))
control_data$TREATMENT <- "Untreated"
colnames(control_data) <- c("CHR","POS","SAMPLE_ID","PRE_TREATMENT","POST_TREATMENT","TREATMENT")


# bring datasets together
data <- dplyr::bind_rows(control_data, bz_data)

# change the labels
data <- data %>%
  mutate(POS = str_replace(POS, c("7029569","7029790"), c("Phe167Tyr","Phe200Tyr")))



# make the plot
plot <- ggplot(data) +
     geom_segment(aes(x="1.PRE", xend="2.POST", y=PRE_TREATMENT, yend=POST_TREATMENT,col=factor(SAMPLE_ID),group=POS), size=1) +
     labs(title="B",x="Sampling time-point",y="Resistant allele frequency",col="Replicate") +
     ylim(0,1)+
     facet_grid(TREATMENT~POS)+
     theme_bw()+theme(text = element_text(size=10))

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

bz_stat.test$TREATMENT <- "Benzimidazole"
bz_stat.test <- bz_stat.test %>%
  mutate(POS = str_replace(POS, c("7029569","7029790"), c("Phe167Tyr","Phe200Tyr")))

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


control_stat.test$TREATMENT <- "Untreated"
control_stat.test <- control_stat.test %>%
  mutate(POS = str_replace(POS, c("7029569","7029790"), c("Phe167Tyr","Phe200Tyr")))

p.data <- dplyr::bind_rows(control_stat.test, bz_stat.test)


# make new plot with p values annotated on it
plot_b <- plot +
     geom_text(data=p.data, aes(x=1.5, y=0.95, group=POS, label = paste('P = ',p.adj)),size=3)
```

#-----------------------------------------------------------------------------------------

# isotype 2 in US farm data

```shell
working dir:
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/US_FIELD/VCF

echo "hcontortus_chr2_Celeg_TT_arrow_pilon 13435823" > btub.positions

cat bam.list > samples.list

vcftools --vcf 2.hcontortus_chr2_Celeg_TT_arrow_pilon.snpeff.vcf --keep samples.list --positions btub.positions --extract-FORMAT-info AD --out us_farms_btub2

for i in `ls *AD.FORMAT`; do
      grep "^hcon" ${i} | awk -F '[\t,]' '{print $1,$2,$4/($3+$4),$6/($5+$6),$8/($7+$8),$10/($9+$10),$12/($11+$12),$14/($13+$14),$16/($15+$16),$18/($17+$18),$20/($19+$20),$22/($21+$22)}' OFS="\t" > ${i%.AD.FORMAT}.ADfreq;
done
```
```R
R
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(rstatix)
library(ggrepel)

us_btub2 <- read.table("us_farms_btub2.ADfreq")
colnames(us_btub2) <- c("CHR","POS","Farm 1","Farm 2","Farm 3","Farm 4","Farm 5","Farm 6","Farm 7","Farm 8","Farm 9","Farm 10")
us_btub2 <- melt(us_btub2, id = c("CHR", "POS"), variable.name = "SAMPLE_ID")

bz_conc <- c(1,19.93,15,7.71,15,61.57,29.6,NA,9.67,29.6)

us_btub2$BZ_CONCENTRATION <- bz_conc
colnames(us_btub2) <- c("CHR","POS","SAMPLE_ID","ALLELE_FREQ","BZ_CONCENTRATION")

us_btub2 <- us_btub2 %>%
  mutate(POS = str_replace(POS, c("13435823"), c("Glu198Val")))

# calculate correlation coefficient between alllele frequency and concentration
af_bz_cor <- cor.test(us_btub2$ALLELE_FREQ, us_btub2$BZ_CONCENTRATION, method = "pearson", use = "complete.obs")

# make the plot
plot_c <- ggplot(us_btub2)+
     geom_smooth(aes(BZ_CONCENTRATION,ALLELE_FREQ),method='lm',col='grey')+
     geom_jitter(aes(BZ_CONCENTRATION,ALLELE_FREQ,col=SAMPLE_ID),size=2)+
     geom_text(aes(10,0.95,label=paste('r = ',signif(af_bz_cor$estimate,3),'\n','P = ',signif(af_bz_cor$p.value,3))),size=3)+
     geom_text_repel(aes(BZ_CONCENTRATION,ALLELE_FREQ,label=SAMPLE_ID,col=SAMPLE_ID),size=3.5)+
     labs(title="C",y="Variant Allele Frequency",x="Benzimidazole concentration",col="US farm ID") +
     ylim(-0.05,1)+
     facet_grid(.~POS)+
     theme_bw()+
     theme(legend.position = "none",text = element_text(size=10))
```

```R
library(patchwork)
plot_a / (plot_b | plot_c)
ggsave("Figure_benzimidazole.pdf", useDingbats=FALSE,width=170,height=140,units="mm")
```


#-----------------------------------------------------------------------------------------

# supplement - isotype 1 data from US farms
working dir:
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/US_FIELD/VCF

cp ../../XQTL_BZ/btub.positions btub1.positions


vcftools --vcf 1.hcontortus_chr1_Celeg_TT_arrow_pilon.snpeff.vcf --keep samples.list --positions btub1.positions --extract-FORMAT-info AD --out us_farms_btub1

grep "^hcon" us_farms_btub1.AD.FORMAT | awk -F '[\t,]' '{print $1,$2,$4/($3+$4),$6/($5+$6),$8/($7+$8),$10/($9+$10),$12/($11+$12),$14/($13+$14),$16/($15+$16),$18/($17+$18),$20/($19+$20),$22/($21+$22)}' OFS="\t" > us_farms_btub1.ADfreq

```R
R
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(rstatix)

us_btub1<- read.table("us_farms_btub1.ADfreq")
colnames(us_btub1) <- c("CHR","POS","Farm 1","Farm 2","Farm 3","Farm 4","Farm 5","Farm 6","Farm 7","Farm 8","Farm 9","Farm 10")
us_btub1 <- melt(us_btub1, id = c("CHR", "POS"), variable.name = "SAMPLE_ID")
colnames(us_btub1) <- c("CHR","POS","SAMPLE_ID","ALLELE_FREQ")
us_btub1 <- us_btub1 %>%
  mutate(POS = str_replace(POS, c("7029569","7029790"), c("Phe167Tyr","Phe200Tyr")))

ggplot(us_btub1,aes(x=SAMPLE_ID,y=ALLELE_FREQ,fill=factor(POS)))+
     geom_bar(position="dodge", stat="identity")+
     labs(title="A", x="Sampling location", y="Resistant allele frequency", fill="Variant")+
     theme_bw()+theme(text = element_text(size=10))

ggsave("FigureSX_USfarm_btub1.pdf", useDingbats=FALSE,width=170,height=100,units="mm")
```




#-----------------------------------------------------------------------------------------

# supplement - ivermectin selection on btub P200


working dir:
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_IVM


grep "pre" bam.list  > ivm_pretreatment_samples.list
grep "post" bam.list  > ivm_posttreatment_samples.list


cp ../XQTL_BZ/btub.positions .

vcftools --vcf XQTL_IVM.raw.snpeff.vcf --keep ivm_pretreatment_samples.list --positions btub.positions --extract-FORMAT-info AD --out ivm_btub_pretreatment
vcftools --vcf XQTL_IVM.raw.snpeff.vcf --keep ivm_posttreatment_samples.list --positions btub.positions --extract-FORMAT-info AD --out ivm_btub_posttreatment




for i in `ls ivm*AD.FORMAT`; do
      grep "^hcon" ${i} | awk -F '[\t,]' '{print $1,$2,$4/($3+$4),$6/($5+$6),$8/($7+$8),$10/($9+$10)}' OFS="\t" > ${i%.AD.FORMAT}.ADfreq;
done

```R
R
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(rstatix)

pre <- read.table("ivm_btub_pretreatment.ADfreq")
colnames(pre) <- c("CHR","POS","R1","R1.2","R2","R3")
pre <- melt(pre, id = c("CHR", "POS"), variable.name = "SAMPLE_ID")

post <- read.table("ivm_btub_posttreatment.ADfreq")
colnames(post) <- c("CHR","POS","R1","R1.2","R2","R3")
post <- melt(post, id = c("CHR", "POS"), variable.name = "SAMPLE_ID")

data <- dplyr::full_join(pre, post, by = c("CHR","POS","SAMPLE_ID"))
data$TREATMENT <- "Ivermectin"
colnames(data) <- c("CHR","POS","SAMPLE_ID","PRE_TREATMENT","POST_TREATMENT","TREATMENT")



# change the labels
data <- data %>%
  mutate(POS = str_replace(POS, c("7029569","7029790"), c("Phe167Tyr","Phe200Tyr")))



# make the plot
plot_a <- ggplot(data) +
     geom_segment(aes(x="1.PRE", xend="2.POST", y=PRE_TREATMENT, yend=POST_TREATMENT,col=factor(SAMPLE_ID),group=POS), size=1) +
     labs(title="A",x="Sampling time-point",y="Resistant allele frequency",col="Replicate") +
     ylim(-0.05,1.05)+
     facet_grid(TREATMENT~POS)+
     theme_bw()+theme(text = element_text(size=10))



# perform pairwise t tests between pre/post for each SNP on BZ treated samples
data_stats <- data %>%
  gather(key = "TREATMENT", value = "FREQ", PRE_TREATMENT, POST_TREATMENT)

stat.test <- data_stats %>%
    group_by(POS) %>%
    pairwise_t_test(
      FREQ ~ TREATMENT, paired = TRUE,
      p.adjust.method = "bonferroni"
      ) %>%
    select(-df, -statistic, -p) # Remove details

stat.test$TREATMENT <- "Ivermectin"

p.data <- stat.test


# make new plot with p values annotated on it
plot_a <- plot_a +
     geom_text(data=p.data, aes(x=1.5, y=0.95, group=POS, label = paste('P = ',p.adj)),size=3)
```


#-------------------------------------------------------------------------------

# correlation between btubulin isotype 1 and ivermectin concentration

working dir:
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/US_FIELD/VCF

cp ../../XQTL_BZ/btub.positions btub1.positions

cat bam.list > samples.list

vcftools --vcf 1.hcontortus_chr1_Celeg_TT_arrow_pilon.snpeff.vcf --keep samples.list --positions btub1.positions --extract-FORMAT-info AD --out us_farms_btub1

grep "^hcon" us_farms_btub1.AD.FORMAT | awk -F '[\t,]' '{print $1,$2,$4/($3+$4),$6/($5+$6),$8/($7+$8),$10/($9+$10),$12/($11+$12),$14/($13+$14),$16/($15+$16),$18/($17+$18),$20/($19+$20),$22/($21+$22)}' OFS="\t" > us_farms_btub1.ADfreq

```R
R
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(rstatix)
library(ggrepel)

us_btub1 <- read.table("us_farms_btub1.ADfreq")
colnames(us_btub1) <- c("CHR","POS","Farm 1","Farm 2","Farm 3","Farm 4","Farm 5","Farm 6","Farm 7","Farm 8","Farm 9","Farm 10")


us_btub1 <- melt(us_btub1, id = c("CHR", "POS"), variable.name = "SAMPLE_ID")

ivm_conc <- c(1.51,1.51,12.04,12.04,9.15,9.15,11.27,11.27,297.9,297.9,619,619,312.1,312.1,NA,NA,0.977,0.977,259.4,259.4)



us_btub1$IVM_CONCENTRATION <- ivm_conc
colnames(us_btub1) <- c("CHR","POS","SAMPLE_ID","ALLELE_FREQ","IVM_CONCENTRATION")

us_btub1 <- us_btub1 %>%
  mutate(POS = str_replace(POS, c("7029569","7029790"), c("Phe167Tyr","Phe200Tyr")))

# calculate correlation coefficient between alllele frequency and concentration
af_ivm_cor <- cor.test(us_btub1$ALLELE_FREQ, us_btub1$IVM_CONCENTRATION, method = "pearson", use = "complete.obs")

P167 <- us_btub1[us_btub1$POS!="Phe167Tyr",]
P200 <- us_btub1[us_btub1$POS!="Phe200Tyr",]

P167cor <- cor.test(P167$ALLELE_FREQ, P167$IVM_CONCENTRATION, method = "pearson", use = "complete.obs")
P200cor <- cor.test(P200$ALLELE_FREQ, P200$IVM_CONCENTRATION, method = "pearson", use = "complete.obs")

P167cor.data <- as.data.frame(P167cor$estimate)
P167cor.data$pvalue <- P167cor$p.value
colnames(P167cor.data) <- c("COR","PVALUE")

P200cor.data <- as.data.frame(P200cor$estimate)
P200cor.data$pvalue <- P200cor$p.value
colnames(P200cor.data) <- c("COR","PVALUE")

cor.data <- dplyr::bind_rows(P167cor.data, P200cor.data)
cor.data$CHR <- "hcontortus_chr1_Celeg_TT_arrow_pilon"
cor.data$POS <-  c("Phe167Tyr","Phe200Tyr")

colnames(cor.data) <- c("COR","PVALUE","CHR","POS")

# make the plot
plot_b <- ggplot(us_btub1)+
     geom_smooth(aes(IVM_CONCENTRATION,ALLELE_FREQ),method='lm',col='grey')+
     geom_jitter(aes(IVM_CONCENTRATION,ALLELE_FREQ,col=SAMPLE_ID),size=2)+
     geom_text_repel(aes(IVM_CONCENTRATION,ALLELE_FREQ,label=SAMPLE_ID,col=SAMPLE_ID),size=2)+
     labs(title="B",y="Resistant Allele Frequency",x="Ivermectin concentration",col="US farm ID") +
     ylim(-0.05,1.05)+
     facet_grid(.~POS)+
     theme_bw()+theme(legend.position = "none",text = element_text(size=10))


plot_b <- plot_b + geom_text(data=cor.data, aes(x=500, y=1, group=POS, label = paste('r = ',signif(COR,3),'\n','P = ',signif(PVALUE,3))),size=3)

library(patchwork)
plot_a + plot_b + plot_layout(ncol=1)

ggsave("FigureSX_USfarm_btub1vsIVM.pdf", useDingbats=FALSE, width=170, height=150, units="mm")
```


#-------------------------------------------------------------------------------
# Figure 4 - levamisole analyses

working dir:
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/05_ANALYSIS/LEV

library(ggplot2)
library(dplyr)




lev <- read.table("XQTL_LEV.merged.fst",header=F)
lev <- lev[lev$V1!="hcontortus_chr_mtDNA_arrow_pilon",]
lev <- dplyr::select(lev, V1, V2, V13)
colnames(lev) <- c("CHR","POS","FST")
lev$LABEL <- "3. F2: Levamisole"
lev$ROW_ID <- 1:nrow(lev)


lev_chr4 <- lev[lev$CHR=="hcontortus_chr4_Celeg_TT_arrow_pilon",]
lev_chr5 <- lev[lev$CHR=="hcontortus_chr5_Celeg_TT_arrow_pilon",]

data <- dplyr::bind_rows(lev_chr4,lev_chr5)

# make the XQTL plot
chr_colours<-c("cornflowerblue")

# genome wide signficance per sample
data_gws <- data %>%
    group_by(CHR) %>%
    summarise(GWS = mean(FST)+3*sd(FST))


peaks <- read.table("XQTL_LEV.peak_windows",header=T)

genes <- read.table("candidategenes.data", header=T)

plot_a <- ggplot(data)+
     geom_hline(data = data_gws, aes(yintercept=GWS), linetype="dashed",col="black")+
     geom_vline(data = genes, aes(xintercept=POS),linetype="dotted",col="darkgrey",size=1)+
     geom_point(aes(POS, FST, group=CHR),size=0.5,col="cornflowerblue")+
     geom_point(data = subset(data,(data$POS >= peaks$PEAK_START_COORD[1]) & (data$POS <= peaks$PEAK_END_COORD[1]) & (FST > data_gws$GWS[1])),aes(POS,FST),col="red",size=1)+
     geom_point(data = subset(data,(data$POS >= peaks$PEAK_START_COORD[2]) & (data$POS <= peaks$PEAK_END_COORD[2]) & (FST > data_gws$GWS[2])),aes(POS,FST),col="red",size=1)+
     geom_point(data = subset(data,(data$POS >= peaks$PEAK_START_COORD[3]) & (data$POS <= peaks$PEAK_END_COORD[3]) & (FST > data_gws$GWS[2])),aes(POS,FST),col="red",size=1)+
     geom_point(data = subset(data,(data$POS >= peaks$PEAK_START_COORD[4]) & (data$POS <= peaks$PEAK_END_COORD[4]) & (FST > data_gws$GWS[2])),aes(POS,FST),col="red",size=1)+
     ylim(0,0.1)+xlim(0,50e6)+
     labs(title="A",x="Chromosome position (5 kbp window)", y="Genetic differentiation (Fst)")+
     theme_bw()+theme(legend.position="none",text = element_text(size=10))+
     facet_grid(CHR~.)



#------


## 03 - Gene model plotter <a name="gene_model_plotter"></a>

Working environment and data

```shell
working dir:
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/05_ANALYSIS/LEV

#gff
ln -fs /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/TRANSCRIPTOME_CURATION/20200407/UPDATED_annotation.gff3 ANNOTATION.gff


```



```R
#install.packages("data.table")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ggbio", version = "3.8")

library(data.table)
library(ggplot2)
library(dplyr)

# load data
# --- genome annotation
gff<-fread(cmd="grep ^hcontortus ANNOTATION.gff")
colnames(gff) <- c("chr","source","feature","start","end","point1","strand","frame","info")


# select gene ID - ACR-8
gene1='HCON_00151270'
gene2='HCON_00151260'



# filter data to select chromosome and mRNA
mrna1_data <- gff[grep(gene1, gff$info), ]
mrna1_data <- mrna1_data[mrna1_data$feature=='mRNA',]
mrna1_data <- cbind(mrna1_data, read.table(text = as.character(mrna1_data$info), sep = ";"))
mrna1_data <- cbind(mrna1_data, read.table(text = as.character(mrna1_data$V3), sep = "=",col.names=c("ID","unique_ID")))
mrna1_id <- head(data.frame(mrna1_data$unique_ID),1)  # gives 1st isoform if multiple
chromosome <- mrna1_data[1,1]
colnames(chromosome)<-c("chromosome_ID")

data1 <- gff[grep(mrna1_id$mrna1_data.unique_ID, gff$info), ]

# filter by feature type
cds1 <- data1[data1$feature=="CDS",]
mrna1 <- data1[data1$feature=="mRNA",]
#gene <- data[data$feature=="gene",]


intron1<-data.frame(head(cds1$end,-1),tail(cds1$start,-1),(tail(cds1$start,-1)-head(cds1$end,-1))/2)
colnames(intron1)<-c("start","end","midpoint")

utr5<-data.frame(head(mrna1$start,1),head(sort(cds1$start),1))
colnames(utr5)<-c("start","end")
utr3<-data.frame(head(mrna1$end,1),tail(sort(cds1$end),1))
colnames(utr3)<-c("start","end")
utr1<-rbind(utr5,utr3)


mrna2_data <- gff[grep(gene2, gff$info), ]
mrna2_data <- mrna2_data[mrna2_data$feature=='mRNA',]
mrna2_data <- cbind(mrna2_data, read.table(text = as.character(mrna2_data$info), sep = ";"))
mrna2_data <- cbind(mrna2_data, read.table(text = as.character(mrna2_data$V3), sep = "=",col.names=c("ID","unique_ID")))
mrna2_id <- head(data.frame(mrna2_data$unique_ID),1)  # gives 1st isoform if multiple
chromosome <- mrna2_data[1,1]
colnames(chromosome)<-c("chromosome_ID")

data2 <- gff[grep(mrna2_id$mrna2_data.unique_ID, gff$info), ]

# filter by feature type
cds2 <- data2[data2$feature=="CDS",]
mrna2 <- data2[data2$feature=="mRNA",]
#gene <- data[data$feature=="gene",]


intron2<-data.frame(head(cds2$end,-1),tail(cds2$start,-1),(tail(cds2$start,-1)-head(cds2$end,-1))/2)
colnames(intron2)<-c("start","end","midpoint")

utr52<-data.frame(head(mrna2$start,1),head(sort(cds2$start),1))
colnames(utr52)<-c("start","end")
utr32<-data.frame(head(mrna2$end,1),tail(sort(cds2$end),1))
colnames(utr32)<-c("start","end")
utr2<-rbind(utr52,utr32)


# deletion coordinates (size)
# 31527022 to 31527119 (97 bp) / 31527121 (99 bp)


# make plot
plot_b <- ggplot()+
#gene1
  geom_rect(data=utr1,aes(xmin=utr1$start,ymin=0.5,xmax=utr1$end,ymax=1.5),fill=NA,col="grey",size=0.4)+
  geom_segment(data=intron1,aes(x=intron1$start,xend=intron1$start+intron1$midpoint,y=1,yend=1.5),size=0.5)+
  geom_segment(data=intron1,aes(x=intron1$start+intron1$midpoint,xend=intron1$end,y=1.5,yend=1),size=0.5)+
  geom_rect(data=cds1,aes(xmin=cds1$start,ymin=0.5,xmax=cds1$end,ymax=1.5),fill="grey",col=NA)+
  geom_text(aes(x=mrna1$end+(0.15*(mrna1$end-mrna1$start)), y=1, label = gene1))+
  geom_segment(aes(x=31521884,xend=31521884,y=0.5,yend=1.5),size=1,col="orange")+  # Ser168Thr
  geom_text(aes(x=31521884, y=0.4),label = "Ser168Thr")+
#gene2
  geom_rect(data=utr2,aes(xmin=utr2$start,ymin=2,xmax=utr2$end,ymax=3),fill=NA,col="grey",size=0.4)+
  geom_segment(data=intron2,aes(x=intron2$start,xend=intron2$start+intron2$midpoint,y=2.5,yend=3),size=0.5)+
  geom_segment(data=intron2,aes(x=intron2$start+intron2$midpoint,xend=intron2$end,y=3,yend=2.5),size=0.5)+
  geom_rect(data=cds2,aes(xmin=cds2$start,ymin=2,xmax=cds2$end,ymax=3),fill="grey",col=NA)+
  geom_text(aes(x=mrna1$end+(0.15*(mrna1$end-mrna1$start)), y=2.5, label = gene2))+
  # acr-8 indel
  geom_segment(aes(x=31527022,xend=31527121,y=0.5,yend=3),size=1,col="red")+
  geom_text(aes(x=31527022, y=0.4),label = "indel")+
  #geom_rect(aes(xmin=31527022,ymin=0,xmax=31527121,ymax=3),fill=NA,col="grey",size=0.4)+
  # plot layout
  theme_classic()+
  #xlab("Genome position (bp)")+
  labs(title="B", x =paste("Chromosome: ",chromosome," position (bp)"))+
  xlim(mrna1$start-(0.1*(mrna1$end-mrna1$start)),mrna1$end+(0.25*(mrna1$end-mrna1$start)))+
  scale_y_reverse(lim=c(3,0.4))+ scale_fill_discrete(guide=FALSE)+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())



#-----

#extracting allele frequencies of acr-8 pos 31521884 to make figure

```shell
#working dir:
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_LEV


grep "pre" bam.list  > lev_pretreatment_samples.list
grep "post" bam.list  > lev_posttreatment_samples.list

grep "pre" ../XQTL_CONTROL/bam.list  > control_pretreatment_samples.list
grep "post" ../XQTL_CONTROL/bam.list  > control_posttreatment_samples.list

vcftools --vcf XQTL_LEV.raw.snpeff.vcf --keep lev_pretreatment_samples.list --positions acr8.positions --extract-FORMAT-info AD --out lev_pretreatment
vcftools --vcf XQTL_LEV.raw.snpeff.vcf --keep lev_posttreatment_samples.list --positions acr8.positions --extract-FORMAT-info AD --out lev_posttreatment

vcftools --gzvcf ../XQTL_CONTROL/5.hcontortus_chr5_Celeg_TT_arrow_pilon.tmp.vcf.gz --keep control_pretreatment_samples.list --positions acr8.positions --extract-FORMAT-info AD --out control_pretreatment
vcftools --gzvcf ../XQTL_CONTROL/5.hcontortus_chr5_Celeg_TT_arrow_pilon.tmp.vcf.gz --keep control_posttreatment_samples.list --positions acr8.positions --extract-FORMAT-info AD --out control_posttreatment

#where "acr8.positions" contains:
#P 31521884
#hcontortus_chr5_Celeg_TT_arrow_pilon 31521884




for i in `ls lev*AD.FORMAT`; do
      grep "^hcon" ${i} | awk -F '[\t,]' '{print $1,$2,$4/($3+$4),$6/($5+$6),$8/($7+$8),$10/($9+$10)}' OFS="\t" > ${i%.AD.FORMAT}.ADfreq;
done

for i in `ls control*AD.FORMAT`; do
      grep "^hcon" ${i} | awk -F '[\t,]' '{print $1,$2,$4/($3+$4),$6/($5+$6),$8/($7+$8)}' OFS="\t" > ${i%.AD.FORMAT}.ADfreq;
done



cd ~/lustre118_link/hc/XQTL/05_ANALYSIS/LEV

ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_LEV/lev_pretreatment.ADfreq
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_LEV/lev_posttreatment.ADfreq
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_LEV/control_pretreatment.ADfreq
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_LEV/control_posttreatment.ADfreq


R
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(rstatix)

lev_pre <- read.table("lev_pretreatment.ADfreq")
colnames(lev_pre) <- c("CHR","POS","R1","R1.2","R2","R3")
lev_pre <- melt(lev_pre, id = c("CHR", "POS"), variable.name = "SAMPLE_ID")

lev_post <- read.table("lev_posttreatment.ADfreq")
colnames(lev_post) <- c("CHR","POS","R1","R1.2","R2","R3")
lev_post <- melt(lev_post, id = c("CHR", "POS"), variable.name = "SAMPLE_ID")

lev_data <- dplyr::full_join(lev_pre, lev_post, by = c("CHR","POS","SAMPLE_ID"))
lev_data$TREATMENT <- "Levamisole"
colnames(lev_data) <- c("CHR","POS","SAMPLE_ID","PRE_TREATMENT","POST_TREATMENT","TREATMENT")


control_pre <- read.table("control_pretreatment.ADfreq")
colnames(control_pre) <- c("CHR","POS","R1","R2","R3")
control_pre <- melt(control_pre, id = c("CHR", "POS"), variable.name = "SAMPLE_ID")

control_post <- read.table("control_posttreatment.ADfreq")
colnames(control_post) <- c("CHR","POS","R1","R2","R3")
control_post <- melt(control_post, id = c("CHR", "POS"), variable.name = "SAMPLE_ID")

control_data <- dplyr::full_join(control_pre, control_post, by = c("CHR","POS","SAMPLE_ID"))
control_data$TREATMENT <- "Untreated"
colnames(control_data) <- c("CHR","POS","SAMPLE_ID","PRE_TREATMENT","POST_TREATMENT","TREATMENT")


# bring datasets together
data <- dplyr::bind_rows(control_data, lev_data)

# change the labels
data <- data %>%
  mutate(POS = str_replace(POS, c("31521884"), c("Ser168Thr")))



# make the plot
plot <- ggplot(data) +
     geom_segment(aes(x="1.PRE", xend="2.POST", y=PRE_TREATMENT, yend=POST_TREATMENT,col=factor(TREATMENT),group=POS), size=1) +
     labs(title="C",x="Sampling time-point",y="Variant allele frequency",col="Replicate") +
     ylim(0,1)+
     theme_bw()+theme(text = element_text(size=10))
     #facet_grid(TREATMENT~POS)+

# perform pairwise t tests between pre/post for each SNP on BZ treated samples
lev_data_stats <- lev_data %>%
  gather(key = "TREATMENT", value = "FREQ", PRE_TREATMENT, POST_TREATMENT)

lev_stat.test <- lev_data_stats %>%
    group_by(POS) %>%
    pairwise_t_test(
      FREQ ~ TREATMENT, paired = TRUE,
      p.adjust.method = "bonferroni"
      ) %>%
    select(-df, -statistic, -p) # Remove details

lev_stat.test$TREATMENT <- "Levamisole"
lev_stat.test <- lev_stat.test %>%
  mutate(POS = str_replace(POS, c("31521884"), c("Ser168Thr")))

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


control_stat.test$TREATMENT <- "Untreated"
control_stat.test <- control_stat.test %>%
  mutate(POS = str_replace(POS, c("31521884"), c("Ser168Thr")))

p.data <- dplyr::bind_rows(control_stat.test, lev_stat.test)


# make new plot with p values annotated on it
plot_c <- plot +
     geom_text(data=p.data, aes(x=1.5, y=0.90, label = paste('P = ',p.adj[2])),size=3, col="red")+
     geom_text(data=p.data, aes(x=1.5, y=0.15, label = paste('P = ',p.adj[1])),size=3, col="blue")

plot_d <- ggplot()+geom_blank()
plot_e <- ggplot()+geom_blank()

library(patchwork)
plot_a + plot_b + ((plot_c / plot_d) | plot_e) + plot_layout(ncol=1, height=c(3,1,3))






ln -s ../../03_MAPPING/XQTL/XQTL_F3_L3_n200_LEV_post_01_21395_3_2.merged.sorted.marked.realigned.bam post.bam
ln -s ../../03_MAPPING/XQTL/XQTL_F3_L3_n200_LEV_post_01_21395_3_2.merged.sorted.marked.realigned.bam.bai post.bam.bai



#-------------------------------------------------------------------------------
# acr-8 deletion plot

library(GenomicAlignments)
library(viridis)
library(stringr)

data_pre <- as.data.frame(readGAlignmentPairs("pre.bam", use.names=TRUE, param=ScanBamParam(which=GRanges("hcontortus_chr5_Celeg_TT_arrow_pilon", IRanges(31525841, 31529149)))))
data_pre$number <- rep(1:100,length.out=nrow(data_pre))
data_pre$treatment <- "3.pre"

data_post <- as.data.frame(readGAlignmentPairs("post.bam", use.names=TRUE, param=ScanBamParam(which=GRanges("hcontortus_chr5_Celeg_TT_arrow_pilon", IRanges(31525841, 31529149)))))
data_post$number <- rep(1:100,length.out=nrow(data_post))
data_post$treatment <- "4.post"



data_ise <- as.data.frame(readGAlignmentPairs("ISE.bam", use.names=TRUE, param=ScanBamParam(which=GRanges("hcontortus_chr5_Celeg_TT_arrow_pilon", IRanges(31525841, 31529149)))))
data_ise$number <- rep(1:100,length.out=nrow(data_ise))
data_ise$treatment <- "1.ISE"

data_uga <- as.data.frame(readGAlignmentPairs("UGA.bam", use.names=TRUE, param=ScanBamParam(which=GRanges("hcontortus_chr5_Celeg_TT_arrow_pilon", IRanges(31525841, 31529149)))))
data_uga$number <- rep(1:100,length.out=nrow(data_uga))
data_uga$treatment <- "2.UGA"

#data <- dplyr::bind_rows(data_ise,data_uga,data_pre,data_post)
data <- dplyr::bind_rows(data_pre,data_post)


softclip_match <- "([:digit:]{1,})S"

cigars1 <- data$cigar.first

softclip1_data <- str_extract(cigars1,softclip_match)
softclip1_data <- gsub("S","",softclip1_data)
softclip1_data[is.na(softclip1_data)] <- 0
softclip1_data <- as.numeric(unlist(softclip1_data))

data$softclip1_length <- softclip1_data


cigars2 <- data$cigar.last
softclip2_match <- "([:digit:]{1,})S"
softclip2_data <- str_extract(cigars2,softclip_match)
softclip2_data <- gsub("S","",softclip2_data)
softclip2_data[is.na(softclip2_data)] <- 0
softclip2_data <- as.numeric(unlist(softclip2_data))

data$softclip2_length <- softclip2_data



plot_e <- ggplot(data)+
          geom_rect(aes(xmin=start.first,ymin=number-0.4,xmax=end.first,ymax=number+0.4,fill=softclip1_length))+
          geom_rect(aes(xmin=start.last,ymin=number-0.4,xmax=end.last,ymax=number+0.4,fill=softclip2_length))+
          #geom_curve(aes(x=start.first,xend=end.last+1,y=number,yend=number+0.1),curvature = -0.05,col="grey")+
          #xlim(31526800, 31527250)+
          coord_cartesian(xlim=c(31526850, 31527250))+
          theme_bw()+
          labs(x="Chromosome position (bp)")+
          theme(legend.position="none",text = element_text(size=10),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
          facet_grid(treatment~.)+
          scale_fill_viridis(option="plasma")


library(patchwork)
plot_a + plot_b + ((plot_c / plot_d) | plot_e) + plot_layout(ncol=1, height=c(3,1,3))


# deletion coordinates (size)
31527022 to 31527119 (97 bp) / 31527121 (99 bp)
