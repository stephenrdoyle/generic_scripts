# playtime

extracting allele frequencies of beta tubulin P167,P198, P200 to make figure


working dir:
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_BZ


vcftools --vcf XQTL_BZ.raw.snpeff.vcf --positions btub.positions --extract-FORMAT-info AD

where "btub.positions" contains:
#P167
hcontortus_chr1_Celeg_TT_arrow_pilon 7029569
#P198
hcontortus_chr1_Celeg_TT_arrow_pilon 7029788
# P200
hcontortus_chr1_Celeg_TT_arrow_pilon 7029790



grep "^hcon" out.AD.FORMAT | awk -F '[\t,]' '{print $1,$2,$4/($3+$4),$6/($5+$6),$8/($7+$8),$10/($9+$10),$12/($11+$12),$14/($13+$14),$16/($15+$16),$18/($17+$18)}' OFS="\t" > bz.freq

hcontortus_chr1_Celeg_TT_arrow_pilon	7029569	0.125	0.197368	0.223301	0.264463	0.15	0.128571	0.131148	0.155172
hcontortus_chr1_Celeg_TT_arrow_pilon	7029790	0.721311	0.805556	0.866071	0.926316	0.534884	0.477612	0.446281	0.347368



awk '{print $1,$2,$7,$8,$9,$10,"bz_pre"}' OFS="\t" bz.freq > bz.freq2
awk '{print $1,$2,$3,$4,$5,$6,"bz_post"}' OFS="\t" bz.freq >> bz.freq2


grep "^hcon" out.AD.FORMAT | awk -F '[\t,]' '{for(i=3;i<=NF;i+=2) print $1,$2,(i+=1)/($i+(i+=1))}'


for ((i=3,j=4;i<=j;i+=2,j+=2)); do \
     grep "^hcon" out.AD.FORMAT | \
     awk -v i=$i -v j=$j -F '[\t,]' '{if(i<NF) print $1,$2,$j/($i+$j)}'
done
