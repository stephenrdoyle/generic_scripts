
# playtime


#install.packages("ggmsa")
#install.packages("seqmagick")
## ## loading the package
library("ggmsa")
library(seqmagick)


sequences <- system.file("extdata", "wbp_acr8.aln",)
print(sequences)


ggmsa( "test", 100, 200, font = NULL, color = "Chemistry_AA")


WBGene00036931|cabrigprjna10731   Caenorhabditis_briggsae
Cni-acr-8|canigoprjna384657  Caenorhabditis_nigoni
WBGene00058468|caremaprjna53967   Caenorhabditis_remanei
Csp11.Scaffold629.g13743|catropprjna53597   Caenorhabditis_tropicalis
ACAC_0000398801|ancantprjeb493    Angiostrongylus_cantonensis
ANCDUO_01121|anduodprjna72581     Ancylostoma_duodenale
ACOC_0000247401|ancostprjeb494    Angiostrongylus_costaricensis
Cang_2012_03_13_00065.g3190|caangaprjna51225     Caenorhabditis_angaria
ANCCAN_17789|ancaniprjna72585     Ancylostoma_caninum
Csp5_scaffold_00005.g324|casiniprjna194557  Caenorhabditis_sinica
CGOC_0000275601|cygoldprjeb498    Cylicostephanus_goldi
arcanus-mkr-S_1590-0.0-mRNA-1|prarcaprjeb27334   Pristionchus_arcanus
Acey_s0072.g653|anceylprjna231479 Ancylostoma_ceylanicum
WBGene00139702|cabrenprjna20035   Caenorhabditis_brenneri
HPOL_0000787701|hepolyprjeb15396  Heligmosomoides_polygyrus
DCO_008622|dicoroprjdb3143   Diploscapter_coronatus
FL83_01328|calateprjna248912 Caenorhabditis_latens
HPLM_0000547301|haplacprjeb509    Haemonchus_placei
fissidentatus-mkr-S40-3.15-mRNA-1|prfissprjeb27334    Pristionchus_fissidentatus
DICVIV_07636|diviviprjna72587     Dictyocaulus_viviparus
entomophagus-mkr-S570-0.22-mRNA-1|prentoprjeb27334    Pristionchus_entomophagus
maxplancki-mkr-S100-2.104-mRNA-1|prmaxpprjeb27334     Pristionchus_maxplancki
HCON_00151270|hacontprjeb506 Haemonchus_contortus
exspectatus-mkr-S_1791-0.14-mRNA-1|prexspprjeb24288   Pristionchus_exspectatus
japonicus-mkr-S267-0.26-mRNA-1|prjapoprjeb27334  Pristionchus_japonicus
OESDEN_01101|oedentprjna72579     Oesophagostomum_dentatum
MicoRS5524-ag_msk-S79-2.36-mRNA-1|mijapoprjeb27334    Micoletzkya_japonica
WBGene00000047|caelegprjna13758   Caenorhabditis_elegans
WBGene00092568|prpaciprjna12644   Pristionchus_pacificus
SVUK_0001037501|stvulgprjeb531    Strongylus_vulgaris
WR25_23289|dipachprjna280107 Diploscapter_pachys
Sp34_X0071520|casp34prjdb5687     Caenorhabditis_inopinata
NECAME_16202|neamerprjna72135     Necator_americanus
mbelari.g26242|mebelaprjeb30104   Mesorhabditis_belari
OTIPU.nOt.2.0.1.g01642|ostipuprjeb15512     Oscheius_tipulae
nAv.1.0.1.g09901|acviteprjeb1697  Acanthocheilonema_viteae
mayeri-mkr-S55-0.24-mRNA-1|prmayeprjeb27334 Pristionchus_mayeri
NBR_0000568801|nibrasprjeb511     Nippostrongylus_brasiliensis



# substitute names
while read old new; do sed -i "s/${old}/${new}/" wbp_acr8.facopy; done < rename

# extract fasta sequences based on the new names to a new fasta
while read old new; do samtools faidx wbp_acr8.facopy ${new}; done < rename > wb_cladeV_acr8.fa


module load mafft/7.407=1

mafft --globalpair --maxiterate 16 --reorder "wb_cladeV_acr8.fa" > "wb_cladeV_acr8.aln"



# removed
entomophagus-mkr-S570-0.22-mRNA-1|prentoprjeb27334    Pristionchus_entomophagus - large insertion

Pristionchus_exspectatus
Cylicostephanus_goldi
Dictyocaulus_viviparus
Ancylostoma_caninum
Caenorhabditis_sinica
Caenorhabditis_brenneri


library(ggmsa)
library(ggplot2)

data <- tidy_msa("wb_cladeV_acr8.aln", start = 225, end = 275)
data[ data == "-" ] <- NA


colours <- c("A"="#CD2127","G"="#CD2127","I"="#CD2127","L"="#CD2127","P"="#CD2127","V"="#CD2127","F"="#AAC64B","W"="#AAC64B","Y"="#AAC64B","D"="#ED6823","E"="#ED6823","R"="#759CD1","H"="#759CD1","K"="#759CD1","S"="#F09AC1","T"="#F09AC1","C"="#EDAB20","M"="#EDAB20","N"="#283983","Q"="#283983","NA"="white")


ggplot(data,aes(y=name, x=position)) +
     geom_tile(aes(fill=character,alpha=0.5),col="grey",na.rm = TRUE) +
     geom_text(aes(label=character),col="black")+
     theme_minimal()+theme(legend.position = "none",text = element_text(size=10))+
     labs(y="Species")+
     scale_fill_manual(values = colours)+
     coord_fixed()






Aliphatic #CD2127   A
Aliphatic #CD2127   G
Aliphatic #CD2127   I
Aliphatic #CD2127   L
Aliphatic #CD2127   P
Aliphatic #CD2127   V
Aromatic  #AAC64B    F
Aromatic  #AAC64B    W
Aromatic  #AAC64B    Y
Acidic    #ED6823   D
Acidic    #ED6823   E
Basic     #759CD1   R
Basic     #759CD1   H
Basic     #759CD1   K
Hydroxylic     #F09AC1   S
Hydroxylic     #F09AC1   T
Sulfur-containing   #EDAB20   C
Sulfur-containing   #EDAB20   M
Amidic    #283983   N
Amidic    #283983   Q
colours <- c("A"="#CD2127","G"="#CD2127","I"="#CD2127","L"="#CD2127","P"="#CD2127","V"="#CD2127","F"="#AAC64B","W"="#AAC64B","Y"="#AAC64B","D"="#ED6823","E"="#ED6823","R"="#759CD1","H"="#759CD1","K"="#759CD1","S"="#F09AC1","T"="#F09AC1","C"="#EDAB20","M"="#EDAB20","N"="#283983","Q"="#283983")
