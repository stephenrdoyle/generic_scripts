
# playtime

# downloaded protein sequences from WBP frmo ortholog set of C. elegans acr-8

WBGene00036931|cabrigprjna10731   Caenorhabditis_briggsae
Cni-acr-8|canigoprjna384657  Caenorhabditis_nigoni
WBGene00058468|caremaprjna53967   Caenorhabditis_remanei
Csp11.Scaffold629.g13743|catropprjna53597   Caenorhabditis_tropicalis
ACAC_0000398801|ancantprjeb493    Angiostrongylus_cantonensis
ANCDUO_01121|anduodprjna72581     Ancylostoma_duodenale
ACOC_0000247401|ancostprjeb494    Angiostrongylus_costaricensis
Cang_2012_03_13_00065.g3190|caangaprjna51225     Caenorhabditis_angaria
arcanus-mkr-S_1590-0.0-mRNA-1|prarcaprjeb27334   Pristionchus_arcanus
Acey_s0072.g653|anceylprjna231479 Ancylostoma_ceylanicum
HPOL_0000787701|hepolyprjeb15396  Heligmosomoides_polygyrus
DCO_008622|dicoroprjdb3143   Diploscapter_coronatus
FL83_01328|calateprjna248912 Caenorhabditis_latens
HPLM_0000547301|haplacprjeb509    Haemonchus_placei
fissidentatus-mkr-S40-3.15-mRNA-1|prfissprjeb27334    Pristionchus_fissidentatus
maxplancki-mkr-S100-2.104-mRNA-1|prmaxpprjeb27334     Pristionchus_maxplancki
HCON_00151270|hacontprjeb506 Haemonchus_contortus
japonicus-mkr-S267-0.26-mRNA-1|prjapoprjeb27334  Pristionchus_japonicus
OESDEN_01101|oedentprjna72579     Oesophagostomum_dentatum
MicoRS5524-ag_msk-S79-2.36-mRNA-1|mijapoprjeb27334    Micoletzkya_japonica
WBGene00000047|caelegprjna13758   Caenorhabditis_elegans
WBGene00092568|prpaciprjna12644   Pristionchus_pacificus
SVUK_0001037501|stvulgprjeb531    Strongylus_vulgaris
WR25_23289|dipachprjna280107 Diploscapter_pachys
Sp34_X0071520|casp34prjdb5687     Caenorhabditis_inopinata
NECAME_16202|neamerprjna72135     Necator_americanus
OTIPU.nOt.2.0.1.g01642|ostipuprjeb15512     Oscheius_tipulae
nAv.1.0.1.g09901|acviteprjeb1697  Acanthocheilonema_viteae
mayeri-mkr-S55-0.24-mRNA-1|prmayeprjeb27334 Pristionchus_mayeri
NBR_0000568801|nibrasprjeb511     Nippostrongylus_brasiliensis

# removed - poor / missing alignment in region of interest
entomophagus-mkr-S570-0.22-mRNA-1|prentoprjeb27334    Pristionchus_entomophagus
exspectatus-mkr-S_1791-0.14-mRNA-1|prexspprjeb24288   Pristionchus_exspectatus
CGOC_0000275601|cygoldprjeb498    Cylicostephanus_goldi
DICVIV_07636|diviviprjna72587     Dictyocaulus_viviparus
ANCCAN_17789|ancaniprjna72585     Ancylostoma_caninum
Csp5_scaffold_00005.g324|casiniprjna194557  Caenorhabditis_sinica
WBGene00139702|cabrenprjna20035   Caenorhabditis_brenneri
mbelari.g26242|mebelaprjeb30104   Mesorhabditis_belari


# substitute names from transcript ID to species name
while read old new; do sed -i "s/${old}/${new}/" wbp_acr8.fa; done < rename

# extract fasta sequences based on the new names to a new fasta
while read old new; do samtools faidx wbp_acr8.fa ${new}; done < rename > wb_cladeV_acr8.fa

# make an alignment using mafft
module load mafft/7.407=1

mafft --localpair --maxiterate 16 --reorder "wb_cladeV_acr8.fa" > "wb_cladeV_acr8.aln"



# make a plot
R
library(ggmsa)
library(ggplot2)

# define the alignment window - chosen to flank the Ser167Thr position in H. contortus
ALIGNMENT_START=225
ALIGNMENT_END=275

# load data and clean up for plotting
data <- tidy_msa("wb_cladeV_acr8.aln", ALIGNMENT_START, ALIGNMENT_END)


# colour scheme - based on amino acid properties
<!-- Aliphatic #CD2127   A
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
Amidic    #283983   Q -->

colours <- c("A"="#CD2127","G"="#CD2127","I"="#CD2127","L"="#CD2127","P"="#CD2127","V"="#CD2127","F"="#AAC64B","W"="#AAC64B","Y"="#AAC64B","D"="#ED6823","E"="#ED6823","R"="#759CD1","H"="#759CD1","K"="#759CD1","S"="#F09AC1","T"="#F09AC1","C"="#EDAB20","M"="#EDAB20","N"="#283983","Q"="#283983","-"="white")

# make the plot
ggplot(data,aes(y=name, x=position)) +
     geom_tile(aes(fill=character),col="grey",na.rm = TRUE) +
     geom_text(aes(label=character),col="black",size=2)+
     geom_text(aes(x=253,y=31),label="*")+
     geom_text(aes(x=253,y=29.5),label="H. contortus\nSer167Thr",size=3)+
     theme_minimal()+theme(legend.position = "none")+
     labs(y="Species", x="Alignment position",text = element_text(size=10))+
     scale_fill_manual(values = colours)

ggsave("acr8_multiple_sequence_alignment.pdf", width=170,height=180,units="mm")
