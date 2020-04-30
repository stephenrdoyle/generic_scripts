# run_bionano_stitch

#######################################
#--- Input information
prefix=Tc_test
genome=canu1.9_purge_redundans.fa
enzyme=BspQI
cmap=/nfs/users/nfs_s/sd21/lustre118_link/tc/raw/bionano/exp_refineFinal1_contigs.cmap
########################################
export PATH=$PATH:/nfs/users/nfs_s/sd21/lustre118_link/software/GENOME_SCAFFOLDING/Irys-scaffolding/KSU_bioinfo_lab/stitch


##--- Make CMAP file for reference
#perl ~/sd21_lustre_link/software/Irys-scaffolding/KSU_bioinfo_lab/assemble_XeonPhi/third-party/fa2cmap_multi.pl -v -i Sp34.v4.1.genome_1.fa -e BspQI

##--- run "sewing machine / stitch"
#sewing_machine.pl -o test_out2 -g exp_refineFinal1_contigs.cmap -p test -e BspQI -f Sp34.v4.1.genome_1.fa -r Sp34.v4.1.genome_1_BspQI.cmap

#--- generate a final report
#perl ~/sd21_lustre_link/software/Irys-scaffolding/KSU_bioinfo_lab/assemble_XeonPhi/write_report.pl -o test_out2 -g exp_refineFinal1_contigs.cmap -p test -e BspQI -f Sp34.v4.1.genome_1.fa
 -r Sp34.v4.1.genome_1_BspQI.cmap






#--- Make CMAP file for reference
perl /nfs/users/nfs_s/sd21/lustre118_link/software/GENOME_SCAFFOLDING/Irys-scaffolding/KSU_bioinfo_lab/assemble_XeonPhi/third-party/fa2cmap_multi.pl -v -i ${genome} -e ${enzyme} -M 10

#--- run "sewing machine / stitch"
sewing_machine.pl -o ${prefix}_out -g ${cmap} -p ${prefix} -e ${enzyme} -f ${genome} -r ${genome%.fa*}_${enzyme}.cmap

#--- generate a final report
perl /nfs/users/nfs_s/sd21/lustre118_link/software/GENOME_SCAFFOLDING/Irys-scaffolding/KSU_bioinfo_lab/assemble_XeonPhi/write_report.pl -o ${prefix}_out -g ${cmap} -p ${prefix} -e ${enzyme} -f ${genome}  -r ${genome%.fa*}_${
enzyme}.cmap






# bionano hybrid scaffold
module load ISG/perl/5.30-0

perl /nfs/users/nfs_s/sd21/lustre118_link/software/GENOME_SCAFFOLDING/bionano/scripts/scripts_fresh/HybridScaffold/hybridScaffold.pl -n canu1.9_purge_redundans.fa -b /nfs/users/nfs_s/sd21/lustre118_link/tc/bionano/exp_refineFinal1_contigs.cmap -c ~sd21/bash_scripts/hybridScaffold_config.xml -o hybridscaffold_out
