# run_augustus_split_by_contigs
- augustus can be slow to run
- it is, however, possible to run it in chunks of X size, annotate each, and then merge them back together again.
- this is a rough workflow of how it was done, way back when I ran augustus (independently of Braker)

```bash
# prepare reference
ln -s ../../REF/HAEM_V4_final.chr.fa

# split reference sequnece into individual fasta files
/software/pathogen/external/apps/usr/bin/splitMfasta.pl HAEM_V4_final.chr.fa --outputpath=./
for f in *.split.*; do NAME=`grep ">" $f`; mv $f ${NAME#>}.fa; done

# prepare a necessary file with base summaries
/software/pathogen/external/apps/usr/bin/summarizeACGTcontent.pl HAEM_V4_final.chr.fa > basesummary.out


# prepare augustus run files
grep "bases" basesummary.out | awk -v PWD=$PWD -v HINTS=hints.gff '{print PWD"/"$3".fa",PWD"/"HINTS,"1",$1}' OFS="\t" > sequences.list

# make the inidivudal augustrus jobs
/software/pathogen/external/apps/usr/bin/createAugustusJoblist.pl --sequences=sequences.list --wrap="#" --overlap=50000 --chunksize=3000000 --outputdir=augustus_split_out --joblist=jobs.lst --jobprefix=augsplit --command \
"/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/augustus_v3.3/bin/augustus \
--species=Hc_V4_chr --strand=both --genemodel=partial --protein=on --introns=on --start=on --stop=on --cds=on --codingseq=on --UTR=off --nc=off --gff3=on --alternatives-from-evidence=on \
--extrinsicCfgFile=/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/augustus_v3.3/config/species/Hc_V4_chr/extrinsic.Hc_V4_chr.modified.cfg \
--AUGUSTUS_CONFIG_PATH=/nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/augustus_v3.3/config/"


# run augustus
mkdir augustus_split_out

for i in augsplit*; do
     echo -e "bsub.py 4 augsplit_log ./${i}" >> run_augsplit;
done        #Max memory used was 1.2Gb

chmod a+x run_augsplit
bsub.py --queue yesterday 1 run_splitter ./run_augsplit

mkdir augsplit_run_dir
mv augsplit* augsplit_run_dir

# merge the individual gffs together
cat augustus_split_out/*gff | /software/pathogen/external/apps/usr/bin/join_aug_pred.pl > HC_V4_augustus_merge.gff

# make a nice GFF
echo -e "##gff-version 3" > HC_V4_augustus_merge.filtered.gff; grep "AUGUSTUS" HC_V4_augustus_merge.gff >> HC_V4_augustus_merge.filtered.gff
```
