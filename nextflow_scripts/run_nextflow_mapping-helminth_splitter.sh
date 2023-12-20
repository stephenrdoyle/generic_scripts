

conda deactivate
module load mapping-helminth/v1.0.8


manifest=${1}
reference=${2}

manifest=qtl.manifest
reference=/nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/QTL/01_REFERENCE/HAEM_V4_final.chr.fa

cat qtl.manifest | sed 1,1d | split -l 100 --suffix-length=1 --numeric-suffixes=1 - qtl.manifest_

for i in ${manifest}_*; do
    sed -i -e '1i\'$'\n''ID,R1,R2' ${i};
done


# setup job conditions
JOBS=$( ls -1 ${manifest}_* | wc -l )

#submit job array to call variants put scaffold / contig
bsub -q long -R'span[hosts=1] select[mem>10000] rusage[mem=10000]' -M10000 -J "mapping_[1-$JOBS]%1" -e mapping.e -o mapping.o "mapping-helminth --input qtl.manifest_\${LSB_JOBINDEX} --reference ${reference} --outdir QTL_MAPPING"
