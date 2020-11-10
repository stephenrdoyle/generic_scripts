

conda env create -f ~/.nextflow/assets/nf-core/eager/environment.yml

conda activate nf-core-eager-2.1.0



fasta=trichuris_trichiura.fa
R1=sample2.R1.fastq.gz
R2=sample2.R2.fastq.gz
prefix=${R1%.R1.fastq.gz}

cpus=7



# Step 1. PREPROCESSING
#--- Create BWA indices if they are not present

bwa index $fasta


#--- Index Fasta file if not specified on CLI

samtools faidx $fasta


#--- Create Sequence Dictionary for FastA

samtools dict $fasta > ${fasta%.*}.dict



#--- FastQC

fastqc ${R1} ${R2}


#--- FastP
* Optional poly-G complexity filtering step before read merging/adapter clipping etc
* Note: Clipping, Merging, Quality Trimning are turned off here - we leave this to adapter removal itself!

complexity_filter_poly_g=false
complexity_filter_poly_g_min=10

#--- single end
#fastp --in1 ${R1} --out1 "${reads[0].baseName}.pG.fq.gz" -A -g --poly_g_min_len "${complexity_filter_poly_g_min}" -Q -L -w ${task.cpus} --json "${reads[0].baseName}"_fastp.json

#--- paired end
fastp --in1=${R1} --in2=${R2} --out1=${R1%.fastq.gz}.pG.fq.gz --out2=${R2%.fastq.gz}.pG.fq.gz -A -g --poly_g_min_len=${complexity_filter_poly_g_min} -Q -L -w ${cpus} --json "${R1%.R1.fastq.gz}"_fastp.json



# Step 2. Adapter Clipping / Read Merging

# adapter removal
clip_forward_adaptor=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
clip_reverse_adaptor=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
clip_readlength=30
clip_min_read_quality=20
min_adap_overlap=1


mkdir -p output

#--- single end

#AdapterRemoval --file1 ${reads[0]} --basename ${base} --gzip --threads ${task.cpus} --trimns --trimqualities --adapter1 ${clip_forward_adaptor} --adapter2 ${clip_reverse_adaptor} --minlength ${clip_readlength} --minquality ${preserve5p}

#--- paired end

AdapterRemoval --file1 ${R1%.fastq.gz}.pG.fq.gz --file2 ${R2%.fastq.gz}.pG.fq.gz  --basename ${R1%.R1.fastq.gz} --trimns --trimqualities --adapter1 ${clip_forward_adaptor} --adapter2 ${clip_reverse_adaptor} --minlength ${clip_readlength} --minquality ${clip_min_read_quality} --gzip --threads ${cpus}

#cat *.collapsed.gz *.collapsed.truncated.gz *.singleton.truncated.gz *.pair1.truncated.gz *.pair2.truncated.gz > output/${base}.combined.fq.gz

#mv *.settings output/



# Step 3 - Mapping with BWA, SAM to BAM, Sort BAM

#--- PE data

bwa mem -t ${cpus} $fasta ${R1%.R1.fastq.gz}.pair1.truncated.gz ${R2%.R2.fastq.gz}.pair2.truncated.gz -R "@RG\\tID:ILLUMINA-${prefix}\\tSM:${prefix}\\tPL:illumina" | samtools sort -O bam - > ${prefix}.mapped.bam
samtools index ${prefix}.mapped.bam

#--- collapsed or SE data
#bwa mem -t ${task.cpus} $fasta $reads -R "@RG\\tID:ILLUMINA-${prefix}\\tSM:${prefix}\\tPL:illumina" | samtools sort -@ ${task.cpus} -O bam - > "${prefix}".mapped.bam
#samtools index "${size}" -@ ${task.cpus} "${prefix}".mapped.bam

#--- flagstat
samtools flagstat ${prefix}.mapped.bam > ${prefix}_flagstat.stats



# Step 4 - Keep unmapped/remove unmapped reads
bam_mapping_quality_threshold=0

samtools view -h -b ${prefix}.mapped.bam -@ ${cpus} -q ${bam_mapping_quality_threshold} -o ${prefix}.filtered.bam
samtools index ${prefix}.filtered.bam



# Step 5a: DeDup / mark duplicates

bamUtils dedup --in ${prefix}.filtered.bam --out ${prefix}_rmdup.bam.tmp
samtools view -h -b ${prefix}_rmdup.bam.tmp | samtools sort -@ ${cpus} - -o ${prefix}_rmdup.bam
samtools index ${prefix}_rmdup.bam
rm *tmp



# Step 6: Preseq
preseq_step_size=1000

preseq c_curve -s ${preseq_step_size} -o ${prefix}.ccurve -B ${prefix}_rmdup.bam


# Step 7a: DMG Assessment
<!-- damageprofiler_length=100
damageprofiler_threshold=15
damageprofiler_yaxis=0.30

damageprofiler=/nfs/users/nfs_s/sd21/lustre118_link/software/ANCIENT/DamageProfiler/jars/DamageProfiler-1.0-java8.jar

bsub.py 10 dm "java -jar ${damageprofiler} -i ${prefix}_rmdup.bam -r ${fasta} -l ${damageprofiler_length} -t ${damageprofiler_threshold} -o . -yaxis_damageplot ${damageprofiler_yaxis}" -->


# Step 8: Qualimap

qualimap bamqc -bam ${prefix}_rmdup.bam -nt ${cpus} -outdir . -outformat "HTML"



<!-- # Step 9: Bedtools

bedtools coverage -a ${anno_file} -b $bam | pigz -p ${task.cpus} > "${bam.baseName}".breadth.gz
bedtools coverage -a ${anno_file} -b $bam -mean | pigz -p ${task.cpus} > "${bam.baseName}".depth.gz -->


# Step 10: PMDtools
pmdtools_range=10
pmdtools_threshold=3
pmdtools_max_reads=10000

#Run Filtering step
samtools calmd -b ${prefix}_rmdup.bam ${fasta} | samtools view -h - | pmdtools --threshold ${pmdtools_threshold} --header | samtools view -@ ${cpus} -Sb - > ${prefix}.pmd.bam

#Run Calc Range step
samtools calmd -b ${prefix}.pmd.bam ${fasta} | samtools view -h - | pmdtools --deamination --range ${pmdtools_range} $treatment -n ${pmdtools_max_reads} > ${prefix}.cpg.range."${pmdtools_range}".txt
samtools index ${prefix}.pmd.bam


# Step 11 - BAM Trimming step using bamUtils
#--- Can be used for UDGhalf protocols to clip off -n bases of each read
bamutils_clip_left=1
bamutils_clip_right=1


bamUtils trimBam ${prefix}.pmd.bam tmp.bam -L ${bamutils_clip_left} -R ${bamutils_clip_right}
samtools sort -@ ${cpus} tmp.bam -o ${prefix}.trimmed.bam
samtools index ${prefix}.trimmed.bam



# Step 12b: Genotyping - UG
 NB: GATK 3.5 is the last release with VCF output in "old" VCF format, not breaking MVA. Therefore we need it (for now at least until downstream tools can read proper 4.2 VCFs... )

 genotyping_tool = ''
 genotyping_source = "raw"
 gatk_ug_jar = ''
 gatk_ug_genotype_model = 'SNP'
 gatk_hc_emitrefconf = 'GVCF'
 gatk_call_conf = '30'
 gatk_ploidy = '2'
 gatk_downsample = '250'
 gatk_ug_out_mode = 'EMIT_VARIANTS_ONLY'
 gatk_hc_out_mode = 'EMIT_VARIANTS_ONLY'
 gatk_dbsnp = ''
 gatk_ug_defaultbasequalities = ''
 freebayes_C = 1
 freebayes_g = 0
 freebayes_p = 2



java -Xmx${task.memory.toGiga()}g -jar ${jar} -T RealignerTargetCreator -R ${fasta} -I ${bam} -nt ${task.cpus} -o ${bam}.intervals ${defaultbasequalities}
java -Xmx${task.memory.toGiga()}g -jar ${jar} -T IndelRealigner -R ${fasta} -I ${bam} -targetIntervals ${bam}.intervals -o ${bam}.realign.bam ${defaultbasequalities}
java -Xmx${task.memory.toGiga()}g -jar ${jar} -T UnifiedGenotyper -R ${fasta} -I ${bam}.realign.bam -o ${bam}.unifiedgenotyper.vcf -nt ${task.cpus} --genotype_likelihoods_model ${gatk_ug_genotype_model} -stand_call_conf ${params.gatk_call_conf} --sample_ploidy ${params.gatk_ploidy} -dcov ${gatk_downsample} --output_mode ${gatk_ug_out_mode} ${defaultbasequalities}
pigz -p ${task.cpus} ${bam}.unifiedgenotyper.vcf




gatk HaplotypeCaller -R ${fasta} -I ${bam} -O ${bam}.haplotypecaller.vcf -stand-call-conf ${params.gatk_call_conf} --sample-ploidy ${params.gatk_ploidy} --output-mode ${params.gatk_hc_out_mode} --emit-ref-confidence ${params.gatk_hc_emitrefconf}

pigz -p ${task.cpus} ${bam}.haplotypecaller.vcf



#  Step 12c: FreeBayes genotyping, should probably add in some options for users to set
freebayes -f ${fasta} -p ${freebayes_p} -C ${freebayes_C} ${skip_coverage} ${bam} > ${bam.baseName}.freebayes.vcf
pigz -p ${task.cpus} ${bam.baseName}.freebayes.vcf



# Step 13: VCF2Genome
run_vcf2genome = false
vcf2genome_outfile = ''
vcf2genome_header = ''
vcf2genome_minc = 5
vcf2genome_minq = 30
vcf2genome_minfreq = 0.8

//MultiVCFAnalyzer Options
run_multivcfanalyzer = false
write_allele_frequencies = false
min_genotype_quality = 30
min_base_coverage = 5
min_allele_freq_hom = 0.9
min_allele_freq_het = 0.9
additional_vcf_files = ''
reference_gff_annotations = 'NA'
reference_gff_exclude = 'NA'
snp_eff_results = 'NA'



pigz -f -d -p ${task.cpus} *.vcf.gz
vcf2genome -draft ${out}.fasta -draftname "${fasta_head}" -in ${vcf.baseName} -minc ${vcf2genome_minc} -minfreq ${vcf2genome_minfreq} -minq ${vcf2genome_minq} -ref ${fasta} -refMod ${out}_refmod.fasta -uncertain ${out}_uncertainy.fasta
pigz -p ${task.cpus} *.fasta
pigz -p ${task.cpus} *.vcf




# Step 13: SNP Table Generation

gunzip -f *.vcf.gz
multivcfanalyzer ${params.snp_eff_results} ${fasta} ${params.reference_gff_annotations} . ${write_freqs} ${params.min_genotype_quality} ${params.min_base_coverage} ${params.min_allele_freq_hom} ${params.min_allele_freq_het} ${params.reference_gff_exclude} *.vcf
pigz -p ${task.cpus} *.tsv *.txt snpAlignment.fasta snpAlignmentIncludingRefGenome.fasta fullAlignment.fasta
rm *.vcf


# Step 14 Mitochondrial to Nuclear Ratio

mtnucratio ${bam} "${params.mtnucratio_header}"



#  * Step 17-B: Metagenomic screening of unmapped reads: Kraken2


kraken2 --db ${krakendb} --threads ${task.cpus} --output $out --report $kreport $fastq

kraken_parse.py -c ${params.metagenomic_min_support_reads} -o $out $kraken_r

merge_kraken_res.py -o $out
