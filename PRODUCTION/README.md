# Mapping, variant calling and quality control

## Table of contents
1. [Raw data](#raw)
2. [Mapping](#mapping)
3. [Variant calling](#variantcalling)
4. [Quality control](#qc)
5. [SnpEff annotation](#snpeff)
6. [Smoove SV calling](#smoove)

## 01 - Raw data <a name="raw"></a>
### Reference genome
```
# Create indexes and a sequence dictionary for the reference genome
bwa index SM_V9.fa
gatk CreateSequenceDictionary --REFERENCE SM_V9.fa
```
## 02 - Mapping <a name="mapping"></a>
### Map sequence reads to reference genome
```
# Trim reads (adapters.fa is list of standard illumina adapters)
bbduk.sh -Xmx30g -in=SAMPLE1_1.fastq.gz -in2=SAMPLE1_2.fastq -out=SAMPLE1_1.fastq.trimmed.gz  -out2=SAMPLE1_2.fastq.trimmed.gz ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo

# Repeat as needed for each read set, for example:
bwa mem -M -t 6 SM_V9.fa SAMPLE1_1.fastq.trimmed.gz SAMPLE1_1.fastq.trimmed.gz | samtools sort -@6 -o SAMPLE1.bam -
```
### Mark PCR duplicates
```
# Mark duplicates
gatk MarkDuplicates --INPUT SAMPLE1.bam --OUTPUT SAMPLE1.markdup.bam --METRICS_FILE SAMPLE1.metrics.txt

# Merge BAM files (from samples where there are multiple sets of FASTQs)
samtools merge -@ 6 SAMPLE1.markdup.merged.bam SAMPLE1.markdup.bam SAMPLE1b.markdup.bam

# Index all BAMs
samtools index SAMPLE1.markdup.merged.bam
```
### Calculate coverage
```
# Create makewindows input
cut -f1,2 SM_V9.fa.fai > SM_V9.chrom.txt

# Create 500 bp windows
bedtools makewindows -g SM_V9.chrom.txt -w 500 > SM_V9.chrom.500bp.bed

# Calculate per-sample coverage
bedtools coverage -sorted -g SM_V9.fa.fai -d -a SM_V9.chrom.500bp.bed -b SAMPLE1.markdup.merged.bam \| datamash -g1,2,3 median 5 mean 5 sstdev 5 > SAMPLE1.cov
awk '{print $1,$2,$3,$4,$5,$6,FILENAME}' SAMPLE1.cov | sed 's/.cov//g' > SAMPLE1.recov
cat *.recov > all.recov.txt
```
## 03 - Variant calling <a name="variantcalling"></a>

### Per-sampling variant calling
```
samtools index SAMPLE1.markdup.merged.bam
gatk HaplotypeCaller --emit-ref-confidence GVCF -I SAMPLE1.markdup.merged.bam -R SM_V9.fa -O SAMPLE1.g.vcf
```
### Combine all samples into a single gVCF and genotype
```
# Combine gVCFs
ls | grep 'g.vcf' > argument.list
gatk CombineGVCFs --arguments_file argument.list --reference SM_V9.fa --output merged_all_samples.g.vcf

# Genotype
gatk GenotypeGVCFs --reference SM_V9.fa --variant merged_all_samples.g.vcf --output merged_all_samples.vcf
```
## 04 - Quality control <a name="qc"></a>
### Calculate quality scores for all variant sites
```
# Produce a table of quality scores for each variant site
gatk VariantsToTable --variant merged_all_samples.vcf -F CHROM -F POS -F TYPE -F QD -F FS -F MQ -F MQRankSum -F ReadPosRankSum -F SOR -F InbreedingCoeff -R SM_V9.fa --output cohort.genotyped.tbl
```
### Separate and filter SNPs
```
# Select SNPs
gatk SelectVariants -R SM_V9.fa --variant merged_all_samples.vcf --select-type-to-include SNP --output merged_all_samples.SNPs.vcf

# Tag low-quality SNPs
gatk VariantFiltration \
--filter-expression "ReadPosRankSum < -8.0" --filter-name "RPRS8" \
--filter-expression "QD < 2.0" --filter-name "QD2" \
--filter-expression "FS > 60.0" --filter-name "FS60" \
--filter-expression "MQ < 40.0" --filter-name "MQ40" \
--filter-expression "MQRankSum < -12.5" --filter-name "MQ12.5" \
--filter-expression "SOR > 3.0" --filter-name "SOR3" \
--variant merged_all_samples.SNPs.vcf \
-R SM_V9.fa  \
--output merged_all_samples.SNPs.tagged.vcf

# Remove low-quality sites
gatk SelectVariants -R SM_V9.fa --variant merged_all_samples.SNPs.tagged.vcf --exclude-filtered --output merged_all_samples.SNPs.filtered.vcf
```
### Separate and filter indels and mixed sites
```
# Select indels and mixed sites
gatk SelectVariants -R SM_V9.fa --variant merged_all_samples.vcf --select-type-to-include SNP --output merged_all_samples.indels_mixed.vcf

# Tag low-quality indels and mixed sites
/lustre/scratch118/infgen/team133/db22/software/gatk-4.1.0.0/gatk VariantFiltration \
--filter-expression "QD < 2.0" --filter-name "QD2" \
--filter-expression "FS > 200.0" --filter-name "FS200" \
--filter-expression "ReadPosRankSum < -20.0" --filter-name "RPRS20" \
--filter-expression "SOR > 10.0" --filter-name "SOR10" \
--variant merged_all_samples.indels_mixed.vcf \
-R SM_V9.fa  \
--output merged_all_samples.indels_mixed.tagged.vcf

# Remove low-quality sites
gatk SelectVariants -R SM_V9.fa --variant merged_all_samples.indels_mixed.tagged.vcf --exclude-filtered --output merged_all_samples.indels_mixed.filtered.vcf
```
### Recombine filtered variants
```
gatk MergeVcfs --INPUT merged_all_samples.SNPs.filtered.vcf --INPUT merged_all_samples.indels_mixed.filtered.vcf --OUTPUT merged_all_samples.filtered.vcf
```
### Remove low-quality samples and variants 
```
# Calculate per-individual missingness rate 
vcftools --vcf merged_all_samples.filtered.vcf.FL1.vcf --missing-indv --out missing_indv

# Filter out individuals with high rates of missing variant calls
awk '$6<=0.10' missing_indv.imiss | grep -v "MISS" | cut -f1  > retain.samples.list
vcftools --vcf merged_all_samples.filtered.vcf.FL1.vcf --keep retain.samples.list --recode-INFO-all --recode --out merged_all_samples.filtered.vcf.FL1.vcf  

# Calculate per-site missingness rate 
vcftools --vcf merged_all_samples.filtered.vcf --missing-site --out missing_site

# Filter out sites with high rates of missing variant calls
awk '$6<=0.05' missing_site.lmiss | grep -v "MISS" | cut -f1,2 > retain.variants.list
vcftools --vcf merged_all_samples.filtered.vcf.FL1.vcf --postions retain.variants.list --recode-INFO-all --recode --out merged_all_samples.filtered.vcf.FL2.vcf 
```
### Move final versions of VCFs to the analysis folder
```
mv merged_all_samples.filtered.vcf.FL2.vcf FREEZE.FULLFILTER.vcf
```
### Accessory VCF 1 - Input for PIXY
```
# Produce an allsites VCF
gatk GenotypeGVCFs --reference SM_V9.fa --variant merged_all_samples.g.vcf --output merged_all_samples.IV.vcf --include-non-variant-sites

# Hard filter
java -Xmx25G -jar gatk-package-4.2.0.0-local.jar VariantFiltration --filter-expression "ReadPosRankSum < -8.0" --filter-name RPRS8 --filter-expression "QD < 2.0" --filter-name QD2 --filter-expression "FS > 60.0" --filter-name FS60 --filter-expression "MQ < 40.0" --filter-name MQ40 --filter-expression "MQRankSum < -12.5" --filter-name MQ12.5 --filter-expression "SOR > 3.0" --filter-name SOR3 --variant merged_all_samples.IV.vcf -R SM_V9.fa --output allsites.tagged.vcf

# Second hard filter (pass.list = list of accessions retained during production of non-accessory VCF)
vcftools --vcf allsites.tagged.vcf --recode --remove-filtered-all --out allsites.filt1.vcf --remove-indels --max-missing 0.8 --min-meanDP 5 --max-meanDP 500 --hwe 0.001 --mac 1 --keep pass.list
```
### Accessory VCF 2 - Inclusion of S. japonicum isolates
```
# Trim reads, map, variant call and merge calls of S. japonicum accessions (as above) independently of all other accessions
# Calculate per-site and per-sample missingness
vcftools --vcf all.SJ.vcf --missing-indv
vcftools --vcf all.SJ.vcf --missing-site

# Subset only to variant sites retained in the primary VCF
bcftools query  -f '%CHROM\t%POS\n' FREEZE.FULLFILTER.vcf > keep.sites
less out.lmiss | awk '$6<1' | awk '$6>0.75' | cut -f1,2 > highmiss.list

# Compress and index
bgzip -@ 12 -c FREEZE.FULLFILTER.vcf > mans.vcf.gz
tabix mans.vcf.gz

bgzip -@ 12 -c all.SJ.F2.vcf > SJ.vcf.gz
tabix SJ.vcf.gz

# Merge
bcftools merge --threads 6 -o merged.vcf mans.vcf.gz SJ.vcf.gz
bcftools view --threads 4 -T keep.sites -o all.SJ.F1.vcf merged.vcf
bcftools view --threads 4 -t 1,2,3,4,5,6,7 -T^highmiss.list -o all.SJ.F2.vcf all.SJ.F1.vcf
```
## 04 - SnpEff annotation <a name="snpeff"></a>
### 
```
# Normalize variants 
bcftools norm -f SM_V9.fa --threads 12 -m - -o ALL.normed.vcf FREEZE.FULLFILTER.vcf

# Build database (assuming you have files in the right directories)
java -jar snpEff.jar build -c snpEff.config -gtf22 SM_V9

# Annotate
java -jar snpEff.jar SM_V9 ALL.normed.vcf > all.snpeff.vcf
```
## 04 - Smoove SV calling <a name="smoove"></a>
### 
```
# Call structural variants per-sample (e.g accession MK0037)
smoove call --outdir results-smoove/ --name MK0037 --fasta SM_V9.fa -p 1 --genotype MK0037.renamed.bam

# Get the union of sites across all samples 
smoove merge --name merged -f SM_V9.fa --outdir OUT results-smoove/*.genotyped.vcf.gz

# Genotype those sites
smoove genotype -d -x -p 1 --name sample-joint --outdir results-genotyped/ --fasta SM_V9.fa --vcf merged.sites.vcf.gz MK0037.renamed.bam

# Merge all single-sample VCFs
smoove paste --name cohort results-genotyped/*.vcf.gz

# Annotate SVs
smoove annotate --gff SM_V9.gff cohort.smoove.square.vcf.gz | bgzip -c > cohort.smoove.square.anno.vcf.gz

# Filter Smoove SV calls
slivar expr \
  --info "variant.call_rate > 0.25 && ((INFO.SVTYPE == 'DEL') ||
(INFO.SVTYPE == 'DUP'))" \
  --sample-expr \
    "LQHET:sample.alts == 1 && (((sample.DHFFC > 0.75) && (INFO.SVTYPE
== 'DEL')) || ((sample.DHFFC < 1.25) && INFO.SVTYPE == 'DUP'))" \
  --sample-expr \
    "HQHET:sample.alts == 1 && (((sample.DHFFC < 0.75) && (INFO.SVTYPE
== 'DEL')) || ((sample.DHFFC > 1.25) && INFO.SVTYPE == 'DUP'))" \
  --sample-expr \
    "HQHA:sample.alts == 2 && (((sample.DHFFC < 0.5) && (INFO.SVTYPE
== 'DEL')) || ((sample.DHFFC > 1.5) && INFO.SVTYPE == 'DUP'))" \
  -o cohort.smoove.square.anno.filt.vcf \
  --vcf cohort.smoove.square.anno.vcf.gz
```
 
