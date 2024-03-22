# Analysis

## Table of contents
1. [Population Structure](#pca)
2. [Admixture statistics](#admix)
3. [Sample relatedness](#kinship)
4. [Smoove SV analyses](#smoove)
5. [Variant within Smp_246790](#TRP)

## 01 - Population structure <a name="pca"></a>
### 
```
# Convert vcf to plink .bed format, select autosomes. 
plink2 --vcf FREEZE.FULLFILTER.vcf --chr SM_V9_1, SM_V9_2, SM_V9_3, SM_V9_4, SM_V9_5, SM_V9_6, SM_V9_7 --make-bed --allow-extra-chr --set-all-var-ids @_# --out autosomes_unfiltered

# Remove variants in strong linkage disequilibrium
plink2 --bfile autosomes_unfiltered --allow-extra-chr --set-all-var-ids @_# --indep-pairwise 50 10 0.2
plink2 --bfile autosomes_unfiltered --allow-extra-chr --set-all-var-ids @_# --extract autosomes_unfiltered.prune.in --out prunedData --make-bed
```
### Principal component analysis
```
plink2 --bfile prunedData --allow-extra-chr --set-all-var-ids @_# --pca
```
### Maximum likelihood phylogeny
```
# Subset VCF to remove sites in LD 
cat autosomes_unfiltered.prune.in | sed 's/_/ /g' > keep.sites
bcftools view -T plink2.prune.in --types snps -o FREEZE.FULLFILTER.pruned.vcf FREEZE.FULLFILTER.vcf

# Convert to phylip format
vcf2phylip.py -i FREEZE.FULLFILTER.pruned.vcf
ascbias.py -p FREEZE.FULLFILTER.pruned.min4.phy

# Run iqtree
iqtree -s out.phy -m MFP+ASC --safe -T 30 -B 1000 --seqtype DNA
```
### Admixture
```
#Fix scaffold names in bim file (ADMIXTURE accepts numerical scaffold names only)
sed -i 's/SM_V9_//g' prunedData.bim

# Produce random list of seeds
shuf -i 0-10000 | head -10 > seed.list

# Run ADMIXTURE of values of K:1-20, using 10 randomly generated seeds.
# There is no way of renaming admixture output based on seed value, so to avoid overwriting output files for each seed replicate, run in their own directory or batch run each seed one at a time. 
parallel --dry-run "admixture -j2 --seed={1} -B1000 prunedData.bed {2} --cv=10" :::: seed.list ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 > run_admix.sh

# For example:
mkdir 6127
cd 6127
admixture -j2 --seed=6127 -B1000 ../prunedData.bed 1 --cv=10
admixture -j2 --seed=6127 -B1000 ../prunedData.bed 2 --cv=10
# ... up to 20

# Add K values to ADMIXTURE output files
parallel --dry-run "sed -e 's/^/{1} /g' autosomes.{1}.Q" ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
cat *.Q > admixture_all.txt

# Get a table of CV scores, found in the stdout files (in our case *.o files)
cat *.o | grep CV | cut -f2 -d "=" | sed 's/)://g' | tr ' ' '\t' > cv_scores.txt
```
### Population diversity and differentiation
```
# Subset the allsites VCF for each chromosome (e.g Chr 1)
bcftools view -t 1 -o allsites.filt1.1.vcf allsites.filt1.vcf

# Run PIXY (pop.list = list of populations to be evaluated/compared, repeat for each CHR)
pixy --stats fst dxy pi 
--vcf allsites.filt1.1.vcf
--zarr_path zarr/ 
--window_size 5000
--populations pop.list
--variant_filter_expression 'DP>=10,GQ>=20,RGQ>=20' 
--invariant_filter_expression 'DP>=10,RGQ>=20' 
--outfile_prefix output/pixy.5000.pop.CHR
```
## 02 - Sample relatedness <a name="kinship"></a>
### ngsRelate
```
# Subset to only unlinked variants (method shown above)
bcftools view -T unlinked.sites -t SM_V9_1,SM_V9_2,SM_V9_3,SM_V9_4,SM_V9_5,SM_V9_6,SM_V9_7 -o FREEZE.FULLFILTER.pruned.autosomes.vcf FREEZE.FULLFILTER.vcf

# Run ngsRelate
ngsRelate -h FREEZE.FULLFILTER.pruned.autosomes.vcf -O all_samples -p 10 

# Repeat using subsets for each population, manually compare and use the subset results. 
```
### Sequoia
```
# Subset to only unlinked variants (method shown above)
bcftools view -T unlinked.sites -t SM_V9_1,SM_V9_2,SM_V9_3,SM_V9_4,SM_V9_5,SM_V9_6,SM_V9_7 -o FREEZE.FULLFILTER.pruned.autosomes.vcf FREEZE.FULLFILTER.vcf

# Convert to sequoia friendly input format
plink2 --vcf FREEZE.FULLFILTER.pruned.autosomes.vcf --make-bed --allow-extra-chr --set-all-var-ids @_# --out autosomes_pruned
plink --bfile autosomes_pruned --allow-extra-chr --make-bed --recodeA --out inputfile_for_sequoia

# Move to R and run Sequoia. 
```
## 03 - Smoove SV analyses <a name="smoove"></a>
### Smoove
```
# Filter SV's for deletions, filter further on columns 1&2 to get deletions within specific regions.
bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\t%SVLEN\t%END\t%SHQ\t%MSHQ\t%DHFFC\n]' cohort.smoove.square.anno.filt.vcf |awk '$6!="0/0"' | awk '$5!="./."' | grep -v INV | grep DEL | awk '$10>3' | awk '$11<0.7' | awk '$7>-600000' | awk '$6!="./."' > deletions.pass.txt
```
## 04 - Variation within Smp_246790 <a name="TRP"></a>
### 
```
# Find mutations within Smp_246790
cat FREEZE.FULLFILTER.pruned.autosomes.vcf | vcfEffOnePerLine.pl | java -jar SnpSift.jar extractFields - CHROM POS REF ALT AF "EFF[*].EFFECT" "ANN[*].GENEID" "EFF[*].IMPACT" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].RANK" | grep "Smp_246790." | grep -e HIGH -e MODERATE -e LOW > FI.mutations.txt
```
