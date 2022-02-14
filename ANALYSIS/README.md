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
## 02 - Admixture statistics <a name="admix"></a>
### Identify fixed S. rodhaini alleles
```
# Get a list of non-admixed S. mansoni samples (showing zero/near-zero S. rodhaini admixture) and create lists for each population S. japonicum, S. rodhaini etc.

# Calculate allele frequencies for 3 population: non-admixed Ugandan samples, non-admixed non-Ugandan samples, S. rodhaini samples). 
plink2 --vcf FREEZE.FULLFILTER.vcf --chr SM_V9_1, SM_V9_2, SM_V9_3, SM_V9_4, SM_V9_5, SM_V9_6, SM_V9_7, SM_V9_PAR1, SM_V9_PAR2 --make-bed --allow-extra-chr --set-all-var-ids @_# --out CHR_all
plink2 --bfile CHR_all --keep not_admixed.list --freq --out na_1.test --set-all-var-ids @_#
plink2 --bfile CHR_all --keep not_admixed_outgroup.list --freq --out na_2.test --set-all-var-ids @_#
plink2 --bfile CHR_all --keep rodhaini.list --freq --out na_3.test --set-all-var-ids @_#

# Identify sites fixed in S. rodhaini, and at low/zero frequency in non-admixed populations
less na_1.test.afreq | awk '$5<0.1' | cut -f2 > X1.list
less na_2.test.afreq | awk '$5==0' | cut -f2 > X2.list
less na_3.test.afreq | awk '$5==1' | cut -f2 > X3.list
cat X1.list X2.list X3.list | sort | uniq -c | awk '$1==3' | tr -s ' ' | sed 's/^ //g' | tr ' ' '\t' | cut -f2 | tr '_' '\t' > pSR.alleles.txt
```
### *f*3 
```
# Subset VCF to only include fixed S. rodhaini alleles
bcftools view -T pSR.alleles.txt -o merged.sub.vcf FREEZE.FULLFILTER.vcf

# Make zarr file
python make_zarr.py

# Run calculate f3
python f3.py
```
### Patterson's D
```
# Subset VCF file to include only autosomes
bcftools view -t SM_V9_1,SM_V9_2,SM_V9_3,SM_V9_4,SM_V9_5,SM_V9_6,SM_V9_7 -o FREEZE.FULLFILTER.autosomes.vcf FREEZE.FULLFILTER.vcf

# Calculate Pattersons D for combinations of trios, (for each population, admixed/non-admixed, unrelated samples)
Dsuite Dtrios FREEZE.FULLFILTER.autosomes.vcf ds.list

# Calculate D and f4 statistics for specified trios
Dsuite Dinvestigate -w 10,1 FREEZE.FULLFILTER.autosomes.vcf input.list trios3.list
Dsuite Dinvestigate -w 5,1 FREEZE.FULLFILTER.autosomes.vcf input.list trios3.list
Dsuite Dinvestigate -w 50,25 FREEZE.FULLFILTER.autosomes.vcf input.list trios3.list
```
### ancestry_hmm
```
# Create a VCF with no missing sites, keeping only biallelic sutes
bcftools view -m2 -M2  -i 'F_MISSING=0' -o ALL.nomiss.vcf FREEZE.FULLFILTER.autosomes.vcf

# Create a list of sites not in linkage disequilibrium (see plink method above)

# Subset VCF to only include unlinked sites, repeat for each population
bcftools view -T unlinked.sites -s pop.list -o ALL.nomiss.pruned.POP1.vcf ALL.nomiss.vcf

# Convert to ahmm format (*.master files describe target and reference populations, unrelated samples), repeat for each subpopulation. 
python vcf2ahmm.py -v ALL.nomiss.pruned.POP1.vcf -s POP1.master -g 1 --min_total 1 -r 4.0953e-8 -o POP1 > POP1.input

# Run a_hmm (repeat for each parameter and subpopulation)
ancestry_hmm -i POP1.input -s POP1 -a 2 0.01 0.99 -p 0 -30 0.10 -p 1 200000 1 -b 1000 1000 --tmin 0 --tmax 100000 --ne 65000 > POP1.log
```
### Sprime
```
# Subset for each chromosome
bcftools view -t SM_V9_1 -o ALL.nomiss.1.vcf ALL.nomiss.vcf

# Phase variants use genetic map generated using uniform per-chromosome recombination rate
java -jar beagle.28Sep18.793.jar gt=ALL.nomiss.1.vcf out=ALL.nomiss.1.beagle map=1_fixed.gmap nthreads=4 iterations=30 burnin=10 ne=65000

# Run Sprime (outgroup.txt=list of non-admixed samples, excl.pop=list of samples not to be analysed, unrelated samples)
java -jar sprime.jar gt=ALL.nomiss.1.beagle.vcf outgroup=outgroup.txt map=plink.gmap out=POP1_1_sprime excludesamples=excl.pop chrom=1 mu=8.1e-9

# Create blocks for each group (Chr 1 shown as example)
cat POP1_1_sprime.score | cut -f1,2,6,7,8 | sort -k3,3 -k1,1 -k2,2nr | awk -F'\t' '!seen[$3]++' > 1.a.test
cat POP1_1_sprime.score | cut -f1,2,6,7,8 | sort -k3,3 -k1,1 -k2,2n | awk -F'\t' '!seen[$3]++' >> 1.a.test
cat 1.a.test | sort -k3,3 | awk '{print $3,$1,$2,$4,$5}' | tr ' ' '\t' | awk '{if(a!=$1) {a=$1; printf "\n%s%s",$0,FS} else {a=$1;$1="";printf $0 }} END {printf "\n" }' > 1.a.b.temp
cat 1.a.b.temp | tr -s ' ' | tr -s '\t' | tr ' ' '\t' | sed 's/$/ gn/g' > POP1.summary

# Calculate genetic distance (Chr 1 shown as example)
less POP1.summary | cut -f2,3,7,8,9,10 | awk '{print $1,$2,$3,$3-$2,$4,$5,$6}' | tr ' ' '\t' | sort -k4,4 -gr | sed 's/-//g' | grerp "^1" | awk '{print $1,(($4*4.095)/1000000),$2,$3,$4,$5,$6,$7}' | sed 's/ / /g' > chr1.segs.txt

# Repeat for all chromosomes, merge into one file
cat chr1.segs.txt chr2.segs.txt ... > all.segs.txt
```
### Proportions of S. rodhaini alleles per sample
```
# Using AD0158 as an example
# Subsample phased VCF to contain only putative S. rodhaini specific alleles. Using merged beagle variants file. E.g:
bcftools query -f [%CHROM\\t%POS\\t%SAMPLE\\t%GT\\n] -o AD0158.query.derived -s AD0158 -T pSR.alleles.txt ALL.nomiss.beagle.vcf
cat AD0158.query.derived | awk '{print \$1,\$2-1,\$2,\$3,\$4}' | sed 's/ /  /g' > AD0158.proc

# Make 2kb windows
bedtools makewindows -g SM_V9.fa.fai -w 2000 | grep -v MITO | grep -v ZSR | grep -v WSR > 2kb_windows.txt

# Intersect variant counts with windows
bedtools intersect -wo -a 2kb_windows.txt -b AD0158.proc |cut -f1,2,3,7,8 | sed 's/$/ 1/g' | sort -k1,1 -k2,2n -k4,4 | datamash -g1,2,3,4,5 count 6 > AD0158.proc.merge

# Get counts of variant sites per window
less AD0158.proc.merge | cut -f1,2,3,6 | datamash -g1,2,3 sum 4 > counts.txt

# Merge counts to get proportions of each genotype per window
join <(sort counts.txt) <(cat AD0158.proc.merge | sed 's/ /_/1' | sed 's/ /_/1' | sort ) | sed 's/_/  /g' > AD0158.final.counts
```
### Allele frequencies of S. rodhaini alleles in admixed populations
```
# Get allele freqencies across admixed populations (unrelated samples)
plink2 --bfile CHR_all --keep all.admixed.samples.list --freq --out all.admixed.samples --set-all-var-ids @_# --extract <(cat pSR.alleles.txt| tr '\t' '_')
```
### Calculate Tajima's D
```
# Calculate Tajima's D for each subpopualation (unrelated samples only)
vcftools --vcf FREEZE.FULLFILTER.vcf --keep POP1.list --TajimaD 5000 --out {2}_TAJIMA_D_POP1_5000
vcftools --vcf FREEZE.FULLFILTER.vcf --keep POP.list --TajimaD 2500 --out {2}_TAJIMA_D_POP1_2500
```
### Neighbour-joining phylogenies (various plots)
```
# Subset and prune vcf file (as above) then calculate distance matrix
plink2 --bfile prunedData --allow-extra-chr --set-all-var-ids @_# --distance square 1-ibs
paste <( cut -f2 prunedData_tree.mdist.id) prunedData_tree.mdist | cat <(cut -f2 prunedData_tree.mdist.id | tr '\n' '\t' | sed -e '1s/^/\t/g') - > autosomes.mdist 
```
## 03 - Sample relatedness <a name="kinship"></a>
### ngsRelate
```
# Subset to only unlinked variants (method shown above)
bcftools view -T unlinked.sites -t SM_V9_1,SM_V8_2,SM_V8_3,SM_V9_4,SM_V9_5,SM_V9_6,SM_V9_7 -o FREEZE.FULLFILTER.pruned.autosomes.vcf FREEZE.FULLFILTER.vcf

# Run ngsRelate
ngsRelate -h FREEZE.FULLFILTER.pruned.autosomes.vcf -O all_samples -p 10 

# Repeat using subsets for each population, manually compare and use the subset results. 
```
### Sequoia
```
# Subset to only unlinked variants (method shown above)
bcftools view -T unlinked.sites -t SM_V9_1,SM_V8_2,SM_V8_3,SM_V9_4,SM_V9_5,SM_V9_6,SM_V9_7 -o FREEZE.FULLFILTER.pruned.autosomes.vcf FREEZE.FULLFILTER.vcf

# Convert to sequoia friendly input format
plink2 --vcf FREEZE.FULLFILTER.pruned.autosomes.vcf --make-bed --allow-extra-chr --set-all-var-ids @_# --out autosomes_pruned
plink --bfile autosomes_pruned --allow-extra-chr --make-bed --recodeA --out inputfile_for_sequoia

# Move to R and run Sequoia. 
```
## 03 - Smoove SV analyses <a name="smoove"></a>
### Smoove
```
```
## 03 - Variantion within Smp_246790 <a name="TRP"></a>
### 
```
```
