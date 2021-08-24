# Rattlesnake venom population genomics

![Rattlesnake venom population genomics](/cover-image.png "cover image")

This repository contains details on the data processing and analysis steps taken in analyses of population genomic variation and signatures of selection in major rattlesnake venom gene regions. This workflow is a companion to the methods description in Schield et al. (in review). Analysis of recombination rates and statistics requiring phased variants are based on [recombination maps](https://figshare.com/articles/dataset/Rattlesnake_Recombination_Maps/11283224) and phased data from [Schield et al. _MBE_ 37: 1272-1294](https://academic.oup.com/mbe/advance-article-abstract/doi/10.1093/molbev/msaa003/5700722).

Lists and reference files (i.e., BED, GFF, etc.) are in the `resources` directory. Shell and Python scripts are in respective `shell` and `python` directories. R scripts are in the `R` directory. Note that you may need to adjust the organization of file locations to suit your environment. This workflow assumes that you are returning to the main working directory after each major section (e.g., mapping, variant calling).

If you have any questions, you can email me at drew.schield[at]colorado.edu.

## Contents

* [Software and dependencies](#software-and-dependencies)
* [General resources](#general-resources)
* [Read filtering](#read-filtering)
* [Read mapping](#read-mapping)
* [Variant calling](#variant-calling)
* [Variant filtering](#variant-filtering)
* [Analysis of copy-number variation](#analysis-of-copy-number-variation)
* [Population structure analysis](#population-structure-analysis)
* [Demographic analysis](#demographic-analysis)
* [Population genetic diversity and differentiation](#population-genetic-diversity-and-differentiation)
	* [Diversity appendix 1: CNV-masking](#diversity-appendix-1-cnv-masking)
	* [Diversity appendix 2: Comparison between venom genic and intergenic regions](#diversity-appendix-1-comparison-between-venom-genic-and-intergenic-regions)
* [Signatures of selection](#signatures-of-selection)
	* [1. Tajima's D](#1-tajimas-d)
	* [2. Fixed differences](#2-fixed-differences)
	* [3. iHS](#3-ihs)
	* [4. ß](#4-ß)
	* [Selection appendix 1: CNV-masking](#selection-appendix-1-cnv-masking)
	* [Selection appendix 2: Venom gene point estimates](#selection-appendix-2-venom-gene-point-estimates)
	* [Selection appendix 3: Estimates for non-venom homologs](#selection-appendix-3-estimates-for-non-venom-homologs)
	* [Selection appendix 4: Point estimates for all genes](#selection-appendix-4-point-estimates-for-all-genes)
	* [Selection appendix 5: Point estimates for 'other' venom genes](#selection-appendix-5-point-estimates-for-other-venom-genes)
	* [Selection appendix 6: Proportion heterozygotes at high ß sites](#selection-appendix-6-proportion-heterozygotes-at-high-ß-sites)
	* [Selection appendix 7: Trans-species polymorphism](#selection-appendix-7-trans-species-polymorphism)
* Recombination rate variation and linkage disequilibrium analysis
* Analysis in R
* [Appendix 1: Mapping statistics](#appendix-1-mapping-statistics)
* [Appendix 2: Genotype quality](#appendix-2-genotype-quality)

## Software and dependencies

The steps described below use the following software and assume that dependencies are on the user path:

* [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [bwa](http://bio-bwa.sourceforge.net/)
* [GATK](https://gatk.broadinstitute.org/hc/en-us) (v3.8-1-0 and v4.0.8.1)
* [htslib](http://www.htslib.org/)
* [samtools](http://www.htslib.org/)
* [bcftools](http://www.htslib.org/)
* [bgzip](http://www.htslib.org/)
* [tabix](http://www.htslib.org/)
* [plink](https://www.cog-genomics.org/plink/)
* [vcftools](https://vcftools.github.io/index.html)
* [bedtools](https://bedtools.readthedocs.io/en/latest/)
* [CNV-seq](https://github.com/hliang/cnv-seq)
* [ADMIXTURE](https://dalexander.github.io/admixture/publications.html)
* [PSMC](https://github.com/lh3/psmc)
* [pixy](https://pixy.readthedocs.io/en/latest/)
* [Glactools](https://github.com/grenaud/glactools)
* [betascan](https://github.com/ksiewert/BetaScan)
* [rehh](https://cran.r-project.org/web/packages/rehh/vignettes/rehh.html) (R package)
* [R](https://cran.r-project.org/)

Note: I installed a number of these dependencies using [conda](https://docs.conda.io/en/latest/).

## General resources

### Processing files

Several files will come up repeatedly throughout this workflow, namely annotation files with coordinates for venom genes and non-venom paralogs across the genome. These files are each located in the `resources` directory.

* `VenomGene_regions_full_updated.gff` - contains venom gene start/stop coordinates.
* `VenomGene_main3vgFams_FULL_annotation.gff` - contains full annotations (UTR, exon, etc.) for SVMP, SVSP, and PLA2 genes.
* `main3vgFams_paralogs_01.12.21.gff` - contains coordinates of non-venom paralogs for comparative analyses.
* `region_SVMP_scaffold-mi1.bed` - contains coordinates of SVMP region of chromosome 9.
* `region_SVSP_scaffold-mi2.bed` - contains coordinates of SVSP region of chromosome 10.
* `region_PLA2_scaffold-mi7.bed` - contains coordinates of PLA2 region of chromosome 15.
* `centromere_mask.bed` - contains approximate coordinates of centromeres on macrochromosomes.

To generate BED files for each of the three main venom families:

```
$grep 'SVMP' VenomGene_regions_full_updated.gff | awk 'BEGIN{OFS="\t"}{print $1, $4-1, $5}' | bedtools sort -i - > gene_SVMP.bed
$grep 'SVSP' VenomGene_regions_full_updated.gff | awk 'BEGIN{OFS="\t"}{print $1, $4-1, $5}' | bedtools sort -i - > gene_SVSP.bed
$grep 'PLA2' VenomGene_regions_full_updated.gff | awk 'BEGIN{OFS="\t"}{print $1, $4-1, $5}' | bedtools sort -i - > gene_PLA2.bed
```

To generate BED files for non-venom paralogs for each of the three main venom families:

```
$grep 'SVMP' main3vgFams_paralogs_01.12.21.gff | grep -v 'scaffold-un' | awk 'BEGIN{OFS="\t"}{print $1, $4-1, $5}' | bedtools sort -i - > non-venom_paralogs_SVMP.bed
$grep 'SVSP' main3vgFams_paralogs_01.12.21.gff | grep -v 'scaffold-un' | awk 'BEGIN{OFS="\t"}{print $1, $4-1, $5}' | bedtools sort -i - > non-venom_paralogs_SVSP.bed
$grep 'PLA2' main3vgFams_paralogs_01.12.21.gff | grep -v 'scaffold-un' | awk 'BEGIN{OFS="\t"}{print $1, $4-1, $5}' | bedtools sort -i - > non-venom_paralogs_PLA2.bed
```

### Reference genome

Use the [Prairie Rattlesnake (_Crotalus viridis_) genome assembly](https://figshare.com/articles/dataset/Prairie_rattlesnake_Crotalus_viridis_genome_assembly/9030782) as the reference in this workflow.

Prepare necessary indexes for downstream analysis.

```
bwa index CroVir_genome_L77pg_16Aug2017.final_rename.fasta
samtools faidx CroVir_genome_L77pg_16Aug2017.final_rename.fasta
./gatk-4.0.8.1/gatk CreateSequenceDictionary -R CroVir_genome_L77pg_16Aug2017.final_rename.fasta
```

## Read filtering

Quality trim and filter raw whole genome resequencing reads using trimmomatic using these settings:

* Remove 5' end bases if quality is below 20
* Remove 3' end bases if quality is below 20
* Minimum read length = 32
* Remove reads if average quality is < 30

### Set up environment

Get raw fastq data into `fastq` directory. <br /> Make a `fastq_filtered` directory for output.

```
mkdir fastq
mkdir fastq_filtered
```

### Run trimmomatic on raw reads

The script below will run trimmomatic on the raw data for each sample in `resources/sample.list`.

__*trimmomatic.sh*__

```
list=$1
for line in `cat $list`; do
	name=$line
	echo processing and filtering ${name}.
	trimmomatic PE -phred33 -threads 16 ./fastq/${name}_R1_001.fastq.gz ./fastq/${name}_R2_001.fastq.gz ./fastq_filtered/${name}_R1_P.trim.fastq.gz ./fastq_filtered/${name}_R1_U.trim.fastq.gz ./fastq_filtered/${name}_R2_P.trim.fastq.gz ./fastq_filtered/${name}_R2_U.trim.fastq.gz LEADING:20 TRAILING:20 MINLEN:32 AVGQUAL:30
done
```

Run the script.

`sh trimmomatic.sh sample.list`

## Read mapping

Use bwa 'mem' to map our filtered reads to the reference genome.

### Set up environment

`mkdir bam`

### Map reads with bwa and sort with samtools

The script below will run bwa mem on the paired, filtered reads per sample and sort the output bam file.

__*bwa_mem.sh*__

```
list=$1
for line in `cat $list`; do
	name=$line
	echo "Mapping filtered $name data to reference"
	bwa mem -t 16 -R "@RG\tID:$name\tLB:CVOS\tPL:illumina\tPU:NovaSeq6000\tSM:$name" CroVir_genome_L77pg_16Aug2017.final_rename.fasta ./fastq_filtered/${name}_R1_P.trim.fastq.gz ./fastq_filtered/${name}_R2_P.trim.fastq.gz | samtools sort -@ 16 -O bam -T temp -o ./bam/$name.bam -
done
```

Run the script.

`sh bwa_mem.sh sample.list`

## Variant calling

Use GATK for individual variant discovery and variant calling among the cohort of samples. This is a two-step process, first using HaplotypeCaller to generate individual genomic VCFs (GVCFs), then using GenotypeGVCFs to call variants among samples and generate an all-sites VCF.

### Set up environment

```
mkdir gvcf
mkdir vcf
```

### Call individual variable sites using HaplotypeCaller

The script below will run GATK HaplotypeCaller on each sample in the `sample.list`. It will also zip and index the GVCF output.

*Note: GATK HaplotypeCaller will take a small eternity to run on all of these samples one-by-one. Consider breaking up the job into smaller lists of samples and running jobs in parallel.*

__*GATK_HaplotypeCaller.sh*__

```
list=$1
for i in `cat $list`; do
	./gatk-4.0.8.1/gatk HaplotypeCaller -R CroVir_genome_L77pg_16Aug2017.final_rename.fasta --ERC GVCF -I ./bam/$i.bam -O ./gvcf/$i.raw.snps.indels.g.vcf
	bgzip ./gvcf/$i.raw.snps.indels.g.vcf
	tabix -p vcf ./gvcf/$i.raw.snps.indels.g.vcf.gz
done
```

Run the script.

`sh GATK_HaplotypeCaller.sh sample.list`

### Call cohort variant sites and generate an 'all-sites' VCF using GenotypeGVCFs

Format a file with paths to the GVCF files to call variants from (this is in `resources/sample.gvcf.list`).

```
java -jar ../gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R CroVir_genome_L77pg_16Aug2017.final_rename.fasta -V sample.gvcf.list -allSites -o ./vcf/cvco+outgroup.raw.vcf.gz
```

## Variant filtering

Now apply the following filters to the raw variants:

1. Mask repeats
2. Set low quality/read depth genotypes, indels, and masked repeats as missing genotypes
3. Set sites with extremely high read depth as missing genotypes and restrict output to chromosome-assigned scaffolds

### 1. Mask repeats

This step will set the format field of the VCF to 'REP' for all repeat element intervals annotated in the reference genome. The repeat annotation GFF is [here](https://figshare.com/articles/dataset/Prairie_rattlesnake_repeat_element_annotation/9031481).

First, convert the GFF to a sorted BED file.

```
gunzip CroVir_genome_L77pg_16Aug2017.repeat.masked.final.out.gff3.gz
grep -v # CroVir_genome_L77pg_16Aug2017.repeat.masked.final.out.gff3 | awk 'BEGIN{OFS="\t"}{print $1,$4-1,$5}' | bedtools sort -i - > CroVir_genome_L77pg_16Aug2017.repeat.masked.final.sort.bed
```

Index the new repeat BED file.

```
gatk-4.0.8.1/gatk IndexFeatureFile -F CroVir_genome_L77pg_16Aug2017.repeat.masked.final.sort.bed
```

Use GATK VariantFiltration to mask repeat bases.

```
java -jar ../gatk-3.8-1-0/GenomeAnalysisTK.jar -T VariantFiltration -R CroVir_genome_L77pg_16Aug2017.final_rename.fasta --mask CroVir_genome_L77pg_16Aug2017.repeat.masked.final.sort.bed --maskName REP --setFilteredGtToNocall --variant ./vcf/cvco+outgroup.raw.vcf.gz --out ./vcf/cvco+outgroup.mask.vcf.gz
```

### 2. Set repeats, indels, and low quality genotypes as missing

Use bcftools to filter and tabix to index the output.

```
bcftools filter --threads 24 -e 'FORMAT/DP<5 | FORMAT/GQ<30 || TYPE="indel" || FILTER="REP"' --set-GTs . -O z -o ./vcf/cvco+outgroup.mask.HardFilter.vcf.gz ./vcf/cvco+outgroup.mask.vcf.gz
tabix -p vcf ./vcf/cvco+outgroup.mask.HardFilter.vcf.gz
```

### 3. Set high coverage sites as missing genotypes and remove unassigned scaffolds

A survey of the read depths at variant sites in the output above found this distribution:

* Mean = 23.9
* Median = 23.87
* 2.5% = 4.9
* 97.5% = 36.24

To avoid potential errors from paralogous mappings, set sites with mean depth > 36.24 as missing genotypes. Also use `resources/chrom.bed` to only keep data from chromosome-assigned scaffolds.

```
bcftools view --threads 16 -R chrom.bed ./vcf/cvco+outgroup.mask.HardFilter.vcf.gz | bcftools filter --threads 16 -e 'MEAN(FORMAT/DP)>36.24' --set-GTs . -O z -o ./vcf/cvco+outgroup.mask.HardFilter.depth.chrom.vcf.gz
tabix -p vcf ./vcf/cvco+outgroup.mask.HardFilter.depth.chrom.vcf.gz
```

### Parse chromosome-specific VCFs

The script below will extract a VCF per chromosome using `resources/chrom.list`.

`mkdir ./vcf/vcf_chrom-specific_cvco+outgroup`

__*parse_chrom_vcf.sh*__

```
chromlist=$1
for chrom in `cat $chromlist`; do
	echo parsing $chrom VCF
	bcftools view --threads 16 -r $chrom -O z -o ./vcf/vcf_chrom-specific_cvco+outgroup/cvco+outgroup.mask.HardFilter.depth.chrom.$chrom.vcf.gz ./vcf/cvco+outgroup.mask.HardFilter.depth.chrom.vcf.gz
	tabix -p vcf ./vcf/vcf_chrom-specific_cvco+outgroup/cvco+outgroup.mask.HardFilter.depth.chrom.$chrom.vcf.gz
done
```

Run the script.

`sh parse_chrom_vcf.sh`

## Analysis of copy-number variation

Compare read depths of pairs of lineages mapped to the Prairie rattlesnake reference to identify regions with evidence of copy number variation.

### Set up environment

```
mkdir cnv
cd cnv
```
### Install CNV-seq dependency

For the Perl script above to work, install the R package cnv in the cnv-seq directory after downloading CNV-seq from GitHub.

The steps detailed on Hilary Parker's R blog [here](https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/) are really helpful for making an R package from what was already present in the directory.

Making a package from the source .R script requires the dependencies:
* devtools
* roxygen2

Also need to install a system dependency `libgit2-dev`:

```
sudo apt-get update -y
sudo apt-get install -y libgit2-dev
```

Open R and install dependencies.

```
R
install.packages("devtools")
install.packages("roxygen2")
library(devtools)
library(roxygen2)
setwd("/data3/venom_population_genomics/1_cnv/cnv-seq")
install("cnv")
library(cnv)
```

If the package does not register "Loading required package: ggplot2", something is wrong.

### Perform analysis

CNV-seq relies on 'hits' derived from samtools, so first obtain hits per sample, then analyze pairs of samples using cnv-seq.pl.

```
mkdir hits
mkdir cnv-seq_output
samtools view -F 4 ../bam/CV0009.bam | perl -lane 'print "$F[2]\t$F[3]"' > hits/Cviridis_CV0009.hits
samtools view -F 4 ../bam/CV0087.bam | perl -lane 'print "$F[2]\t$F[3]"' > hits/Coreganus_CV0087.hits
```

### Prune unassigned scaffold hits from results

```
for i in *.hits; do grep -v 'scaffold-un' $i > $i.chrom; done
rm *.hits
```
 
### Perform CNV analysis using cnv-seq.pl

The total length of chromosome-assigned scaffolds is 1,296,980,472 bp.

Use the following settings:

* log2 = 0.6
* p = 0.001
* window-size = 10000 (overwrites log2 and p commands!)
* bigger-window = 1.5
* annotate
* minimum-windows = 4

```
cd hits
../cnv-seq/cnv-seq.pl --test Coreganus_CV0087.hits.chrom --ref Cviridis_CV0009.hits.chrom --genome-size 1296980472 --log2 0.6 --p 0.001 --window-size 10000 --bigger-window 1.5 --annotate --minimum-windows 4
mv *.cnv ../cnv-seq_output/
mv *.count ../cnv-seq_output/
```

### Extract chromosome-specific hits for venom-linked microchromosomes

```
for i in *.chrom; do grep -w 'scaffold-mi1' $i > $i.mi1; done
for i in *.chrom; do grep 'scaffold-mi2' $i > $i.mi2; done
for i in *.chrom; do grep 'scaffold-mi7' $i > $i.mi7; done
```

### Perform chromosome-specific analyses

The lengths of microchromosomes 1, 2, and 7 are:

* scaffold-mi1 (chromosome 9) = 22521304
* scaffold-mi2 (chromosome 10)= 19978503
* scaffold-mi7 (chromosome 15)= 12380205
	
Used the following settings:

* log2 = 0.6
* p = 0.001
* bigger-window = 1.5
* annotate
* minimum-windows = 4

```
../cnv-seq/cnv-seq.pl --test Coreganus_CV0087.hits.chrom.mi1 --ref Cviridis_CV0009.hits.chrom.mi1 --genome-size 22521304 --log2 0.6 --p 0.001 --bigger-window 1.5 --annotate --minimum-windows 4
../cnv-seq/cnv-seq.pl --test Coreganus_CV0087.hits.chrom.mi2 --ref Cviridis_CV0009.hits.chrom.mi2 --genome-size 19978503 --log2 0.6 --p 0.001 --bigger-window 1.5 --annotate --minimum-windows 4
../cnv-seq/cnv-seq.pl --test Coreganus_CV0087.hits.chrom.mi7 --ref Cviridis_CV0009.hits.chrom.mi7 --genome-size 12380205 --log2 0.6 --p 0.001 --bigger-window 1.5 --annotate --minimum-windows 4
mv *.cnv ../cnv-seq_output
mv *.count ../cnv-seq_output
```

### Perform higher-resolution analysis on chromosome 15

```
../cnv-seq/cnv-seq.pl --test Coreganus_CV0087.hits.chrom.mi7 --ref Cviridis_CV0009.hits.chrom.mi7 --genome-size 12380205 --log2 0.6 --p 0.001 --window-size 500 --bigger-window 1.5 --annotate --minimum-windows 4
mv *.cnv ../cnv-seq_output/
mv *.count ../cnv-seq_output/
```

### Format CNV intervals for masking

Use awk to generate BED files with significant CNV coordinates per microchromosome.

```
tail -n +2 cnv-seq_output/Coreganus_CV0087.hits.chrom.mi1-vs-Cviridis_CV0009.hits.chrom.mi1.log2-0.6.pvalue-0.001.minw-4.cnv | awk '{if($9>0) print $0}' | sed 's/"//g' | awk 'BEGIN{OFS="\t"}{print $1, $2-1, $3, $11, $12}' > cnv.mi1.CV-CO.bed
tail -n +2 cnv-seq_output/Coreganus_CV0087.hits.chrom.mi2-vs-Cviridis_CV0009.hits.chrom.mi2.log2-0.6.pvalue-0.001.minw-4.cnv | awk '{if($9>0) print $0}' | sed 's/"//g' | awk 'BEGIN{OFS="\t"}{print $1, $2-1, $3, $11, $12}' > cnv.mi2.CV-CO.bed
tail -n +2 cnv-seq_output/Coreganus_CV0087.hits.chrom.mi7-vs-Cviridis_CV0009.hits.chrom.mi7.window-500.minw-4.cnv | awk '{if($9>0) print $0}' | sed 's/"//g' | awk 'BEGIN{OFS="\t"}{print $1, $2-1, $3, $11, $12}' > cnv.mi7.CV-CO.bed
```

Use bedtools to intersect CNVs with venom gene regions for masking.

```
bedtools intersect -a cnv.mi1.CV-CO.bed -b /data3/venom_population_genomics/venom_annotations/region_SVMP_scaffold-mi1.bed > cnv.mi1.CV-CO.SVMP_region.bed
bedtools intersect -a cnv.mi2.CV-CO.bed -b /data3/venom_population_genomics/venom_annotations/region_SVSP_scaffold-mi2.bed > cnv.mi2.CV-CO.SVSP_region.bed
bedtools intersect -a cnv.mi7.CV-CO.bed -b /data3/venom_population_genomics/venom_annotations/region_PLA2_scaffold-mi7.bed > cnv.mi7.CV-CO.PLA2_region.bed
```

Also, make a concatenated version.

```
cat cnv.mi1.CV-CO.SVMP_region.bed cnv.mi2.CV-CO.SVSP_region.bed cnv.mi7.CV-CO.PLA2_region.bed > cnv.CV-CO.venom_regions.bed
```

## Population structure analysis

Use ADMIXTURE to infer the most likely K genetic clusters in the data.

### Set up environment

```
mkdir population_structure
cd population_structure
mkdir input
```

### Make input SNP VCF

Use the all-sites VCF generated above to extract biallelic ingroup SNPs with:
* minor-allele frequencies > 0.05
* Thinned by 1 kb to reduce effects of linkage
* At least 60% of individuals with present genotypes

This will refer to the ingroup sample list in `resources/sample.cvco.list`.

```
vcftools --gzvcf ../vcf/cvco+outgroup.mask.HardFilter.vcf.gz --recode --stdout --keep sample.cvco.list --bed chrom.bed --min-alleles 2 --max-alleles 2 --maf 0.05 --thin 1000 --max-missing 0.4 | bgzip -c > ../vcf/cvco.ingroup.mask.HardFilter.chrom.snp.maf05.thin1kb.miss04.vcf.gz
```

### Convert VCF to ADMIXTURE input

Use plink to convert the VCF to .ped input format read by ADMIXTURE.

```
plink --vcf ../vcf/cvco.ingroup.mask.HardFilter.chrom.snp.maf05.thin1kb.miss04.vcf.gz --make-bed --out ./input/cvco.ingroup.mask.HardFilter.chrom.snp.maf05.thin1kb.miss04 --allow-extra-chr --recode12
```

### Run ADMIXTURE over a series of K values

The script below will perform ADMIXTURE analyses for K 1-16.

__*run_admixture.sh*__

```
ped=$1
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; do
	admixture --cv $ped $K | tee log${K}.out
done
```

Run the script.

```
sh run_admixture ./input/cvco.ingroup.mask.HardFilter.chrom.snp.maf05.thin10kb.miss04.ped
```

### Evaluate the results

Look at the CV error of each K value to determine the best-supported number of genetic clusters.

```
grep -h CV log*.out
```

## Demographic analysis

Use the pairwise Markovian sequential coalescent (PSMC) to estimate effective population size through time.

The steps below follow [Al Ludington's PSMC workflow](https://github.com/a-lud/snakes-demographic-history) fairly closely.

### Set up environment

```
mkdir psmc_analysis
cd psmc_analysis
mkdir log
mkdir input
mkdir results
```

### Install local implementation of PSMC

```
git clone https://github.com/lh3/psmc.git
cd psmc
make
cd utils
make
```

### Choose samples and prepare inputs

PSMC takes diploid sequences from an individual as input. These can be generated using mappings (i.e., bam files) and individual variant calls using samtools/bcftools.

Individual inputs will be generated for a single male from each population:
* CV0632 (CV1) 
* CV0860 (CV2)
* CV0151 (CO1)
* CV0781 (CO2)

The list of samples is in `resources/sample.psmc.list`.

#### 1. Call consensus variants per individual and index output

This will also mask sites overlapping repeat annotations.

```
bcftools mpileup -C 50 -q 30 -Q 25 -Ou -f ../CroVir_genome_L77pg_16Aug2017.final_rename.fasta ../bam/CV0632.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=5 & DP<=50" | bcftools view --exclude-types indels -T ^CroVir_genome_L77pg_16Aug2017.repeat.masked.final.sort.bed | bcftools sort --temp-dir ./input/temp_CV0632 -Oz -o ./input/CV1_CV0632.vcf.gz
bcftools mpileup -C 50 -q 30 -Q 25 -Ou -f ../CroVir_genome_L77pg_16Aug2017.final_rename.fasta ../bam/CV0860.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=5 & DP<=50" | bcftools view --exclude-types indels -T ^CroVir_genome_L77pg_16Aug2017.repeat.masked.final.sort.bed | bcftools sort --temp-dir ./input/temp_CV0860 -Oz -o ./input/CV2_CV0860.vcf.gz
bcftools mpileup -C 50 -q 30 -Q 25 -Ou -f ../CroVir_genome_L77pg_16Aug2017.final_rename.fasta ../bam/CV0151.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=5 & DP<=50" | bcftools view --exclude-types indels -T ^CroVir_genome_L77pg_16Aug2017.repeat.masked.final.sort.bed | bcftools sort --temp-dir ./input/temp_CV0151 -Oz -o ./input/CO1_CV0151.vcf.gz
bcftools mpileup -C 50 -q 30 -Q 25 -Ou -f ../CroVir_genome_L77pg_16Aug2017.final_rename.fasta ../bam/CV0781.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=5 & DP<=50" | bcftools view --exclude-types indels -T ^CroVir_genome_L77pg_16Aug2017.repeat.masked.final.sort.bed | bcftools sort --temp-dir ./input/temp_CV0781 -Oz -o ./input/CO2_CV0781.vcf.gz

tabix -p vcf ./input/CV1_CV0632.vcf.gz
tabix -p vcf ./input/CV2_CV0860.vcf.gz
tabix -p vcf ./input/CO1_CV0151.vcf.gz
tabix -p vcf ./input/CO2_CV0781.vcf.gz
```

#### 2. Get diploid sequence

```
bcftools view ./input/CV1_CV0632.vcf.gz | vcfutils.pl vcf2fq | gzip > ./input/CV1_CV0632.fastq.gz
bcftools view ./input/CV2_CV0860.vcf.gz | vcfutils.pl vcf2fq | gzip > ./input/CV2_CV0860.fastq.gz
bcftools view ./input/CO1_CV0151.vcf.gz | vcfutils.pl vcf2fq | gzip > ./input/CO1_CV0151.fastq.gz
bcftools view ./input/CO2_CV0781.vcf.gz | vcfutils.pl vcf2fq | gzip > ./input/CO2_CV0781.fastq.gz
```

#### 3. Convert to psmcfa format

```
./psmc/utils/fq2psmcfa -q 20 ./input/CV1_CV0632.fastq.gz > ./input/CV1_CV0632.psmcfa
./psmc/utils/fq2psmcfa -q 20 ./input/CV2_CV0860.fastq.gz > ./input/CV2_CV0860.psmcfa
./psmc/utils/fq2psmcfa -q 20 ./input/CO1_CV0151.fastq.gz > ./input/CO1_CV0151.psmcfa
./psmc/utils/fq2psmcfa -q 20 ./input/CO2_CV0781.fastq.gz > ./input/CO2_CV0781.psmcfa
```

### Choice of time segment patterns

These are specified by the '-p' flag in PSMC analysis.

Run analyses with these patterns:
* 4+5*3+4 (default)
* 4+10*3+6+8
* 4+25*2+4+6
* 4+30*2+4+6+10

### Run main analyses

```
./psmc/psmc -N25 -t12 -r5 -p "4+30*2+4+6+10" -o ./results/pattern4.CV1_CV0632.diploid.psmc ./input/CV1_CV0632.psmcfa
./psmc/psmc -N25 -t12 -r5 -p "4+30*2+4+6+10" -o ./results/pattern4.CV2_CV0860.diploid.psmc ./input/CV2_CV0860.psmcfa
./psmc/psmc -N25 -t12 -r5 -p "4+30*2+4+6+10" -o ./results/pattern4.CO1_CV0151.diploid.psmc ./input/CO1_CV0151.psmcfa
./psmc/psmc -N25 -t12 -r5 -p "4+30*2+4+6+10" -o ./results/pattern4.CO2_CV0781.diploid.psmc ./input/CO2_CV0781.psmcfa
```

Repeat the commands above with the different time segment patterns, if desired.

### Bootstrap analysis

Perform a series of bootstraps to vet main analysis results.

#### 1. Generate 'split' psmcfa inputs for bootstrapping

```
./psmc/utils/splitfa ./input/CV1_CV0632.psmcfa > ./input/CV1_CV0632-split.psmcfa
./psmc/utils/splitfa ./input/CV2_CV0860.psmcfa > ./input/CV2_CV0860-split.psmcfa
./psmc/utils/splitfa ./input/CO1_CV0151.psmcfa > ./input/CO1_CV0151-split.psmcfa
./psmc/utils/splitfa ./input/CO2_CV0781.psmcfa > ./input/CO2_CV0781-split.psmcfa
```

#### 2. Perform bootstrapping analysis

This will perform 100 bootstrap replicates per sample.

```
mkdir ./results/bootstrap
seq 100 | xargs -i echo ./psmc/psmc -N25 -t15 -r5 -b -p "4+30*2+4+6+10" -o ./results/bootstrap/CV1_CV0632.round-{}.psmc ./input/CV1_CV0632-split.psmcfa | sh
seq 100 | xargs -i echo ./psmc/psmc -N25 -t15 -r5 -b -p "4+30*2+4+6+10" -o ./results/bootstrap/CV2_CV0860.round-{}.psmc ./input/CV2_CV0860-split.psmcfa | sh
seq 100 | xargs -i echo ./psmc/psmc -N25 -t15 -r5 -b -p "4+30*2+4+6+10" -o ./results/bootstrap/CO1_CV0151.round-{}.psmc ./input/CO1_CV0151-split.psmcfa | sh
seq 100 | xargs -i echo ./psmc/psmc -N25 -t15 -r5 -b -p "4+30*2+4+6+10" -o ./results/bootstrap/CO2_CV0781.round-{}.psmc ./input/CO2_CV0781-split.psmcfa | sh
```

### Combine main and bootstrap results and plot

Concatenate main and bootstrap results files.

```
cat ./results/pattern4.CV1_CV0632.diploid.psmc ./results/bootstrap/CV1_CV0632.round-*.psmc > ./results/pattern4_combined.CV1_CV0632.psmc
cat ./results/pattern4.CV2_CV0860.diploid.psmc ./results/bootstrap/CV2_CV0860.round-*.psmc > ./results/pattern4_combined.CV2_CV0860.psmc
cat ./results/pattern4.CO1_CV0151.diploid.psmc ./results/bootstrap/CO1_CV0151.round-*.psmc > ./results/pattern4_combined.CO1_CV0151.psmc
cat ./results/pattern4.CO2_CV0781.diploid.psmc ./results/bootstrap/CO2_CV0781.round-*.psmc > ./results/pattern4_combined.CO2_CV0781.psmc
```

Plot combined results.

```
./psmc/utils/psmc2history.pl ./results/pattern4_combined.CV1_CV0632.psmc | ./psmc/utils/history2ms.pl > ms-cmd.sh
./psmc/utils/psmc_plot.pl -u 0.2e-08 -g 3 ./results/pattern4_combined.CV1_CV0632 ./results/pattern4_combined.CV1_CV0632.psmc
./psmc/utils/psmc2history.pl ./results/pattern4_combined.CV2_CV0860.psmc | ./psmc/utils/history2ms.pl > ms-cmd.sh
./psmc/utils/psmc_plot.pl -u 0.2e-08 -g 3 ./results/pattern4_combined.CV2_CV0860 ./results/pattern4_combined.CV2_CV0860.psmc
./psmc/utils/psmc2history.pl ./results/pattern4_combined.CO1_CV0151.psmc | ./psmc/utils/history2ms.pl > ms-cmd.sh
./psmc/utils/psmc_plot.pl -u 0.2e-08 -g 3 ./results/pattern4_combined.CO1_CV0151 ./results/pattern4_combined.CO1_CV0151.psmc
./psmc/utils/psmc2history.pl ./results/pattern4_combined.CO2_CV0781.psmc | ./psmc/utils/history2ms.pl > ms-cmd.sh
./psmc/utils/psmc_plot.pl -u 0.2e-08 -g 3 ./results/pattern4_combined.CO2_CV0781 ./results/pattern4_combined.CO2_CV0781.psmc
```

Convert EPS output to PDF.

```
cd ./results
epstopdf.pl pattern4_combined.CV1_CV0632.eps
epstopdf.pl pattern4_combined.CV2_CV0860.eps
epstopdf.pl pattern4_combined.CO1_CV0151.eps
epstopdf.pl pattern4_combined.CO2_CV0781.eps
```

## Population genetic diversity and differentiation

Estimate population genetic diversity and differentiation parameters (π, dxy, Fst) across the genome in sliding windows using [pixy](https://pixy.readthedocs.io/en/latest/).

### Set up environment

```
mkdir pixy
cd pixy
mkdir pixy_results
mkdir pixy_zarr
```

Pixy analyses are guided by a two-column, tab-delimited population map with population IDs for each sample. The sample names should match exactly the names in the VCF header. The population map is in `resources/pixy.popmap`.

### Run pixy analysis

*Note: I installed pixy in its own conda environment based in Python 3.6.*

#### 1. Activate pixy conda environment

```
conda deactivate
conda activate pixy
```

#### 2. Run analysis on each chromosome VCF

This script will execute pixy on each of the chromosome-specific VCFs parsed above. It takes the chromosome list in `resources/chrom.list`, a window size (integer), window abbreviation (e.g., 1 kb) as input, and also wants to know if you want to use existing zarr files from a previous run (yes/no).

__*pixyloop.sh*__

```
list=$1
window=$2
abbrev=$3
ans=$4
for chrom in `cat $list`; do
	pixy --stats pi fst dxy --vcf ../vcf/vcf_chrom-specific_cvco+outgroup/cvco+outgroup.mask.HardFilter.depth.chrom.$chrom.vcf.gz --zarr_path ./pixy_zarr --reuse_zarr $ans --window_size $window --populations pixy.popmap --variant_filter_expression 'DP>=5' --invariant_filter_expression 'DP>=5' --outfile_prefix ./pixy_results/pixy.$chrom.$abbrev
done
```

Run the script with various window sizes.

```
sh pixyloop.sh chrom.list 100000 100kb no
sh pixyloop.sh chrom.list 10000 10kb yes
sh pixyloop.sh chrom.list 1000 1kb yes
```

#### 3. Run higher-resolution analysis on chromosome 15

```
pixy --stats pi fst dxy --vcf ../vcf/vcf_chrom-specific_cvco+outgroup/cvco+outgroup.mask.HardFilter.depth.chrom.scaffold-mi7.vcf.gz --zarr_path ./pixy_zarr --reuse_zarr yes --window_size 250 --populations pixy.popmap --variant_filter_expression 'DP>=5' --invariant_filter_expression 'DP>=5' --outfile_prefix ./pixy_results/pixy.scaffold-mi7.250bp
```

### Concatenate results from different chromosomes

The script below will combine chromosome-specific results into a single file per statistic.

__*concatenatePixyResults.sh*__

```
abbrev=$1
# concatenate pi output
head -n 1 ./pixy_results/pixy.scaffold-ma1.${abbrev}_pi.txt > ./pixy_results/pixy.all.${abbrev}_pi.txt
for chrom in `cat chrom.list`; do
	tail -n +2 ./pixy_results/pixy.${chrom}.${abbrev}_pi.txt >> ./pixy_results/pixy.all.${abbrev}_pi.txt
done
# concatenate dxy output
head -n 1 ./pixy_results/pixy.scaffold-ma1.${abbrev}_dxy.txt > ./pixy_results/pixy.all.${abbrev}_dxy.txt
for chrom in `cat chrom.list`; do
	tail -n +2 ./pixy_results/pixy.${chrom}.${abbrev}_dxy.txt >> ./pixy_results/pixy.all.${abbrev}_dxy.txt
done
# concatenate fst output
head -n 1 ./pixy_results/pixy.scaffold-ma1.${abbrev}_fst.txt > ./pixy_results/pixy.all.${abbrev}_fst.txt
for chrom in `cat chrom.list`; do
	tail -n +2 ./pixy_results/pixy.${chrom}.${abbrev}_fst.txt >> ./pixy_results/pixy.all.${abbrev}_fst.txt
done
```

Run the script to concatenate results at different window sizes.

```
sh concatenatePixyResults.sh 100kb
sh concatenatePixyResults.sh 10kb
sh concatenatePixyResults.sh 1kb
```

### Diversity appendix 1: CNV-masking

Mask results in genomic windows overlapping with significant copy-number variation between C. viridis and C. oreganus.

#### Set up environment

```
mkdir cnv_masked_results
cd cnv_masked_results
mkdir mask_bed
```

#### Use bedtools intersect to get 10 kb windows to mask

This uses the CNV BED file in `resources/cnv.CV-CO.venom_regions.bed`.

```
for chrom in scaffold-mi1 scaffold-mi2 scaffold-mi7; do tail -n +2 ../pixy/pixy_results/pixy.$chrom.10kb_pi.txt | grep 'CO1' | awk 'BEGIN{OFS="\t"}{print $2,$3,$4}' | bedtools intersect -u -a - -b cnv.CV-CO.venom_regions.bed > mask_bed/cnv_filter.10kb.$chrom.pixy.bed; done
```

#### Use bedtools intersect to get 1 kb and 250 bp windows to mask

```
for chrom in scaffold-mi1 scaffold-mi2 scaffold-mi7; do tail -n +2 ../pixy/pixy_results/pixy.$chrom.1kb_pi.txt | grep 'CO1' | awk 'BEGIN{OFS="\t"}{print $2,$3,$4}' | bedtools intersect -u -a - -b mask_bed/cnv_filter.10kb.$chrom.pixy.bed > mask_bed/cnv_filter.1kb.$chrom.pixy.bed; done
tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi7.250bp_pi.txt | grep 'CO1' | awk 'BEGIN{OFS="\t"}{print $2,$3,$4}' | bedtools intersect -u -a - -b mask_bed/cnv_filter.10kb.scaffold-mi7.pixy.bed > mask_bed/cnv_filter.250bp.scaffold-mi7.pixy.bed
```

#### Mask results including CO1 and CO2 in CNV windows as 'NA's

The script `cnvMaskSelection.py` has code blocks to recognize different results files and mask windows based on BED files in `resources`.

The script is in the `python` directory. Its third argument asks whether to print the input file header (yes/no).

Mask pixy results with script.

```
for chrom in scaffold-mi1 scaffold-mi2 scaffold-mi7; do for window in 10kb 1kb; do head -n 1 ../pixy/pixy_results/pixy.$chrom.${window}_pi.txt > pi.$chrom.$window.co1.cnvMask.txt; python cnvMaskSelection.py /data3/venom_population_genomics/3_selection/pixy/pixy_results/pixy.$chrom.${window}_pi.txt ./mask_bed/cnv_filter.$window.$chrom.pixy.bed no | grep 'CO1' >> pi.$chrom.$window.co1.cnvMask.txt; done; done
for chrom in scaffold-mi1 scaffold-mi2 scaffold-mi7; do for window in 10kb 1kb; do head -n 1 ../pixy/pixy_results/pixy.$chrom.${window}_pi.txt > pi.$chrom.$window.co2.cnvMask.txt; python cnvMaskSelection.py /data3/venom_population_genomics/3_selection/pixy/pixy_results/pixy.$chrom.${window}_pi.txt ./mask_bed/cnv_filter.$window.$chrom.pixy.bed no | grep 'CO2' >> pi.$chrom.$window.co2.cnvMask.txt; done; done
head -n 1 ../pixy/pixy_results/pixy.scaffold-mi7.250bp_pi.txt > pi.scaffold-mi7.250bp.co1.cnvMask.txt; python cnvMaskSelection.py ../pixy/pixy_results/pixy.scaffold-mi7.250bp_pi.txt ./mask_bed/cnv_filter.250bp.scaffold-mi7.pixy.bed no | grep 'CO1' >> pi.scaffold-mi7.250bp.co1.cnvMask.txt
head -n 1 ../pixy/pixy_results/pixy.scaffold-mi7.250bp_pi.txt > pi.scaffold-mi7.250bp.co2.cnvMask.txt; python cnvMaskSelection.py ../pixy/pixy_results/pixy.scaffold-mi7.250bp_pi.txt ./mask_bed/cnv_filter.250bp.scaffold-mi7.pixy.bed no | grep 'CO2' >> pi.scaffold-mi7.250bp.co2.cnvMask.txt

for chrom in scaffold-mi1 scaffold-mi2 scaffold-mi7; do for window in 10kb 1kb; do head -n 1 ../pixy/pixy_results/pixy.$chrom.${window}_dxy.txt > dxy.$chrom.$window.cv1co1.cnvMask.txt; python cnvMaskSelection.py ../pixy/pixy_results/pixy.$chrom.${window}_dxy.txt ./mask_bed/cnv_filter.$window.$chrom.pixy.bed no | grep -P "CV1\tCO1" >> dxy.$chrom.$window.cv1co1.cnvMask.txt; done; done
for chrom in scaffold-mi1 scaffold-mi2 scaffold-mi7; do for window in 10kb 1kb; do head -n 1 ../pixy/pixy_results/pixy.$chrom.${window}_dxy.txt > dxy.$chrom.$window.co1co2.cnvMask.txt; python cnvMaskSelection.py ../pixy/pixy_results/pixy.$chrom.${window}_dxy.txt ./mask_bed/cnv_filter.$window.$chrom.pixy.bed no | grep -P "CO1\tCO2" >> dxy.$chrom.$window.co1co2.cnvMask.txt; done; done
head -n 1 ../pixy/pixy_results/pixy.scaffold-mi7.250bp_dxy.txt > dxy.scaffold-mi7.250bp.cv1co1.cnvMask.txt; python cnvMaskSelection.py ../pixy/pixy_results/pixy.scaffold-mi7.250bp_dxy.txt ./mask_bed/cnv_filter.250bp.scaffold-mi7.pixy.bed no | grep -P "CV1\tCO1" >> dxy.scaffold-mi7.250bp.cv1co1.cnvMask.txt
head -n 1 ../pixy/pixy_results/pixy.scaffold-mi7.250bp_dxy.txt > dxy.scaffold-mi7.250bp.co1co2.cnvMask.txt; python cnvMaskSelection.py ../pixy/pixy_results/pixy.scaffold-mi7.250bp_dxy.txt ./mask_bed/cnv_filter.250bp.scaffold-mi7.pixy.bed no | grep -P "CO1\tCO2" >> dxy.scaffold-mi7.250bp.co1co2.cnvMask.txt

for chrom in scaffold-mi1 scaffold-mi2 scaffold-mi7; do for window in 10kb 1kb; do head -n 1 ../pixy/pixy_results/pixy.$chrom.${window}_fst.txt > fst.$chrom.$window.cv1co1.cnvMask.txt; python cnvMaskSelection.py ../pixy/pixy_results/pixy.$chrom.${window}_fst.txt ./mask_bed/cnv_filter.$window.$chrom.pixy.bed no | grep -P "CV1\tCO1" >> fst.$chrom.$window.cv1co1.cnvMask.txt; done; done
for chrom in scaffold-mi1 scaffold-mi2 scaffold-mi7; do for window in 10kb 1kb; do head -n 1 ../pixy/pixy_results/pixy.$chrom.${window}_fst.txt > fst.$chrom.$window.co1co2.cnvMask.txt; python cnvMaskSelection.py ../pixy/pixy_results/pixy.$chrom.${window}_fst.txt ./mask_bed/cnv_filter.$window.$chrom.pixy.bed no | grep -P "CO1\tCO2" >> fst.$chrom.$window.co1co2.cnvMask.txt; done; done
head -n 1 ../pixy/pixy_results/pixy.scaffold-mi7.250bp_fst.txt > fst.scaffold-mi7.250bp.cv1co1.cnvMask.txt; python cnvMaskSelection.py ../pixy/pixy_results/pixy.scaffold-mi7.250bp_fst.txt ./mask_bed/cnv_filter.250bp.scaffold-mi7.pixy.bed no | grep -P "CV1\tCO1" >> fst.scaffold-mi7.250bp.cv1co1.cnvMask.txt
head -n 1 ../pixy/pixy_results/pixy.scaffold-mi7.250bp_fst.txt > fst.scaffold-mi7.250bp.co1co2.cnvMask.txt; python cnvMaskSelection.py ../pixy/pixy_results/pixy.scaffold-mi7.250bp_fst.txt ./mask_bed/cnv_filter.250bp.scaffold-mi7.pixy.bed no | grep -P "CO1\tCO2" >> fst.scaffold-mi7.250bp.co1co2.cnvMask.txt
```

### Diversity appendix 2: Comparison between venom genic and intergenic regions

Compare π, dxy, and Fst between genic and intergenic intervals of venom gene regions.

#### Set up environment

```
mkdir gene_vs_intergenic
cd gene_vs_intergenic
```

#### Use bedtools to parse intergenic and genic windows in the SVMP region

Intergenic:

```
echo -e "chrom\tstart\tend\tpi_mean" > svmp_intergenic.pi_mean.cv1.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi1.1kb_pi.txt | grep 'CV1' | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -wb -a region_SVMP_scaffold-mi1.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_SVMP.bed >> svmp_intergenic.pi_mean.cv1.txt
echo -e "chrom\tstart\tend\tpi_mean" > svmp_intergenic.pi_mean.cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi1.1kb_pi.txt | grep 'CV2' | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -wb -a region_SVMP_scaffold-mi1.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_SVMP.bed >> svmp_intergenic.pi_mean.cv2.txt
echo -e "chrom\tstart\tend\tpi_mean" > svmp_intergenic.pi_mean.co1.txt; tail -n +2 ../cnv_masked_results/pi.scaffold-mi1.1kb.co1.cnvMask.txt | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -wb -a region_SVMP_scaffold-mi1.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_SVMP.bed >> svmp_intergenic.pi_mean.co1.txt
echo -e "chrom\tstart\tend\tpi_mean" > svmp_intergenic.pi_mean.co2.txt; tail -n +2 ../cnv_masked_results/pi.scaffold-mi1.1kb.co2.cnvMask.txt | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -wb -a region_SVMP_scaffold-mi1.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_SVMP.bed >> svmp_intergenic.pi_mean.co2.txt

echo -e "chrom\tstart\tend\tdxy_mean" > svmp_intergenic.dxy_mean.cv1co1.txt; tail -n +2 ../cnv_masked_results/dxy.scaffold-mi1.1kb.cv1co1.cnvMask.txt | grep -P 'CV1\tCO1' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -wb -a region_SVMP_scaffold-mi1.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_SVMP.bed >> svmp_intergenic.dxy_mean.cv1co1.txt
echo -e "chrom\tstart\tend\tdxy_mean" > svmp_intergenic.dxy_mean.cv1cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi1.1kb_dxy.txt | grep -P 'CV1\tCV2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -wb -a region_SVMP_scaffold-mi1.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_SVMP.bed >> svmp_intergenic.dxy_mean.cv1cv2.txt
echo -e "chrom\tstart\tend\tdxy_mean" > svmp_intergenic.dxy_mean.co1co2.txt; tail -n +2 ../cnv_masked_results/dxy.scaffold-mi1.1kb.co1co2.cnvMask.txt | grep -P 'CO1\tCO2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -wb -a region_SVMP_scaffold-mi1.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_SVMP.bed >> svmp_intergenic.dxy_mean.co1co2.txt

echo -e "chrom\tstart\tend\tfst_mean" > svmp_intergenic.fst_mean.cv1co1.txt; tail -n +2 ../cnv_masked_results/fst.scaffold-mi1.1kb.cv1co1.cnvMask.txt | grep -v 'nan' | grep -P 'CV1\tCO1' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -wb -a region_SVMP_scaffold-mi1.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_SVMP.bed >> svmp_intergenic.fst_mean.cv1co1.txt
echo -e "chrom\tstart\tend\tfst_mean" > svmp_intergenic.fst_mean.cv1cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi1.1kb_fst.txt | grep -v 'nan' | grep -P 'CV1\tCV2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -wb -a region_SVMP_scaffold-mi1.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_SVMP.bed >> svmp_intergenic.fst_mean.cv1cv2.txt
echo -e "chrom\tstart\tend\tfst_mean" > svmp_intergenic.fst_mean.co1co2.txt; tail -n +2 ../cnv_masked_results/fst.scaffold-mi1.1kb.co1co2.cnvMask.txt | grep -v 'nan' | grep -P 'CO1\tCO2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -wb -a region_SVMP_scaffold-mi1.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_SVMP.bed >> svmp_intergenic.fst_mean.co1co2.txt
```

Genic:

```
echo -e "chrom\tstart\tend\tpi_mean" > svmp_genic.pi_mean.cv1.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi1.1kb_pi.txt | grep 'CV1' | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -a - -b gene_SVMP.bed >> svmp_genic.pi_mean.cv1.txt
echo -e "chrom\tstart\tend\tpi_mean" > svmp_genic.pi_mean.cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi1.1kb_pi.txt | grep 'CV2' | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -a - -b gene_SVMP.bed >> svmp_genic.pi_mean.cv2.txt
echo -e "chrom\tstart\tend\tpi_mean" > svmp_genic.pi_mean.co1.txt; tail -n +2 ../cnv_masked_results/pi.scaffold-mi1.1kb.co1.cnvMask.txt | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -a - -b gene_SVMP.bed >> svmp_genic.pi_mean.co1.txt
echo -e "chrom\tstart\tend\tpi_mean" > svmp_genic.pi_mean.co2.txt; tail -n +2 ../cnv_masked_results/pi.scaffold-mi1.1kb.co2.cnvMask.txt | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -a - -b gene_SVMP.bed >> svmp_genic.pi_mean.co2.txt

echo -e "chrom\tstart\tend\tdxy_mean" > svmp_genic.dxy_mean.cv1co1.txt; tail -n +2 ../cnv_masked_results/dxy.scaffold-mi1.1kb.cv1co1.cnvMask.txt | grep -P 'CV1\tCO1' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b gene_SVMP.bed >> svmp_genic.dxy_mean.cv1co1.txt
echo -e "chrom\tstart\tend\tdxy_mean" > svmp_genic.dxy_mean.cv1cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi1.1kb_dxy.txt | grep -P 'CV1\tCV2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b gene_SVMP.bed >> svmp_genic.dxy_mean.cv1cv2.txt
echo -e "chrom\tstart\tend\tdxy_mean" > svmp_genic.dxy_mean.co1co2.txt; tail -n +2 ../cnv_masked_results/dxy.scaffold-mi1.1kb.co1co2.cnvMask.txt | grep -P 'CO1\tCO2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b gene_SVMP.bed >> svmp_genic.dxy_mean.co1co2.txt

echo -e "chrom\tstart\tend\tfst_mean" > svmp_genic.fst_mean.cv1co1.txt; tail -n +2 ../cnv_masked_results/fst.scaffold-mi1.1kb.cv1co1.cnvMask.txt | grep -v 'nan' | grep -P 'CV1\tCO1' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b gene_SVMP.bed >> svmp_genic.fst_mean.cv1co1.txt
echo -e "chrom\tstart\tend\tfst_mean" > svmp_genic.fst_mean.cv1cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi1.1kb_fst.txt | grep -v 'nan' | grep -P 'CV1\tCV2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b gene_SVMP.bed >> svmp_genic.fst_mean.cv1cv2.txt
echo -e "chrom\tstart\tend\tfst_mean" > svmp_genic.fst_mean.co1co2.txt; tail -n +2 ../cnv_masked_results/fst.scaffold-mi1.1kb.co1co2.cnvMask.txt | grep -v 'nan' | grep -P 'CO1\tCO2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b gene_SVMP.bed >> svmp_genic.fst_mean.co1co2.txt
```

#### Use bedtools to parse intergenic and genic windows in the SVSP region

Intergenic:

```
echo -e "chrom\tstart\tend\tpi_mean" > svsp_intergenic.pi_mean.cv1.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi2.1kb_pi.txt | grep 'CV1' | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -wb -a region_SVSP_scaffold-mi2.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_SVSP.bed >> svsp_intergenic.pi_mean.cv1.txt
echo -e "chrom\tstart\tend\tpi_mean" > svsp_intergenic.pi_mean.cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi2.1kb_pi.txt | grep 'CV2' | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -wb -a region_SVSP_scaffold-mi2.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_SVSP.bed >> svsp_intergenic.pi_mean.cv2.txt
echo -e "chrom\tstart\tend\tpi_mean" > svsp_intergenic.pi_mean.co1.txt; tail -n +2 ../cnv_masked_results/pi.scaffold-mi2.1kb.co1.cnvMask.txt | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -wb -a region_SVSP_scaffold-mi2.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_SVSP.bed >> svsp_intergenic.pi_mean.co1.txt
echo -e "chrom\tstart\tend\tpi_mean" > svsp_intergenic.pi_mean.co2.txt; tail -n +2 ../cnv_masked_results/pi.scaffold-mi2.1kb.co2.cnvMask.txt | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -wb -a region_SVSP_scaffold-mi2.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_SVSP.bed >> svsp_intergenic.pi_mean.co2.txt

echo -e "chrom\tstart\tend\tdxy_mean" > svsp_intergenic.dxy_mean.cv1co1.txt; tail -n +2 ../cnv_masked_results/dxy.scaffold-mi2.1kb.cv1co1.cnvMask.txt | grep -P 'CV1\tCO1' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -wb -a region_SVSP_scaffold-mi2.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_SVSP.bed >> svsp_intergenic.dxy_mean.cv1co1.txt
echo -e "chrom\tstart\tend\tdxy_mean" > svsp_intergenic.dxy_mean.cv1cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi2.1kb_dxy.txt | grep -P 'CV1\tCV2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -wb -a region_SVSP_scaffold-mi2.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_SVSP.bed >> svsp_intergenic.dxy_mean.cv1cv2.txt
echo -e "chrom\tstart\tend\tdxy_mean" > svsp_intergenic.dxy_mean.co1co2.txt; tail -n +2 ../cnv_masked_results/dxy.scaffold-mi2.1kb.co1co2.cnvMask.txt | grep -P 'CO1\tCO2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -wb -a region_SVSP_scaffold-mi2.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_SVSP.bed >> svsp_intergenic.dxy_mean.co1co2.txt

echo -e "chrom\tstart\tend\tfst_mean" > svsp_intergenic.fst_mean.cv1co1.txt; tail -n +2 ../cnv_masked_results/fst.scaffold-mi2.1kb.cv1co1.cnvMask.txt | grep -v 'nan' | grep -P 'CV1\tCO1' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -wb -a region_SVSP_scaffold-mi2.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_SVSP.bed >> svsp_intergenic.fst_mean.cv1co1.txt
echo -e "chrom\tstart\tend\tfst_mean" > svsp_intergenic.fst_mean.cv1cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi2.1kb_fst.txt | grep -v 'nan' | grep -P 'CV1\tCV2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -wb -a region_SVSP_scaffold-mi2.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_SVSP.bed >> svsp_intergenic.fst_mean.cv1cv2.txt
echo -e "chrom\tstart\tend\tfst_mean" > svsp_intergenic.fst_mean.co1co2.txt; tail -n +2 ../cnv_masked_results/fst.scaffold-mi2.1kb.co1co2.cnvMask.txt | grep -v 'nan' | grep -P 'CO1\tCO2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -wb -a region_SVSP_scaffold-mi2.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_SVSP.bed >> svsp_intergenic.fst_mean.co1co2.txt
```

Genic:

```
echo -e "chrom\tstart\tend\tpi_mean" > svsp_genic.pi_mean.cv1.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi2.1kb_pi.txt | grep 'CV1' | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -a - -b gene_SVSP.bed >> svsp_genic.pi_mean.cv1.txt
echo -e "chrom\tstart\tend\tpi_mean" > svsp_genic.pi_mean.cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi2.1kb_pi.txt | grep 'CV2' | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -a - -b gene_SVSP.bed >> svsp_genic.pi_mean.cv2.txt
echo -e "chrom\tstart\tend\tpi_mean" > svsp_genic.pi_mean.co1.txt; tail -n +2 ../cnv_masked_results/pi.scaffold-mi2.1kb.co1.cnvMask.txt | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -a - -b gene_SVSP.bed >> svsp_genic.pi_mean.co1.txt
echo -e "chrom\tstart\tend\tpi_mean" > svsp_genic.pi_mean.co2.txt; tail -n +2 ../cnv_masked_results/pi.scaffold-mi2.1kb.co2.cnvMask.txt | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -a - -b gene_SVSP.bed >> svsp_genic.pi_mean.co2.txt

echo -e "chrom\tstart\tend\tdxy_mean" > svsp_genic.dxy_mean.cv1co1.txt; tail -n +2 ../cnv_masked_results/dxy.scaffold-mi2.1kb.cv1co1.cnvMask.txt | grep -P 'CV1\tCO1' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b gene_SVSP.bed >> svsp_genic.dxy_mean.cv1co1.txt
echo -e "chrom\tstart\tend\tdxy_mean" > svsp_genic.dxy_mean.cv1cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi2.1kb_dxy.txt | grep -P 'CV1\tCV2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b gene_SVSP.bed >> svsp_genic.dxy_mean.cv1cv2.txt
echo -e "chrom\tstart\tend\tdxy_mean" > svsp_genic.dxy_mean.co1co2.txt; tail -n +2 ../cnv_masked_results/dxy.scaffold-mi2.1kb.co1co2.cnvMask.txt | grep -P 'CO1\tCO2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b gene_SVSP.bed >> svsp_genic.dxy_mean.co1co2.txt

echo -e "chrom\tstart\tend\tfst_mean" > svsp_genic.fst_mean.cv1co1.txt; tail -n +2 ../cnv_masked_results/fst.scaffold-mi2.1kb.cv1co1.cnvMask.txt | grep -v 'nan' | grep -P 'CV1\tCO1' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b gene_SVSP.bed >> svsp_genic.fst_mean.cv1co1.txt
echo -e "chrom\tstart\tend\tfst_mean" > svsp_genic.fst_mean.cv1cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi2.1kb_fst.txt | grep -v 'nan' | grep -P 'CV1\tCV2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b gene_SVSP.bed >> svsp_genic.fst_mean.cv1cv2.txt
echo -e "chrom\tstart\tend\tfst_mean" > svsp_genic.fst_mean.co1co2.txt; tail -n +2 ../cnv_masked_results/fst.scaffold-mi2.1kb.co1co2.cnvMask.txt | grep -v 'nan' | grep -P 'CO1\tCO2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b gene_SVSP.bed >> svsp_genic.fst_mean.co1co2.txt
```

#### Use bedtools to parse intergenic and genic windows in the PLA2 region

Intergenic:

```
echo -e "chrom\tstart\tend\tpi_mean" > pla2_intergenic.pi_mean.cv1.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi7.1kb_pi.txt | grep 'CV1' | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -wb -a region_PLA2_scaffold-mi7.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_PLA2.bed >> pla2_intergenic.pi_mean.cv1.txt
echo -e "chrom\tstart\tend\tpi_mean" > pla2_intergenic.pi_mean.cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi7.1kb_pi.txt | grep 'CV2' | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -wb -a region_PLA2_scaffold-mi7.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_PLA2.bed >> pla2_intergenic.pi_mean.cv2.txt
echo -e "chrom\tstart\tend\tpi_mean" > pla2_intergenic.pi_mean.co1.txt; tail -n +2 ../cnv_masked_results/pi.scaffold-mi7.1kb.co1.cnvMask.txt | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -wb -a region_PLA2_scaffold-mi7.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_PLA2.bed >> pla2_intergenic.pi_mean.co1.txt
echo -e "chrom\tstart\tend\tpi_mean" > pla2_intergenic.pi_mean.co2.txt; tail -n +2 ../cnv_masked_results/pi.scaffold-mi7.1kb.co2.cnvMask.txt | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -wb -a region_PLA2_scaffold-mi7.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_PLA2.bed >> pla2_intergenic.pi_mean.co2.txt

echo -e "chrom\tstart\tend\tdxy_mean" > pla2_intergenic.dxy_mean.cv1co1.txt; tail -n +2 ../cnv_masked_results/dxy.scaffold-mi7.1kb.cv1co1.cnvMask.txt | grep -P 'CV1\tCO1' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -wb -a region_PLA2_scaffold-mi7.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_PLA2.bed >> pla2_intergenic.dxy_mean.cv1co1.txt
echo -e "chrom\tstart\tend\tdxy_mean" > pla2_intergenic.dxy_mean.cv1cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi7.1kb_dxy.txt | grep -P 'CV1\tCV2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -wb -a region_PLA2_scaffold-mi7.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_PLA2.bed >> pla2_intergenic.dxy_mean.cv1cv2.txt
echo -e "chrom\tstart\tend\tdxy_mean" > pla2_intergenic.dxy_mean.co1co2.txt; tail -n +2 ../cnv_masked_results/dxy.scaffold-mi7.1kb.co1co2.cnvMask.txt | grep -P 'CO1\tCO2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -wb -a region_PLA2_scaffold-mi7.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_PLA2.bed >> pla2_intergenic.dxy_mean.co1co2.txt

echo -e "chrom\tstart\tend\tfst_mean" > pla2_intergenic.fst_mean.cv1co1.txt; tail -n +2 ../cnv_masked_results/fst.scaffold-mi7.1kb.cv1co1.cnvMask.txt | grep -v 'nan' | grep -P 'CV1\tCO1' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -wb -a region_PLA2_scaffold-mi7.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_PLA2.bed >> pla2_intergenic.fst_mean.cv1co1.txt
echo -e "chrom\tstart\tend\tfst_mean" > pla2_intergenic.fst_mean.cv1cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi7.1kb_fst.txt | grep -v 'nan' | grep -P 'CV1\tCV2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -wb -a region_PLA2_scaffold-mi7.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_PLA2.bed >> pla2_intergenic.fst_mean.cv1cv2.txt
echo -e "chrom\tstart\tend\tfst_mean" > pla2_intergenic.fst_mean.co1co2.txt; tail -n +2 ../cnv_masked_results/fst.scaffold-mi7.1kb.co1co2.cnvMask.txt | grep -v 'nan' | grep -P 'CO1\tCO2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -wb -a region_PLA2_scaffold-mi7.bed -b - | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' | bedtools intersect -v -a - -b gene_PLA2.bed >> pla2_intergenic.fst_mean.co1co2.txt
```

Genic:

```
echo -e "chrom\tstart\tend\tpi_mean" > pla2_genic.pi_mean.cv1.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi7.1kb_pi.txt | grep 'CV1' | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -a - -b gene_PLA2.bed >> pla2_genic.pi_mean.cv1.txt
echo -e "chrom\tstart\tend\tpi_mean" > pla2_genic.pi_mean.cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi7.1kb_pi.txt | grep 'CV2' | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -a - -b gene_PLA2.bed >> pla2_genic.pi_mean.cv2.txt
echo -e "chrom\tstart\tend\tpi_mean" > pla2_genic.pi_mean.co1.txt; tail -n +2 ../cnv_masked_results/pi.scaffold-mi7.1kb.co1.cnvMask.txt | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -a - -b gene_PLA2.bed >> pla2_genic.pi_mean.co1.txt
echo -e "chrom\tstart\tend\tpi_mean" > pla2_genic.pi_mean.co2.txt; tail -n +2 ../cnv_masked_results/pi.scaffold-mi7.1kb.co2.cnvMask.txt | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -a - -b gene_PLA2.bed >> pla2_genic.pi_mean.co2.txt

echo -e "chrom\tstart\tend\tdxy_mean" > pla2_genic.dxy_mean.cv1co1.txt; tail -n +2 ../cnv_masked_results/dxy.scaffold-mi7.1kb.cv1co1.cnvMask.txt | grep -P 'CV1\tCO1' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b gene_PLA2.bed >> pla2_genic.dxy_mean.cv1co1.txt
echo -e "chrom\tstart\tend\tdxy_mean" > pla2_genic.dxy_mean.cv1cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi7.1kb_dxy.txt | grep -P 'CV1\tCV2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b gene_PLA2.bed >> pla2_genic.dxy_mean.cv1cv2.txt
echo -e "chrom\tstart\tend\tdxy_mean" > pla2_genic.dxy_mean.co1co2.txt; tail -n +2 ../cnv_masked_results/dxy.scaffold-mi7.1kb.co1co2.cnvMask.txt | grep -P 'CO1\tCO2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b gene_PLA2.bed >> pla2_genic.dxy_mean.co1co2.txt

echo -e "chrom\tstart\tend\tfst_mean" > pla2_genic.fst_mean.cv1co1.txt; tail -n +2 ../cnv_masked_results/fst.scaffold-mi7.1kb.cv1co1.cnvMask.txt | grep -v 'nan' | grep -P 'CV1\tCO1' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b gene_PLA2.bed >> pla2_genic.fst_mean.cv1co1.txt
echo -e "chrom\tstart\tend\tfst_mean" > pla2_genic.fst_mean.cv1cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi7.1kb_fst.txt | grep -v 'nan' | grep -P 'CV1\tCV2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b gene_PLA2.bed >> pla2_genic.fst_mean.cv1cv2.txt
echo -e "chrom\tstart\tend\tfst_mean" > pla2_genic.fst_mean.co1co2.txt; tail -n +2 ../cnv_masked_results/fst.scaffold-mi7.1kb.co1co2.cnvMask.txt | grep -v 'nan' | grep -P 'CO1\tCO2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b gene_PLA2.bed >> pla2_genic.fst_mean.co1co2.txt
```

## Signatures of selection

Quantify a suite of population genetic estimators designed to test various predictions of neutrality versus directional and balancing selection.

The statistics focused on are:
1. Tajima's D - tests for deviations in the site-frequency spectrum from neutrality
2. Fixed differences (df) - tests for signatures of selective sweeps
3. iHS - tests for extended haplotype homozygosity surrounding target(s) of selection
4. ß - tests for allele frequency correlation surrounding balanced polymorphism.

### 1. Tajima's D

Perform sliding window analyses of Tajima's D across the genome.

#### Set up environment

```
mkdir tajima_d
cd tajima_d
mkdir results
```

This step requires sample lists per population. These are in `/resources`:

* pop.list.cv.colorado
* pop.list.cv.montana
* pop.list.co.california
* pop.list.co.idaho

#### Estimate Tajima's D in sliding windows in each population

The script below will run sliding window analysis on each population. It expects a window size (integer) and window size abbreviation (e.g., 1 kb) as arguments.

__*windowTajimaD.sh*__
```
window=$1
abbrev=$2
for chrom in `cat chrom.list`; do
	vcftools --gzvcf ../vcf/vcf_chrom-specific_cvco+outgroup/cvco+outgroup.mask.HardFilter.depth.chrom.$chrom.vcf.gz --keep pop.list.cv.colorado --TajimaD $window --out ./results/cv.colorado.$chrom.$abbrev
	vcftools --gzvcf ../vcf/vcf_chrom-specific_cvco+outgroup/cvco+outgroup.mask.HardFilter.depth.chrom.$chrom.vcf.gz --keep pop.list.cv.montana --TajimaD $window --out ./results/cv.montana.$chrom.$abbrev
	vcftools --gzvcf ../vcf/vcf_chrom-specific_cvco+outgroup/cvco+outgroup.mask.HardFilter.depth.chrom.$chrom.vcf.gz --keep pop.list.co.california --TajimaD $window --out ./results/co.california.$chrom.$abbrev
	vcftools --gzvcf ../vcf/vcf_chrom-specific_cvco+outgroup/cvco+outgroup.mask.HardFilter.depth.chrom.$chrom.vcf.gz --keep pop.list.co.idaho --TajimaD $window --out ./results/co.idaho.$chrom.$abbrev
done
```

Run the script.

```
sh window_tajima_d.sh 100000 100kb
sh window_tajima_d.sh 10000 10kb
sh window_tajima_d.sh 1000 1kb
```

#### Run higher-resolution analyses on chromosome 15

```
vcftools --gzvcf /media/drewschield/DuskBucket/crotalus/vcf/vcf_chrom-specific_cvco+outgroup/cvco+outgroup.mask.HardFilter.depth.chrom.scaffold-mi7.vcf.gz --keep pop.list.cv.colorado --TajimaD 250 --out cv.colorado.scaffold-mi7.250bp
vcftools --gzvcf /media/drewschield/DuskBucket/crotalus/vcf/vcf_chrom-specific_cvco+outgroup/cvco+outgroup.mask.HardFilter.depth.chrom.scaffold-mi7.vcf.gz --keep pop.list.cv.montana --TajimaD 250 --out cv.montana.scaffold-mi7.250bp
vcftools --gzvcf /media/drewschield/DuskBucket/crotalus/vcf/vcf_chrom-specific_cvco+outgroup/cvco+outgroup.mask.HardFilter.depth.chrom.scaffold-mi7.vcf.gz --keep pop.list.co.california --TajimaD 250 --out co.california.scaffold-mi7.250bp
vcftools --gzvcf /media/drewschield/DuskBucket/crotalus/vcf/vcf_chrom-specific_cvco+outgroup/cvco+outgroup.mask.HardFilter.depth.chrom.scaffold-mi7.vcf.gz --keep pop.list.co.idaho --TajimaD 250 --out co.idaho.scaffold-mi7.250bp
```

#### Concatenate results

```
for pop in cv.colorado cv.montana co.california co.idaho; do for window in 100kb 10kb 1kb; do head -n 1 ./results/$pop.scaffold-ma1.$window.Tajima.D > $pop.all.$window.Tajima.D; for chrom in `cat chrom.list`; do tail -n +2 ./results/$pop.$chrom.$window.Tajima.D >> $pop.all.$window.Tajima.D; done; done; done
```

### 2. Fixed differences

Regions under positive selection will have a relative abundance of fixed differences between populations. In contrast, regions subject to balancing selection may have fewer fixed differences relative to neutral expectations due to the maintenance of multiple alleles. Here, a strategy to examine fixed differences (df) isto calculate site-based Fst, then calculate the frequency of Fst = 1 SNPs in sliding windows.

#### Set up environment

```
mkdir df
cd df
mkdir results_fst
mkdir results_df
```

#### Perform per-site Fst analysis

This script will perform pairwise Fst analysis per site between populations. It uses the same population lists as in Tajima's D analysis above, found in `./resources`.

__*calcFst.sh*__
```
list1=$1
list2=$2
pop1=$3
pop2=$4
for chrom in `cat chrom.list`; do
	vcftools --gzvcf ../vcf/vcf_chrom-specific_cvco+outgroup/cvco+outgroup.mask.HardFilter.depth.chrom.$chrom.vcf.gz --min-alleles 2 --max-alleles 2 --maf 0.05 --weir-fst-pop $list1 --weir-fst-pop $list2 --out ./results_fst/fst.$chrom.$pop1.$pop2
done
```

Run the script.

```
sh calc_fst.sh pop.list.cv.colorado pop.list.cv.montana cv1 cv2
sh calc_fst.sh pop.list.cv.colorado pop.list.co.california cv1 co1
sh calc_fst.sh pop.list.co.california pop.list.co.idaho co1 co2
```

#### Format results to calculate frequency of fixed differences (df) in sliding windows

Parse all SNP positions, and positions with Fst = 1.

```
cd ./results_fst
for f in fst.*.fst; do chrom=`echo $f | cut -d. -f2`; pop1=`echo $f | cut -d. -f3`; pop2=`echo $f | cut -d. -f4`; tail -n +2 $f | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2}' > snp.$chrom.$pop1.$pop2.weir.bed; done
for f in fst.*.fst; do chrom=`echo $f | cut -d. -f2`; pop1=`echo $f | cut -d. -f3`; pop2=`echo $f | cut -d. -f4`; awk 'BEGIN{OFS="\t"}{if($3==1) print $1,$2-1,$2,$3}' $f > df.$chrom.$pop1.$pop2.weir.bed; done
cd ..
```

Use bedtools intersect to count snp positions in sliding windows. This step uses windowed BED files in `./resources/CroVir_genome_{100kb,10kb,1kb,250bp}_window.bed`.

```
cd ./results_fst
for chrom in `cat ../chrom.list`; do for pair in cv1.cv2 cv1.co1 co1.co2; do grep -w $chrom CroVir_genome_100kb_window.bed | bedtools intersect -a - -b snp.$chrom.$pair.weir.bed -c >> window.100kb.snp.$chrom.$pair.txt; done; done
for chrom in `cat ../chrom.list`; do for pair in cv1.cv2 cv1.co1 co1.co2; do grep -w $chrom CroVir_genome_10kb_window.bed | bedtools intersect -a - -b snp.$chrom.$pair.weir.bed -c >> window.10kb.snp.$chrom.$pair.txt; done; done
for chrom in `cat ../chrom.list`; do for pair in cv1.cv2 cv1.co1 co1.co2; do grep -w $chrom CroVir_genome_1kb_window.bed | bedtools intersect -a - -b snp.$chrom.$pair.weir.bed -c >> window.1kb.snp.$chrom.$pair.txt; done; done

grep -w scaffold-mi7 CroVir_genome_250bp_window.bed | bedtools intersect -a - -b snp.scaffold-mi7.cv1.co1.weir.bed -c >> window.250bp.snp.scaffold-mi7.cv1.co1.txt
grep -w scaffold-mi7 CroVir_genome_250bp_window.bed | bedtools intersect -a - -b snp.scaffold-mi7.cv1.cv2.weir.bed -c >> window.250bp.snp.scaffold-mi7.cv1.cv2.txt
grep -w scaffold-mi7 CroVir_genome_250bp_window.bed | bedtools intersect -a - -b snp.scaffold-mi7.co1.co2.weir.bed -c >> window.250bp.snp.scaffold-mi7.co1.co2.txt
```

Use bedtools intersect to get counts of fixed differences in windows.

```
for chrom in `cat chrom.list`; do for pair in cv1.cv2 cv1.co1 co1.co2; do grep -w $chrom CroVir_genome_100kb_window.bed | bedtools intersect -a - -b df.$chrom.$pair.weir.bed -c >> window.100kb.df.$chrom.$pair.txt; done; done
for chrom in `cat chrom.list`; do for pair in cv1.cv2 cv1.co1 co1.co2; do grep -w $chrom CroVir_genome_10kb_window.bed | bedtools intersect -a - -b df.$chrom.$pair.weir.bed -c >> window.10kb.df.$chrom.$pair.txt; done; done
for chrom in `cat chrom.list`; do for pair in cv1.cv2 cv1.co1 co1.co2; do grep -w $chrom CroVir_genome_1kb_window.bed | bedtools intersect -a - -b df.$chrom.$pair.weir.bed -c >> window.1kb.df.$chrom.$pair.txt; done; done

grep -w scaffold-mi7 CroVir_genome_250bp_window.bed | bedtools intersect -a - -b df.scaffold-mi7.cv1.co1.weir.bed -c >> window.250bp.df.scaffold-mi7.cv1.co1.txt
grep -w scaffold-mi7 CroVir_genome_250bp_window.bed | bedtools intersect -a - -b df.scaffold-mi7.cv1.cv2.weir.bed -c >> window.250bp.df.scaffold-mi7.cv1.cv2.txt
grep -w scaffold-mi7 CroVir_genome_250bp_window.bed | bedtools intersect -a - -b df.scaffold-mi7.co1.co2.weir.bed -c >> window.250bp.df.scaffold-mi7.co1.co2.txt
```

### 3. iHS

Use iHS, a measure of haplotype diversity ([Voight et al. 2006](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0040072)), for scans of selection. This section uses the R package rehh to calculate iHS from phased VCF input. There is a great tutorial on using rehh for this type of analysis from Joanna Meier and Mark Ravinet [here](https://speciationgenomics.github.io/haplotypes/).

#### Set up environment

```
mkdir rehh
cd rehh
mkdir data
mkdir chrom_results
mkdir significance
```

#### Perform haplotype diversity analysis in rehh R package

The details of this analysis are in `./R/rehhCalculations.R`. The phased VCF input data can be retrieved [here](https://figshare.com/articles/dataset/Phased_VCFs/16415556).

#### Format rehh output

Convert iHS tables to BED format.

```
cd ./chrom_results
for i in *.txt; do file=`echo $i | sed 's/.txt/.bed/g'`; awk 'BEGIN{OFS="\t"} NR>1 {print $1,$2-1,$2,$3,$4}' $i > $file; done
```

Concatenate chromosome-specific results.

```
cat cv.ma1_ihs.bed cv.ma2_ihs.bed cv.ma3_ihs.bed cv.ma4_ihs.bed cv.ma5_ihs.bed cv.ma6_ihs.bed cv.ma7_ihs.bed cv.Z_ihs.bed cv.mi1_ihs.bed cv.mi2_ihs.bed cv.mi3_ihs.bed cv.mi4_ihs.bed cv.mi5_ihs.bed cv.mi6_ihs.bed cv.mi7_ihs.bed cv.mi8_ihs.bed cv.mi9_ihs.bed cv.mi10_ihs.bed > cv.all_ihs.bed
cat co.ma1_ihs.bed co.ma2_ihs.bed co.ma3_ihs.bed co.ma4_ihs.bed co.ma5_ihs.bed co.ma6_ihs.bed co.ma7_ihs.bed co.Z_ihs.bed co.mi1_ihs.bed co.mi2_ihs.bed co.mi3_ihs.bed co.mi4_ihs.bed co.mi5_ihs.bed co.mi6_ihs.bed co.mi7_ihs.bed co.mi8_ihs.bed co.mi9_ihs.bed co.mi10_ihs.bed > co.all_ihs.bed
```

Use bedtools map to calculate mean iHS in sliding windows.

```
for window in 100kb 10kb 1kb; do for pop in cv co; do echo -e "chrom\tstart\tend\tiHS\tp-value" > ../$pop.all_ihs.$window.txt | bedtools map -a /data3/venom_population_genomics/general/CroVir_genome_${window}_window.bed -b $pop.all_ihs.bed -c 4,5 -o mean >> ../$pop.all_ihs.$window.txt; done; done
echo -e "chrom\tstart\tend\tiHS\tp-value" > ../co.mi7_ihs.250bp.txt | grep -w 'scaffold-mi7' co.all_ihs.bed | bedtools map -a /data3/venom_population_genomics/general/CroVir_genome_250bp_window.bed -b - -c 4,5 -o mean >> ../co.mi7_ihs.250bp.txt
echo -e "chrom\tstart\tend\tiHS\tp-value" > ../cv.mi7_ihs.250bp.txt | grep -w 'scaffold-mi7' cv.all_ihs.bed | bedtools map -a /data3/venom_population_genomics/general/CroVir_genome_250bp_window.bed -b - -c 4,5 -o mean >> ../cv.mi7_ihs.250bp.txt
cd ..
```

#### Significance testing

Use random resampling of genome-wide iHS values to compare with mean values for venom gene regions and to quantify how often the same or greater values are observed.

Extract 10 kb windows for SVMP and SVSP regions, 1 kb windows for PLA2 region. This uses region BED files in `./resources`.

```
cd significance
tail -n +2 ../cv.all_ihs.10kb.txt | bedtools intersect -wa -a - -b region_SVMP_scaffold-mi1.bed > region_SVMP_scaffold-mi1.10kb.bed
tail -n +2 ../cv.all_ihs.10kb.txt | bedtools intersect -wa -a - -b region_SVSP_scaffold-mi2.bed > region_SVSP_scaffold-mi2.10kb.bed
tail -n +2 ../cv.all_ihs.1kb.txt | bedtools intersect -wa -a - -b region_PLA2_scaffold-mi7.bed > region_PLA2_scaffold-mi7.1kb.bed
```

Perform blocked permutations using the script below. These can be compared to observed means in venom gene regions.

__*permutationsIHS.sh*__
```
vbed=$1
data=$2
lines=`wc -l $vbed | cut -d' ' -f 1`
lines_fix=`echo "$(($lines-1))"`
for perm in $(seq 1 10000); do
	mean=`tail -n +2 $data | shuf -n 1 | grep -f - -A $lines_fix $data | grep -v -P "\.\t\." | awk '{print $4}' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }'`
	echo "${perm}\t${mean}"
done
```

Run the script.

```
echo -e "permutation\tiHS" > cv1.svmp.iHS.10kb_permutations.txt; sh permutationsIHS.sh region_SVMP_scaffold-mi1.10kb.bed ../cv.all_ihs.10kb.txt >> cv1.svmp.iHS.10kb_permutations.txt
echo -e "permutation\tiHS" > co1.svmp.iHS.10kb_permutations.txt; sh permutationsIHS.sh region_SVMP_scaffold-mi1.10kb.bed ../co.all_ihs.10kb.txt >> co1.svmp.iHS.10kb_permutations.txt
echo -e "permutation\tiHS" > cv1.svsp.iHS.10kb_permutations.txt; sh permutationsIHS.sh region_SVSP_scaffold-mi2.10kb.bed ../cv.all_ihs.10kb.txt >> cv1.svsp.iHS.10kb_permutations.txt
echo -e "permutation\tiHS" > co1.svsp.iHS.10kb_permutations.txt; sh permutationsIHS.sh region_SVSP_scaffold-mi2.10kb.bed ../co.all_ihs.10kb.txt >> co1.svsp.iHS.10kb_permutations.txt
echo -e "permutation\tiHS" > cv1.pla2.iHS.1kb_permutations.txt; sh permutationsIHS.sh region_PLA2_scaffold-mi7.1kb.bed ../cv.all_ihs.1kb.txt >> cv1.pla2.iHS.1kb_permutations.txt
echo -e "permutation\tiHS" > co1.pla2.iHS.1kb_permutations.txt; sh permutationsIHS.sh region_PLA2_scaffold-mi7.1kb.bed ../co.all_ihs.1kb.txt >> co1.pla2.iHS.1kb_permutations.txt
```

Calculate p-values by querying the proportion of values that exceed venom region means.

Mean iHS for CV1 and CO1 populations:

| Population | SVMP      | SVSP      | PLA2      |
|------------|-----------|-----------|-----------|
| CV1        | 0.977623  | 0.8214266 | 0.7144463 |
| CO1        | 0.7008945 | 0.7098776 | 0.3739584 |

```
$tail -n +2 cv1.svmp.iHS.10kb_permutations.txt | awk '$2>0.977623' | wc -l
$tail -n +2 co1.svmp.iHS.10kb_permutations.txt | awk '$2>0.7008945' | wc -l
$tail -n +2 cv1.svsp.iHS.10kb_permutations.txt | awk '$2>0.8214266' | wc -l
$tail -n +2 co1.svsp.iHS.10kb_permutations.txt | awk '$2>0.7098776' | wc -l
$tail -n +2 cv1.pla2.iHS.1kb_permutations.txt | awk '$2>0.7144463' | wc -l
$tail -n +2 co1.pla2.iHS.1kb_permutations.txt | awk '$2>0.3739584' | wc -l
```

### 4. ß

Balancing selection can produce clusters of intermediate-frequency alleles surrounding a balanced polymorphism. [BetaScan](https://github.com/ksiewert/BetaScan) is designed to identify this signature in a genome scan approach. The method is described in [Siewert and Voight 2017](https://academic.oup.com/mbe/article/34/11/2996/3988103?login=true).

#### Set up environment

```
mkdir beta
cd beta
mkdir data
mkdir input
mkdir results
mkdir parameter_permutations
mkdir significance
```

*Note: BetaScan uses ACF format, which can be converted from VCF format using [glactools](https://github.com/grenaud/glactools), so make sure it's installed.*

#### Convert VCF data to ACF format

These analyses will use the same phased variants as in iHS analysis above, which can be downloaded [here](https://figshare.com/articles/dataset/Phased_VCFs/16415556).

Perform format conversion with glactools for each scaffold in `./resources/chrom.list`.

```
for chrom in `cat chrom.list`; do glactools vcfm2acf --onlyGT --fai ../CroVir_genome_L77pg_16Aug2017.final_rename.fasta.fai ./data/viridis.phased.$chrom.vcf > ./input/cv1.phased.$chrom.acf.gz; done
for chrom in `cat chrom.list`; do glactools vcfm2acf --onlyGT --fai ../CroVir_genome_L77pg_16Aug2017.final_rename.fasta.fai ./data/oreganus.phased.$chrom.vcf > ./input/co1.phased.$chrom.acf.gz; done
```

#### Convert ACF to BetaScan input format

```
for pop in cv1 co1; do for chrom in `cat chrom.list`; do glactools acf2betascan --fold ./input/$pop.phased.$chrom.acf.gz | gzip > ./input/$pop.phased.$chrom.beta.txt.gz; done; done
```

#### Run BetaScan

```
for pop in cv1 co1; do for chrom in `cat chrom.list`; do python BetaScan.py -i ./input/$pop.phased.$chrom.beta.txt.gz -fold -o ./results/$pop.phased.$chrom.betascores.txt; done; done
```

#### Convert results to BED format, concatenate, and calculate mean ß in sliding windows

```
cd ./results
for beta in *.betascores.txt; do file=`echo $beta | sed 's/.txt/.bed/g'`; chrom=`echo $beta | cut -d. -f3`; tail -n +2 $beta | awk -v chromvar="$chrom" 'BEGIN{OFS="\t"}{print chromvar,$1-1,$1,$2}' > $file; done

cat cv1.phased.scaffold-ma1.betascores.bed cv1.phased.scaffold-ma2.betascores.bed cv1.phased.scaffold-ma3.betascores.bed cv1.phased.scaffold-ma4.betascores.bed cv1.phased.scaffold-ma5.betascores.bed cv1.phased.scaffold-ma6.betascores.bed cv1.phased.scaffold-ma7.betascores.bed cv1.phased.scaffold-Z.betascores.bed cv1.phased.scaffold-mi1.betascores.bed cv1.phased.scaffold-mi2.betascores.bed cv1.phased.scaffold-mi3.betascores.bed cv1.phased.scaffold-mi4.betascores.bed cv1.phased.scaffold-mi5.betascores.bed cv1.phased.scaffold-mi6.betascores.bed cv1.phased.scaffold-mi7.betascores.bed cv1.phased.scaffold-mi8.betascores.bed cv1.phased.scaffold-mi9.betascores.bed cv1.phased.scaffold-mi10.betascores.bed > cv1.phased.all.betascores.bed
cat co1.phased.scaffold-ma1.betascores.bed co1.phased.scaffold-ma2.betascores.bed co1.phased.scaffold-ma3.betascores.bed co1.phased.scaffold-ma4.betascores.bed co1.phased.scaffold-ma5.betascores.bed co1.phased.scaffold-ma6.betascores.bed co1.phased.scaffold-ma7.betascores.bed co1.phased.scaffold-Z.betascores.bed co1.phased.scaffold-mi1.betascores.bed co1.phased.scaffold-mi2.betascores.bed co1.phased.scaffold-mi3.betascores.bed co1.phased.scaffold-mi4.betascores.bed co1.phased.scaffold-mi5.betascores.bed co1.phased.scaffold-mi6.betascores.bed co1.phased.scaffold-mi7.betascores.bed co1.phased.scaffold-mi8.betascores.bed co1.phased.scaffold-mi9.betascores.bed co1.phased.scaffold-mi10.betascores.bed > co1.phased.all.betascores.bed

echo -e "chrom\tstart\tend\tBeta1*" > cv1.phased.all.betascores.100kb.txt; bedtools map -a /data3/venom_population_genomics/general/CroVir_genome_100kb_window.bed -b cv1.phased.all.betascores.bed -c 4 -o mean >> cv1.phased.all.betascores.100kb.txt
echo -e "chrom\tstart\tend\tBeta1*" > cv1.phased.all.betascores.10kb.txt; bedtools map -a /data3/venom_population_genomics/general/CroVir_genome_10kb_window.bed -b cv1.phased.all.betascores.bed -c 4 -o mean >> cv1.phased.all.betascores.10kb.txt
echo -e "chrom\tstart\tend\tBeta1*" > cv1.phased.all.betascores.1kb.txt; bedtools map -a /data3/venom_population_genomics/general/CroVir_genome_1kb_window.bed -b cv1.phased.all.betascores.bed -c 4 -o mean >> cv1.phased.all.betascores.1kb.txt
echo -e "chrom\tstart\tend\tBeta1*" > co1.phased.all.betascores.100kb.txt; bedtools map -a /data3/venom_population_genomics/general/CroVir_genome_100kb_window.bed -b co1.phased.all.betascores.bed -c 4 -o mean >> co1.phased.all.betascores.100kb.txt
echo -e "chrom\tstart\tend\tBeta1*" > co1.phased.all.betascores.10kb.txt; bedtools map -a /data3/venom_population_genomics/general/CroVir_genome_10kb_window.bed -b co1.phased.all.betascores.bed -c 4 -o mean >> co1.phased.all.betascores.10kb.txt
echo -e "chrom\tstart\tend\tBeta1*" > co1.phased.all.betascores.1kb.txt; bedtools map -a /data3/venom_population_genomics/general/CroVir_genome_1kb_window.bed -b co1.phased.all.betascores.bed -c 4 -o mean >> co1.phased.all.betascores.1kb.txt
echo -e "chrom\tstart\tend\tBeta1*" > cv1.phased.all.betascores.250bp.txt; grep -w 'scaffold-mi7' /data3/venom_population_genomics/general/CroVir_genome_250bp_window.bed | bedtools map -a - -b cv1.phased.scaffold-mi7.betascores.bed -c 4 -o mean >> cv1.phased.scaffold-mi7.betascores.250bp.txt
echo -e "chrom\tstart\tend\tBeta1*" > co1.phased.all.betascores.250bp.txt; grep -w 'scaffold-mi7' /data3/venom_population_genomics/general/CroVir_genome_250bp_window.bed | bedtools map -a - -b co1.phased.scaffold-mi7.betascores.bed -c 4 -o mean >> co1.phased.scaffold-mi7.betascores.250bp.txt

cd ..
```

#### Perform analysis with permuted settings

Run BetaScan with a series of different window size (-w) and scaling constant (-p) parameters.

```
for p in 2 5 10 20; do for w in 500 1000 2000; do python BetaScan.py -i ./input/cv1.phased.scaffold-mi1.beta.txt.gz -fold -p $p -w $w -o ./permutations/cv1.phased.scaffold-mi1.p$p.w$w.betascores.txt; done; done
```

#### Significance testing

BetaScan does not output a measure of significance by default. Use random resampling of the genome-wide dataset compared to values observed in venom gene regions to quantify how often the same or greater values are observed.

Extract 10 kb windows for SVMP and SVSP regions, 1 kb windows for PLA2.

```
cd ./significance
tail -n +2 ../results/cv1.phased.all.betascores.10kb.txt | bedtools intersect -wa -a - -b /data3/venom_population_genomics/general/region_SVMP_scaffold-mi1.bed > region_SVMP_scaffold-mi1.10kb.bed
tail -n +2 ../results/cv1.phased.all.betascores.10kb.txt | bedtools intersect -wa -a - -b /data3/venom_population_genomics/general/region_SVSP_scaffold-mi2.bed > region_SVSP_scaffold-mi2.10kb.bed
tail -n +2 ../results/cv1.phased.all.betascores.1kb.txt | bedtools intersect -wa -a - -b /data3/venom_population_genomics/general/region_PLA2_scaffold-mi7.bed > region_PLA2_scaffold-mi7.1kb.bed
```

Perform permutations that can be compared to observed means in venom regions with this script.

__*permutationsBeta.sh*__
```
vbed=$1
data=$2
lines=`wc -l $vbed | cut -d' ' -f 1`
lines_fix=`echo "$(($lines-1))"`
for perm in $(seq 1 10000); do
	mean=`tail -n +2 $data | shuf -n 1 | grep -f - -A $lines_fix $data | awk '{print $4}' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }'`
	echo "${perm}\t${mean}"
done
```

Run the script.

```
echo -e "permutation\tbeta" > cv1.svmp.beta.10kb_permutations.txt; sh permutations.sh region_SVMP_scaffold-mi1.10kb.bed ../results/cv1.phased.all.betascores.10kb.txt >> cv1.svmp.beta.10kb_permutations.txt
echo -e "permutation\tbeta" > co1.svmp.beta.10kb_permutations.txt; sh permutations.sh region_SVMP_scaffold-mi1.10kb.bed ../results/co1.phased.all.betascores.10kb.txt >> co1.svmp.beta.10kb_permutations.txt
echo -e "permutation\tbeta" > cv1.svsp.beta.10kb_permutations.txt; sh permutations.sh region_SVSP_scaffold-mi2.10kb.bed ../results/cv1.phased.all.betascores.10kb.txt >> cv1.svsp.beta.10kb_permutations.txt
echo -e "permutation\tbeta" > co1.svsp.beta.10kb_permutations.txt; sh permutations.sh region_SVSP_scaffold-mi2.10kb.bed ../results/co1.phased.all.betascores.10kb.txt >> co1.svsp.beta.10kb_permutations.txt
echo -e "permutation\tbeta" > cv1.pla2.beta.1kb_permutations.txt; sh permutations.sh region_PLA2_scaffold-mi7.1kb.bed ../results/cv1.phased.all.betascores.1kb.txt >> cv1.pla2.beta.1kb_permutations.txt
echo -e "permutation\tbeta" > co1.pla2.beta.1kb_permutations.txt; sh permutations.sh region_PLA2_scaffold-mi7.1kb.bed ../results/co1.phased.all.betascores.1kb.txt >> co1.pla2.beta.1kb_permutations.txt
```

Calculate p-values by querying the proportion of values that exceed venom region means.

Mean Beta for CV1 and CO1 populations:

| Population | SVMP      | SVSP      | PLA2      |
|------------|-----------|-----------|-----------|
| CV1        | 2.292093  | 2.073221  | 1.084576  |
| CO1        | 2.137636  | 1.425567  | 1.605109  |


```
tail -n +2 cv1.svmp.beta.10kb_permutations.txt | awk '$2>2.292093' | wc -l
tail -n +2 co1.svmp.beta.10kb_permutations.txt | awk '$2>2.137636' | wc -l
tail -n +2 cv1.svsp.beta.10kb_permutations.txt | awk '$2>2.073221' | wc -l
tail -n +2 co1.svsp.beta.10kb_permutations.txt | awk '$2>1.425567' | wc -l
tail -n +2 cv1.pla2.beta.1kb_permutations.txt | awk '$2>1.084576' | wc -l
tail -n +2 co1.pla2.beta.1kb_permutations.txt | awk '$2>1.605109' | wc -l
```

### Selection appendix 1: CNV-masking

Use `cnvMaskSelection.py` in the `python` directory to mask results overlapping with copy-number variants in CO1 and CO2. See [Diversity appendix 1: CNV-masking](#diversity-appendix-1-cnv-masking) for more details on the script.

#### Set up environment

```
cd cnv_masked_results
```

Make sure that the appropriate BED files from `resources` are placed in the `mask_bed` subdirectory.

#### Mask Tajima's D results

```
for pop in co.california co.idaho; do head -n 1 ../tajima_d/$pop.all.10kb.Tajima.D > $pop.all.10kb.cnvMask.Tajima.D; python cnvMaskSelection.py ../tajima_d/$pop.all.10kb.Tajima.D ./mask_bed/cnv_filter.10kb.bed no >> $pop.all.10kb.cnvMask.Tajima.D; done
for pop in co.california co.idaho; do head -n 1 ../tajima_d/$pop.all.1kb.Tajima.D > $pop.all.1kb.cnvMask.Tajima.D; python cnvMaskSelection.py ../tajima_d/$pop.all.1kb.Tajima.D ./mask_bed/cnv_filter.1kb.bed no >> $pop.all.1kb.cnvMask.Tajima.D; done
for pop in co.california co.idaho; do head -n 1 ../tajima_d/$pop.scaffold-mi7.250bp.Tajima.D > $pop.scaffold-mi7.250bp.cnvMask.Tajima.D; python cnvMaskSelection.py ../tajima_d/$pop.scaffold-mi7.250bp.Tajima.D ./mask_bed/cnv_filter.250bp.bed no >> $pop.scaffold-mi7.250bp.cnvMask.Tajima.D; done
```

#### Mask df results

```
for pop in cv1.co1 co1.co2; do head -n 1 ../df/results_df/window.10kb.df_prop.all.$pop.txt > window.10kb.df_prop.all.$pop.cnvMask.txt; python cnvMaskSelection.py ../df/results_df/window.10kb.df_prop.all.$pop.txt ./mask_bed/cnv_filter.10kb.bed no >> window.10kb.df_prop.all.$pop.cnvMask.txt; done
for pop in cv1.co1 co1.co2; do head -n 1 ../df/results_df/window.1kb.df_prop.all.$pop.txt > window.1kb.df_prop.all.$pop.cnvMask.txt; python cnvMaskSelection.py ../df/results_df/window.1kb.df_prop.all.$pop.txt ./mask_bed/cnv_filter.1kb.bed no >> window.1kb.df_prop.all.$pop.cnvMask.txt; done
for pop in cv1.co1 co1.co2; do head -n 1 ../df/results_df/window.250bp.df_prop.scaffold-mi7.$pop.txt > window.250bp.df_prop.scaffold-mi7.$pop.cnvMask.txt; python cnvMaskSelection.py ../df/results_df/window.250bp.df_prop.scaffold-mi7.$pop.txt ./mask_bed/cnv_filter.250bp.bed no >> window.250bp.df_prop.scaffold-mi7.$pop.cnvMask.txt; done
```

#### Mask iHS results

```
$head -n 1 ../rehh/co.all_ihs.10kb.txt > co.all_ihs.10kb.cnvMask.txt; python cnvMaskSelection.py ../rehh/co.all_ihs.10kb.txt ./mask_bed/cnv_filter.10kb.bed no >> co.all_ihs.10kb.cnvMask.txt
$head -n 1 ../rehh/co.all_ihs.1kb.txt > co.all_ihs.1kb.cnvMask.txt; python cnvMaskSelection.py ../rehh/co.all_ihs.1kb.txt ./mask_bed/cnv_filter.1kb.bed no >> co.all_ihs.1kb.cnvMask.txt
$head -n 1 ../rehh/co.mi7_ihs.250bp.txt > co.mi7_ihs.250bp.cnvMask.txt; python cnvMaskSelection.py ../rehh/co.all_ihs.1kb.txt ./mask_bed/cnv_filter.250bp.bed no >> co.mi7_ihs.250bp.cnvMask.txt
```

#### Mask ß results

```
head -n 1 /data3/venom_population_genomics/3_selection/beta/results/co1.phased.all.betascores.10kb.txt > co1.phased.all.betascores.10kb.cnvMask.txt; python cnvMaskSelection.py /data3/venom_population_genomics/3_selection/beta/results/co1.phased.all.betascores.10kb.txt ./mask_bed/cnv_filter.10kb.bed no >> co1.phased.all.betascores.10kb.cnvMask.txt
head -n 1 /data3/venom_population_genomics/3_selection/beta/results/co1.phased.all.betascores.1kb.txt > co1.phased.all.betascores.1kb.cnvMask.txt; python cnvMaskSelection.py /data3/venom_population_genomics/3_selection/beta/results/co1.phased.all.betascores.1kb.txt ./mask_bed/cnv_filter.1kb.bed no >> co1.phased.all.betascores.1kb.cnvMask.txt
echo -e "chrom\tstart\tend\tBeta1*" > co1.phased.scaffold-mi7.betascores.250bp.cnvMask.txt; python cnvMaskSelection.py /data3/venom_population_genomics/3_selection/beta/results/co1.phased.scaffold-mi7.betascores.250bp.txt ./mask_bed/cnv_filter.250bp.bed no >> co1.phased.scaffold-mi7.betascores.250bp.cnvMask.txt
```

### Selection appendix 2: Venom gene point estimates

It will be useful to quantify gene-specific mean values for population genetic estimators for the major venom genes.

#### Set up environment

```
mkdir gene_means
cd gene_means
```

These analyses will use the coordinates in `resources/gene_venom.bed`.

#### Calculate mean π per venom gene

```
echo -e "chrom\tstart\tend\tpi_mean" > svmp.pi_mean.cv1.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi1.1kb_pi.txt | grep 'CV1' | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools map -a gene_SVMP.bed -b - -c 4 -o mean >> svmp.pi_mean.cv1.txt
echo -e "chrom\tstart\tend\tpi_mean" > svmp.pi_mean.cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi1.1kb_pi.txt | grep 'CV2' | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools map -a gene_SVMP.bed -b - -c 4 -o mean >> svmp.pi_mean.cv2.txt
echo -e "chrom\tstart\tend\tpi_mean" > svmp.pi_mean.co1.txt; tail -n +2 ../cnv_masked_results/pi.scaffold-mi1.1kb.co1.cnvMask.txt | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools map -a gene_SVMP.bed -b - -c 4 -o mean >> svmp.pi_mean.co1.txt
echo -e "chrom\tstart\tend\tpi_mean" > svmp.pi_mean.co2.txt; tail -n +2 ../cnv_masked_results/pi.scaffold-mi1.1kb.co2.cnvMask.txt | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools map -a gene_SVMP.bed -b - -c 4 -o mean >> svmp.pi_mean.co2.txt

echo -e "chrom\tstart\tend\tpi_mean" > svsp.pi_mean.cv1.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi2.1kb_pi.txt | grep 'CV1' | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools map -a gene_SVSP.bed -b - -c 4 -o mean >> svsp.pi_mean.cv1.txt
echo -e "chrom\tstart\tend\tpi_mean" > svsp.pi_mean.cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi2.1kb_pi.txt | grep 'CV2' | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools map -a gene_SVSP.bed -b - -c 4 -o mean >> svsp.pi_mean.cv2.txt
echo -e "chrom\tstart\tend\tpi_mean" > svsp.pi_mean.co1.txt; tail -n +2 ../cnv_masked_results/pi.scaffold-mi2.1kb.co1.cnvMask.txt | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools map -a gene_SVSP.bed -b - -c 4 -o mean >> svsp.pi_mean.co1.txt
echo -e "chrom\tstart\tend\tpi_mean" > svsp.pi_mean.co2.txt; tail -n +2 ../cnv_masked_results/pi.scaffold-mi2.1kb.co2.cnvMask.txt | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools map -a gene_SVSP.bed -b - -c 4 -o mean >> svsp.pi_mean.co2.txt

echo -e "chrom\tstart\tend\tpi_mean" > pla2.pi_mean.cv1.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi7.1kb_pi.txt | grep 'CV1' | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools map -a gene_PLA2.bed -b - -c 4 -o mean >> pla2.pi_mean.cv1.txt
echo -e "chrom\tstart\tend\tpi_mean" > pla2.pi_mean.cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi7.1kb_pi.txt | grep 'CV2' | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools map -a gene_PLA2.bed -b - -c 4 -o mean >> pla2.pi_mean.cv2.txt
echo -e "chrom\tstart\tend\tpi_mean" > pla2.pi_mean.co1.txt; tail -n +2 ../cnv_masked_results/pi.scaffold-mi7.1kb.co1.cnvMask.txt | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools map -a gene_PLA2.bed -b - -c 4 -o mean >> pla2.pi_mean.co1.txt
echo -e "chrom\tstart\tend\tpi_mean" > pla2.pi_mean.co2.txt; tail -n +2 ../cnv_masked_results/pi.scaffold-mi7.1kb.co2.cnvMask.txt | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools map -a gene_PLA2.bed -b - -c 4 -o mean >> pla2.pi_mean.co2.txt
```

#### Calculate mean dxy per venom gene

```
echo -e "chrom\tstart\tend\tdxy_mean" > svmp.dxy_mean.cv1co1.txt; tail -n +2 ../cnv_masked_results/dxy.scaffold-mi1.1kb.cv1co1.cnvMask.txt | grep -P 'CV1\tCO1' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools map -a gene_SVMP.bed -b - -c 4 -o mean >> svmp.dxy_mean.cv1co1.txt
echo -e "chrom\tstart\tend\tdxy_mean" > svmp.dxy_mean.cv1cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi1.1kb_dxy.txt | grep -P 'CV1\tCV2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools map -a gene_SVMP.bed -b - -c 4 -o mean >> svmp.dxy_mean.cv1cv2.txt
echo -e "chrom\tstart\tend\tdxy_mean" > svmp.dxy_mean.co1co2.txt; tail -n +2 ../cnv_masked_results/dxy.scaffold-mi1.1kb.co1co2.cnvMask.txt | grep -P 'CO1\tCO2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools map -a gene_SVMP.bed -b - -c 4 -o mean >> svmp.dxy_mean.co1co2.txt

echo -e "chrom\tstart\tend\tdxy_mean" > svsp.dxy_mean.cv1co1.txt; tail -n +2 ../cnv_masked_results/dxy.scaffold-mi2.1kb.cv1co1.cnvMask.txt | grep -P 'CV1\tCO1' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools map -a gene_SVSP.bed -b - -c 4 -o mean >> svsp.dxy_mean.cv1co1.txt
echo -e "chrom\tstart\tend\tdxy_mean" > svsp.dxy_mean.cv1cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi2.1kb_dxy.txt | grep -P 'CV1\tCV2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools map -a gene_SVSP.bed -b - -c 4 -o mean >> svsp.dxy_mean.cv1cv2.txt
echo -e "chrom\tstart\tend\tdxy_mean" > svsp.dxy_mean.co1co2.txt; tail -n +2 ../cnv_masked_results/dxy.scaffold-mi2.1kb.co1co2.cnvMask.txt | grep -P 'CO1\tCO2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools map -a gene_SVSP.bed -b - -c 4 -o mean >> svsp.dxy_mean.co1co2.txt

echo -e "chrom\tstart\tend\tdxy_mean" > pla2.dxy_mean.cv1co1.txt; tail -n +2 ../cnv_masked_results/dxy.scaffold-mi7.1kb.cv1co1.cnvMask.txt | grep -P 'CV1\tCO1' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools map -a gene_PLA2.bed -b - -c 4 -o mean >> pla2.dxy_mean.cv1co1.txt
echo -e "chrom\tstart\tend\tdxy_mean" > pla2.dxy_mean.cv1cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi7.1kb_dxy.txt | grep -P 'CV1\tCV2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools map -a gene_PLA2.bed -b - -c 4 -o mean >> pla2.dxy_mean.cv1cv2.txt
echo -e "chrom\tstart\tend\tdxy_mean" > pla2.dxy_mean.co1co2.txt; tail -n +2 ../cnv_masked_results/dxy.scaffold-mi7.1kb.co1co2.cnvMask.txt | grep -P 'CO1\tCO2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools map -a gene_PLA2.bed -b - -c 4 -o mean >> pla2.dxy_mean.co1co2.txt
```

#### Calculate mean Fst per venom gene

```
echo -e "chrom\tstart\tend\tfst_mean" > svmp.fst_mean.cv1co1.txt; tail -n +2 ../cnv_masked_results/fst.scaffold-mi1.1kb.cv1co1.cnvMask.txt | grep -v 'nan' | grep -P 'CV1\tCO1' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools map -a gene_SVMP.bed -b - -c 4 -o mean >> svmp.fst_mean.cv1co1.txt
echo -e "chrom\tstart\tend\tfst_mean" > svmp.fst_mean.cv1cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi1.1kb_fst.txt | grep -v 'nan' | grep -P 'CV1\tCV2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools map -a gene_SVMP.bed -b - -c 4 -o mean >> svmp.fst_mean.cv1cv2.txt
echo -e "chrom\tstart\tend\tfst_mean" > svmp.fst_mean.co1co2.txt; tail -n +2 ../cnv_masked_results/fst.scaffold-mi1.1kb.co1co2.cnvMask.txt | grep -v 'nan' | grep -P 'CO1\tCO2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools map -a gene_SVMP.bed -b - -c 4 -o mean >> svmp.fst_mean.co1co2.txt

echo -e "chrom\tstart\tend\tfst_mean" > svsp.fst_mean.cv1co1.txt; tail -n +2 ../cnv_masked_results/fst.scaffold-mi2.1kb.cv1co1.cnvMask.txt | grep -v 'nan' | grep -P 'CV1\tCO1' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools map -a gene_SVSP.bed -b - -c 4 -o mean >> svsp.fst_mean.cv1co1.txt
echo -e "chrom\tstart\tend\tfst_mean" > svsp.fst_mean.cv1cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi2.1kb_fst.txt | grep -v 'nan' | grep -P 'CV1\tCV2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools map -a gene_SVSP.bed -b - -c 4 -o mean >> svsp.fst_mean.cv1cv2.txt
echo -e "chrom\tstart\tend\tfst_mean" > svsp.fst_mean.co1co2.txt; tail -n +2 ../cnv_masked_results/fst.scaffold-mi2.1kb.co1co2.cnvMask.txt | grep -v 'nan' | grep -P 'CO1\tCO2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools map -a gene_SVSP.bed -b - -c 4 -o mean >> svsp.fst_mean.co1co2.txt

echo -e "chrom\tstart\tend\tfst_mean" > pla2.fst_mean.cv1co1.txt; tail -n +2 ../cnv_masked_results/fst.scaffold-mi7.1kb.cv1co1.cnvMask.txt | grep -v 'nan' | grep -P 'CV1\tCO1' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools map -a gene_PLA2.bed -b - -c 4 -o mean >> pla2.fst_mean.cv1co1.txt
echo -e "chrom\tstart\tend\tfst_mean" > pla2.fst_mean.cv1cv2.txt; tail -n +2 ../pixy/pixy_results/pixy.scaffold-mi7.1kb_fst.txt | grep -v 'nan' | grep -P 'CV1\tCV2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools map -a gene_PLA2.bed -b - -c 4 -o mean >> pla2.fst_mean.cv1cv2.txt
echo -e "chrom\tstart\tend\tfst_mean" > pla2.fst_mean.co1co2.txt; tail -n +2 ../cnv_masked_results/fst.scaffold-mi7.1kb.co1co2.cnvMask.txt | grep -v 'nan' | grep -P 'CO1\tCO2' | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools map -a gene_PLA2.bed -b - -c 4 -o mean >> pla2.fst_mean.co1co2.txt
```

#### Calculate mean Tajima's D per venom gene

```
echo -e "chrom\tstart\tend\tTajimaD" > svmp.TajimaD_mean.cv1.txt; tail -n +2 ../tajima_d/cv.colorado.all.1kb.Tajima.D | grep -v 'nan' | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+1000,$4}' | bedtools sort -i - | bedtools map -a gene_SVMP.bed -b - -c 4 -o mean >> svmp.TajimaD_mean.cv1.txt
echo -e "chrom\tstart\tend\tTajimaD" > svsp.TajimaD_mean.cv1.txt; tail -n +2 ../tajima_d/cv.colorado.all.1kb.Tajima.D | grep -v 'nan' | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+1000,$4}' | bedtools sort -i - | bedtools map -a gene_SVSP.bed -b - -c 4 -o mean >> svsp.TajimaD_mean.cv1.txt
echo -e "chrom\tstart\tend\tTajimaD" > pla2.TajimaD_mean.cv1.txt; tail -n +2 ../tajima_d/cv.colorado.scaffold-mi7.250bp.Tajima.D | grep -v 'nan' | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+250,$4}' | bedtools sort -i - | bedtools map -a gene_PLA2.bed -b - -c 4 -o mean >> pla2.TajimaD_mean.cv1.txt

echo -e "chrom\tstart\tend\tTajimaD" > svmp.TajimaD_mean.co1.txt; tail -n +2 ../cnv_masked_results/co.california.all.1kb.cnvMask.Tajima.D | grep -v 'nan' | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+1000,$4}' | bedtools sort -i - | bedtools map -a gene_SVMP.bed -b - -c 4 -o mean >> svmp.TajimaD_mean.co1.txt
echo -e "chrom\tstart\tend\tTajimaD" > svsp.TajimaD_mean.co1.txt; tail -n +2 ../cnv_masked_results/co.california.all.1kb.cnvMask.Tajima.D | grep -v 'nan' | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+1000,$4}' | bedtools sort -i - | bedtools map -a gene_SVSP.bed -b - -c 4 -o mean >> svsp.TajimaD_mean.co1.txt
echo -e "chrom\tstart\tend\tTajimaD" > pla2.TajimaD_mean.co1.txt; tail -n +2 ../cnv_masked_results/co.california.scaffold-mi7.250bp.cnvMask.Tajima.D | grep -v 'nan' | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+250,$4}' | bedtools sort -i - | bedtools map -a gene_PLA2.bed -b - -c 4 -o mean >> pla2.TajimaD_mean.co1.txt

echo -e "chrom\tstart\tend\tTajimaD" > svmp.TajimaD_mean.cv2.txt; tail -n +2 ../tajima_d/cv.montana.all.1kb.Tajima.D | grep -v 'nan' | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+1000,$4}' | bedtools sort -i - | bedtools map -a gene_SVMP.bed -b - -c 4 -o mean >> svmp.TajimaD_mean.cv2.txt
echo -e "chrom\tstart\tend\tTajimaD" > svsp.TajimaD_mean.cv2.txt; tail -n +2 ../tajima_d/cv.montana.all.1kb.Tajima.D | grep -v 'nan' | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+1000,$4}' | bedtools sort -i - | bedtools map -a gene_SVSP.bed -b - -c 4 -o mean >> svsp.TajimaD_mean.cv2.txt
echo -e "chrom\tstart\tend\tTajimaD" > pla2.TajimaD_mean.cv2.txt; tail -n +2 ../tajima_d/cv.montana.scaffold-mi7.250bp.Tajima.D | grep -v 'nan' | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+250,$4}' | bedtools sort -i - | bedtools map -a gene_PLA2.bed -b - -c 4 -o mean >> pla2.TajimaD_mean.cv2.txt

echo -e "chrom\tstart\tend\tTajimaD" > svmp.TajimaD_mean.co2.txt; tail -n +2 ../cnv_masked_results/co.idaho.all.1kb.cnvMask.Tajima.D | grep -v 'nan' | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+1000,$4}' | bedtools sort -i - | bedtools map -a gene_SVMP.bed -b - -c 4 -o mean >> svmp.TajimaD_mean.co2.txt
echo -e "chrom\tstart\tend\tTajimaD" > svsp.TajimaD_mean.co2.txt; tail -n +2 ../cnv_masked_results/co.idaho.all.1kb.cnvMask.Tajima.D | grep -v 'nan' | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+1000,$4}' | bedtools sort -i - | bedtools map -a gene_SVSP.bed -b - -c 4 -o mean >> svsp.TajimaD_mean.co2.txt
echo -e "chrom\tstart\tend\tTajimaD" > pla2.TajimaD_mean.co2.txt; tail -n +2 ../cnv_masked_results/co.idaho.scaffold-mi7.250bp.cnvMask.Tajima.D | grep -v 'nan' | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+250,$4}' | bedtools sort -i - | bedtools map -a gene_PLA2.bed -b - -c 4 -o mean >> pla2.TajimaD_mean.co2.txt
```

#### Calculate mean df per venom gene

```
echo -e "chrom\tstart\tend\tdf_mean" > svmp.df_mean.cv1co1.txt; tail -n +2 ../cnv_masked_results/window.1kb.df_prop.all.cv1.co1.cnvMask.txt | grep -v 'nan' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | bedtools sort -i - | bedtools map -a gene_SVMP.bed -b - -c 4 -o mean >> svmp.df_mean.cv1co1.txt
echo -e "chrom\tstart\tend\tdf_mean" > svsp.df_mean.cv1co1.txt; tail -n +2 ../cnv_masked_results/window.1kb.df_prop.all.cv1.co1.cnvMask.txt | grep -v 'nan' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | bedtools sort -i - | bedtools map -a gene_SVSP.bed -b - -c 4 -o mean >> svsp.df_mean.cv1co1.txt
echo -e "chrom\tstart\tend\tdf_mean" > pla2.df_mean.cv1co1.txt; tail -n +2 ../cnv_masked_results/window.250bp.df_prop.scaffold-mi7.cv1.co1.cnvMask.txt | grep -v 'nan' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | bedtools sort -i - | bedtools map -a gene_PLA2.bed -b - -c 4 -o mean >> pla2.df_mean.cv1co1.txt
```

#### Calculate mean iHS per venom gene

*Note: The slick little first awk command checks for and removes missing data lines. They're lame!*

```
echo -e "chrom\tstart\tend\tiHS" > svmp.iHS_mean.cv1.txt; tail -n +2 ../rehh/cv.all_ihs.1kb.txt | awk '{ if ( $4 != "." ) { print $0; } }' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' | bedtools sort -i - | bedtools map -a gene_SVMP.bed -b - -c 4 -o mean >> svmp.iHS_mean.cv1.txt
echo -e "chrom\tstart\tend\tiHS" > svsp.iHS_mean.cv1.txt; tail -n +2 ../rehh/cv.all_ihs.1kb.txt | awk '{ if ( $4 != "." ) { print $0; } }' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' | bedtools sort -i - | bedtools map -a gene_SVSP.bed -b - -c 4 -o mean >> svsp.iHS_mean.cv1.txt
echo -e "chrom\tstart\tend\tiHS" > pla2.iHS_mean.cv1.txt; tail -n +2 ../rehh/cv.mi7_ihs.250bp.txt | awk '{ if ( $4 != "." ) { print $0; } }' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' | bedtools sort -i - | bedtools map -a gene_PLA2.bed -b - -c 4 -o mean >> pla2.iHS_mean.cv1.txt

echo -e "chrom\tstart\tend\tiHS" > svmp.iHS_mean.co1.txt; tail -n +2 ../cnv_masked_results/co.all_ihs.1kb.cnvMask.txt | awk '{ if ( $4 != "." ) { print $0; } }' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' | bedtools sort -i - | bedtools map -a gene_SVMP.bed -b - -c 4 -o mean >> svmp.iHS_mean.co1.txt
echo -e "chrom\tstart\tend\tiHS" > svsp.iHS_mean.co1.txt; tail -n +2 ../cnv_masked_results/co.all_ihs.1kb.cnvMask.txt | awk '{ if ( $4 != "." ) { print $0; } }' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' | bedtools sort -i - | bedtools map -a gene_SVSP.bed -b - -c 4 -o mean >> svsp.iHS_mean.co1.txt
echo -e "chrom\tstart\tend\tiHS" > pla2.iHS_mean.co1.txt; tail -n +2 ../cnv_masked_results/co.mi7_ihs.250bp.cnvMask.txt | awk '{ if ( $4 != "." ) { print $0; } }' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' | bedtools sort -i - | bedtools map -a gene_PLA2.bed -b - -c 4 -o mean >> pla2.iHS_mean.co1.txt
```

#### Calculate mean ß per venom gene

```
echo -e "chrom\tstart\tend\tiHS" > svmp.beta_mean.cv1.txt; tail -n +2 ../beta/results/cv1.phased.all.betascores.1kb.txt | awk '{ if ( $4 != "." ) { print $0; } }' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' | bedtools sort -i - | bedtools map -a gene_SVMP.bed -b - -c 4 -o mean >> svmp.beta_mean.cv1.txt
echo -e "chrom\tstart\tend\tiHS" > svsp.beta_mean.cv1.txt; tail -n +2 ../beta/results/cv1.phased.all.betascores.1kb.txt | awk '{ if ( $4 != "." ) { print $0; } }' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' | bedtools sort -i - | bedtools map -a gene_SVSP.bed -b - -c 4 -o mean >> svsp.beta_mean.cv1.txt
echo -e "chrom\tstart\tend\tiHS" > pla2.beta_mean.cv1.txt; tail -n +2 ../beta/results/cv1.phased.scaffold-mi7.betascores.250bp.txt | awk '{ if ( $4 != "." ) { print $0; } }' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' | bedtools sort -i - | bedtools map -a gene_PLA2.bed -b - -c 4 -o mean >> pla2.beta_mean.cv1.txt

echo -e "chrom\tstart\tend\tiHS" > svmp.beta_mean.co1.txt; tail -n +2 ../cnv_masked_results/co1.phased.all.betascores.1kb.cnvMask.txt | awk '{ if ( $4 != "." ) { print $0; } }' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' | bedtools sort -i - | bedtools map -a gene_SVMP.bed -b - -c 4 -o mean >> svmp.beta_mean.co1.txt
echo -e "chrom\tstart\tend\tiHS" > svsp.beta_mean.co1.txt; tail -n +2 ../cnv_masked_results/co1.phased.all.betascores.1kb.cnvMask.txt | awk '{ if ( $4 != "." ) { print $0; } }' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' | bedtools sort -i - | bedtools map -a gene_SVSP.bed -b - -c 4 -o mean >> svsp.beta_mean.co1.txt
echo -e "chrom\tstart\tend\tiHS" > pla2.beta_mean.co1.txt; tail -n +2 ../cnv_masked_results/co1.phased.scaffold-mi7.betascores.250bp.cnvMask.txt | awk '{ if ( $4 != "." ) { print $0; } }' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' | bedtools sort -i - | bedtools map -a gene_PLA2.bed -b - -c 4 -o mean >> pla2.beta_mean.co1.txt
```

### Selection appendix 3: Estimates for non-venom homologs

Parse results files for variation in non-venom homologs, for comparison to venom gene results.

These steps use the non-venom homolog BED coordinates in `resources/non-venom_paralogs_[SVMP,SVSP,PLA2].bed`

#### Extract π, dxy, and Fst results for non-venom homologs

```
cd pixy
mkdir non-venom

for pop in CV1 CV2 CO1 CO2; do for venom in SVMP SVSP PLA2; do echo -e "chromosome\twindow_pos_1\twindow_pos_2\tavg_pi" > ./non-venom/$pop.non-venom_${venom}.10kb_pi.txt; tail -n +2 ./pixy_results/pixy.all.10kb_pi.txt | grep $pop | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -a - -b non-venom_paralogs_${venom}.bed >> ./non-venom/$pop.non-venom_${venom}.10kb_pi.txt; done; done
for pop in CV1 CV2 CO1 CO2; do for venom in SVMP SVSP PLA2; do echo -e "chromosome\twindow_pos_1\twindow_pos_2\tavg_pi" > ./non-venom/$pop.non-venom_${venom}.1kb_pi.txt; cat ./pixy_results/pixy.*.1kb_pi.txt | grep $pop | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' | bedtools intersect -a - -b non-venom_paralogs_${venom}.bed >> ./non-venom/$pop.non-venom_${venom}.1kb_pi.txt; done; done

for pop in cv1co1 cv1cv2 co1co2; do for venom in SVMP SVSP PLA2; do echo -e "chromosome\twindow_pos_1\twindow_pos_2\tavg_dxy" > ./non-venom/$pop.non-venom_${venom}.10kb_dxy.txt; tail -n +2 ./pixy_results/dxy.10kb.$pop.txt | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b non-venom_paralogs_${venom}.bed >> ./non-venom/$pop.non-venom_${venom}.10kb_dxy.txt; done; done
for venom in SVMP SVSP PLA2; do echo -e "chromosome\twindow_pos_1\twindow_pos_2\tavg_dxy" > ./non-venom/cv1co1.non-venom_${venom}.1kb_dxy.txt; cat ./pixy_results/pixy.*.1kb_dxy.txt | grep -P "CV1\tCO1" | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b non-venom_paralogs_${venom}.bed >> ./non-venom/cv1co1.non-venom_${venom}.1kb_dxy.txt; done
for venom in SVMP SVSP PLA2; do echo -e "chromosome\twindow_pos_1\twindow_pos_2\tavg_dxy" > ./non-venom/cv1cv2.non-venom_${venom}.1kb_dxy.txt; cat ./pixy_results/pixy.*.1kb_dxy.txt | grep -P "CV1\tCV2" | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b non-venom_paralogs_${venom}.bed >> ./non-venom/cv1cv2.non-venom_${venom}.1kb_dxy.txt; done
for venom in SVMP SVSP PLA2; do echo -e "chromosome\twindow_pos_1\twindow_pos_2\tavg_dxy" > ./non-venom/co1co2.non-venom_${venom}.1kb_dxy.txt; cat ./pixy_results/pixy.*.1kb_dxy.txt | grep -P "CO1\tCO2" | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b non-venom_paralogs_${venom}.bed >> ./non-venom/co1co2.non-venom_${venom}.1kb_dxy.txt; done

for pop in cv1co1 cv1cv2 co1co2; do for venom in SVMP SVSP PLA2; do echo -e "chromosome\twindow_pos_1\twindow_pos_2\tavg_wc_fst" > ./non-venom/$pop.non-venom_${venom}.10kb_fst.txt; tail -n +2 ./pixy_results/fst.10kb.$pop.txt | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b non-venom_paralogs_${venom}.bed >> ./non-venom/$pop.non-venom_${venom}.10kb_fst.txt; done; done
for venom in SVMP SVSP PLA2; do echo -e "chromosome\twindow_pos_1\twindow_pos_2\tavg_wc_fst" > ./non-venom/cv1co1.non-venom_${venom}.1kb_fst.txt; cat ./pixy_results/pixy.*.1kb_fst.txt | grep -P "CV1\tCO1" | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b non-venom_paralogs_${venom}.bed >> ./non-venom/cv1co1.non-venom_${venom}.1kb_fst.txt; done
for venom in SVMP SVSP PLA2; do echo -e "chromosome\twindow_pos_1\twindow_pos_2\tavg_wc_fst" > ./non-venom/cv1cv2.non-venom_${venom}.1kb_fst.txt; cat ./pixy_results/pixy.*.1kb_fst.txt | grep -P "CV1\tCV2" | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b non-venom_paralogs_${venom}.bed >> ./non-venom/cv1cv2.non-venom_${venom}.1kb_fst.txt; done
for venom in SVMP SVSP PLA2; do echo -e "chromosome\twindow_pos_1\twindow_pos_2\tavg_wc_fst" > ./non-venom/co1co2.non-venom_${venom}.1kb_fst.txt; cat ./pixy_results/pixy.*.1kb_fst.txt | grep -P "CO1\tCO2" | awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6}' | bedtools intersect -a - -b non-venom_paralogs_${venom}.bed >> ./non-venom/co1co2.non-venom_${venom}.1kb_fst.txt; done
```

#### Extract Tajima's D results for non-venom homologs

```
cd tajima_d

for pop in cv.colorado cv.montana co.california co.idaho; do for venom in SVMP SVSP PLA2; do head -n 1 $pop.all.10kb.Tajima.D > $pop.non-venom_${venom}.10kb.Tajima.D; tail -n +2 $pop.all.10kb.Tajima.D | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+10000,$3,$4}' | bedtools intersect -a - -b non-venom_paralogs_${venom}.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$5}' >> $pop.non-venom_${venom}.10kb.Tajima.D; done; done
for pop in cv.colorado cv.montana co.california co.idaho; do for venom in SVMP SVSP PLA2; do head -n 1 $pop.all.1kb.Tajima.D > $pop.non-venom_${venom}.1kb.Tajima.D; tail -n +2 $pop.all.1kb.Tajima.D | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+10000,$3,$4}' | bedtools intersect -a - -b non-venom_paralogs_${venom}.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$5}' >> $pop.non-venom_${venom}.1kb.Tajima.D; done; done
```

#### Extract df results for non-venom homologs

```
cd df/results_df

for pop in cv1.co1 cv1.cv2 co1.co2; do for venom in SVMP SVSP PLA2; do head -n 1 window.10kb.df_prop.all.$pop.txt > window.10kb.df_prop.non-venom_${venom}.$pop.txt; tail -n +2 window.10kb.df_prop.all.$pop.txt | bedtools intersect -a - -b non-venom_paralogs_${venom}.bed >> window.10kb.df_prop.non-venom_${venom}.$pop.txt; done; done
for pop in cv1.co1 cv1.cv2 co1.co2; do for venom in SVMP SVSP PLA2; do head -n 1 window.1kb.df_prop.all.$pop.txt > window.1kb.df_prop.non-venom_${venom}.$pop.txt; tail -n +2 window.1kb.df_prop.all.$pop.txt | bedtools intersect -a - -b non-venom_paralogs_${venom}.bed >> window.1kb.df_prop.non-venom_${venom}.$pop.txt; done; done
```

#### Extract iHS results for non-venom homologs

```
cd rehh

for pop in cv co; do for venom in SVMP SVSP PLA2; do head -n 1 $pop.all_ihs.10kb.txt > $pop.non-venom_${venom}_ihs.10kb.txt; tail -n +2 $pop.all_ihs.10kb.txt | bedtools intersect -a - -b non-venom_paralogs_${venom}.bed >> $pop.non-venom_${venom}_ihs.10kb.txt; done; done
for pop in cv co; do for venom in SVMP SVSP PLA2; do head -n 1 $pop.all_ihs.1kb.txt > $pop.non-venom_${venom}_ihs.1kb.txt; tail -n +2 $pop.all_ihs.1kb.txt | bedtools intersect -a - -b non-venom_paralogs_${venom}.bed >> $pop.non-venom_${venom}_ihs.1kb.txt; done; done
```

#### Extract ß results for non-venom homologs

```
cd beta/results

for pop in cv1 co1; do for venom in SVMP SVSP PLA2; do head -n 1 $pop.phased.all.betascores.10kb.txt > $pop.phased.non-venom_${venom}.betascores.10kb.txt; tail -n +2 $pop.phased.all.betascores.10kb.txt | bedtools intersect -a - -b non-venom_paralogs_${venom}.bed >> $pop.phased.non-venom_${venom}.betascores.10kb.txt; done; done
for pop in cv1 co1; do for venom in SVMP SVSP PLA2; do head -n 1 $pop.phased.all.betascores.1kb.txt > $pop.phased.non-venom_${venom}.betascores.1kb.txt; tail -n +2 $pop.phased.all.betascores.1kb.txt | bedtools intersect -a - -b non-venom_paralogs_${venom}.bed >> $pop.phased.non-venom_${venom}.betascores.1kb.txt; done; done
```

### Selection appendix 4: Point estimates for all genes

Calculate mean iHS and ß for all genes in the genome annotation.

#### Set up environment

```
cd gene_means
```

These commands will use the gene coordinates in `resources/gene.all.bed`.

#### Calculate mean iHS and ß per gene

```
echo -e "chrom\tstart\tend\tiHS" > gene.all.iHS_mean.cv1.txt; tail -n +2 ../rehh/cv.all_ihs.1kb.txt | awk '{ if ( $4 != "." ) { print $0; } }' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' | bedtools sort -i - | bedtools map -a gene.all.bed -b - -c 4 -o mean >> gene.all.iHS_mean.cv1.txt
echo -e "chrom\tstart\tend\tiHS" > gene.all.iHS_mean.co1.txt; tail -n +2 ../cnv_masked_results/co.all_ihs.1kb.cnvMask.txt | awk '{ if ( $4 != "." ) { print $0; } }' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' | bedtools sort -i - | bedtools map -a gene.all.bed -b - -c 4 -o mean >> gene.all.iHS_mean.co1.txt

echo -e "chrom\tstart\tend\tbeta" > gene.all.beta_mean.cv1.txt; tail -n +2 ../beta/results/cv1.phased.all.betascores.1kb.txt | awk '{ if ( $4 != "." ) { print $0; } }' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' | bedtools sort -i - | bedtools map -a gene.all.bed -b - -c 4 -o mean >> gene.all.beta_mean.cv1.txt
echo -e "chrom\tstart\tend\tbeta" > gene.all.beta_mean.co1.txt; tail -n +2 ../cnv_masked_results/co1.phased.all.betascores.1kb.cnvMask.txt | awk '{ if ( $4 != "." ) { print $0; } }' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' | bedtools sort -i - | bedtools map -a gene.all.bed -b - -c 4 -o mean >> gene.all.beta_mean.co1.txt
```

### Selection appendix 5: Point estimates for 'other' venom genes

#### Set up environment

```
mkdir otherVG
cd otherVG
```

These commands will use the coordinates in `resources/otherVG.bed`.

#### Calculate mean Tajima's D for 'other' venom genes

```
echo -e "chrom\tstart\tend\tgene\tTajimaD" > otherVG.TajimaD_mean.cv1.txt; tail -n +2 ../tajima_d/cv.colorado.all.1kb.Tajima.D | grep -v 'nan' | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+1000,$4}' | bedtools sort -i - | bedtools map -a otherVG.sort.bed -b - -c 4 -o mean >> otherVG.TajimaD_mean.cv1.txt
echo -e "chrom\tstart\tend\tgene\tTajimaD" > otherVG.TajimaD_mean.co1.txt; tail -n +2 ../cnv_masked_results/co.california.all.1kb.cnvMask.Tajima.D | grep -v 'nan' | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+1000,$4}' | bedtools sort -i - | bedtools map -a otherVG.sort.bed -b - -c 4 -o mean >> otherVG.TajimaD_mean.co1.txt
```

#### Calculate mean iHS for 'other' venom genes

```
echo -e "chrom\tstart\tend\tgene\tiHS" > otherVG.iHS_mean.cv1.txt; tail -n +2 ../rehh/cv.all_ihs.1kb.txt | awk '{ if ( $4 != "." ) { print $0; } }' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' | bedtools sort -i - | bedtools map -a otherVG.sort.bed -b - -c 4 -o mean >> otherVG.iHS_mean.cv1.txt
echo -e "chrom\tstart\tend\tgene\tiHS" > otherVG.iHS_mean.co1.txt; tail -n +2 ../cnv_masked_results/co.all_ihs.1kb.cnvMask.txt | awk '{ if ( $4 != "." ) { print $0; } }' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' | bedtools sort -i - | bedtools map -a otherVG.sort.bed -b - -c 4 -o mean >> otherVG.iHS_mean.co1.txt
```

#### Calculate mean ß for 'other' venom genes

```
echo -e "chrom\tstart\tend\tgene\tbeta" > otherVG.beta_mean.cv1.txt; tail -n +2 ../beta/results/cv1.phased.all.betascores.1kb.txt | awk '{ if ( $4 != "." ) { print $0; } }' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' | bedtools sort -i - | bedtools map -a otherVG.sort.bed -b - -c 4 -o mean >> otherVG.beta_mean.cv1.txt
echo -e "chrom\tstart\tend\tgene\tbeta" > otherVG.beta_mean.co1.txt; tail -n +2 ../cnv_masked_results/co1.phased.all.betascores.1kb.cnvMask.txt | awk '{ if ( $4 != "." ) { print $0; } }' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' | bedtools sort -i - | bedtools map -a otherVG.sort.bed -b - -c 4 -o mean >> otherVG.beta_mean.co1.txt
```

### Selection appendix 6: Proportion heterozygotes at high ß sites

Quantify the proportion heterozygotes within populations at variable sites, then focus on SNPs with high ß.

#### Set up environment

```
mkdir heterozygosity
cd heterozygosity
```

These commands will use the population lists and sliding window BED files in `resources`.

#### Quantify the number of heterozygotes per site using vcftools (HWE output)

These analyses will be performed on phased variants in CV1 and CO1 populations from [here](https://figshare.com/articles/dataset/Phased_VCFs/16415556).

```
for chrom in scaffold-mi1 scaffold-mi2 scaffold-mi7; do vcftools --vcf viridis.phased.$chrom.vcf --keep pop.list.cv.colorado --hardy --out cv1.$chrom; done
for chrom in scaffold-mi1 scaffold-mi2 scaffold-mi7; do vcftools --vcf oreganus.phased.$chrom.vcf --keep pop.list.co.california --hardy --out co1.$chrom; done
```

#### Format outputs so they are useful

The script `formatHeterozygosity.py` in the `python` subdirectory is designed to remove sites with no data and quantify the proportion of heterozygotes at each SNP. It also removes sites with HWE p-values < 0.0001, and allows the user to specify the minimum number of genotyped individuals at a SNP to be included (e.g., 10).

```
python formatHeterozygosity.py cv1.scaffold-mi1.hwe 18 > het.cv1.scaffold-mi1.txt
python formatHeterozygosity.py cv1.scaffold-mi2.hwe 18 > het.cv1.scaffold-mi2.txt
python formatHeterozygosity.py cv1.scaffold-mi7.hwe 18 > het.cv1.scaffold-mi7.txt
python formatHeterozygosity.py co1.scaffold-mi1.hwe 17 > het.co1.scaffold-mi1.txt
python formatHeterozygosity.py co1.scaffold-mi2.hwe 17 > het.co1.scaffold-mi2.txt
python formatHeterozygosity.py co1.scaffold-mi7.hwe 17 > het.co1.scaffold-mi7.txt
```

#### Use bedtools to quantify mean proportion of heterozygotes in sliding windows

```
for chrom in scaffold-mi1 scaffold-mi2 scaffold-mi7; do for window in 1kb 10kb 100kb; do echo -e "chrom\tstart\tend\thom1\thet\thom2\tprop_het" > het.cv1.$chrom.$window.txt; tail -n +2 het.cv1.$chrom.txt | bedtools map -a CroVir_genome_${window}_window.bed -b - -c 4,5,6,8 -o mean | grep -w $chrom >> het.cv1.$chrom.$window.txt; done; done
for chrom in scaffold-mi1 scaffold-mi2 scaffold-mi7; do for window in 1kb 10kb 100kb; do echo -e "chrom\tstart\tend\thom1\thet\thom2\tprop_het" > het.co1.$chrom.$window.txt; tail -n +2 het.co1.$chrom.txt | bedtools map -a CroVir_genome_${window}_window.bed -b - -c 4,5,6,8 -o mean | grep -w $chrom >> het.co1.$chrom.$window.txt; done; done
```

#### Mask CNVs

```
cd cnv_mask_results
head -n 1 ../heterozygosity/het.co1.scaffold-mi1.10kb.txt > het.co1.scaffold-mi1.10kb.cnvMask.txt; python cnvMaskSelection.py ../heterozygosity/het.co1.scaffold-mi1.10kb.txt ./mask_bed/cnv_filter.10kb.bed no >> het.co1.scaffold-mi1.10kb.cnvMask.txt
head -n 1 ../heterozygosity/het.co1.scaffold-mi1.1kb.txt > het.co1.scaffold-mi1.1kb.cnvMask.txt; python cnvMaskSelection.py ../heterozygosity/het.co1.scaffold-mi1.1kb.txt ./mask_bed/cnv_filter.1kb.bed no >> het.co1.scaffold-mi1.1kb.cnvMask.txt
head -n 1 ../heterozygosity/het.co1.scaffold-mi2.10kb.txt > het.co1.scaffold-mi2.10kb.cnvMask.txt; python cnvMaskSelection.py ../heterozygosity/het.co1.scaffold-mi2.10kb.txt ./mask_bed/cnv_filter.10kb.bed no >> het.co1.scaffold-mi2.10kb.cnvMask.txt
head -n 1 ../heterozygosity/het.co1.scaffold-mi2.1kb.txt > het.co1.scaffold-mi2.1kb.cnvMask.txt; python cnvMaskSelection.py ../heterozygosity/het.co1.scaffold-mi2.1kb.txt ./mask_bed/cnv_filter.1kb.bed no >> het.co1.scaffold-mi2.1kb.cnvMask.txt
head -n 1 ../heterozygosity/het.co1.scaffold-mi7.10kb.txt > het.co1.scaffold-mi7.10kb.cnvMask.txt; python cnvMaskSelection.py ../heterozygosity/het.co1.scaffold-mi7.10kb.txt ./mask_bed/cnv_filter.10kb.bed no >> het.co1.scaffold-mi7.10kb.cnvMask.txt
head -n 1 ../heterozygosity/het.co1.scaffold-mi7.1kb.txt > het.co1.scaffold-mi7.1kb.cnvMask.txt; python cnvMaskSelection.py ../heterozygosity/het.co1.scaffold-mi7.1kb.txt ./mask_bed/cnv_filter.1kb.bed no >> het.co1.scaffold-mi7.1kb.cnvMask.txt
```

#### Extract sites with ß scores > 95th quantile

CV1 95th quantile = 3.345
CO1 95th quantile = 3.841

```
awk '$4 > 3.344' ../beta/results/cv1.phased.scaffold-mi1.betascores.bed > beta.outliers.cv1.scaffold-mi1.bed
awk '$4 > 3.344' ../beta/results/cv1.phased.scaffold-mi2.betascores.bed > beta.outliers.cv1.scaffold-mi2.bed
awk '$4 > 3.344' ../beta/results/cv1.phased.scaffold-mi7.betascores.bed > beta.outliers.cv1.scaffold-mi7.bed
awk '$4 > 3.84' ../beta/results/co1.phased.scaffold-mi1.betascores.bed > beta.outliers.co1.scaffold-mi1.bed
awk '$4 > 3.84' ../beta/results/co1.phased.scaffold-mi2.betascores.bed > beta.outliers.co1.scaffold-mi2.bed
awk '$4 > 3.84' ../beta/results/co1.phased.scaffold-mi7.betascores.bed > beta.outliers.co1.scaffold-mi7.bed
```

#### Determine overlap between formatted HWE data and ß outliers

```
echo -e "chrom\tstart\tend\thom1\thet\thom2\tprop_het" > het.cv1.scaffold-mi1.beta.txt; tail -n +2 het.cv1.scaffold-mi1.txt | bedtools intersect -a - -b beta.outliers.cv1.scaffold-mi1.bed >> het.cv1.scaffold-mi1.beta.txt
echo -e "chrom\tstart\tend\thom1\thet\thom2\tprop_het" > het.cv1.scaffold-mi2.beta.txt; tail -n +2 het.cv1.scaffold-mi2.txt | bedtools intersect -a - -b beta.outliers.cv1.scaffold-mi2.bed >> het.cv1.scaffold-mi2.beta.txt
echo -e "chrom\tstart\tend\thom1\thet\thom2\tprop_het" > het.cv1.scaffold-mi7.beta.txt; tail -n +2 het.cv1.scaffold-mi7.txt | bedtools intersect -a - -b beta.outliers.cv1.scaffold-mi7.bed >> het.cv1.scaffold-mi7.beta.txt
echo -e "chrom\tstart\tend\thom1\thet\thom2\tprop_het" > het.co1.scaffold-mi1.beta.txt; tail -n +2 het.co1.scaffold-mi1.txt | bedtools intersect -a - -b beta.outliers.co1.scaffold-mi1.bed >> het.co1.scaffold-mi1.beta.txt
echo -e "chrom\tstart\tend\thom1\thet\thom2\tprop_het" > het.co1.scaffold-mi2.beta.txt; tail -n +2 het.co1.scaffold-mi2.txt | bedtools intersect -a - -b beta.outliers.co1.scaffold-mi2.bed >> het.co1.scaffold-mi2.beta.txt
echo -e "chrom\tstart\tend\thom1\thet\thom2\tprop_het" > het.co1.scaffold-mi7.beta.txt; tail -n +2 het.co1.scaffold-mi7.txt | bedtools intersect -a - -b beta.outliers.co1.scaffold-mi7.bed >> het.co1.scaffold-mi7.beta.txt
```

#### Extract per site values within venom genes

```
echo -e "chrom\tstart\tend\thom1\thet\thom2\tprop_het" > het.beta.cv1.svmp.txt; tail -n +2 het.cv1.scaffold-mi1.beta.txt | bedtools intersect -a - -b gene.SVMP.bed >> het.beta.cv1.svmp.txt
echo -e "chrom\tstart\tend\thom1\thet\thom2\tprop_het" > het.beta.co1.svmp.txt; tail -n +2 het.co1.scaffold-mi1.beta.txt | bedtools intersect -a - -b gene.SVMP.co.bed >> het.beta.co1.svmp.txt
echo -e "chrom\tstart\tend\thom1\thet\thom2\tprop_het" > het.beta.cv1.svsp.txt; tail -n +2 het.cv1.scaffold-mi2.beta.txt | bedtools intersect -a - -b gene.SVSP.bed >> het.beta.cv1.svsp.txt
echo -e "chrom\tstart\tend\thom1\thet\thom2\tprop_het" > het.beta.co1.svsp.txt; tail -n +2 het.co1.scaffold-mi2.beta.txt | bedtools intersect -a - -b gene.SVSP.co.bed >> het.beta.co1.svsp.txt
echo -e "chrom\tstart\tend\thom1\thet\thom2\tprop_het" > het.beta.cv1.pla2.txt; tail -n +2 het.cv1.scaffold-mi7.beta.txt | bedtools intersect -a - -b gene.PLA2.bed >> het.beta.cv1.pla2.txt
echo -e "chrom\tstart\tend\thom1\thet\thom2\tprop_het" > het.beta.co1.pla2.txt; tail -n +2 het.co1.scaffold-mi7.beta.txt | bedtools intersect -a - -b gene.PLA2.co.bed >> het.beta.co1.pla2.txt
```

#### Extract per site values within all genes

```
echo -e "chrom\tstart\tend\thom1\thet\thom2\tprop_het" > het.beta.cv1.scaffold-mi1.gene.txt; tail -n +2 het.cv1.scaffold-mi1.beta.txt | bedtools intersect -a - -b gene.all.scaffold-mi1.bed | bedtools intersect -v -a - -b gene.SVMP.bed >> het.beta.cv1.scaffold-mi1.gene.txt
echo -e "chrom\tstart\tend\thom1\thet\thom2\tprop_het" > het.beta.co1.scaffold-mi1.gene.txt; tail -n +2 het.co1.scaffold-mi1.beta.txt | bedtools intersect -a - -b gene.all.scaffold-mi1.bed | bedtools intersect -v -a - -b gene.SVMP.co.bed >> het.beta.co1.scaffold-mi1.gene.txt
echo -e "chrom\tstart\tend\thom1\thet\thom2\tprop_het" > het.beta.cv1.scaffold-mi2.gene.txt; tail -n +2 het.cv1.scaffold-mi2.beta.txt | bedtools intersect -a - -b gene.all.scaffold-mi2.bed | bedtools intersect -v -a - -b gene.SVSP.bed >> het.beta.cv1.scaffold-mi2.gene.txt
echo -e "chrom\tstart\tend\thom1\thet\thom2\tprop_het" > het.beta.co1.scaffold-mi2.gene.txt; tail -n +2 het.co1.scaffold-mi2.beta.txt | bedtools intersect -a - -b gene.all.scaffold-mi2.bed | bedtools intersect -v -a - -b gene.SVSP.co.bed >> het.beta.co1.scaffold-mi2.gene.txt
echo -e "chrom\tstart\tend\thom1\thet\thom2\tprop_het" > het.beta.cv1.scaffold-mi7.gene.txt; tail -n +2 het.cv1.scaffold-mi7.beta.txt | bedtools intersect -a - -b gene.all.scaffold-mi7.bed | bedtools intersect -v -a - -b gene.PLA2.bed >> het.beta.cv1.scaffold-mi7.gene.txt
echo -e "chrom\tstart\tend\thom1\thet\thom2\tprop_het" > het.beta.co1.scaffold-mi7.gene.txt; tail -n +2 het.co1.scaffold-mi7.beta.txt | bedtools intersect -a - -b gene.all.scaffold-mi7.bed | bedtools intersect -v -a - -b gene.PLA2.co.bed >> het.beta.co1.scaffold-mi7.gene.txt
```

### Selection appendix 7: Trans-species polymorphism

Produce phased alignments to query for trans-species polymorphisms and private polymorphisms. The alignments will be generated from the phased VCF data for CV1 and CO1 populations, and outgroup C. atrox data from the all-sites VCFs.

#### Set up environment

```
mkdir alignments
cd alignments
mkdir input
```

#### Format input VCFs

The phased VCF files need to be bgzipped and tabix indexed.

```
bgzip -c viridis.phased.scaffold-mi1.vcf > ./input/viridis.phased.scaffold-mi1.vcf.gz
tabix -p vcf ./input/viridis.phased.scaffold-mi1.vcf.gz

bgzip -c oreganus.phased.scaffold-mi1.vcf > ./input/oreganus.phased.scaffold-mi1.vcf.gz
tabix -p vcf ./input/oreganus.phased.scaffold-mi1.vcf.gz

bgzip -c viridis.phased.scaffold-mi2.vcf > ./input/viridis.phased.scaffold-mi2.vcf.gz
tabix -p vcf ./input/viridis.phased.scaffold-mi2.vcf.gz

bgzip -c oreganus.phased.scaffold-mi2.vcf > ./input/oreganus.phased.scaffold-mi2.vcf.gz
tabix -p vcf ./input/oreganus.phased.scaffold-mi2.vcf.gz

bgzip -c viridis.phased.scaffold-mi7.vcf > ./input/viridis.phased.scaffold-mi7.vcf.gz
tabix -p vcf ./input/viridis.phased.scaffold-mi7.vcf.gz

bgzip -c oreganus.phased.scaffold-mi7.vcf > ./input/oreganus.phased.scaffold-mi7.vcf.gz
tabix -p vcf ./input/oreganus.phased.scaffold-mi7.vcf.gz
```

#### Format fasta files for venom-linked chromosomes with both haplotypes per sample

This script will extract data from the VCF files and produce 'consensus' sequences per haplotype based on these data and the reference genome. It also requires sample lists for CV1 and CO1 in `resources`. The last line produces a 'sequential' fasta output without sequence text wrapping - comment this line out or delete the sequential output if you don't care about this.

__*fastaMake.sh*__
```
for chrom in scaffold-mi1 scaffold-mi2 scaffold-mi7; do
	touch venom.$chrom.fasta
	samtools faidx ../CroVir_genome_L77pg_16Aug2017.final_rename.fasta $chrom | bcftools consensus ../vcf/vcf_chrom-specific_cvco+outgroup/cvco+outgroup.mask.HardFilter.depth.chrom.$chrom.vcf.gz -s CA0346 -H 1 -M N -i 'type="snp" | alt="."' >> venom.$chrom.fasta
	sed -i.bak "s/$chrom/atrox_CA0346_A/" venom.$chrom.fasta
	samtools faidx ../CroVir_genome_L77pg_16Aug2017.final_rename.fasta $chrom | bcftools consensus ../vcf/vcf_chrom-specific_cvco+outgroup/cvco+outgroup.mask.HardFilter.depth.chrom.$chrom.vcf.gz -s CA0346 -H 2 -M N -i 'type="snp" | alt="."' >> venom.$chrom.fasta
	sed -i.bak "s/$chrom/atrox_CA0346_B/" venom.$chrom.fasta
	for indv in `cat pop.list.cv1`; do
		samtools faidx ../CroVir_genome_L77pg_16Aug2017.final_rename.fasta $chrom | bcftools consensus ./input/viridis.phased.$chrom.vcf.gz -s $indv -H 1 >> venom.$chrom.fasta
		sed -i.bak "s/$chrom/cv1_${indv}_A/" venom.$chrom.fasta
		samtools faidx ../CroVir_genome_L77pg_16Aug2017.final_rename.fasta $chrom | bcftools consensus ./input/viridis.phased.$chrom.vcf.gz -s $indv -H 2 >> venom.$chrom.fasta
		sed -i.bak "s/$chrom/cv1_${indv}_B/" venom.$chrom.fasta
	done
	for indv in `cat pop.list.co1`; do
		samtools faidx ../CroVir_genome_L77pg_16Aug2017.final_rename.fasta $chrom | bcftools consensus ./input/oreganus.phased.$chrom.vcf.gz -s $indv -H 1 >> venom.$chrom.fasta
		sed -i.bak "s/$chrom/co1_${indv}_A/" venom.$chrom.fasta
		samtools faidx ../CroVir_genome_L77pg_16Aug2017.final_rename.fasta $chrom | bcftools consensus ./input/oreganus.phased.$chrom.vcf.gz -s $indv -H 2 >> venom.$chrom.fasta
		sed -i.bak "s/$chrom/co1_${indv}_B/" venom.$chrom.fasta
	done
	sed -i.bak "s/*/N/g" venom.$chrom.fasta
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < venom.$chrom.fasta | tail -n +2 > venom.$chrom.sequential.fasta
done
```

Run the script.

```
sh fastaMake.sh
```

## Recombination rate variation and linkage disequilibrium analysis

## Analysis in R

## Appendix 1: Mapping statistics

## Appendix 2: Genotype quality



































