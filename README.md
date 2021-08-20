# Rattlesnake venom population genomics

![Rattlesnake venom population genomics](/cover-image.png "cover image")

This repository contains details on the data processing and analysis steps taken in analyses of population genomic variation and signatures of selection in major rattlesnake venom gene regions. This workflow is a companion to the methods description in Schield et al. (in review). Analysis of recombination rates and statistics requiring phased variants are based on [recombination maps](https://figshare.com/articles/dataset/Rattlesnake_Recombination_Maps/11283224) and phased data from [Schield et al. _MBE_ 37: 1272-1294](https://academic.oup.com/mbe/advance-article-abstract/doi/10.1093/molbev/msaa003/5700722).

Lists and reference files (i.e., BED, GFF, etc.) are in the `resources` directory. Shell and Python scripts are in respective `shell` and `python` directories. R scripts are in the `R` directory. Note that you may need to adjust the organization of file locations to suit your environment.

## Contents

* [Software and dependencies](#software-and-dependencies)
* [General resources](#general-resources)
* [Read filtering](#read-filtering)
* [Read mapping](#read-mapping)
* Variant calling
* Variant filtering
* Analysis of copy-number variation
* Population structure analysis
* Demographic analysis
* Population genetic diversity and differentiation
* Signatures of selection
* Recombination rate variation and linkage disequilibrium
* Analysis in R

### Software and dependencies

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
* [betascan](https://github.com/ksiewert/BetaScan)
* [rehh](https://cran.r-project.org/web/packages/rehh/vignettes/rehh.html) - R package.
* [R](https://cran.r-project.org/)

Note: I installed a number of these dependencies using [conda](https://docs.conda.io/en/latest/).

### General resources

#### Processing files

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

#### Reference genome

Use the [Prairie Rattlesnake (_Crotalus viridis_) genome assembly](https://figshare.com/articles/dataset/Prairie_rattlesnake_Crotalus_viridis_genome_assembly/9030782) as the reference in this workflow.

Prepare necessary indexes for downstream analysis.

```
bwa index CroVir_genome_L77pg_16Aug2017.final_rename.fasta
samtools faidx CroVir_genome_L77pg_16Aug2017.final_rename.fasta
./gatk-4.0.8.1/gatk CreateSequenceDictionary -R CroVir_genome_L77pg_16Aug2017.final_rename.fasta
```

### Read filtering

Quality trim and filter raw whole genome resequencing reads using trimmomatic using these settings:

* Remove 5' end bases if quality is below 20
* Remove 3' end bases if quality is below 20
* Minimum read length = 32
* Remove reads if average quality is < 30

#### Set up environment

Get raw fastq data into `fastq` directory. <br /> Make a `fastq_filtered` directory for output.

```
mkdir fastq
mkdir fastq_filtered
```

#### Run trimmomatic on raw reads

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

### Read mapping

Use bwa 'mem' to map our filtered reads to the reference genome.

#### Set up environment

`mkdir bam`

#### Map reads with bwa and sort with samtools

The script below will run bwa mem on the paired, filtered reads per sample and sort the output bam file.

__*bwa_mem.sh*__

```
list=$1
for line in `cat $list`; do
	name=$line
	echo "Mapping filtered $name data to reference"
	bwa mem -t 16 -R "@RG\tID:$name\tLB:CVOS\tPL:illumina\tPU:NovaSeq6000\tSM:$name" ../genome_crotalus/CroVir_genome_L77pg_16Aug2017.final_rename.fasta ./fastq_filtered/${name}_R1_P.trim.fastq.gz ./fastq_filtered/${name}_R2_P.trim.fastq.gz | samtools sort -@ 16 -O bam -T temp -o ./2ref_male/bam/$name.bam -
done
```

Run the script.

`sh bwa_mem.sh sample.list`

### Variant calling

Use GATK for individual variant discovery and variant calling among the cohort of samples. This is a two-step process, first using HaplotypeCaller to generate individual genomic VCFs (gVCFs), then using GenotypeGVCFs to call variants among samples and generate an all-sites VCF.

#### Set up environment

```
mkdir gvcf
mkdir vcf
```

#### Call individual variable sites using HaplotypeCaller

The script below will run GATK HaplotypeCaller on each sample in the `sample.list`. It will also zip and index the gVCF output.

*Note: GATK HaplotypeCaller will take a small eternity to run on all of these samples one-by-one. Consider breaking up the job into smaller lists of samples and running jobs in parallel.*

GATK_HaplotypeCaller.sh

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

#### Call cohort variant sites and generate an 'all-sites' VCF using GenotypeGVCFs



