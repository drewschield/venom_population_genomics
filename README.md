# Rattlesnake venom population genomics

![Rattlesnake venom population genomics](/cover-image.png "cover image")

This repository contains details on the data processing and analysis steps taken in analyses of population genomic variation and signatures of selection in major rattlesnake venom gene regions. This workflow is a companion to the methods description in Schield et al. (in review). Analysis of recombination rates and statistics requiring phased variants are based on [recombination maps](https://figshare.com/articles/dataset/Rattlesnake_Recombination_Maps/11283224) and phased data from [Schield et al. _MBE_ 37: 1272-1294](https://academic.oup.com/mbe/advance-article-abstract/doi/10.1093/molbev/msaa003/5700722).

Lists and reference files (i.e., BED, GFF, etc.) are in the `resources` directory. Shell and Python scripts are in respective `shell` and `python` directories. R scripts are in the `R` directory. Note that you may need to adjust the organization of file locations to suit your environment.

## Contents

* [Software and dependencies](#software-and-dependencies)
* [General resources](#general-resources)
* [Read filtering](#read-filtering)
* Read mapping
* Quantifying mapping results
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

### Read filtering

We'll quality trim and filter raw whole genome resequencing reads using `trimmomatic` using these settings:

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
