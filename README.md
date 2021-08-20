# Rattlesnake venom population genomics

This repository contains details on the data processing and analysis steps taken in analyses of population genomic variation and signatures of selection in major rattlesnake venom gene regions. This workflow is a companion to the methods description in Schield et al. (in review).

The steps described below rely on the following software, and assume that dependencies are on the user path:

* trimmomatic
* bwa
* GATK (v3.8-1-0 and v4.0.8.1)
* htslib
* samtools
* bcftools
* bgzip
* tabix
* vcftools
* bedtools
* CNV-seq
* ADMIXTURE
* PSMC
* [pixy](https://pixy.readthedocs.io/en/latest/)
* betascan
* R

Note: I installed and maintain a number of these dependencies within a [conda](https://docs.conda.io/en/latest/) environment.

