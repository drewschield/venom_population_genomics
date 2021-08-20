# Rattlesnake venom population genomics

This repository contains details on the data processing and analysis steps taken in analyses of population genomic variation and signatures of selection in major rattlesnake venom gene regions. This workflow is a companion to the methods description in Schield et al. (in review). Analysis of recombination rates and statistics requiring phased variants are based on [recombination maps](https://figshare.com/articles/dataset/Rattlesnake_Recombination_Maps/11283224) and phased data from [Schield et al. _MBE_ 37: 1272-1294](https://academic.oup.com/mbe/advance-article-abstract/doi/10.1093/molbev/msaa003/5700722).

The steps described below use on the following software, and assume that dependencies are on the user path:

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

Note: I installed a number of these dependencies using [conda](https://docs.conda.io/en/latest/).

Lists and miscellaneous files are in the `processing_files` directory.
Shell and Python scripts are in respective `shell` and `python` directories.
R scripts are in the `R` directory.


