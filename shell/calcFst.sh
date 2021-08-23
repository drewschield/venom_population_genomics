list1=$1
list2=$2
pop1=$3
pop2=$4
for chrom in `cat chrom.list`; do
	vcftools --gzvcf ../vcf/vcf_chrom-specific_cvco+outgroup/cvco+outgroup.mask.HardFilter.depth.chrom.$chrom.vcf.gz --min-alleles 2 --max-alleles 2 --maf 0.05 --weir-fst-pop $list1 --weir-fst-pop $list2 --out ./results_fst/fst.$chrom.$pop1.$pop2
done
