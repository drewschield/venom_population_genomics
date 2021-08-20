chromlist=$1
for chrom in `cat $chromlist`; do
	echo parsing $chrom VCF
	bcftools view --threads 16 -r $chrom -O z -o ./vcf/vcf_chrom-specific_cvco+outgroup/cvco+outgroup.mask.HardFilter.depth.chrom.$chrom.vcf.gz ./vcf/cvco+outgroup.mask.HardFilter.depth.chrom.vcf.gz
	tabix -p vcf ./vcf/vcf_chrom-specific_cvco+outgroup/cvco+outgroup.mask.HardFilter.depth.chrom.$chrom.vcf.gz
done
