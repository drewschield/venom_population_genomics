list=$1
window=$2
abbrev=$3
ans=$4
for chrom in `cat $list`; do
	pixy --stats pi fst dxy --vcf ../vcf/vcf_chrom-specific_cvco+outgroup/cvco+outgroup.mask.HardFilter.depth.chrom.$chrom.vcf.gz --zarr_path ./pixy_zarr --reuse_zarr $ans --window_size $window --populations pixy.popmap --variant_filter_expression 'DP>=5' --invariant_filter_expression 'DP>=5' --outfile_prefix ./pixy_results/pixy.$chrom.$abbrev
done
