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
