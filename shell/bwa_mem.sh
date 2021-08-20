list=$1
for line in `cat $list`; do
	name=$line
	echo "Mapping filtered $name data to reference"
	bwa mem -t 16 -R "@RG\tID:$name\tLB:CVOS\tPL:illumina\tPU:NovaSeq6000\tSM:$name" CroVir_genome_L77pg_16Aug2017.final_rename.fasta ./fastq_filtered/${name}_R1_P.trim.fastq.gz ./fastq_filtered/${name}_R2_P.trim.fastq.gz | samtools sort -@ 16 -O bam -T temp -o ./bam/$name.bam -
done
