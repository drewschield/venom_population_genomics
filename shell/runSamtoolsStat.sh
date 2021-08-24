list=$1
for indv in `cat $list`; do
	echo calculating mapping statistics for $indv
	samtools stat -@ 8 ../bam/$indv.bam > $indv.bam.stat.txt
done