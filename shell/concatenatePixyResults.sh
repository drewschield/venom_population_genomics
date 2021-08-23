abbrev=$1
# concatenate pi output
head -n 1 ./pixy_results/pixy.scaffold-ma1.${abbrev}_pi.txt > ./pixy_results/pixy.all.${abbrev}_pi.txt
for chrom in `cat chrom.list`; do
	tail -n +2 ./pixy_results/pixy.${chrom}.${abbrev}_pi.txt >> ./pixy_results/pixy.all.${abbrev}_pi.txt
done
# concatenate dxy output
head -n 1 ./pixy_results/pixy.scaffold-ma1.${abbrev}_dxy.txt > ./pixy_results/pixy.all.${abbrev}_dxy.txt
for chrom in `cat chrom.list`; do
	tail -n +2 ./pixy_results/pixy.${chrom}.${abbrev}_dxy.txt >> ./pixy_results/pixy.all.${abbrev}_dxy.txt
done
# concatenate fst output
head -n 1 ./pixy_results/pixy.scaffold-ma1.${abbrev}_fst.txt > ./pixy_results/pixy.all.${abbrev}_fst.txt
for chrom in `cat chrom.list`; do
	tail -n +2 ./pixy_results/pixy.${chrom}.${abbrev}_fst.txt >> ./pixy_results/pixy.all.${abbrev}_fst.txt
done
