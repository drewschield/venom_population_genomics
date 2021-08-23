vbed=$1
data=$2
lines=`wc -l $vbed | cut -d' ' -f 1`
lines_fix=`echo "$(($lines-1))"`
for perm in $(seq 1 10000); do
	mean=`tail -n +2 $data | shuf -n 1 | grep -f - -A $lines_fix $data | awk '{print $4}' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }'`
	echo "${perm}\t${mean}"
done
