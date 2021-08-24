import sys

print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % ('chrom', 'start', 'end', 'mean', 'p0.500', 'p0.025', 'p0.975')

with open(sys.argv[1], 'r') as rho:
	next(rho)
	for line in rho:
		chrom = line.split()[0]
		start = line.split()[1]
		end = line.split()[2]
		if 'scaffold-mi' in str(chrom):
			print line.strip()
		else:
			if chrom == 'scaffold-ma1':
				if int(start) >= 156000000 and int(end) <= 167000000:
					print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chrom, start, end, 'NA', 'NA', 'NA', 'NA')
				else:
					print line.strip()
			if chrom == 'scaffold-ma2':
				if int(start) >= 135000000 and int(end) <= 145000000:
					print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chrom, start, end, 'NA', 'NA', 'NA', 'NA')
				else:
					print line.strip()
			if chrom == 'scaffold-ma3':
				if int(start) >= 80000000 and int(end) <= 94000000:
					print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chrom, start, end, 'NA', 'NA', 'NA', 'NA')
				else:
					print line.strip()
			if chrom == 'scaffold-ma4':
				if int(start) >= 34000000 and int(end) <= 48000000:
					print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chrom, start, end, 'NA', 'NA', 'NA', 'NA')
				else:
					print line.strip()
			if chrom == 'scaffold-ma5':
				if int(start) >= 15000000 and int(end) <= 28000000:
					print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chrom, start, end, 'NA', 'NA', 'NA', 'NA')
				else:
					print line.strip()
			if chrom == 'scaffold-ma6':
				if int(start) >= 19000000 and int(end) <= 36000000:
					print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chrom, start, end, 'NA', 'NA', 'NA', 'NA')
				else:
					print line.strip()
			if chrom == 'scaffold-ma7':
				if int(start) >= 38000000 and int(end) <= 50000000:
					print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chrom, start, end, 'NA', 'NA', 'NA', 'NA')
				else:
					print line.strip()
			if chrom == 'scaffold-Z':
				if int(start) >= 60000000 and int(end) <= 73000000:
					print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chrom, start, end, 'NA', 'NA', 'NA', 'NA')
				else:
					print line.strip()
			