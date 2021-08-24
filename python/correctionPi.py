import sys

print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % ('chrom', 'start', 'end', 'mean_corr', 'p0.500_corr', 'p0.025_corr', 'p0.975_corr')

with open(sys.argv[1],'r') as rhopi:
	next(rhopi)
	for line in rhopi:
		chrom = line.split()[0]
		start = line.split()[1]
		end = line.split()[2]
		mean = line.split()[3]
		p50 = line.split()[4]
		p025 = line.split()[5]
		p975 = line.split()[6]
		pi = line.split()[7]

		if float(pi) > float(0):
			if mean != 'NA':
				mean_corr = float(mean)/float(pi)
				p50_corr = float(p50)/float(pi)
				p025_corr = float(p025)/float(pi)
				p975_corr = float(p975)/float(pi)
			else:
				mean_corr = 'NA'
				p50_corr = 'NA'
				p025_corr = 'NA'
				p975_corr = 'NA'
		else:
			mean_corr = mean
			p50_corr = p50
			p025_corr = p025
			p975_corr = p975

		print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chrom, start, end, mean_corr, p50_corr, p025_corr, p975_corr)
