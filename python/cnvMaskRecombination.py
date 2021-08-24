"""
usage: python cnv_mask.py <cnv-overlap.bed> <input.txt[assumes header present!]>
"""
import sys

cnvs = []

for cnv in open(sys.argv[1],'r'):
	chrom = cnv.split()[0]
	start = str(cnv.split()[1])
	end = str(cnv.split()[2])
	dat = chrom+'_'+start+'_'+end
	cnvs.append(dat)

with open(sys.argv[2],'r') as rho:
	header = next(rho)
	print header.strip()
	for line in rho:
		chrom = line.split()[0]
		start = line.split()[1]
		end = line.split()[2]
		mean = line.split()[3]
		med = line.split()[4]
		low = line.split()[5]
		high = line.split()[6]

		info = chrom+'_'+str(start)+'_'+str(end)
		
		if info in cnvs:
			print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chrom,start,end,'NA','NA','NA','NA')
		else:
			print line.strip()
