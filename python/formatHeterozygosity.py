import sys

print "chrom"+"\t"+"start"+"\t"+"end"+"\t"+"hom1"+"\t"+"het"+"\t"+"hom2"+"\t"+"pval"+"\t"+"prop_het"

with open(sys.argv[1],'r') as input:
	next(input)
	for line in input:
		chrom = line.split()[0]
		start = int(line.split()[1]) - 1
		end = line.split()[1]
		pval = line.split()[5]
		var = line.split()[2]
		if var != '0/0/0':
			if float(pval) > 0.0001:
				hom1 = var.split('/')[0]
				het = var.split('/')[1]
				hom2 = var.split('/')[2]
				total = int(hom1) + int(het) + int(hom2)
				prophet = float(het)/float(total)
				if total >= int(sys.argv[2]):
					print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chrom,start,end,hom1,het,hom2,pval,prophet)
		
		
