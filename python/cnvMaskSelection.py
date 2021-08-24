import sys

# Get file name - splits path to take the name after the last '/' occurence.
file = sys.argv[1]
filename = str(file.rsplit('/',1)[1])

head = str(sys.argv[3])

cnvs = []
for line in open(sys.argv[2],'r'):
	chrom = line.split()[0]
	start = line.split()[1]
	end = line.split()[2]
	cnv = chrom+'_'+start+'_'+end
	cnvs.append(cnv)

# Parse file and check if windows overlap with CNVs
if '_pi' in filename:
	with open(file, 'r') as input:
		header = input.readlines()[0].strip()
		if head == 'yes':
			print header
	with open(file, 'r') as data:
		next(data)
		for line in data:
			pop = line.split()[0]
			chrom = line.split()[1]
			start = line.split()[2]
			end = line.split()[3]
			info = chrom+'_'+start+'_'+end
			if info in cnvs:
				print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (pop,chrom,start,end,'NA','NA','NA','NA','NA')
			else:
				print line.strip()

if '_dxy' in file:
	with open(file, 'r') as input:
		header = input.readlines()[0].strip()
		if head == 'yes':
			print header
	with open(file, 'r') as data:
		next(data)
		for line in data:
			pop1 = line.split()[0]
			pop2 = line.split()[1]
			chrom = line.split()[2]
			start = line.split()[3]
			end = line.split()[4]
			info = chrom+'_'+start+'_'+end
			if info in cnvs:
				print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (pop1,pop2,chrom,start,end,'NA','NA','NA','NA','NA')
			else:
				print line.strip()

if '_fst' in file:
	with open(file, 'r') as input:
		header = input.readlines()[0].strip()
		if head == 'yes':
			print header
	with open(file, 'r') as data:
		next(data)
		for line in data:
			pop1 = line.split()[0]
			pop2 = line.split()[1]
			chrom = line.split()[2]
			start = line.split()[3]
			end = line.split()[4]
			info = chrom+'_'+start+'_'+end
			if info in cnvs:
				print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (pop1,pop2,chrom,start,end,'NA','NA')
			else:
				print line.strip()

if 'Tajima.D' in filename:
	with open(file, 'r') as input:
		header = input.readlines()[0].strip()
		if head == 'yes':
			print header
	with open(file, 'r') as data:
		next(data)
		for line in data:
			chrom = line.split()[0]
			start = line.split()[1]
			if '10kb' in filename:
				end = int(start) + 10000
			if '1kb' in filename:
				end = int(start) + 1000
			if '250bp' in filename:
				end = int(start) + 250
			info = chrom+'_'+str(start)+'_'+str(end)
			if info in cnvs:
				print "%s\t%s\t%s\t%s" % (chrom,start,'NA','NA')
			else:
				print line.strip()

if 'df_prop' in filename:
	with open(file, 'r') as input:
		header = input.readlines()[0].strip()
		if head == 'yes':
			print header
	with open(file, 'r') as data:
		next(data)
		for line in data:
			chrom = line.split()[0]
			start = line.split()[1]
			end = line.split()[2]
			info = chrom+'_'+str(start)+'_'+str(end)
			if info in cnvs:
				print "%s\t%s\t%s\t%s\t%s\t%s" % (chrom,start,end,'NA','NA','NA')
			else:
				print line.strip()

if '_ihs' in filename:
	with open(file, 'r') as input:
		header = input.readlines()[0].strip()
		if head == 'yes':
			print header
	with open(file, 'r') as data:
		next(data)
		for line in data:
			chrom = line.split()[0]
			start = line.split()[1]
			end = line.split()[2]
			info = chrom+'_'+str(start)+'_'+str(end)
			if info in cnvs:
				print "%s\t%s\t%s\t%s\t%s" % (chrom,start,end,'NA','NA')
			else:
				print line.strip()

if 'beta' in filename:
	with open(file, 'r') as input:
		header = input.readlines()[0].strip()
		if head == 'yes':
			print header
	with open(file, 'r') as data:
		next(data)
		for line in data:
			chrom = line.split()[0]
			start = line.split()[1]
			end = line.split()[2]
			info = chrom+'_'+str(start)+'_'+str(end)
			if info in cnvs:
				print "%s\t%s\t%s\t%s" % (chrom,start,end,'NA')
			else:
				print line.strip()

if 'het.' in filename:
	with open(file, 'r') as input:
		header = input.readlines()[0].strip()
		if head == 'yes':
			print header
	with open(file, 'r') as data:
		next(data)
		for line in data:
			chrom = line.split()[0]
			start = line.split()[1]
			end = line.split()[2]
			info = chrom+'_'+str(start)+'_'+str(end)
			if info in cnvs:
				print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chrom,start,end,'NA','NA','NA','NA')
			else:
				print line.strip()
