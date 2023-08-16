import sys

genome_index = sys.argv[1]

chromosomes_hg38 = ['chr' + str(i) for i in range(1,23)] + ['chrX', 'chrY']
chromosomes_hg37 = [str(i) for i in range(1,23)] + ['X', 'Y']

chrom_to_length = {}
for line in open(genome_index, 'r'):
	fields = line.split()
	if fields[0] in chromosomes_hg38 or fields[0] in chromosomes_hg37:
		chrom_to_length[fields[0]] = fields[1]

counter = 0

no_contigs = True
for line in sys.stdin:
	if line.startswith('##'):
		if line.startswith('##contig=<'):
			continue
		print(line.strip())
		continue
	if line.startswith('#'):
		if no_contigs:
			for chrom,length in chrom_to_length.items():
				print("##contig=<ID=" + chrom + ",length=" + length + ">")
		print(line.strip())
		continue
	fields = line.split()
	if fields[0] not in chromosomes_hg38 and fields[0] not in chromosomes_hg37:
		continue
	print(line.strip())
