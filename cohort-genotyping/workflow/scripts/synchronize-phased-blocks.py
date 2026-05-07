import sys
from collections import defaultdict
import gzip

vcf = sys.argv[1]
block_to_stats = defaultdict(lambda: [0,0])

# read VCF once and record the number of 0s and 1s on the first haplotype
for line in gzip.open(vcf, 'rt'):
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	if not 'PS' in fields[8]:
		continue
	gt_index = fields[8].split(':').index('GT')
	ps_index = fields[8].split(':').index('PS')
	gt_string = fields[9].split(':')[gt_index]
	ps_string = fields[9].split(':')[ps_index]
	if gt_string == "0|1":
		block_to_stats[ps_string][0] += 1
	elif gt_string == "1|0":
		block_to_stats[ps_string][1] += 1
	else:
		assert False


# read the VCF again and synchronize / filter phased blocks
for line in gzip.open(vcf, 'rt'):
	if line.startswith('#'):
		print(line.strip())
		continue
	fields = line.strip().split()
	if not 'PS' in fields[8]:
		sys.stderr.write("Skipping " + fields[0] + ":" + fields[1] + " since PS tag is missing.\n")
		print(line.strip())
		continue	
	gt_index = fields[8].split(':').index('GT')
	ps_index = fields[8].split(':').index('PS')
	gt_string = fields[9].split(':')[gt_index]
	ps_string = fields[9].split(':')[ps_index]

	# check whether to consider the phased block
	if block_to_stats[ps_string][0] > block_to_stats[ps_string][1]:
		# don't flip genotypes
		print(line.strip())
	elif block_to_stats[ps_string][1] > block_to_stats[ps_string][0]:
		# flip genotypes
		assert gt_string in ["1|0", "0|1"]
		new_genotypes = fields[9].split(':')
		new_genotypes[gt_index] = "1|0" if gt_string == "0|1" else "0|1"
		fields[9] = ':'.join(new_genotypes)
		print('\t'.join(fields))
	else:
		sys.stderr.write("Skipping " + fields[0] + ":" + fields[1] + " since it's unclear which haplotype is the dominant one.\n")

