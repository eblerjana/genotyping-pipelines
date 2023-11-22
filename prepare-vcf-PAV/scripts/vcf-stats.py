import sys
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(prog='vcf-stats.py', description=__doc__)
parser.add_argument('--sample-info', metavar='SAMPLEINFO', help="gender information for vcf samples.")
args = parser.parse_args()


class Stats:
	def __init__(self, name):
		self.name = name
		self.number_of_records = 0
		self.number_of_ids = 0
		self.sample_to_ac = defaultdict(lambda: 0)

	def print(self, samples, info = None):
		print('Statistics for ' + self.name)
		print('Number of Records:\t' + str(self.number_of_records))
		print('Number of IDs:\t' + str(self.number_of_ids))
		print('Number of non-0/0 genotypes per sample:')
		
		
		if info:
			sample_to_info = {}
			for line in open(info, 'r'):
				fields = line.strip().split()
				sample_to_info[fields[0]] = fields[1]
		for sample in samples:
			if info:
				print(sample + ' (' + sample_to_info[sample] + ')' + ':\t' + str(self.sample_to_ac[sample]))
			else:
				print(sample + ':\t' + str(self.sample_to_ac[sample]))


current_chrom = ''
stats = {'genome': Stats('genome')}
samples = []

for line in sys.stdin:
	if line.startswith('##'):
		continue
	if line.startswith('#'):
		fields = line.strip().split()
		samples = fields[9:]
		continue
	fields = line.strip().split()
	if fields[0] != current_chrom:
		current_chrom = fields[0]
		stats[current_chrom] = Stats(current_chrom)
	
	stats[current_chrom].number_of_records += 1
	stats['genome'].number_of_records += 1
	
	info_field = {f.split('=')[0] : f.split('=')[1] for f in fields[7].split(';') if '=' in f}
	assert 'ID' in info_field

	ids = set([])
	for allele in info_field['ID'].split(','):
		for id in allele.split(':'):
			ids.add(id)
	stats[current_chrom].number_of_ids += len(ids)
	stats['genome'].number_of_ids += len(ids)
	
	gt_index = fields[8].split(':').index('GT')
	for sample, geno in zip(samples, fields[9:]):
		g = geno.split(':')[gt_index]
		assert '|' in g
		genotype = g.split('|')
		if (genotype[0] not in ['0', '.']) or (genotype[1] not in ['0', '.']):
			stats[current_chrom].sample_to_ac[sample] += 1
			stats['genome'].sample_to_ac[sample] += 1

# print out statistics
for r in sorted(stats.keys()):
	stats[r].print(samples,args.sample_info)
	print('')
	
