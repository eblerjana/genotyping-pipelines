import sys
import argparse
from collections import defaultdict
import gzip

class VariantStatistics:
	def __init__(self):
		self.confusion_matrix = [[0 for x in range(4)] for y in range(4)]
		self.varid = None
		self.correct = 0
		self.wrong = 0
		self.not_typed = 0
		self.absent_in_truth = 0
		self.total = 0

	def to_str(self):
		return '\t'.join([
			self.varid,
			str(self.total),
			str((self.correct / max(self.correct + self.wrong, 1)) * 100.0),
			str((self.wrong / max(self.correct + self.wrong, 1)) * 100.0),
			str((self.not_typed / max(self.correct + self.wrong + self.not_typed, 1)) * 100.0),
			str(self.correct),
			str(self.wrong),
			str(self.not_typed),
			str(self.absent_in_truth)		
		] + [str(self.confusion_matrix[i][j]) for i in range(3) for j in range(4)])
		


def convert_to_genotype(gt_string):
	if '.' in gt_string:
		return 3
	if '|' in gt_string:
		alleles = gt_string.split('|')
		return int(alleles[0]) + int(alleles[1])
	if '/' in gt_string:
		alleles = gt_string.split('/')
		return int(alleles[0]) + int(alleles[1])
	


def check_genotyping_accuracy(baseline_vcf, callset_vcf, outname, samples, column_prefix):
	"""
	Check genotyping accuracy across samples variant-wise.
	"""

	header_line = '\t'.join(['variant_id',
		column_prefix + '_considered_samples',
		column_prefix + '_correct [%]',
		column_prefix + '_wrong [%]',
		column_prefix + '_not_typed [%]',
		column_prefix + '_correct',
		column_prefix + '_wrong',
		column_prefix + '_not_typed',
		column_prefix + '_absent_in_truth',
		column_prefix + '_0/0_typed_0/0',
		column_prefix + '_0/0_typed_0/1',
		column_prefix + '_0/0_typed_1/1',
		column_prefix + '_0/0_typed_./.',
		column_prefix + '_0/1_typed_0/0',
		column_prefix + '_0/1_typed_0/1',
		column_prefix + '_0/1_typed_1/1',
		column_prefix + '_0/1_typed_./.',
		column_prefix + '_1/1_typed_0/0',
		column_prefix + '_1/1_typed_0/1',
		column_prefix + '_1/1_typed_1/1',
		column_prefix + '_1/1_typed_./.'
		]) + '\n'

	with open(outname, 'w') as outfile:

		outfile.write(header_line)
		baseline_samples = {}
		baseline_genotypes = defaultdict(list)

		# parse the baseline file
		for line in gzip.open(baseline_vcf, 'rt'):
			if line.startswith('##'):
				# header
				continue
			fields = line.strip().split()
			if line.startswith('#'):
				for i,s in enumerate(fields[9:]):
					baseline_samples[s] = i
				# check if all samples are present in baseline
				for s in samples:
					if not s in baseline_samples:
						raise RuntimeError("Requested sample " + s + " is not contained in baseline VCF.")
						sys.exit(1)
				continue
			# parse the genotypes
			genotypes = [None] * len(fields[9:])
			genotype_index = fields[8].split(':').index('GT') 
			for i,genotype in enumerate(fields[9:]):
				genotypes[i] = convert_to_genotype(genotype.split(':')[genotype_index])
			info_fields = { i.split('=')[0] : i.split('=')[1].strip() for i in fields[7].split(';') if '=' in i}
			assert 'ID' in info_fields
			var_id = info_fields['ID']
			baseline_genotypes[var_id] = genotypes


		# parse the callset file
		callset_samples = []
		for line in gzip.open(callset_vcf, 'rt'):
			if line.startswith('##'):
				# header
				continue
			fields = line.strip().split()
			if line.startswith('#'):
				for s in fields[9:]:
					callset_samples.append(s)
				# check if all samples are present in callset
				for s in samples:
					if not s in callset_samples:
						raise RuntimeError("Requested sample " + s + " is not contained in callset VCF.")
						sys.exit(1)
				continue
			info_fields = { i.split('=')[0] : i.split('=')[1].strip() for i in fields[7].split(';') if '=' in i}
			assert 'ID' in info_fields
			var_id = info_fields['ID']
			if var_id not in baseline_genotypes:
				raise RuntimeError("Callset VCF contains variant not present in the baseline set.")
				sys.exit(1)

			# compare the genotypes to baseline genotypes
			genotype_index = fields[8].split(':').index('GT')
			variant_stats = VariantStatistics()
			variant_stats.varid = var_id
			for sample,gt_string in zip(callset_samples, fields[9:]):
				if sample in samples:
					genotype = convert_to_genotype(gt_string.split(':')[genotype_index])
					truth_genotype = baseline_genotypes[var_id][baseline_samples[sample]]
					if truth_genotype != 3:
						variant_stats.confusion_matrix[truth_genotype][genotype] += 1
						if genotype == 3:
							variant_stats.not_typed += 1
						elif genotype == truth_genotype:
							variant_stats.correct += 1
						else:
							variant_stats.wrong += 1
					else:
						variant_stats.absent_in_truth += 1
					variant_stats.total += 1

			# print stats for variant to file
			outfile.write(variant_stats.to_str() + '\n')		
	


if __name__ == "__main__":

	# baseline: multisample VCF containing variants and panel of true genotypes
	# callset: multisample VCF containing genotype predictions per sample
	parser = argparse.ArgumentParser(prog='evaluate-self-genotyping.py', description=__doc__)
	parser.add_argument('baseline', metavar='BASELINE', help='multisample baseline VCF (ground truth genotypes).')
	parser.add_argument('callset', metavar='CALLSET', help='multisample callset VCF (genotyped variants).')
	parser.add_argument('outfile', metavar='OUTFILE', help='output file name.')
	parser.add_argument('samples', metavar='SAMPLES', help='TSV file with list of samples to consider.')
	parser.add_argument('column_prefix', metavar='PREFIX', help='prefix of the output column names.')
	args = parser.parse_args()

	samples = [line.strip() for line in open(args.samples, 'r')]
	check_genotyping_accuracy(args.baseline, args.callset, args.outfile, samples, args.column_prefix)

