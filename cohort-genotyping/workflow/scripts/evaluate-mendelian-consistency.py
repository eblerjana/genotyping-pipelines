import sys
import argparse
import gzip
from variantclassifier import VariantType, determine_variant_from_line
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np 

def check_mendelian_consistency(child_gt, parent1_gt, parent2_gt):
	child = child_gt.get_alleles()
	parent1 = parent1_gt.get_alleles()
	parent2 = parent2_gt.get_alleles()
	if child[0] in parent1 and child[1] in parent2:
		return True
	if child[0] in parent2 and child[1] in parent1:
		return True
	return False

class Genotype:
	def __init__(self, alleles, is_phased):
		self._alleles = alleles
		self._phased = is_phased

	def __str__(self):
		if self._phased:
			return '|'.join([str(i) for i in self._alleles])
		else:
			return '/'.join([str(i) for i in self._alleles])

	def get_alleles(self):
		return self._alleles

	def get_ploidy(self):
		return len(self._alleles)

	def is_phased(self):
		return self._phased

	def __eq__(self, other):
		return sorted(self._alleles) == sorted(other._alleles)

	def is_hom_ref(self):
		return all([int(a) == 0 for a in self._alleles])

	def is_none(self):
		return self._alleles == []

def genotype_from_string(gt_string):
	is_phased = False
	alleles = []
	if '.' in gt_string:
		# untyped
		return Genotype(alleles, is_phased)

	if '|' in gt_string:
		is_phased = True
		alleles = [int(allele) for allele in gt_string.split('|')]
	elif '/' in gt_string:
		is_phased = False
		alleles = [int(allele) for allele in gt_string.split('/')]
	else:
		assert False
	return Genotype(alleles, is_phased)

def parse_trios(ped_file, samples):
	samples_to_include = set()
	with open(samples, 'r') as listing:
		for line in listing:
			samples_to_include.add(line.strip())
	trios = {}
	for line in open(ped_file, 'r'):
		if line.startswith('#'):
			continue
		fields = line.split()
		# map trio_name -> [child, parent1, parent2, nr consistent, nr_inconsistent, nr_untyped, nr total]
		if any([s not in samples_to_include for s in fields[1:4]]):
			continue
		trios[(fields[1], fields[2], fields[3])] = [fields[1], fields[2], fields[3], 0, 0, 0, 0]
	return trios


def remove_samples_from_header(header_line, samples):
	fields = header_line.split()
	result = fields[:9]
	for s in fields[9:]:
		if s not in samples:
			result.append(s)
	return '\t'.join(result)


def run_statistics(vcf_file, ped_file, samples, table_tsv, column_prefix):
	# map child -> [parents]
	trios = parse_trios(ped_file, samples)
	# histograms[vartype][i] = number of variants consistent in i trios
	histograms = { vartype.name : [0]*(len(trios)+1) for vartype in VariantType }
	histograms['multi_var'] = [0] * (len(trios)+1)
	untyped_vars = 0
	multitype_vars = 0

	with open(table_tsv, 'w') as out_tsv:
		header = '\t'.join(['variant_id', column_prefix + '_mendelian_consistent_trios', column_prefix + '_alternative_transmitted', column_prefix + '_considered_trios', column_prefix + '_all_0/0', column_prefix + '_all_0/1', column_prefix + '_all_1/1', column_prefix + '_untyped_alleles_present'])
		out_tsv.write(header + '\n')
		header = None
		sample_to_index = {}
		for record in gzip.open(vcf_file, 'rt'):
			if record.startswith('##'):
				continue
			if record.startswith('#'):
				header = record.strip().split()
				for i,f in enumerate(header):
					sample_to_index[f] = i
				continue
			assert header is not None
			fields = record.split()
			genotype_index = fields[8].split(':').index('GT')
			consistent_trios = 0
			total_consistent_trios = 0
			alt_transmitted = 0
			all_het = 0
			all_abs = 0
			all_present = 0
			assert len(fields[4].split(',')) == 1
			info_fields = { i.split('=')[0] : i.split('=')[1].strip() for i in fields[7].split(';') if '=' in i} 
			assert 'ID' in info_fields
			assert len(info_fields['ID'].split(',')) == 1
			variant_id = info_fields['ID']
			total_trios = 0
			untyped_trios = 0
			if 'X' in fields[0] or 'Y' in fields[0]:
				continue
			vartype = determine_variant_from_line(record)
			assert isinstance(vartype, VariantType), 'Unexpected return value: {}'.format(type(vartype))
			untyped = False
			for name in trios:
				field_child = fields[sample_to_index[trios[name][0]]].split(':')
				field_parent1 = fields[sample_to_index[trios[name][1]]].split(':')
				field_parent2 = fields[sample_to_index[trios[name][2]]].split(':')
				gt_child = genotype_from_string(field_child[genotype_index])
				gt_parent1 = genotype_from_string(field_parent1[genotype_index])
				gt_parent2 = genotype_from_string(field_parent2[genotype_index])
				if any([g.is_none() for g in [gt_child, gt_parent1, gt_parent2]]):
					untyped_trios += 1
					if not untyped:
						untyped = True
						untyped_vars += 1
				elif all([ g == gt_child for g in [gt_parent1, gt_parent2]]):
					# all genotypes same, automatically consistent
					if gt_child == Genotype([0,0], False):
						all_abs += 1
					elif gt_child == Genotype([0,1], False):
						all_het += 1
					elif gt_child == Genotype([1,1], False):
						all_present += 1
					else:
						assert(False)
					total_consistent_trios += 1
				else:
					total_trios += 1
					consistent = check_mendelian_consistency(gt_child, gt_parent1, gt_parent2)
					if consistent:
						consistent_trios += 1
						total_consistent_trios += 1
						# check how often alt allele was transmitted
						alt_transmitted += sum(a!=0 for a in gt_child.get_alleles())
			assert total_trios + all_abs + all_het + all_present + untyped_trios == len(trios)
			if not untyped:
				histograms[vartype.name][total_consistent_trios] += 1
			# write to output
			out_tsv.write('\t'.join([variant_id, str(consistent_trios), str(alt_transmitted), str(total_trios), str(all_abs), str(all_het), str(all_present), str(untyped_trios)]) + '\n')
	# output the results
	print('variants_untyped\t' + str(untyped_vars))
	print('\t'.join(['nr_trios'] + [str(i) for i in range(0,len(trios)+1)]))
	for vartype in VariantType:
		print('\t'.join([vartype.name] + [str(i) for i in histograms[vartype.name]]))
		

def run_plot_statistics(outprefix, tsv):
	type_to_numbers = {}
	nr_trios = -1
	for line in open(tsv, 'r'):
		if line.startswith('variants_untyped') or line.startswith('nr_trios') or line.startswith('variants_multitype'):
			continue
		fields = line.split()
		nr_trios = len(fields)-1
		type_to_numbers[fields[0]] = [int(i) for i in fields[1:]]
	all_types = [0] * nr_trios
	x_values = [i for i in range(nr_trios)]
	for var, numbers in type_to_numbers.items():
		assert len(numbers) == nr_trios
		plt.bar(x_values, numbers)
		plt.title(var)
		plt.ylabel('count')
		plt.yscale('log')
		plt.xlabel('nr of mendelian consistent trios')
#		plt.xticks(x_values)
		plt.savefig(outprefix + '_statistics_' + var + '.pdf')
		plt.close()
		for i,number in enumerate(numbers):
			all_types[i] += number
	plt.bar(x_values, all_types)
	plt.title('all variants')
	plt.ylabel('count')
	plt.yscale('log')
	plt.xlabel('nr of mendelian consistent trios')
	plt.xticks(x_values)
	plt.savefig(outprefix + '_statistics_all.pdf')
	plt.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='mendelian-consistency.py', description=__doc__)
	subparsers = parser.add_subparsers(dest="subparser_name")

	parser_statistics = subparsers.add_parser('statistics', help='print statistics.')
	parser_statistics.add_argument('-vcf', metavar='VCF', required=True, help='Multisample VCF-file with genotypes.')
	parser_statistics.add_argument('-ped', metavar='PED', required=True, help='Trio relationships to be considered.')
	parser_statistics.add_argument('-samples', metavar='SAMPLES', required=True, help='Samples to include')
	parser_statistics.add_argument('-table', metavar='TABLE', required=True, help='Write statistics per variant ID.')
	parser_statistics.add_argument('-column-prefix', metavar='COLUMNPREFIX', required=True, help='Column prefix in output.')

	parser_plot = subparsers.add_parser('plot', help='create plots from outputs of commands statistics and filter.')
	parser_plot.add_argument('-statistics', required=True, metavar='TSV', help='tsv produced by command statistics.')
	parser_plot.add_argument('outprefix', metavar='OUTPREFIX', help='prefix of the output files.')
	args = parser.parse_args()

	if args.subparser_name == 'statistics':
		run_statistics(args.vcf, args.ped, args.samples, args.table, args.column_prefix)
	if args.subparser_name == 'plot':
		run_plot_filter(args.outprefix, args.statistics)
