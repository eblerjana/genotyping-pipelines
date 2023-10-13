#!/usr/bin/python

import sys
import argparse
import re
import pyfaidx
from collections import defaultdict


class VcfRecord:
	"""
	Represents a VCF record.
	"""
	def __init__(self, line, sample_names):
		fields = line.strip().split()
		self._chrom = fields[0]
		self._pos = int(fields[1])
		self._id = fields[2]
		self._ref = fields[3]
		alt_alleles = fields[4].split(',')
		self._alt = alt_alleles
		self._qual = fields[5]
		self._filter = fields[6]
		self._info = fields[7]
		self._format = fields[8]
		if len(sample_names) != len(fields[9:]):
			raise RuntimeError("VcfRecord: number of sample names does not match the number of sample columns in the VCF line.")
		self._sample_to_field = {}
		for name, data in zip(sample_names, fields[9:]):
			self._sample_to_field[name] = data


	def __eq__(self, other):
		if self._chrom != other._chrom:
			return False
		if self._pos != other._pos:
			return False
		if self._id != other._id:
			return False
		if self._ref != other._ref:
			return False
		if self._alt != other._alt:
			return False
		if self._qual  != other._qual:
			return False
		if self._filter != other._filter:
			return False
		if self._info != other._info:
			return False
		if self._format != other._format:
			return False
		if self._sample_to_field != other._sample_to_field:
			return False
		return True

	def __lt__(self, other):
		if self._chrom != other._chrom:
			return self._chrom < other._chrom
		else:
			return self._pos < other._pos

	def chrom(self):
		return self._chrom

	def pos(self):
		return self._pos

	def id(self):
		return self._id

	def ref(self):
		return self._ref

	def alt(self):
		return self._alt

	def same_variant_location(self, other):
		"""
		Return true, if the variant is representing the same
		event as variant "other". This is the case if the
		position, ref allele and alt alleles are identical.
		"""

		if self._chrom != other._chrom:
			return False
		
		if self._pos != other._pos:
			return False

		if self._ref != other._ref:
			return False
		
		alt_alleles = set(self._alt)
		other_alt_alleles = set(other._alt)
		if alt_alleles != other_alt_alleles:
			return False
		return True

	def combine_variants(self, other):
		"""
		Combine two variants into one record. The variants
		need to match in terms of position and reference allele.
		"""
		
		# make sure that the records do not have overlapping samples
		if bool(set([s for s in self._sample_to_field]) & set([s for s in other._sample_to_field])):
			raise RuntimeError('VcfRecord: combine_variants: variants can only be combined if the sample columns are disjoint.')

		if not self.same_variant_location(other):
			raise RuntimeError('VcfRecord: combine_variants: variants can only be combined if they represent the same allele.')
		# determine indices of matching allele(s) between both variants
		new_indices = {}
		for i,b in enumerate(other._alt):
			if b in self._alt:
				index = self._alt.index(b)
				new_indices[i+1] = index + 1
			else:
				self._alt.append(b)
				new_indices[i+1] = len(self._alt)
		# add sample columns and update genotypes
		# index of GT field within samples
		gt_index = other._format.split(':').index('GT')
		for name, data in other._sample_to_field.items():
			genotype = data.split(':')[gt_index]
			updated_genotype = []
			separator = '|' if '|' in genotype else '/'

			for a in genotype.split(separator):
				if int(a) in new_indices:
					updated_genotype.append(str(new_indices[int(a)]))
				else:
					updated_genotype.append(a)
			updated_sample_field = data.split(':')
			updated_sample_field[gt_index] = separator.join(updated_genotype)
			self._sample_to_field[name] = ':'.join(updated_sample_field)


	def str_record(self, sample_names = [], is_phased = False):
		sample_columns = list(set([s for s in self._sample_to_field.keys()] + sample_names))
		outline = [self._chrom, str(self._pos), self._id, self._ref, ','.join(self._alt), self._qual, self._filter, self._info, self._format]
		for sample in sorted(sample_columns):
			if sample in self._sample_to_field:
				outline.append(self._sample_to_field[sample])
			else:
				outline.append("0|0" if is_phased else "0/0")
		return '\t'.join(outline)


def find_overlaps(variants):
	"""
	Finds records sharing position and
	REF allele. These are the ones that
	are to be combined into one record later.
	"""
	assert len(variants) > 1
	current_position = variants[0].pos()
	current_chrom = variants[0].chrom()
	current_cluster = [variants[0]]
	for v in variants[1:]:
		if (v.pos() != current_position) or (v.chrom() != current_chrom) or (v.ref() != current.ref()):
			yield current_cluster
			current_cluster = []
		current_chrom = v.chrom()
		current_cluster.append(v)
		current_pos = v.pos()
	if len(current_cluster) != 0:
		yield current_cluster


def create_combined_cluster(cluster):
	"""
	Given records that share locations and
	REF alleles, combined them into one single
	record.
	"""
	assert len(cluster) > 1
	combined_record = cluster[0]
	for c in cluster[1:]:
		combined_record.combine_variants(c)
	return c


def run_combine_vcfs(vcfs):
	vcf_to_samples = {}
	variants = []
	for vcf in vcfs:
		for line in open(vcf, 'r'):
			if line.startswith('##'):
				continue
			if line.startswith('#'):
				# extract list of variants
				samples = line.strip().split('\t')[9:]
				vcf_to_samples[vcf] = samples
				continue
			assert vcf in vcf_to_samples
			variant = VcfRecord(line, vcf_to_samples[vcf])
			variants.append(variant)
	variants.sort()
	print(variants)
	for cluster in find_overlaps(variants):
		combined_record = create_combined_record(cluster)
		
			
# TODOs:
# write header of VCF properly
# remove information other than the GT field (write header accordingly)
# when and when not to combine variants??
# Plan: 
#	restrict script to biallelic VCFs only --> first check if my merging script is compartible with that
#	require identical alleles for variants to be the same
#	sorting also according to alleles, so that identical records are always ascending
	 
		
	

if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='combine_vcfs.py', description=__doc__)
	parser.add_argument('-vcfs', metavar='VCFs', required=True, nargs='+', help='VCF files to be combined (at least two files required).')
	args = parser.parse_args()
	
	if len(args.vcfs) < 2:
		raise RuntimeError('At least two VCFs are required.')
	run_combine_vcfs(args.vcfs)
