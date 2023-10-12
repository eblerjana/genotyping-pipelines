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
	def __init__(self, line, ):
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
		self._samples = fields[9:]

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
		if self._samples != other._samples:
			return False
		return True

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
		position and ref allele are identical.
		"""

		if self._chrom != other._chrom:
			return False
		
		if self._pos != other._pos:
			return False

		if self._ref != other._ref:
			return False
		
		return True

	def combine_variants(self, other):
		"""
		Combine two variants into one record. The variants
		need to match in terms of position and reference allele.
		"""

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
		for i in range(len(other._samples)):
			genotype = other._samples[i].split(':')[gt_index]
			updated_genotype = []
			separator = '|' if '|' in genotype else '/'

			for a in genotype.split(separator):
				if int(a) in new_indices:
					updated_genotype.append(str(new_indices[int(a)]))
				else:
					updated_genotype.append(a)
			updated_sample_field = other._samples[i].split(':')
			updated_sample_field[gt_index] = separator.join(updated_genotype)
			self._samples.append(':'.join(updated_sample_field))
	
	def add_samples(self, other_samples):
		"""
		Add sample information as additional columns to
		the VCF record.
		"""
		self._samples += other_samples
	
	def __str__(self):
		outline = '\t'.join([self._chrom, str(self._pos), self._id, self._ref, ','.join(self._alt), self._qual, self._filter, self._info, self._format, '\t'.join(self._samples)])
		return outline
		
		
# TODO: in merging, how to know which samples need to be set to 0/0 for missing variants?
# maybe keep track in variant record which genotypes correspond to which samples / indices in final VCF sample columns?
def run_combine_vcfs(vcfs):
	vcf_to_samples = {}
	variants = 
	for vcf in vcfs:
		for line in open(vcf, 'r'):
			if line.startswith('##'):
				continue
			if line.startswith('#'):
				# extract list of variants
				samples = line.strip().split('\t')[9:]
				vcf_to_samples[vcf] = samples
				continue
			
			 
		
	

if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='combine_vcfs.py', description=__doc__)
	parser.add_argument('-vcfs', metavar='VCFs', required=True, nargs='+', help='VCF files to be combined (at least two files required).')
	args = parser.parse_args()
	
	if len(args.vcfs) < 2:
		raise RuntimeError('At least two VCFs are required.')
	run_combine_vcfs(args.vcfs)
