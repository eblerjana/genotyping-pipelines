#!/usr/bin/python

import sys
import argparse
import re
import pyfaidx
import gzip

def parse_info(fields):
	info_fields = { k.split('=')[0] : k.split('=')[1] for k in fields[7].split(';') if '=' in k }
	return info_fields


def insert_sequence(sequence, allele, offset, start, end):
	"""
	Given a sequence (reference) insert the given allele.
	"""

	included = sequence[0:(start-offset)] + allele + sequence[(end-offset):]
	offset -= len(included) - len(sequence)
	return included, offset

# biallelic VCF
biallelic = sys.argv[1]
# multiallelic VCF
multiallelic = sys.argv[2]
# FASTA with reference sequence
fastafile = sys.argv[3]

# read all biallelic variants and store
# var_id -> [ref_allele, alt_allele]
id_to_sequence = {}

nr_biallelic_ids = 0
biallelic_ids = set([])
nr_multiallelic_ids = 0
multiallelic_ids = set([])
nr_bubbles = 0

fasta = pyfaidx.Fasta(fastafile, as_raw=True, sequence_always_upper=True)

print("Read FASTA sequence.")

for line in gzip.open(biallelic, 'rt'):
	if line.startswith('#'):
		continue
	fields = line.split()
	alt_alleles = fields[4].split(',')
	assert len(alt_alleles) == 1
	info_fields = parse_info(fields)
	assert 'ID' in info_fields
	if info_fields['ID'] in id_to_sequence:
		raise Exception('ID ' + info_fields['ID'] + ' is duplicated in the biallelic VCF.')
	nr_biallelic_ids += 1
	biallelic_ids.add(info_fields['ID'])
	genotypes = [g.split('|') for g in fields[9:]]
	id_to_sequence[info_fields['ID']] = [fields[3], alt_alleles[0], genotypes]

print("Read biallelic VCF.")

# read multiallelic VCF and check whether allele sequences match IDs
for line in gzip.open(multiallelic, 'rt'):
	if line.startswith('#'):
		continue
	fields = line.split()
	info_fields = parse_info(fields)
	assert 'ID' in info_fields
	alt_alleles = [a.upper() for a in fields[4].split(',')]
	ids = info_fields['ID'].split(',')
	assert len(alt_alleles) == len(ids)
	start = int(fields[1])
	end = start + len(fields[3])

	nr_bubbles += 1

	# determine left- and right-most variant coordinates
	# from IDs (due to normalization there might be different
	# coordinates in multi/biallelic VCFs!
	unique_ids = set([])
	left_most = float('inf')
	right_most = 0
	for id in ids:
		for unique in id.split(':'):
			unique_ids.add(unique)
			multiallelic_ids.add(unique)
			id_start = int(unique.split('-')[1])
			left_most = min(left_most, id_start)
			right_most = max(right_most, id_start + len(id_to_sequence[unique][0]))
	nr_multiallelic_ids += len(unique_ids)

	if len(alt_alleles) != len(set(alt_alleles)):
		print('Duplicate ALT allele sequence at ' + fields[0] + fields[1])

	# generate alt allele by inserting IDs into reference
	# and compare to the alt_allele given in the multiallelic VCF.
	# Results should be identical, otherwise there is a conflict.
	for alt_allele, allele_id in zip(alt_alleles, ids):
		single_ids = allele_id.split(':')
		chr_name = fields[0].split('.')[-1]
#		expected_alt = fasta[chr_name][left_most -1 :  right_most-1]
		expected_alt = fasta[chr_name][int(fields[1]) : right_most-1] 
		offset = 0
		for id in single_ids:
			position = int(id.split('-')[1])
			allele_sequence = id_to_sequence[id][1]
			reference_sequence = id_to_sequence[id][0]
			vcf_pos = int(fields[1])
			if vcf_pos > position:
				# handle case where flanking base was added to DEL/INS
				assert vcf_pos-1 == position
				position = vcf_pos
				allele_sequence = id_to_sequence[id][1][1:]
				reference_sequence = id_to_sequence[id][0][1:]

			fasta_alt = fasta[chr_name][(position-1): (position-1 + len(reference_sequence))]
			if vcf_pos > position:
				raise Exception('Position of ID is before bubble position for allele ' + id  + ' (' + str(vcf_pos) + ')' )
			assert vcf_pos <= position
			index = position - vcf_pos
			if not reference_sequence == fasta_alt:
				raise Exception('Reference sequence for allele ' + id + ' does not match reference genome.' + fasta[chr_name][position-10:position+100])
			if id not in id_to_sequence:
				raise Exception('ID ' + id + ' found in multi-allelic VCF is not present in bi-allelic VCF.')
			if not expected_alt[(index-offset): (index + len(reference_sequence) - offset)] == reference_sequence:
				e_len = len(expected_alt[(index-offset): (index + len(reference_sequence) - offset)])
				r_len = len(reference_sequence)
				assert e_len < r_len
				expected = expected_alt[(index-offset): (index + len(reference_sequence) - offset)] + reference_sequence[e_len:r_len]
				if expected != reference_sequence:
					print(expected)
					print(reference_sequence)
					print('ID ' + id + ': reference allele does not match reference sequence at this position.')
			expected_alt, offset = insert_sequence(expected_alt, allele_sequence, offset, index, index + len(reference_sequence) )
		if alt_allele != expected_alt.upper()[:len(alt_allele)]:
			print(alt_allele)
			print(expected_alt)
			print('ID ' + allele_id + ' does not match the allele sequence in the multi-allelic VCF.')
	# check if genotypes in both VCFs match
	multi_genotypes = fields[9:]
	all_ids = set([])
	for allele_id in ids:
		for id in allele_id.split(':'):
			all_ids.add(id)
	all_ids = list(all_ids)
	for i,geno in enumerate(fields[9:]):
		alleles = geno.split('|')
		for j,allele in enumerate(alleles):
			if allele == '.':
				# make sure j-th haplotype is '.' for all these alleles
				for id in all_ids:
					if id_to_sequence[id][2][i][j] != '.':
						raise Exception('Genotype at position ' + fields[0] + ':' + fields[1] + ' for ID ' + id + ' and sample ' + str(i) + 'does not match.' )
				continue
			# look up all IDs on this allele
			allele_ids = ids[int(allele)-1].split(':') if allele != '0' else []
			# check if all alleles at this locus were correctly typed as absent/present
			for id in all_ids:
				expected_allele = '1' if id in allele_ids else '0'
				if id_to_sequence[id][2][i][j] != expected_allele:
					raise Exception('Genotype at position ' + fields[0] + ':' + fields[1] + ' for ID ' + id + ' and sample ' + str(i) + 'does not match.' )


if nr_biallelic_ids != nr_multiallelic_ids:
	raise Exception("Biallelic VCF contains IDs not present in the multiallelic version.")
print('VCFs are valid!')
print('number of variant alleles: ' + str(nr_biallelic_ids))
print('number of bubbles: ' + str(nr_bubbles))

