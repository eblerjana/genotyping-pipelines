import sys
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import seaborn as sns

from variantclassifier import VariantType, determine_type_from_id

type_to_biallelic = defaultdict(lambda: 0)
type_to_multiallelic = defaultdict(lambda: 0)

for line in sys.stdin:
	if line.startswith('#'):
		continue
	fields = line.split()
	info_fields = { i.split('=')[0] : i.split('=')[1].strip() for i in fields[7].split(';') if '=' in i}
	assert 'ID' in info_fields
	var_ids = set([])
	for allele in info_fields['ID'].split(','):
		for varid in allele.split(':'):
			var_ids.add(varid)
	alleles = [fields[3]] + [a for a in fields[4].split(',')]
	nr_alleles = len(alleles)

	for varid in var_ids:
		vartype = determine_type_from_id(varid)
		if nr_alleles > 2:
			type_to_multiallelic[vartype] += 1
		else:
			type_to_biallelic[vartype] += 1
snps = [0,0]
indels = [0,0]
large_ins = [0,0]
large_del = [0,0]
large_compl = [0,0]

for vartype in [VariantType.snp, VariantType.small_insertion, VariantType.small_deletion, VariantType.small_complex, VariantType.midsize_insertion, VariantType.midsize_deletion, VariantType.midsize_complex, VariantType.large_insertion, VariantType.large_deletion, VariantType.large_complex]:
	if vartype == VariantType.snp:
		snps[0] += type_to_biallelic[vartype]
		snps[1] += type_to_multiallelic[vartype]
	elif vartype in [VariantType.small_insertion, VariantType.small_deletion, VariantType.small_complex, VariantType.midsize_insertion, VariantType.midsize_deletion, VariantType.midsize_complex]:
		indels[0] += type_to_biallelic[vartype]
		indels[1] += type_to_multiallelic[vartype]
	elif vartype == VariantType.large_insertion:
		large_ins[0] += type_to_biallelic[vartype]
		large_ins[1] += type_to_multiallelic[vartype]
	elif vartype == VariantType.large_deletion:
		large_del[0] += type_to_biallelic[vartype]
		large_del[1] += type_to_multiallelic[vartype]
	else:
		assert vartype == VariantType.large_complex
		large_compl[0] += type_to_biallelic[vartype]
		large_compl[1] += type_to_multiallelic[vartype]

print('vartype\t#alleles in biallelic regions\t#alleles in multiallelic regions')
print('\t'.join(['SNVs', str(snps[0]), str(snps[1])] ))
print('\t'.join(['indels', str(indels[0]), str(indels[1])] ))
print('\t'.join(['large deletions', str(large_del[0]), str(large_del[1])] ))
print('\t'.join(['large insertions', str(large_ins[0]), str(large_ins[1])] ))
print('\t'.join(['large complex', str(large_compl[0]), str(large_compl[1])] ))
