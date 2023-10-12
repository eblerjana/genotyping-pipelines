from enum import Enum
import vcf

class VariantType(Enum):
	snp = 0
	small_insertion = 1
	small_deletion = 2
	small_complex = 3
	midsize_insertion = 4
	midsize_deletion = 5
	midsize_complex = 6
	large_insertion = 7
	large_deletion = 8
	large_complex = 9


def determine_type_from_ids(ids):
	ids = set(ids)
	results = []
	for var in ids:
		var_type = var.split('-')[2]
		var_len = 1 if var_type == "SNV" else int(var.split('-')[-1])
		if var_type == "SNV":
			results.append(VariantType.snp)
			continue
		if var_len < 20:
			if var_type == "INS":
				results.append(VariantType.small_insertion)
			elif var_type == "DEL":
				results.append(VariantType.small_deletion)
			else:
				results.append(VariantType.small_complex)
			continue

		if var_len >= 20 and var_len < 50:
			if var_type == "INS":
				results.append(VariantType.midsize_insertion)
			elif var_type == "DEL":
				results.append(VariantType.midsize_deletion)
			else:
				results.append(VariantType.midsize_complex)
			continue

		if var_len >= 50:
			if var_type == "INS":
				results.append(VariantType.large_insertion)
			elif var_type == "DEL":
				results.append(VariantType.large_deletion)
			else:
				results.append(VariantType.large_complex)
		
	return results


def determine_type_from_record(record):
	"""
	Determines the variant type.
	"""
	if 'ID' in record.INFO:
		allele_ids = record.INFO['ID']
		# handle merged IDs
		all_ids = []
		for i in record.INFO['ID']:
			for j in i.split(':'):
				all_ids.append(j)

		return determine_type_from_ids(all_ids)
	else:
		alleles = [record.REF] + record.ALT
		varlen = max([len(a) for a in alleles])

		if record.is_snp:
			return [VariantType.snp]

		is_deletion = record.var_subtype == 'del' or record.var_subtype == 'DEL'
		is_insertion = record.var_subtype == 'ins' or record.var_subtype == 'INS'

		if '<' in str(record.ALT[0]):
			is_deletion = str(record.ALT[0])[1:-1] == 'DEL'
			is_insertion = str(record.ALT[0])[1:-1] == 'INS'	

		if varlen < 20:
			if is_insertion:
				return [VariantType.small_insertion]
			if is_deletion:
				return [VariantType.small_deletion]
			return [VariantType.small_complex]

		if varlen >= 20 and varlen <= 50:
			if is_insertion:
				return [VariantType.midsize_insertion]
			if is_deletion:
				return [VariantType.midsize_deletion]
			return [VariantType.midsize_complex]

		if varlen > 50:
			if is_insertion:
				return [VariantType.large_insertion]
			if is_deletion:
				return [VariantType.large_deletion]
			return [VariantType.large_complex]

def determine_type_from_info(info_string):
	info_fields = { i.split('=')[0] : i.split('=')[1] for i in info_string.split(';') if "=" in i}
	variant_ids = []
	for v in info_fields['ID'].split(','):
		for n in v.split(':'):
			variant_ids.append(n)
	return determine_type_from_ids(variant_ids)


def determine_class_from_line(line):
	"""
	Determine variant class (SNV, INS, DEL)
	from VCF line.
	"""
	splitted = line.split()
	alleles = [splitted[3]] + [s for s in splitted[4].split(',')]

	if all([len(a) == 1 for a in alleles]):
		return 'SNV'

	if len(alleles) > 2:
		return 'COMPLEX'

	is_deletion = (len(alleles[0]) > 1) and (len(alleles[1]) == 1)
	is_insertion = (len(alleles[0]) == 1) and (len(alleles[1]) > 1)

	if is_deletion:
		return 'DEL'
	elif is_insertion:
		return 'INS'
	else:
		return 'COMPLEX'
