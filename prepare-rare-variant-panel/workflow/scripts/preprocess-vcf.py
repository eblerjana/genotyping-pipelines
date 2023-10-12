import sys


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


nr_alleles = 0
skipped_alleles = 0

for line in sys.stdin:
	if line.startswith("##"):
		print(line.strip())
		continue
	fields = line.strip().split()
	if line.startswith("#"):
		print('##INFO=<ID=ID,Number=A,Type=String,Description=\"Variant IDs per ALT allele.\">')
		# remove sample information and create pseudo sample carrying all alleles
		print("\t".join(fields[:9] + ["rare-variants"]))
		continue
	alleles = [fields[3]] + [a for a in fields[4].split(',')]
	allele_id = "-".join([fields[0], fields[1], determine_class_from_line(line), "rare" + str(nr_alleles), str(max(1, max([len(a) for a in alleles])-1))])
	format_field = { f.split("=")[0] : f.split("=")[1] for f in fields[7].split(";") if "=" in f}
	
	# make sure the VCF is biallelic
	assert len(alleles) == 2

	# if MC VCF, make sure to only keep LV=0 variants
	if "LV" in format_field:
		if int(format_field["LV"]) > 0:
			skipped_alleles += 1
			continue
	
	# add ID for the allele. Overwrites existing ID fields that
	# might have been present before
	format_field["ID"] = allele_id
	fields[7] = ";".join([k + "=" + v for k,v in format_field.items()])
	print("\t".join(fields[:9] + ["1"]))
	nr_alleles += 1

sys.stderr.write("Wrote " + str(nr_alleles) + " to the output VCF.\n")
sys.stderr.write("Skipped " + str(skipped_alleles) + " from the input VCF.\n")
