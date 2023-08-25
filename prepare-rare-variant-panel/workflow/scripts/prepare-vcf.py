import sys
import argparse


# filter out symbolic records
# filter out multi-allelic records
# filter out overlapping records

def process_vcf():
	total_records = 0
	skipped_records = 0
	current_end = 0
	allele_index = 0
	for line in sys.stdin:
		if line.startswith('##'):
			# header line
			continue
		# skip genotypes
		fields = line.strip().split()[:9]
		if line.startswith('#'):
			print('\t'.join(fields + ['rare-alleles']))
		# make sure that the record is not symbolic
		if ('<' in fields[3]) or ('<' in fields[4]):
			skipped_records += 1
			sys.stderr.write("Skipping record at position " + fields[0] + ':' + fields[1] + " because it is symbolic.\n")
		alleles = [a for a in fields[4].split(',')]
		if len(alleles) > 1:
			# multiallelic variant
			skipped_records += 1
			sys.stderr.write("Skipping record at position " + fields[0] + ':' + fields[1] + " because it is multiallelic.\n")
		start = int(fields[1])
		end = start + len(fields[3]))
		# check if current record overlaps the previous one
		if start < current_end:
			skipped_records += 1
			sys.stderr.write("Skipping record at position " + fields[0] + ':' + fields[1] + " because it overlaps the previous one.\n")
			continue
		current_end = end
		total_records += 1

		# add IDs to all alleles
		info_string = fields[7]
		info_fields = { i.split('=')[0] : i.split('=')[1] for i in info_string.split(';') if '=' in i}
		chrom = fields[0]
		var_len = len(alleles[0])
		var_type = 'rare' + str(allele_index)
		var_id = '-'.join([chrom, start, var_type, str(var_len)])
		allele_index += 1
		info_fields['ID'] = var_id
		updated_info = []
		for k,v in info_fields.items():
			updated_info.append(k + '=' + v)
		fields[7] = ';'.join(updated_info)
		# if there is additional info present in sample field, remove it
		fields[8] = 'GT'		
		print('\t'.join(fields + ['1']))

	sys.stderr("Wrote " + str(total_records - skipped_records) + "/" + str(total_records) + " to output VCF.")


if __name__ == "__main__":
	process_vcf()
