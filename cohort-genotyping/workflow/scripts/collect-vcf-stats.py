import sys
import argparse


def collect_stats(genotype_stats, outname, column_prefix):


	with open(outname, 'w') as outfile:

		header = [
			"variant_id",
			column_prefix + '_allele_freq',
			column_prefix + '_alternative_alleles',
			column_prefix + '_total_alleles',
			column_prefix + '_unknown_alleles'
		]

		if genotype_stats:
			header.extend([
				column_prefix + '_heterozygosity',
				column_prefix + '_heterozygous_genotypes',
				column_prefix + '_total_genotypes',
				column_prefix + '_unique_kmers',
				column_prefix + '_GQ>=200'
			])
		outfile.write('\t'.join(header) + '\n')

		nr_samples = 0
		for line in sys.stdin:
			if line.startswith("##"):
				continue
			fields = line.strip().split()
			if line.startswith("#"):
				nr_samples = len(fields[9:])
				continue
			info_fields = { i.split('=')[0] : i.split('=')[1].strip() for i in fields[7].split(';') if '=' in i}

			allele_freq = info_fields["AF"]
			ac = info_fields["AC"]
			an = info_fields["AN"]
			unknown = (nr_samples * 2) - int(info_fields["AN"])
			variant_id = info_fields['ID']

			outfile.write('\t'.join([
				variant_id,
				allele_freq,
				ac,
				an,
				str(unknown)
			]))

			if genotype_stats:
				total_genotypes = (nr_samples) - int(unknown / 2)
				heterozygosity = int(info_fields["AC_Het"]) / total_genotypes
				heterozygous_genotypes = info_fields["AC_Het"]
				unique_kmers = info_fields["UK"]
				high_qual_gt = info_fields["HGQ"]

				outfile.write('\t' + '\t'.join([
					str(heterozygosity),
					heterozygous_genotypes,
					str(total_genotypes),
					unique_kmers,
					high_qual_gt
				]) + '\n')
			else:
				outfile.write('\n')
 



if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog='collect-vcf-stats.py', description=__doc__)
	parser.add_argument('--genotyping-stats', action='store_true', help='Include genotype specific statistics')
	parser.add_argument('-outname', metavar='OUTNAME', required=True, help='Name of the output file' )
	parser.add_argument('-column-prefix', metavar='COLUMN_PREFIX', required=True, help='column prefix.')
	args = parser.parse_args()

	collect_stats(args.genotyping_stats, args.outname, args.column_prefix)
