import sys
import argparse
import random
import gzip

class Variant:
	"""
	Represents a variant record.
	"""
	def __init__(self, chrom, start, ref_allele, alt_allele, var_index: int, var_id = '.', var_ident = None):
		self._chrom = chrom
		self._start = start
		self._ref_allele = ref_allele
		self._alt_allele = alt_allele
		self._end = start + len(ref_allele)
		self._index = var_index
		self._id = var_id
		self._ident = var_ident
		self._last_in_cluster = False

	def chrom(self):
		return self._chrom

	def start(self):
		return self._start

	def end(self):
		return self._end

	def ref(self):
		return self._ref_allele

	def alt(self):
		return self._alt_allele

	def index(self):
		return self._index

	def id(self):
		return self._id

	def set_last_in_cluster(self):
		self._last_in_cluster = True

	def is_last(self):
		return self._last_in_cluster

	def ident(self):
		return self._ident



def variants_overlap(var1: Variant, var2: Variant, margin = 0):
	"""
	checks whether two variant intervals are overlapping.
	"""
	if var1.chrom() != var2.chrom():
		return False

	start1 = max(var1.start()-margin, 0)
	end1 = var1.end() + margin

	start2 = max(var2.start()-margin, 0)
	end2 = var2.end() + margin

	if start1 < start2:
		if start2 < end1:
			return True
		else:
			return False

	else:
		if start1 < end2:
			return True
		else:
			return False

class Cluster:
	"""
	Represents clusters of Variant objects.
	"""
	def __init__(self):
		self._variants = []

	def add_var(self, variant):
		self._variants.append(variant)

	def merge(self, cluster):
		self._variants = cluster._variants

	def overlaps_last(self, variant: Variant, margin = 0):
		# overlap last variant overlaps
		if not self._variants:
			return False
		else:
			return variants_overlap(variant, self._variants[-1], margin)

	def __iter__(self):
		for v in self._variants:
			yield v
	

def group_variants(variants, margin = 0):
	"""
	Given a list of variant clusters, merge non-overlapping clusters.
	"""

	offset = 0 #variants[0].index()
	clusters = []
	for var in variants:
		added_to_cluster = False
		for c in clusters:
			if not c.overlaps_last(var, margin):
				c.add_var(var)
				added_to_cluster = True
				break
		if not added_to_cluster:
			new_c = Cluster()
			new_c.add_var(var)
			clusters.append(new_c)

	## return sorted list (sorted by index)
	assignments = [None] * len(variants)
	for i,c in enumerate(clusters):
		for v in c:
			assignments[v.index()] = i

	return assignments


def print_vcf(variants, assignments):
	"""
	Given variants and their group assignments,
	write a VCF with one haplotype by group.
	"""

	if len(variants) != len(assignments):
		raise Exception("Error in print_vcf: number of variants does not match the group assignments.")

	# total number of paths
	n_paths = max(assignments)

	# print VCF header
	print("##fileformat=VCFv4.2")
	print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
	print('##INFO=<ID=ID,Number=A,Type=String,Description=\"Variant IDs per ALT allele.\">')
	path_names = ["path" + str(i) for i in range(n_paths + 1)]
	print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + '\t'.join(path_names))

	# print variant records and haplotype paths
	for var, a in zip(variants, assignments):
		vcf_line = [var.chrom(), str(var.start()), var.id(), var.ref(), var.alt(), '.', 'PASS', "ID=" + var.ident(), 'GT']
		paths = ['0' if i != a else '1' for i in range(n_paths + 1)]
		print('\t'.join(vcf_line + paths))


def print_individual_vcfs(variants, assignments, prefix):
	"""
	Given variants and their group assignments,
	write a VCF for each group.
	"""

	if len(variants) != len(assignments):
		raise Exception("Error in print_vcf: number of variants does not match the group assignments.")

	n_paths = max(assignments)
	for group in range(n_paths + 1):
		with open(prefix + '_path' + str(group) + '.vcf', 'w') as outfile:
			outfile.write("##fileformat=VCFv4.2\n")
			outfile.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
			outfile.write('##INFO=<ID=ID,Number=A,Type=String,Description=\"Variant IDs per ALT allele.\">\n')
			header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "path" + str(group) + '\n'
			outfile.write(header_line)
			for var, a in zip(variants, assignments):
				if a == group:
					vcf_line = [var.chrom(), str(var.start()), var.id(), var.ref(), var.alt(), '.', 'PASS', "ID=" + var.ident(), 'GT', '1']
					outfile.write('\t'.join(vcf_line) + '\n')



def run_paths(filename, single, margin):
	# list of all variants in VCF
	all_variants = []
	# assigns a group to each variant
	groups = []
	total_variants = 0
	skipped_variants = 0

	current_cluster = []
	prev_end = 0
	prev_chrom = None
	for line in gzip.open(filename, 'rt'):
		if line.startswith('#'):
			# header line
			continue
		total_variants += 1
		fields = line.strip().split()
		chrom = fields[0]
		start = int(fields[1])
		reference = fields[3]

		# make sure VCF is biallelic
		alt_alleles = fields[4].split(',')
		if len(alt_alleles) > 1:
			raise Exception(filename + ' contains multi-allelic records. Make sure to convert VCF to biallelic (bcftools norm -m-) before running this program.')
		alt = alt_alleles[0]

		# skip symbolic records
		if "<" in alt:
			sys.stderr.write("Skipping variant at position " + chrom + ":" + str(start) + " because it is symbolic.\n")
			skipped_variants += 1
			continue

		# skip break-ends
		if "SVTYPE=BND" in line:
			sys.stderr.write("Skipping variant at position " + chrom + ":" + str(start) + " because it is a BND record.\n")
			skipped_variants += 1
			continue

		# make sure VCF is sorted
		if len(all_variants) > 0:
			if (chrom == all_variants[-1].chrom) and (start < all_variants[-1].start()):
				raise Exception(filename + ' is not sorted.')

		if ((start >= prev_end) or (prev_chrom != chrom)) and current_cluster:
			# process current cluster and group non-overlapping variants
			assignments = group_variants(current_cluster, margin)
			assert len(assignments) == len(current_cluster)
			groups += assignments
			all_variants += current_cluster
			all_variants[-1].set_last_in_cluster()
			current_cluster = []

		var_id = fields[2]
		var_ident = None
		format_field = {f.split("=")[0] : f.split("=")[1] for f in fields[7].split(";") if "=" in f}
		if "ID" in format_field:
			var_ident = format_field["ID"]
		var_index = len(current_cluster)
		variant = Variant(chrom, start, reference, alt, var_index, var_id, var_ident)

		current_cluster.append(variant)
		prev_chrom = chrom
		prev_end = max(variant.end(), prev_end)

	# handle last cluster (in case there is one)
	if current_cluster:
		assignments = group_variants(current_cluster, margin)
		assert len(assignments) == len(current_cluster)
		groups += assignments
		all_variants += current_cluster

	# randomly distiribute groupings to paths
	mapping = [i for i in range(max(groups) + 1)]
	for i in range(len(all_variants)):
		groups[i] = mapping[groups[i]]
		if all_variants[i].is_last():
			random.shuffle(mapping)


	if single:
		print_individual_vcfs(all_variants,groups, single)
	else:
		print_vcf(all_variants, groups)

	# print log
	sys.stderr.write("\n########## Summary ##########\n")
	sys.stderr.write("Number of variants read from " + filename + ": " + str(total_variants) + "\n")
	sys.stderr.write("Number of variants written " + str(len(all_variants)) + "\n")
	sys.stderr.write("Number of variants skipped " + str(skipped_variants) + "\n")



if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='create_paths.py', description=__doc__)
	parser.add_argument('-vcf', metavar='VCF', required=True, help='callset in VCF format.')
	parser.add_argument('-single', metavar='PREFIX', default=None, help='write one VCF per path to <PREFIX>_<pathID>.vcf')
	parser.add_argument('-margin', metavar='MARGIN', type=int, default=0, help='margin to consider around variants.')
	args = parser.parse_args()
	run_paths(args.vcf, args.single, args.margin)
