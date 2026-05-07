import sys
import argparse


class VariantWindow:
	def __init__(self):
		"""
		Empty window.
		"""
		self.genotypes = []
		self.nr_ref = 0
		self.nr_alt = 0

	def add_gt(self, genotype):
		"""
		Add a new genotype to the window.
		"""
		if not genotype in ["0|1", "1|0"]:
			raise RuntimeError("Only heterozygous genotypes are allowed.")
		self.genotypes.append(genotype)
		if genotype == "0|1":
			self.nr_ref += 1
		else:
			self.nr_alt += 1

	def advance(self, genotype):
		"""
		Add a new genotype to the window and
		delete leftmost one.
		"""
		self.add_gt(genotype)
		left_gt = self.genotypes.pop(0)
		if left_gt == "0|1":
			self.nr_ref -= 1
		else:
			self.nr_alt -= 1

	def check_position(self, i):
		"""
		Check whether to correct this position. 
		"""
		if i >= len(self.genotypes):
			raise RuntimeError("VariantWindow: index out of bounds.")
			sys.exit(1)
		if self.genotypes[i] == "0|1":
			if self.nr_alt > self.nr_ref:
				self.genotypes[i] = "1|0"
				self.nr_alt += 1
				self.nr_ref -= 1
				return True
			else:
				return False
		else:
			if self.nr_ref > self.nr_alt:
				self.genotypes[i] = "0|1"
				self.nr_alt -= 1
				self.nr_ref += 1
				return True
			else:
				return False



def find_errors(k, input_stream):
	"""
	Find positions of errors to correct in
	reference haplotype.
	"""

	variant_window = VariantWindow()
	lines_read = 0
	first_window = True
	lines = []
	nr_windows = 0

	for line in input_stream:
		if line.startswith('#'):
			# header line
			yield line.strip()
			continue
		fields = line.strip().split()
		if len(fields[4].split(',')) > 1:
			raise RuntimeError("VCF is not biallelic.")
			sys.exit(1)
		if len(fields) != 10:
			raise RuntimeError("VCF must contain exactly one sample.")
			sys.exit(1)

		genotype_index = fields[8].split(':').index("GT")
		genotype_str = fields[9].strip().split(':')[genotype_index]
		if genotype_str == "1/1":
			genotype_str = "1|1"
		if genotype_str == "0/0":
			genotype_str = "0|0"
		if genotype_str not in ["1|1", "0|1", "1|0", "0|0"]:
			sys.stderr.write("Skipping genotype " + genotype_str + " at position " + fields[0] + ":" + fields[1] + "\n")
			continue
		
		if genotype_str == "1|1":
			# always correct hom genotype
			yield line.strip()
			continue

		if genotype_str == "0|0":
			continue

		if lines_read < 2*k+1:
			variant_window.add_gt(genotype_str)
			if lines_read >= k:
				lines.append(line.strip())
		else:
			nr_windows += 1
			current_line = lines.pop(0)
			if variant_window.check_position(k):
				yield current_line
			lines.append(line.strip())
			variant_window.advance(genotype_str)
		lines_read += 1

	sys.stderr.write("Considered " + str(nr_windows) + " windows.\n")



if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog='find_errors.py', description=__doc__)
	parser.add_argument('-window-width', metavar="WIDTH", type=int, default=15, help="Size of the window to consider.")
	args = parser.parse_args()

	for line in find_errors(args.window_width, sys.stdin):
		print(line)
