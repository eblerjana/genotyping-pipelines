import sys

for line in sys.stdin:
	if line.startswith("#"):
		print(line.strip())
		continue
	fields = line.strip().split()
	gt_index = fields[8].split(':').index('GT')
	gt_field = fields[9].split(':')
	gt = gt_field[gt_index]

	if gt == "1|0":
		gt_field[gt_index] = "0|1"

	fields[9] = ":".join(gt_field)
	print("\t".join(fields))	
