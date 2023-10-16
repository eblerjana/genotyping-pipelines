import sys

for line in sys.stdin:
	if line.startswith("##"):
		print(line.strip())
		continue
	fields = line.strip().split()
	if line.startswith("#"):
		# combine paths to diploid pseudo-samples
		samples = len(fields[9:])
		if samples % 2 == 0:
			new_samples = ["pseudo-sample" + str(i) for i in range(int(samples/2))]
		else:
			new_samples = ["pseudo-sample" + str(i) for i in range(int((samples + 1)/2))]
		print("\t".join(fields[:9] + new_samples))
		continue
	
	current_gt = ""
	current_index = 0
	gt_fields = fields[9:]
	genotypes = []
	
	while current_index < samples:
		if (current_index % 2 != 0) and (current_index > 0):
			current_gt += "|" + gt_fields[current_index]
			genotypes.append(current_gt)
			current_gt = ""
		else:
			current_gt = gt_fields[current_index]
		current_index += 1
	if not("|" in current_gt) and (samples % 2 != 0):
		current_gt += "|" + current_gt
		genotypes.append(current_gt)	
	print("\t".join(fields[:9] + genotypes))
