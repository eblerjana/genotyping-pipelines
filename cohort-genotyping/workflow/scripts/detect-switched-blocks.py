import sys


genome_index_file = sys.argv[1]
threshold = 25

chr_to_len = {}
for line in open(genome_index_file, 'r'):
	fields = line.strip().split()
	chr_to_len[fields[0]] = fields[1]

block_coord = None
block_len = 0
block_open = False

prev_gt = "0|1"
prev_chr = ""

for line in sys.stdin:
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	gt_pos = fields[8].split(':').index('GT')
	cur_gt = fields[9].split(':')[gt_pos]
	cur_chr = fields[0]


	# check if end of chromosome is reached
	if prev_chr != cur_chr:
		if block_len > threshold:
			print('\t'.join(block_coord + [str(chr_to_len[prev_chr]), str(block_len)]))
		block_coord = None
		block_len = 0
		prev_chr = cur_chr
		prev_gt = "0|1"
		continue

	if (prev_gt == "0|1") and (cur_gt == "0|1"):
		prev_gt = cur_gt
		continue

	if (prev_gt == "1|0") and (cur_gt == "1|0"):
		assert block_len > 0
		block_len += 1
		prev_gt = cur_gt
		continue

	if (prev_gt == "0|1") and (cur_gt == "1|0"):
		assert block_len == 0
		block_coord = [cur_chr, fields[1]]
		block_len = 1
		prev_gt = cur_gt
		continue

	if (prev_gt == "1|0") and (cur_gt == "0|1"):
		assert block_len > 0
		if block_len > threshold:
			print('\t'.join(block_coord + [str(int(fields[1]) - 1), str(block_len)]))
		block_coord = None
		block_len = 0
		prev_gt = cur_gt
