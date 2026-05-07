import sys

fai = sys.argv[1]
outprefix = sys.argv[2]

for line in open(fai, 'r'):
	fields = line.strip().split()
	with open(outprefix + fields[0] + '.txt', 'w') as outfile:
		outfile.write(fields[0] + '\n')

