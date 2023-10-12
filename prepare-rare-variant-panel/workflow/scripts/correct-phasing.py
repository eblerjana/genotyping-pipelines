import sys

for line in sys.stdin:
	if line.startswith('#'):
		print(line.strip())
		continue
	fields = line.strip().split()
	for i in range(9, len(fields)):
		if '/' in fields[i]:
			fields[i] = fields[i].replace("/", "|")
	print("\t".join(fields))
