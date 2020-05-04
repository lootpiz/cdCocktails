import os
import sys

infile_name = sys.argv[1]

targets = []
infile = open("target.txt")
for line in infile.readlines():
	target = line.strip()
	targets.append(target)
infile.close()

infile = open(infile_name)
outfile_name = infile_name.replace(".txt", "_out.txt")
outfile = open(outfile_name, 'w')

header = True
for line in infile.readlines():
	if header:
		outfile.write(line)
		header = False
	else:
		items = line.strip().split("\t")
		ensg = items[1]
		if ensg in targets:
			outfile.write(line)
infile.close()
outfile.close()
		
print("Infile     : " + infile_name)
print("Target file: target.txt")
print("Outfile    : " + outfile_name)
