import os
import sys
from collections import defaultdict

drugs = {}
infile = open("mapping_drug.txt")
for line in infile.readlines()[1:]:
	items = line.strip().split("\t")
	nci_id = items[0]
	syn_id = items[1]
	drugs[nci_id] = syn_id
infile.close()

infile_combo = open("combo.txt")
outfile_combo = open("combo_input.txt", 'w')
for line in infile_combo.readlines()[1:]:
	items = line.strip().split("\t")
	drug = drugs[items[0]]
	outfile_combo.write("186\t" + str(drug) + "\t" + "\t".join(str(v) for v in items[1:5]) + "\n")
infile_combo.close()
outfile_combo.close()
