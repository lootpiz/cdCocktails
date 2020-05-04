input_file_name = "ComboDrugGrowth_Nov2017.csv"

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
import os
import csv
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

cells = {}
infile = open("mapping_sample.txt")
for line in infile.readlines()[1:]:
	items = line.strip().split("\t")
	name = items[0]
	syn_id = items[2]
	cells[name] = syn_id
infile.close()

sheet_dic_mono = defaultdict(list)
sheet_dic_combo = defaultdict(list)

with open(input_file_name) as csvfile:
	infile = csv.reader(csvfile, delimiter=',')
	for line in infile:
		if line[0] != "COMBODRUGSEQ": # skip header
			if '\x1a' in line[28]:
				cell_name = line[28][-1]
			else:
				cell_name = line[28]
			cell = cells[cell_name]
			drugA = drugs[line[8]]
			concA = line[11]

			hasDrugB = line[14]
			if hasDrugB != "": # combo
				drugB = drugs[line[14]]
				concB = line[17]

				key = str(cell) + "$" + str(drugA) + "$" + str(concA) + "$" + str(drugB) + "$" + str(concB)

				if line[20] != "":
					viability = float(line[20])/100.0
					sheet_dic_combo[key].append(viability)

			else: # mono
				if line[20] != "":
					viability = float(line[20])/100.0
					key = str(cell) + "$" + str(drugA) + "$" + str(concA)
					sheet_dic_mono[key].append(viability)

outfile_mono = open("mono_input.txt", 'w')
for key, values in sheet_dic_mono.iteritems():
	outfile_mono.write('\t'.join(key.split('$')))
	outfile_mono.write("\t" + '\t'.join(str(v) for v in values))
	outfile_mono.write("\n")
outfile_mono.close()


outfile_combo = open("combo_input.txt", 'w')
for key, values in sheet_dic_combo.iteritems():
	outfile_combo.write('\t'.join(key.split('$')))
	outfile_combo.write("\t" + '\t'.join(str(v) for v in values))
	outfile_combo.write("\n")
outfile_combo.close()
