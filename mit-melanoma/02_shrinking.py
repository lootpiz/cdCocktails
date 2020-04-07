import os
import sys

drugs = {}
infile = open("mapping_drug.txt")
for line in infile.readlines()[1:]:
	items = line.strip().split("\t")
	mit_id = items[0]
	syn_id = items[2]
	drugs[mit_id] = syn_id
infile.close()

cells = {}
infile = open("mapping_sample.txt")
for line in infile.readlines()[1:]:
	items = line.strip().split("\t")
	mit_id = items[0]
	syn_id = items[2]
	cells[mit_id] = syn_id
infile.close()

infile_combo = open("combo.txt")
outfile_combo = open("combo_input.txt", 'w')
for line in infile_combo.readlines()[1:]:
	items = line.strip().split("\t")
	cell = cells[items[0]]
	drugA = drugs[items[1]]
	concA = items[2]
	drugB = drugs[items[3]]
	concB = items[4]	
	viability = items[5]
	outfile_combo.write(str(cell) + "\t" + str(drugA) + "\t" + str(concA) + "\t" + str(drugB) + "\t" + str(concB) +  "\t" + str(viability) + "\n")
infile_combo.close()
outfile_combo.close()

infile_mono = open("mono.txt")
outfile_mono = open("mono_input.txt", 'w')
for line in infile_mono.readlines()[1:]:
	items = line.strip().split("\t")
	cell = cells[items[0]]
	drugA = drugs[items[1]]
	concA = items[2]
	viability = items[3]
	outfile_mono.write(str(cell) + "\t" + str(drugA) + "\t" + str(concA) + "\t" + str(viability) + "\n")
infile_mono.close()
outfile_mono.close()
