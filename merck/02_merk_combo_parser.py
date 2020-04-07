#-*- coding: utf-8 -*-
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
input_file_name = "156849_1_supp_1_w2lrww.xls"
outfile_name_prefix = "combo"
key_indices = [1, 2, 3, 4, 5]
ignore_indices = [0, 6, 11, 12]


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
import os
import sys
from xlrd import open_workbook
from xlrd.sheet import ctype_text  
from django.utils.encoding import smart_str, smart_unicode
from collections import defaultdict

def isnumber(x):
	returnValue = True
	try:
		float(x)
	except ValueError:
		returnValue = False
	return returnValue

## Drug master
master_drugs = {}
drug_infile = open("mapping_drug.txt")
for line in drug_infile.readlines():
        items = line.strip().split("\t")
        name = items[0]
        idDrug = items[2]
        master_drugs[name] = idDrug
drug_infile.close()

master_samples = {}
cell_infile = open("mapping_sample.txt")
for line in cell_infile.readlines():
        items = line.strip().split("\t")
        name = items[0]
        idSample = items[2]
        master_samples[name] = idSample
cell_infile.close()

## File loading
book = open_workbook(input_file_name)

print "=" * 50
print "Loaded      : " + input_file_name
print "No of sheets: " + str(book.nsheets)
for sheet_name in book.sheet_names():
	print " - " + sheet_name
print "=" * 50

# for each sheet(s)
for sheet in book.sheets():
	# sheet info
	sheet_name = sheet.name
	number_of_rows = sheet.nrows
	number_of_cols = sheet.ncols

	print " -" * 25
	print " Processing : " + sheet_name
	print " Rows       : " + str(number_of_rows)
	print " Columns    : " + str(number_of_cols)

	# header and key info
	header = sheet.row(0)
	key_names = []
	for key_idx in key_indices:
		key_names.append("\"" + smart_str(header[key_idx].value) + "\"")
	print " Key         : " + "-".join(key_names)

	# for each row
	sheet_dic = defaultdict(list)
	for row_idx in range(1, number_of_rows):
		key_values = []
		values = []
		for col_idx in range(0, number_of_cols):
			if col_idx not in ignore_indices:
				cell_obj = sheet.cell(row_idx, col_idx)
				cell_value = cell_obj.value

				if col_idx in key_indices:
					key_values.append(cell_value)
				else:
					if not isnumber(cell_value):
						cell_value = "NA"
					else:
						values.append(cell_value)

		sheet_dic['$'.join(str(v) for v in key_values)].append(values)

	# write files
	outfile_value = open(outfile_name_prefix +  ".txt", 'w')
	for key, values in sheet_dic.iteritems():
		items = key.split('$')
		idSample = master_samples[items[0]]
		idDrugA = master_drugs[items[1]]
		concA = items[2]
		idDrugB = master_drugs[items[3]]
		concB = items[4]
		outfile_value.write(str(idSample) + "\t" + str(idDrugA) + "\t" + str(concA) + "\t" + str(idDrugB) + "\t" + str(concB))
		for value in values:
			outfile_value.write("\t" + '\t'.join(str(v) for v in value))
		outfile_value.write("\n")

	outfile_value.close()

	print " -" * 25

print "=" * 50
print "Done"
print "Outfile 2   : " + outfile_name_prefix + ".txt" 
print "=" * 50
