#source = "NCI-ALMANAC"
source = "2"

import os
import sys
import glob
import pymysql

def quotation(element):
	return_value = ""
	if element == "NA":
		return_value = "NULL"
	elif element == "NaN":
		return_value = "NULL"
	else:
		return_value = "\"" + element + "\""
	return return_value

db = pymysql.connect("localhost","root","gmldnjs!!","synergxdb" )
cursor = db.cursor()

value_types = ["raw", "bliss", "loewe", "hsa", "zip"]

for value_type in value_types:
	file_names = sorted(glob.glob("tbl.almanac." + value_type + ".txt"))

	# name	drugA	drugB	concA	concB	experiments
	#OCUBM	ABT-888	Fluorouracil	0	0	1
	#OCUBM	ABT-888	Fluorouracil	0	0.35	NA
	#OCUBM	ABT-888	Fluorouracil	0	1.08	0.97966250606605

	for file_name in file_names:
		updated_cnt = 0
		inserted_cnt = 0

		infile = open(file_name)
		for line in infile.readlines():
			items = line.strip().split('\t')
			idSample = items[0]
			idDrugA = items[1]
			idDrugB = items[2]
			concA = items[3]
			concB = items[4]
			value = items[5]

			query = "select idCombo_Design from Combo_Design where idSample = " + str(idSample) + " and idDrugA = " + str(idDrugA) + " and idDrugB = " + str(idDrugB)
			cursor.execute(query)
			result = cursor.fetchone()
			if type(result) != type(None):
				idCombo_Design = result[0]
			
				query = "select idCombo_Matrix from Combo_Matrix where idCombo_Design = " + str(idCombo_Design) + " and idSource = " + quotation(source) + " and concA = " + str(concA) + " and concB = " + str(concB) 
				cursor.execute(query)
				result = cursor.fetchone()
				if type(result) != type(None):
					idCombo_Matrix = result[0]
					query = "update Combo_Matrix set " + value_type + "_matrix = " + quotation(value) + " where idCombo_Matrix = " + str(idCombo_Matrix)
					cursor.execute(query)
					updated_cnt = updated_cnt + 1
				else:
					query = "insert into Combo_Matrix (idCombo_Design, concA, concB, " + value_type + "_matrix, idSource) values(" + str(idCombo_Design) + ", " + str(concA) + ", " + str(concB) + ", " + quotation(value) + ", " + quotation(source) + ")"
					cursor.execute(query)
					inserted_cnt = inserted_cnt + 1
			else:
				print("Warning:: " + name + "\t" + drugA + " : " + drugB)
				print(query)
		print("    -- The number of updated records : " + str(updated_cnt))
		print("    -- The number of inserted records: " + str(inserted_cnt))

infile.close()
db.commit()
cursor.close()
db.close()
