#source = "NCI-ALMANAC"
idSource = "2"

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

value_types = ["comboscore"]

for value_type in value_types:
	file_names = sorted(glob.glob("tbl.almanac." + value_type + ".txt"))

	for file_name in file_names:
		updated_cnt = 0
		inserted_cnt = 0
		na_cnt = 0

		infile = open(file_name)
		for line in infile.readlines():
			items = line.strip().split('\t')
			if len(items) > 5:
				idSample = items[0]
				idDrugA = items[1]
				idDrugB = items[2]
				concA = items[3]
				concB = items[4]
				comboScore = items[5]

				query = "select idCombo_Design from Combo_Design where idSample = " + str(idSample) + " and idDrugA = " + str(idDrugA) + " and idDrugB = " + str(idDrugB)
				cursor.execute(query)
				result = cursor.fetchone()
				if type(result) != type(None):
					idCombo_Design = result[0]
					query = "select idCombo_Matrix from Combo_Matrix where idCombo_Design = " + str(idCombo_Design) + " and concA = " + quotation(concA) + " and concB = " + quotation(concB) + " and idSource = " + str(idSource)
					cursor.execute(query)
					result = cursor.fetchone()
					if type(result) != type(None):
						idCombo_Matrix = result[0]
						query = "update Combo_Matrix set combo_score = " + quotation(comboScore) + " where idCombo_Matrix = " + str(idCombo_Matrix)
						cursor.execute(query)
						updated_cnt = updated_cnt + 1
				else:
					print("Warning 2:: " + str(idSample) + "\t" + str(idDrugA) + "\t" + str(idDrugB) + "\t" + quotation(concA) + "\t" + str(concB))
					print(query)
			else:
				na_cnt = na_cnt + 1

		print("    -- The number of updated records : " + str(updated_cnt))
		print("    -- The number of inserted records: " + str(inserted_cnt))
		print("    -- The number of NAs             : " + str(na_cnt))

infile.close()
db.commit()
cursor.close()
db.close()
