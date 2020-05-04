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

cnt = 0
cnt2 = 0

infile = open("comboscore.txt")
for line in infile.readlines():
	items = line.strip().split('\t')
	idSample = items[0]
	idDrugA = items[1]
	idDrugB = items[2]
	comboScore = items[3]

	query = "select idCombo_Design from Combo_Design where idSample = " + str(idSample) + " and idDrugA = " + str(idDrugA) + " and idDrugB = " + str(idDrugB)
	cursor.execute(query)
	idCombo_Design = cursor.fetchone()
	if type(idCombo_Design) == type(None):
		query = "insert into Combo_Design (idSample, idDrugA, idDrugB) values(" + str(idSample) + ", " + str(idDrugA) + ", " + str(idDrugB) + ")"
		cursor.execute(query)
		idCombo_Design = db.insert_id()
		cnt = cnt + 1
	else:
		idCombo_Design = idCombo_Design[0]

	query = "update Synergy_Score set comboscore = " + quotation(comboScore) + " where idCombo_Design = " + str(idCombo_Design) + " and idSource = " + str(idSource)
	cursor.execute(query)
	cnt2 = cnt2 + 1

print("The number of combo designs: " + str(cnt))
print("The number of comboscores  : " + str(cnt2))

infile.close()
db.commit()
cursor.close()
db.close()
