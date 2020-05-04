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

cnt = 0
cnt2 = 0

#Cellline	DrugA	DrugB	Bliss	Loewe	HSA	ZIP
#786-0	Melphalan	Dexrazoxane	NA	NA	NA	-9e-04
#786-0	7-Ethyl-10-hydroxycamptothecin	Melphalan	NA	NA	NA	-8e-04
#786-0	Melphalan	Vinorelbine_tartrate	NA	NA	NA	-0.0011
infile = open("tbl.almanac.combo_summary.txt")
for line in infile.readlines():
	items = line.strip().split('\t')
	idSample = items[0]
	idDrugA = items[1]
	idDrugB = items[2]
	bliss = items[3]
	loewe = items[4]
	hsa = items[5]
	zip = items[6]

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

	query = "insert ignore into Synergy_Score (idCombo_Design, bliss, loewe, hsa, zip, idSource) values(" + str(idCombo_Design) + ", " + quotation(bliss) + ", " + quotation(loewe) + ", " + quotation(hsa) + ", " + quotation(zip) + ", " + quotation(source) + ")"
	cursor.execute(query)
	cnt2 = cnt2 + 1

print("The number of combo designs: " + str(cnt))
print("The number of experiments  : " + str(cnt2))

infile.close()
db.commit()
cursor.close()
db.close()
