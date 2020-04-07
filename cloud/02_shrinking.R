input <- read.table("combo.txt", header=T, sep="\t", stringsAsFactor=F)
meta <- read.table("concentrations.txt", header=T, sep="\t", stringsAsFactor=F)

master_drugs <- read.table("mapping_drug.txt", header=T, sep="\t")
idSample <- 160

outfile <- c()
for(row_idx in c(1:nrow(input))){
	row <- input[row_idx, ]

	drug_id_X1 <- row[1]
	drug_name_X1 <- meta$DRUG_NAME[which(meta$CLOUD_ID == as.character(drug_id_X1))]
	drug_X1_conc <- meta$CONCENTRATION[which(meta$CLOUD_ID == as.character(drug_id_X1))]
	drug_X1_y <- unlist(row[3])
	idDrugA <- master_drugs$id[which(master_drugs$name == drug_name_X1)]

	drug_id_X2 <- row[2]
	drug_name_X2 <- meta$DRUG_NAME[which(meta$CLOUD_ID == as.character(drug_id_X2))]
	drug_X2_conc <- meta$CONCENTRATION[which(meta$CLOUD_ID == as.character(drug_id_X2))]
	drug_X2_y <- unlist(row[4])
	idDrugB <- master_drugs$id[which(master_drugs$name == drug_name_X2)]

	drug_combo <- unlist(row[5])

	result <- c(idSample, idDrugA, drug_X1_conc, drug_X1_y, idDrugB, drug_X2_conc, drug_X2_y, drug_combo)
	outfile <- rbind(outfile, result)
}

write.table(outfile, file="combo_input.txt", row.names=F, col.names=F, quote=F, sep="\t")

q("no")

