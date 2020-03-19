library(synergyfinder)

input <- read.table("../data/CLOUD/cloud_input.txt", header=T, sep="\t", stringsAsFactor=F)
meta <- read.table("../data/CLOUD/conc_out.txt", header=T, sep="\t", stringsAsFactor=F)

master_drugs <- read.table("../data/mapping/CLOUD_drug.txt", header=T, sep="\t")
master_cells <- read.table("../data/mapping/CLOUD_sample.txt", header=T, sep="\t")

cell_name <- "KBM7"
cell_name_lbl <- "KBM-7"
idSample <- 160

drug_names <- unique(union(input[,1], input[,2]))

outfile_stats <- c()

raw_matrix <- c()
bliss_matrix <- c()
loewe_matrix <- c()
hsa_matrix <- c()
zip_matrix <- c()

for(row_idx in c(1:nrow(input))){
	print(row_idx)
	row <- input[row_idx, ]

	drug_id_X1 <- row[1]
	drug_name_X1 <- meta$DRUG_NAME[which(meta$CLOUD_ID == as.character(drug_id_X1))]
	drug_X1_conc <- meta$CONCENTRATION[which(meta$CLOUD_ID == as.character(drug_id_X1))]
	drug_X1_y <- row[3]
	idDrugA <- master_drugs$id[which(master_drugs$name == drug_name_X1)]

	drug_id_X2 <- row[2]
	drug_name_X2 <- meta$DRUG_NAME[which(meta$CLOUD_ID == as.character(drug_id_X2))]
	drug_X2_conc <- meta$CONCENTRATION[which(meta$CLOUD_ID == as.character(drug_id_X2))]
	drug_X2_y <- row[4]
	idDrugB <- master_drugs$id[which(master_drugs$name == drug_name_X2)]

	drug_combo <- row[5]

	Z <- array(c(0, as.numeric(drug_X1_y), as.numeric(drug_X2_y), as.numeric(drug_combo)), dim=c(2,2))
	rownames(Z) <- c("0", as.character(drug_X1_conc))
	colnames(Z) <- c("0", as.character(drug_X2_conc))

	expr_info <- data.frame(drug.col = drug_name_X2, drug.row = drug_name_X1, concUnit = "µM", blockIDs = 1)
	data <- list(dose.response.mats = list(Z), drug.pairs = expr_info)

	bliss_score <- NA; bliss_mat <- NA; bliss <- NA
	loewe_score <- NA; loewe_mat <- NA; loewe <- NA
	hsa_score <- NA; hsa_mat <- NA; hsa <- NA
	zip_score <- NA; zip_mat <- NA; zip <- NA

	bliss_exe <- tryCatch({ bliss <- CalculateSynergy(data=data, correction=F, method="Bliss")
		}, error = function(err) { print("Bliss error") ; bliss <- NA
	})
	if(!is.na(bliss_exe)){
		bliss_score <- round(mean(bliss$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)
		bliss_mat <- bliss$scores[[1]]/100.0
	} else {
		bliss_mat <- matrix(rep(NA, length(c(Z))), ncol=ncol(Z), nrow=nrow(Z))
		rownames(bliss_mat) <- c("0", as.character(drug_X1_conc))
		colnames(bliss_mat) <- c("0", as.character(drug_X2_conc))
	}

	hsa_exe <- tryCatch({ hsa <- CalculateSynergy(data=data, correction=F, method="HSA")
		}, error = function(err) { print("hsa error") ; hsa <- NA
	})
	if(!is.na(hsa_exe)){
		hsa_score <- round(mean(hsa$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)
		hsa_mat <- hsa$scores[[1]]/100.0
	} else {
		hsa_mat <- matrix(rep(NA, length(c(Z))), ncol=ncol(Z), nrow=nrow(Z))
		rownames(hsa_mat) <- c("0", as.character(drug_X1_conc))
		colnames(hsa_mat) <- c("0", as.character(drug_X2_conc))
	}

	result <- c(idSample, idDrugA, idDrugB, bliss_score, loewe_score, hsa_score, zip_score)

	outfile_stats <- rbind(outfile_stats, result)

	if(idDrugA < idDrugB) {
		Z_out <- data.frame(
			name = rep(idSample, length(c(Z))),
			drugA = rep(idDrugA, length(c(Z))),
			drugB = rep(idDrugB, length(c(Z))),
			concA = as.numeric(rep(colnames(Z), each=nrow(Z))),
			concB = as.numeric(rep(rownames(Z), ncol(Z))),
			experiments = as.numeric(c(Z))
		)
	} else {
		Z_out <- data.frame(
			name = rep(idSample, length(c(Z))),
			drugA = rep(idDrugB, length(c(Z))),
			drugB = rep(idDrugA, length(c(Z))),
			concA = as.numeric(rep(rownames(Z), each=ncol(Z))),
			concB = as.numeric(rep(colnames(Z), nrow(Z))),
			experiments = as.numeric(c(Z))
		)
	}
	raw_matrix <- rbind(raw_matrix, Z_out)

	if(idDrugA < idDrugB) {
		bliss_out <- data.frame(
			name = rep(idSample, length(c(bliss_mat))),
			drugA = rep(idDrugA, length(c(bliss_mat))),
			drugB = rep(idDrugB, length(c(bliss_mat))),
			concA = as.numeric(rep(c("0", as.character(drug_X2_conc)), each=nrow(bliss_mat))),
			concB = as.numeric(rep(c("0", as.character(drug_X1_conc)), ncol(bliss_mat))),
			bliss = as.numeric(c(bliss_mat))
		)
	} else {
	bliss_out <- data.frame(
			name = rep(idSample, length(c(bliss_mat))),
			drugA = rep(idDrugB, length(c(bliss_mat))),
			drugB = rep(idDrugA, length(c(bliss_mat))),
			concA = as.numeric(rep(c("0", as.character(drug_X1_conc)), ncol(bliss_mat))),
			concB = as.numeric(rep(c("0", as.character(drug_X2_conc)), each=nrow(bliss_mat))),
			bliss = as.numeric(c(bliss_mat))
		)
	}
	bliss_matrix <- rbind(bliss_matrix, bliss_out)

	if(idDrugA < idDrugB) {
		hsa_out <- data.frame(
			name = rep(idSample, length(c(hsa_mat))),
			drugA = rep(idDrugA, length(c(hsa_mat))),
			drugB = rep(idDrugB, length(c(hsa_mat))),
			concA = as.numeric(rep(c("0", as.character(drug_X2_conc)), each=nrow(hsa_mat))),
			concB = as.numeric(rep(c("0", as.character(drug_X1_conc)), ncol(hsa_mat))),
			hsa = as.numeric(c(hsa_mat))
		)
	} else {
	hsa_out <- data.frame(
			name = rep(idSample, length(c(hsa_mat))),
			drugA = rep(idDrugB, length(c(hsa_mat))),
			drugB = rep(idDrugA, length(c(hsa_mat))),
			concA = as.numeric(rep(c("0", as.character(drug_X1_conc)), ncol(hsa_mat))),
			concB = as.numeric(rep(c("0", as.character(drug_X2_conc)), each=nrow(hsa_mat))),
			hsa = as.numeric(c(hsa_mat))
		)
	}
	hsa_matrix <- rbind(hsa_matrix, hsa_out)
}


colnames(outfile_stats) <- c("idSample", "idDrugA", "idDrugB", "Bliss", "Loewe", "HSA", "ZIP")
write.table(outfile_stats, file=paste0("../results/tbl.cloud.combo_summary.txt"), row.names=F, col.names=F, quote=F, sep="\t")

write.table(raw_matrix, file=paste("../results/tbl.cloud.raw.txt",sep =""), row.names=F, col.names=T, quote=F, sep="\t")
write.table(bliss_matrix, file=paste("../results/tbl.cloud.bliss.txt",sep =""), row.names=F, col.names=T, quote=F, sep="\t")
write.table(hsa_matrix, file=paste("../results/tbl.cloud.hsa.txt",sep =""), row.names=F, col.names=T, quote=F, sep="\t")

quit("no")


