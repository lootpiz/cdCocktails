library(synergyfinder)

all_files <- list.files("../data/VISAGE/")
cell_names <- sub('.txt', '', all_files)
drug_names <- c("BI2536", "TH588")

master_drugs <- read.table("../data/resource/master_mapping_drug.txt", header=T, sep="\t", stringsAsFactor=F)
master_cells <- read.table("../data/resource/master_mapping_sample.txt", header=T, sep="\t", stringsAsFactor=F)

outfile_stats <- c()

raw_matrix <- c()
bliss_matrix <- c()
loewe_matrix <- c()
hsa_matrix <- c()
zip_matrix <- c()

idDrugA <- 2015
drug_name_X1_id <- 2015
drug_name_X1 <- master_drugs[2,2]
idDrugB <- 164
drug_name_X2_id <- 164
drug_name_X2 <- master_drugs[1,2]


for(cell_name in cell_names) {
	print(cell_name)
	dat <- read.table(paste0("../data/VISAGE/", cell_name, ".txt"), header=F, sep="\t", stringsAsFactor=F)
	Z <- array(as.numeric(as.matrix(dat[c(2:nrow(dat)), c(2:ncol(dat))])), dim=c(6, 10))
	colnames(Z) <- dat[1, c(2:ncol(dat))]
	rownames(Z) <- dat[c(2:nrow(dat)), 1]

	combo_row_concentration <- rownames(Z)
	combo_col_concentration <- colnames(Z)

	idSample <- master_cells$id[which(master_cells$name == cell_name)]
	sample_name_id <- idSample

	meta <- data.frame(drug.col = drug_name_X1, drug.row = drug_name_X2, concUnit = "µM", blockIDs = 1)
	data <- list(dose.response.mats = list(Z*100), drug.pairs = meta)

	bliss_score <- NA; bliss_mat <- NA; bliss <- NA
	loewe_score <- NA; loewe_mat <- NA; loewe <- NA
	hsa_score <- NA; hsa_mat <- NA; hsa <- NA
	zip_score <- NA; zip_mat <- NA; zip <- NA

	bliss_exe <- tryCatch({ bliss <- CalculateSynergy(data=data, method="Bliss", correction=FALSE)
		}, error = function(err) { print("Bliss error") ; bliss <- NA
	})
	if(!is.na(bliss_exe)){
		bliss_score <- round(mean(bliss$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)
		bliss_mat <- bliss$scores[[1]]/100.0
	} else {
		bliss_mat <- matrix(rep(NA, length(c(Z))), ncol=ncol(Z), nrow=nrow(Z))
		rownames(bliss_mat) <- c(combo_row_concentration)
		colnames(bliss_mat) <- c(combo_col_concentration)
	}

	loewe_exe <- tryCatch({ loewe <- CalculateSynergy(data=data, method="Loewe", correction=FALSE)
		}, error = function(err) { print("Loewe error") ; loewe <- NA
	})
	if(!is.na(loewe_exe)){
		loewe_score <- round(mean(loewe$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)
		loewe_mat <- loewe$scores[[1]]/100.0
	} else {
		loewe_mat <- matrix(rep(NA, length(c(Z))), ncol=ncol(Z), nrow=nrow(Z))
		rownames(loewe_mat) <- c(combo_row_concentration)
		colnames(loewe_mat) <- c(combo_col_concentration)
	}

	hsa_exe <- tryCatch({ 	hsa <- CalculateSynergy(data=data, method="HSA", correction=FALSE)
		}, error = function(err) { print("HSA error"); hsa <- NA
	})
	if(!is.na(hsa_exe)){
		hsa_score <- round(mean(hsa$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)
		hsa_mat <- hsa$scores[[1]]/100.0
	} else {
		hsa_mat <- matrix(rep(NA, length(c(Z))), ncol=ncol(Z), nrow=nrow(Z))
		rownames(hsa_mat) <- c(combo_row_concentration)
		colnames(hsa_mat) <- c(combo_col_concentration)
	} 

	zip_exe <- tryCatch({ zip <- CalculateSynergy(data=data, method="ZIP", correction=FALSE) 
		}, error = function(err) { print("ZIP error"); zip <- NA 
	})
	if(!is.na(zip_exe)){
		zip_score <- round(mean(zip$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)			
		zip_mat <- zip$scores[[1]]/100.0
	} else {
		zip_mat <- matrix(rep(NA, length(c(Z))), ncol=ncol(Z), nrow=nrow(Z))
		rownames(zip_mat) <- c(combo_row_concentration)
		colnames(zip_mat) <- c(combo_col_concentration)
	} 

	result <- c(idSample, idDrugA, idDrugB, bliss_score, loewe_score, hsa_score, zip_score)
	outfile_stats <- rbind(outfile_stats, result)

	Z_out <- data.frame(
			name = rep(sample_name_id, length(c(Z))),
			drugA = rep(drug_name_X2_id, length(c(Z))),
			drugB = rep(drug_name_X1_id, length(c(Z))),
			concA = as.numeric(rep(rownames(Z), each=ncol(Z))),
			concB = as.numeric(rep(colnames(Z), nrow(Z))),
			experiments = as.numeric(c(Z))
	)
	raw_matrix <- rbind(raw_matrix, Z_out)

	bliss_out <- data.frame(
			name = rep(sample_name_id, length(c(bliss_mat))),
			drugA = rep(drug_name_X2_id, length(c(bliss_mat))),
			drugB = rep(drug_name_X1_id, length(c(bliss_mat))),
			concA = as.numeric(rep(rownames(Z), each=ncol(Z))),
			concB = as.numeric(rep(colnames(Z), nrow(Z))),
			bliss = as.numeric(c(bliss_mat))
	)
	bliss_matrix <- rbind(bliss_matrix, bliss_out)

	loewe_out <- data.frame(
			name = rep(sample_name_id, length(c(loewe_mat))),
			drugA = rep(drug_name_X2_id, length(c(loewe_mat))),
			drugB = rep(drug_name_X1_id, length(c(loewe_mat))),
			concA = as.numeric(rep(rownames(Z), each=ncol(Z))),
			concB = as.numeric(rep(colnames(Z), nrow(Z))),
			loewe = as.numeric(c(loewe_mat))
	)
	loewe_matrix <- rbind(loewe_matrix, loewe_out)

	hsa_out <- data.frame(
			name = rep(sample_name_id, length(c(hsa_mat))),
			drugA = rep(drug_name_X2_id, length(c(hsa_mat))),
			drugB = rep(drug_name_X1_id, length(c(hsa_mat))),
			concA = as.numeric(rep(rownames(Z), each=ncol(Z))),
			concB = as.numeric(rep(colnames(Z), nrow(Z))),
			hsa = as.numeric(c(hsa_mat))
	)
	hsa_matrix <- rbind(hsa_matrix, hsa_out)

	zip_out <- data.frame(
			name = rep(sample_name_id, length(c(zip_mat))),
			drugA = rep(drug_name_X2_id, length(c(zip_mat))),
			drugB = rep(drug_name_X1_id, length(c(zip_mat))),
			concA = as.numeric(rep(rownames(Z), each=ncol(Z))),
			concB = as.numeric(rep(colnames(Z), nrow(Z))),
			zip = as.numeric(c(zip_mat))
	)
	zip_matrix <- rbind(zip_matrix, zip_out)
					

} 

colnames(outfile_stats) <- c("idSample", "idDrugA", "idDrugB", "Bliss", "Loewe", "HSA", "ZIP")
write.table(outfile_stats, file="../output/tbl.visage.combo_summary.txt", row.names=F, col.names=F, quote=F, sep="\t")

write.table(raw_matrix, file=paste("../output/tbl.visage.raw.txt",sep =""), row.names=F, col.names=T, quote=F, sep="\t")
write.table(bliss_matrix, file=paste("../output/tbl.visage.bliss.txt",sep =""), row.names=F, col.names=T, quote=F, sep="\t")
write.table(loewe_matrix, file=paste("../output/tbl.visage.loewe.txt",sep =""), row.names=F, col.names=T, quote=F, sep="\t")
write.table(hsa_matrix, file=paste("../output/tbl.visage.hsa.txt",sep =""), row.names=F, col.names=T, quote=F, sep="\t")
write.table(zip_matrix, file=paste("../output/tbl.visage.zip.txt",sep =""), row.names=F, col.names=T, quote=F, sep="\t")

quit("no")
