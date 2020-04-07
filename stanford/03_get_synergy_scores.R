idSample <- 186
idDrugA <- 64
idDrugB <- 92

library(synergyfinder)

combo <- read.table("combo_input.txt", header=F, sep="\t", stringsAsFactor=F)
colnames(combo) <- c("cellline", "drugA", "cell0", "cell1", "cell2", "cell3") 

combo_idx <- intersect(which(combo$cellline == idSample), which(combo$drugA == as.character(idDrugA)))

if(length(combo_idx) > 0) {
	combo_mat <- combo[combo_idx,]
	Z <- array(c(as.numeric(combo_mat[,c(3:6)])), dim=c(2,2))
	colnames(Z) <- c(0, 400)
	rownames(Z) <- c(0, 5)

	combo_col_concentration <- colnames(Z)
	combo_row_concentration <- rownames(Z)

	meta <- data.frame(drug.col = idDrugA, drug.row = idDrugB, concUnit = "microM", blockIDs = 1)
	data <- list(dose.response.mats = list(block=(1-Z)), drug.pairs = meta) # Inhibition rate!
	data

	bliss_score <- NA; bliss_mat <- NA
	loewe_score <- NA; loewe_mat <- NA
	hsa_score <- NA; hsa_mat <- NA
	zip_score <- NA; zip_mat <- NA

	bliss_exe <- tryCatch({ bliss <- CalculateSynergy(data=data, method="Bliss", adjusted=F)
		}, error = function(err) { bliss <- NA
	})
	if(!is.na(bliss_exe)){
		bliss_score <- round(mean(bliss$scores$block[-1,-1], na.rm=TRUE), 4)
		bliss_mat <- bliss$scores$block
	} else {
		bliss_mat <- matrix(rep(NA, length(c(Z))), ncol=ncol(Z), nrow=nrow(Z))
		colnames(bliss_mat) <- c(combo_col_concentration)
		rownames(bliss_mat) <- c(combo_row_concentration)
	}

	loewe_exe <- tryCatch({ loewe <- CalculateSynergy(data=data, method="Loewe", adjusted=F)
		}, error = function(err) { loewe <- NA
	})
	if(!is.na(loewe_exe)){
		loewe_score <- round(mean(loewe$scores$block[-1,-1], na.rm=TRUE), 4)
		loewe_mat <- loewe$scores$block
	} else {
		loewe_mat <- matrix(rep(NA, length(c(Z))), ncol=ncol(Z), nrow=nrow(Z))
		colnames(loewe_mat) <- c(combo_col_concentration)
		rownames(loewe_mat) <- c(combo_row_concentration)
	}

	hsa_exe <- tryCatch({ 	hsa <- CalculateSynergy(data=data, method="HSA", adjusted=F)
		}, error = function(err) { hsa <- NA
	})
	if(!is.na(hsa_exe)){
		hsa_score <- round(mean(hsa$scores$block[-1,-1], na.rm=TRUE), 4)
		hsa_mat <- hsa$scores$block
	} else {
		hsa_mat <- matrix(rep(NA, length(c(Z))), ncol=ncol(Z), nrow=nrow(Z))
		colnames(hsa_mat) <- c(combo_col_concentration)
		rownames(hsa_mat) <- c(combo_row_concentration)
	} 

	zip_exe <- tryCatch({ zip <- CalculateSynergy(data=data, method="ZIP", adjusted=F)
		}, error = function(err) { zip <- NA 
	})
	if(!is.na(zip_exe)){
		zip_score <- round(mean(zip$scores$block[-1,-1], na.rm=TRUE), 4)			
		zip_mat <- zip$scores$block
	} else {
		zip_mat <- matrix(rep(NA, length(c(Z))), ncol=ncol(Z), nrow=nrow(Z))
		colnames(zip_mat) <- c(combo_col_concentration)
		rownames(zip_mat) <- c(combo_row_concentration)
	}

	result <- c(idSample, idDrugA, idDrugB, bliss_score, loewe_score, hsa_score, zip_score)
	names(result) <- c("Cellline", "DrugA", "DrugB", "Bliss", "Loewe", "HSA", "ZIP")


	raw_matrix <- data.frame(
		name = rep(idSample, length(c(Z))),
		drugA = rep(idDrugA, length(c(Z))),
		drugB = rep(idDrugB, length(c(Z))),
		concA = as.numeric(rep(colnames(Z), each=nrow(Z))),
		concB = as.numeric(rep(rownames(Z), ncol(Z))),
		experiments = as.numeric(c(Z))
	)
	bliss_matrix <- data.frame(
		name = rep(idSample, length(c(bliss_mat))),
		drugA = rep(idDrugA, length(c(bliss_mat))),
		drugB = rep(idDrugB, length(c(bliss_mat))),
		concA = as.numeric(rep(colnames(Z), each=nrow(Z))),
		concB = as.numeric(rep(rownames(Z), ncol(Z))),
		bliss = as.numeric(c(bliss_mat))
	)

	loewe_matrix <- data.frame(
		name = rep(idSample, length(c(loewe_mat))),
		drugA = rep(idDrugA, length(c(loewe_mat))),
		drugB = rep(idDrugB, length(c(loewe_mat))),
		concA = as.numeric(rep(colnames(Z), each=nrow(Z))),
		concB = as.numeric(rep(rownames(Z), ncol(Z))),
		loewe = as.numeric(c(loewe_mat))
	)

	hsa_matrix <- data.frame(
		name = rep(idSample, length(c(hsa_mat))),
		drugA = rep(idDrugA, length(c(hsa_mat))),
		drugB = rep(idDrugB, length(c(hsa_mat))),
		concA = as.numeric(rep(colnames(Z), each=nrow(Z))),
		concB = as.numeric(rep(rownames(Z), ncol(Z))),
		hsa = as.numeric(c(hsa_mat))
	)

	zip_matrix <- data.frame(
		name = rep(idSample, length(c(zip_mat))),
		drugA = rep(idDrugA, length(c(zip_mat))),
		drugB = rep(idDrugB, length(c(zip_mat))),
		concA = as.numeric(rep(colnames(Z), each=nrow(Z))),
		concB = as.numeric(rep(rownames(Z), ncol(Z))),
		zip = as.numeric(c(zip_mat))
	)

	write.table(t(result), file=paste0("tbl.stanford.combo_summary_", idSample, "__", idDrugA, "__", idDrugB, ".txt"),
		row.names=F, col.names=F, quote=F, sep="\t")
	write.table(raw_matrix, file=paste0("tbl.stanford.raw_", idSample, "__", idDrugA, "__", idDrugB, ".txt"),
		row.names=F, col.names=F, quote=F, sep="\t")
	write.table(bliss_matrix, file=paste0("tbl.stanford.bliss_", idSample, "__", idDrugA, "__", idDrugB, ".txt"),
		row.names=F, col.names=F, quote=F, sep="\t")
	write.table(loewe_matrix, file=paste0("tbl.stanford.loewe_", idSample, "__", idDrugA, "__", idDrugB, ".txt"),
		row.names=F, col.names=F, quote=F, sep="\t")
	write.table(hsa_matrix, file=paste0("tbl.stanford.hsa_", idSample, "__", idDrugA, "__", idDrugB, ".txt"),
		row.names=F, col.names=F, quote=F, sep="\t")
	write.table(zip_matrix, file=paste0("tbl.stanford.zip_", idSample, "__", idDrugA, "__", idDrugB, ".txt"),
		row.names=F, col.names=F, quote=F, sep="\t")

}	

q("no")
