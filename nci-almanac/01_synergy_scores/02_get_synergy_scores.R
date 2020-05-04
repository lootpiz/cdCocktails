idSample <- 42
idDrugA <- 11
idDrugB <- 97


library(synergyfinder)
library(PharmacoGx)

getHill <- function (x, pars) 
{
    return(pars[[2]] + (1 - pars[[2]])/(1 + (10^x/10^pars[[3]])^pars[[1]]))
}

mono_ncol <- max(count.fields("mono_input.txt"))
mono <- read.table("mono_input.txt", header=F, sep="\t", fill=TRUE, stringsAsFactor=F, col.names=paste0('V', seq_len(mono_ncol)))

combo <- read.table("combo_input.txt", header=F, sep="\t", fill=TRUE, stringsAsFactor=F)
colnames(combo) <- c("cellline", "drugA", "concentrationA", "drugB", "concentrationB", "viability") 

drugA_idx <- intersect(which(mono[,1] == idSample), which(mono[,2] == idDrugA))
if(length(drugA_idx) > 0) {
	drugA_mat <- mono[drugA_idx, c(3:ncol(mono))]
	drugA_mat <- drugA_mat[order(drugA_mat[,1]),] # order by concentration

	drugA_concentration <- rep(drugA_mat[,1], ncol(drugA_mat[,c(2:ncol(drugA_mat))]))
	drugA_viability <- as.vector(unlist(drugA_mat[,c(2:ncol(drugA_mat))]))
	drugA <- data.frame(concentration = drugA_concentration, viability = drugA_viability)
	drugA_NA_idx <- which(is.na(drugA$viability))
	if(length(drugA_NA_idx) > 1) {
		drugA <- drugA[-drugA_NA_idx,]
	}
}

drugB_idx <- intersect(which(mono[,1] == idSample), which(mono[,2] == idDrugB))
if(length(drugB_idx) > 0) {
	drugB_mat <- mono[drugB_idx, c(3:ncol(mono))]
	drugB_mat <- drugA_mat[order(drugB_mat[,1]),]

	drugB_concentration <- rep(drugB_mat[,1], ncol(drugB_mat[,c(2:ncol(drugB_mat))]))
	drugB_viability <- as.vector(unlist(drugB_mat[,c(2:ncol(drugB_mat))]))
	drugB <- data.frame(concentration = drugB_concentration, viability = drugB_viability)
	drugB_NA_idx <- which(is.na(drugB$viability))
	if(length(drugB_NA_idx) > 1) {
		drugB <- drugB[-drugB_NA_idx,]
	}
}

combo_idx <- intersect(which(combo$cellline == idSample), 
	intersect(which(combo$drugA == as.character(idDrugA)), 
		which(combo$drugB == as.character(idDrugB))))
		
if(length(combo_idx) > 0) {
	combo_mat <- combo[combo_idx,]
	combo_A_concentration <- unique(combo_mat$concentrationA)[order(unique(combo_mat$concentrationA))]
	combo_B_concentration <- unique(combo_mat$concentrationB)[order(unique(combo_mat$concentrationB))]
	combo_mat_z <- as.matrix(combo_mat[order(combo_mat$concentrationB, combo_mat$concentrationA), 6])

	Z_tmp <- array(combo_mat_z, dim=c(length(combo_B_concentration), length(combo_A_concentration)))
	colnames(Z_tmp) <- combo_A_concentration
	rownames(Z_tmp) <- combo_B_concentration

	params_1 <- PharmacoGx::logLogisticRegression(conc = drugA$concentration, viability = drugA$viability, viability_as_pct=FALSE)
	params_1[[3]] <- log10(params_1[[3]])
	imputed_viability_drugA <- PharmacoGx::getHill(log10(combo_A_concentration), unlist(params_1))
	imputed_viability_drugA

	params_2 <- logLogisticRegression(conc = drugB$concentration, viability = drugB$viability, viability_as_pct=FALSE)
	params_2[[3]] <- log10(params_2[[3]])	
	imputed_viability_drugB <- getHill(log10(combo_B_concentration), unlist(params_2))
	imputed_viability_drugB

	Z <- rbind(imputed_viability_drugA, Z_tmp)
	Z <- cbind(c(1,imputed_viability_drugB), Z)

	colnames(Z) <- c(0, combo_A_concentration)
	rownames(Z) <- c(0, combo_B_concentration)

	combo_col_concentration <- colnames(Z)
	combo_row_concentration <- rownames(Z)

	meta <- data.frame(drug.col = idDrugA, drug.row = idDrugB, concUnit = "M", blockIDs = 1)
	data <- list(dose.response.mats = list(block=Z), drug.pairs = meta)

	bliss_score <- NA; bliss_mat <- NA
	loewe_score <- NA; loewe_mat <- NA
	hsa_score <- NA; hsa_mat <- NA
	zip_score <- NA; zip_mat <- NA

	bliss_exe <- tryCatch({ bliss <- CalculateSynergy(data=data, method="Bliss", adjusted=F)
		}, error = function(err) { print("Bliss error") ; bliss <- NA
	})
	if(!is.na(bliss_exe)){
		bliss_score <- round(mean(bliss$scores[[1]][-1,-1], na.rm=TRUE), 4)
		bliss_mat <- bliss$scores[[1]]
	} else {
		bliss_mat <- matrix(rep(NA, length(c(Z))), ncol=ncol(Z), nrow=nrow(Z))
		colnames(bliss_mat) <- c(combo_col_concentration)
		rownames(bliss_mat) <- c(combo_row_concentration)
	}

	loewe_exe <- tryCatch({ loewe <- CalculateSynergy(data=data, method="Loewe", adjusted=F)
		}, error = function(err) { print("Loewe error") ; loewe <- NA
	})
	if(!is.na(loewe_exe)){
		loewe_score <- round(mean(loewe$scores[[1]][-1,-1], na.rm=TRUE), 4)
		loewe_mat <- loewe$scores[[1]]
	} else {
		loewe_mat <- matrix(rep(NA, length(c(Z))), ncol=ncol(Z), nrow=nrow(Z))
		colnames(loewe_mat) <- c(combo_col_concentration)
		rownames(loewe_mat) <- c(combo_row_concentration)
	}

	hsa_exe <- tryCatch({ 	hsa <- CalculateSynergy(data=data, method="HSA", adjusted=F)
		}, error = function(err) { print("HSA error"); hsa <- NA
	})
	if(!is.na(hsa_exe)){
		hsa_score <- round(mean(hsa$scores[[1]][-1,-1], na.rm=TRUE), 4)
		hsa_mat <- hsa$scores[[1]]
	} else {
		hsa_mat <- matrix(rep(NA, length(c(Z))), ncol=ncol(Z), nrow=nrow(Z))
		colnames(hsa_mat) <- c(combo_col_concentration)
		rownames(hsa_mat) <- c(combo_row_concentration)
	} 

	zip_exe <- tryCatch({ zip <- CalculateSynergy(data=data, method="ZIP", adjusted=F) 
		}, error = function(err) { print("ZIP error"); zip <- NA 
	})
	if(!is.na(zip_exe)){
		zip_score <- round(mean(zip$scores[[1]][-1,-1], na.rm=TRUE), 4)			
		zip_mat <- zip$scores[[1]]
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

	write.table(t(result), file=paste0("/u/hseo/research/almanac/results/s", idSample, "/tbl.nci.combo_summary_", idSample, "__", idDrugA, "__", idDrugB, ".txt"),
		row.names=F, col.names=F, quote=F, sep="\t")
	write.table(raw_matrix, file=paste0("/u/hseo/research/almanac/results/s", idSample, "/tbl.nci.raw_", idSample, "__", idDrugA, "__", idDrugB, ".txt"),
		row.names=F, col.names=F, quote=F, sep="\t")
	write.table(bliss_matrix, file=paste0("/u/hseo/research/almanac/results/s", idSample, "/tbl.nci.bliss_", idSample, "__", idDrugA, "__", idDrugB, ".txt"),
		row.names=F, col.names=F, quote=F, sep="\t")
	write.table(loewe_matrix, file=paste0("/u/hseo/research/almanac/results/s", idSample, "/tbl.nci.loewe_", idSample, "__", idDrugA, "__", idDrugB, ".txt"),
		row.names=F, col.names=F, quote=F, sep="\t")
	write.table(hsa_matrix, file=paste0("/u/hseo/research/almanac/results/s", idSample, "/tbl.nci.hsa_", idSample, "__", idDrugA, "__", idDrugB, ".txt"),
		row.names=F, col.names=F, quote=F, sep="\t")
	write.table(zip_matrix, file=paste0("/u/hseo/research/almanac/results/s", idSample, "/tbl.nci.zip_", idSample, "__", idDrugA, "__", idDrugB, ".txt"),
		row.names=F, col.names=F, quote=F, sep="\t")

}	

q("no")
