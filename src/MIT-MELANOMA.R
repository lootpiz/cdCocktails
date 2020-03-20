library(synergyfinder)

mono <- read.table("../data/MIT-MELANOMA/mit-melanoma_single.txt", header=T, sep="\t", stringsAsFactor=F)
combo <- read.table("../data/MIT-MELANOMA/mit-melanoma_combo.txt", header=T, sep="\t", stringsAsFactor=F)

master_drugs <- read.table("../data/mapping/MIT-MELANOMA_drug.txt", header=T, sep="\t", stringsAsFactor=F)
master_cells <- read.table("../data/mapping/MIT-MELANOMA_sample.txt", header=T, sep="\t", stringsAsFactor=F)

cell_names <- unique(mono[order(mono[,1], decreasing = T),1])
drug_names <- unique(mono[order(mono[,2]),2])

outfile_stats <- c()

raw_matrix <- c()
bliss_matrix <- c()
loewe_matrix <- c()
hsa_matrix <- c()
zip_matrix <- c()

for(cell_name in cell_names) {
	for(drug_idx in c(1:length(drug_names))){
		drug_name_X1 <- drug_names[drug_idx]
		for(drug_idx_x2 in c(1:length(drug_names))) {
			drug_name_X2 <- drug_names[drug_idx_x2]

			if (drug_name_X1 != drug_name_X2) {
				label <- paste("Cellline: ", cell_name, ", Drug_X1: ", drug_name_X1, ", Drug_X2: ", drug_name_X2, sep="")
				print(label)

				idSample <- master_cells$id[which(master_cells$name == cell_name)]
				idDrugA <- master_drugs$id[which(master_drugs$name == drug_name_X1)]
				idDrugB <- master_drugs$id[which(master_drugs$name == drug_name_X2)]
				sample_name_id <- idSample
				drug_name_X1_id <- idDrugA
				drug_name_X2_id <- idDrugB

				target_idx_X1 <- intersect(which(mono[,1] == cell_name), which(mono[,2] == drug_name_X1))
				if(length(target_idx_X1) > 0 ){
					cell_drug_mat <- mono[target_idx_X1, c(3,4)]
					cell_drug_mat <- cell_drug_mat[order(cell_drug_mat$Dose),]
					one <- data.frame(concentration = cell_drug_mat$Dose, 
						viability = cell_drug_mat$Inhibition)
					one_na_idx <- which(is.na(one$viability))
					if(length(one_na_idx) > 0) {
						one <- one[-one_na_idx,]
					}
					
				}

				target_idx_X2 <- intersect(which(mono[,1] == cell_name), which(mono[,2] == drug_name_X2))
				if(length(target_idx_X2) > 0 ){
					cell_drug_mat <- mono[target_idx_X2, c(3,4)]
					cell_drug_mat <- cell_drug_mat[order(cell_drug_mat$Dose),]
					two <- data.frame(concentration = cell_drug_mat$Dose, 
						viability = cell_drug_mat$Inhibition)
					two_na_idx <- which(is.na(two$viability))
					if(length(two_na_idx) > 0) {
						two <- two[-two_na_idx,]
					}
				}

				target_idx <- intersect(which(combo$Cellline == cell_name), 
							intersect(which(combo$DrugA == as.character(drug_name_X1)), 
								which(combo$DrugB == as.character(drug_name_X2))))
		
				if(length(target_idx) > 0) {
					combo_mat <- combo[target_idx,]
					combo_A_concentration <- unique(combo_mat$DoseA)[order(unique(combo_mat$DoseA))]
					combo_col_concentration <- combo_A_concentration
					combo_B_concentration <- unique(combo_mat$DoseB)[order(unique(combo_mat$DoseB))]
					combo_row_concentration <- combo_B_concentration
					combo_mat_z <- as.matrix(combo_mat[order(combo_mat$DoseA, combo_mat$DoseB), 6])
					Z_tmp <- array(combo_mat_z, dim=c(length(combo_B_concentration), length(combo_A_concentration)))
					if(length(unlist(Z_tmp)) > 1) {
						Z_tmp[1,2] <- NA
						Z_tmp[2,1] <- NA
					}
					colnames(Z_tmp) <- combo_A_concentration
					rownames(Z_tmp) <- combo_B_concentration

					Z <- cbind(unlist(two$viability), Z_tmp)
					Z <- rbind(c(100, one$viability[which(one$concentration == colnames(Z_tmp))]), Z)

					colnames(Z) <- c(0, combo_A_concentration)
					rownames(Z) <- c(0, combo_B_concentration)

					meta <- data.frame(drug.col = drug_name_X1, drug.row = drug_name_X2, concUnit = "µM", blockIDs = 1)
					data <- list(dose.response.mats = list(Z), drug.pairs = meta)

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
						rownames(bliss_mat) <- c(0,combo_row_concentration)
						colnames(bliss_mat) <- c(0,combo_col_concentration)
					}

					hsa_exe <- tryCatch({ hsa <- CalculateSynergy(data=data, correction=F, method="HSA")
						}, error = function(err) { print("HSA error") ; hsa <- NA
					})
					if(!is.na(hsa_exe)){
						hsa_score <- round(mean(hsa$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)
						hsa_mat <- hsa$scores[[1]]/100.0
					} else {
						hsa_mat <- matrix(rep(NA, length(c(Z))), ncol=ncol(Z), nrow=nrow(Z))
						rownames(hsa_mat) <- c(0,combo_row_concentration)
						colnames(hsa_mat) <- c(0,combo_col_concentration)
					}


					result <- c(idSample, idDrugA, idDrugB, bliss_score, loewe_score, hsa_score, zip_score)

					outfile_stats <- rbind(outfile_stats, result)

	if(idDrugA < idDrugB) {
		Z_out <- data.frame(
			name = rep(sample_name_id, length(c(Z))),
			drugA = rep(drug_name_X1_id, length(c(Z))),
			drugB = rep(drug_name_X2_id, length(c(Z))),
			concA = as.numeric(rep(colnames(Z), each=nrow(Z))),
			concB = as.numeric(rep(rownames(Z), ncol(Z))),
			experiments = as.numeric(c(Z))
		)
	} else {
		Z_out <- data.frame(
			name = rep(sample_name_id, length(c(Z))),
			drugA = rep(drug_name_X2_id, length(c(Z))),
			drugB = rep(drug_name_X1_id, length(c(Z))),
			concA = as.numeric(rep(rownames(Z), each=ncol(Z))),
			concB = as.numeric(rep(colnames(Z), nrow(Z))),
			experiments = as.numeric(c(Z))
		)
	}
	raw_matrix <- rbind(raw_matrix, Z_out)

	if(idDrugA < idDrugB) {
		bliss_out <- data.frame(
			name = rep(sample_name_id, length(c(bliss_mat))),
			drugA = rep(drug_name_X1_id, length(c(bliss_mat))),
			drugB = rep(drug_name_X2_id, length(c(bliss_mat))),
			concA = as.numeric(rep(c(0,combo_col_concentration), each=nrow(bliss_mat))),
			concB = as.numeric(rep(c(0,combo_row_concentration), ncol(bliss_mat))),
			bliss = as.numeric(c(bliss_mat))
		)
	} else {
		bliss_out <- data.frame(
			name = rep(sample_name_id, length(c(bliss_mat))),
			drugA = rep(drug_name_X2_id, length(c(bliss_mat))),
			drugB = rep(drug_name_X1_id, length(c(bliss_mat))),
			concA = as.numeric(rep(c(0,combo_row_concentration), ncol(bliss_mat))),
			concB = as.numeric(rep(c(0,combo_col_concentration), each=nrow(bliss_mat))),
			bliss = as.numeric(c(bliss_mat))
		)
	}
	bliss_matrix <- rbind(bliss_matrix, bliss_out)
	
	if(idDrugA < idDrugB) {
		hsa_out <- data.frame(
			name = rep(sample_name_id, length(c(hsa_mat))),
			drugA = rep(drug_name_X1_id, length(c(hsa_mat))),
			drugB = rep(drug_name_X2_id, length(c(hsa_mat))),
			concA = as.numeric(rep(c(0,combo_col_concentration), each=nrow(hsa_mat))),
			concB = as.numeric(rep(c(0,combo_row_concentration), ncol(hsa_mat))),
			hsa = as.numeric(c(hsa_mat))
		)
	} else {
		hsa_out <- data.frame(
			name = rep(sample_name_id, length(c(hsa_mat))),
			drugA = rep(drug_name_X2_id, length(c(hsa_mat))),
			drugB = rep(drug_name_X1_id, length(c(hsa_mat))),
			concA = as.numeric(rep(c(0,combo_row_concentration), ncol(hsa_mat))),
			concB = as.numeric(rep(c(0,combo_col_concentration), each=nrow(hsa_mat))),
			hsa = as.numeric(c(hsa_mat))
		)
	}
	hsa_matrix <- rbind(hsa_matrix, hsa_out)
	


				} 
			}		
		}
	}
	
}

colnames(outfile_stats) <- c("idSample", "idDrugA", "idDrugB", "Bliss", "Loewe", "HSA", "ZIP")
write.table(outfile_stats, file=paste0("../results/tbl.mit-melanoma.combo_summary.txt"), row.names=F, col.names=F, quote=F, sep="\t")

write.table(raw_matrix, file=paste("../results/tbl.mit-melanoma.raw.txt",sep =""), row.names=F, col.names=T, quote=F, sep="\t")
write.table(bliss_matrix, file=paste("../results/tbl.mit-melanoma.bliss.txt",sep =""), row.names=F, col.names=T, quote=F, sep="\t")
write.table(hsa_matrix, file=paste("../results/tbl.mit-melanoma.hsa.txt",sep =""), row.names=F, col.names=T, quote=F, sep="\t")

quit("no")


