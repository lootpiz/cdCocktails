library(synergyfinder)

mono <- read.table("../data/YALE-PDAC/yale-pdac_single.txt", header=T, sep="\t", stringsAsFactor=F)
combo <- read.table("../data/YALE-PDAC/yale-pdac_combo.txt", header=T, sep="\t", stringsAsFactor=F)

master_drugs <- read.table("../data/mapping/YALE-PDAC_drug.txt", header=T, sep="\t")
master_cells <- read.table("../data/mapping/YALE-PDAC_sample.txt", header=T, sep="\t")

cell_names <- unique(combo[order(combo[,2]),2])
drug_names <- unique(union(combo[,4], combo[,8]))

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

				if(idDrugA > idDrugB) {
					idTmp <- idDrugA
					idDrugA <- idDrugB
					idDrugB <- idTmp
				}

				target_idx_X1 <- intersect(which(mono[,2] == cell_name), which(mono[,4] == drug_name_X1))
				if(length(target_idx_X1) > 0 ){
					cell_drug_mat <- mono[target_idx_X1, c(5,6)]
					cell_drug_mat <- cell_drug_mat[order(cell_drug_mat$conc_first_agent),]
					one <- data.frame(concentration = cell_drug_mat$conc_first_agent, 
						viability = cell_drug_mat$first_agent_inhibition)
					one_na_idx <- which(is.na(one$viability))
					if(length(one_na_idx) > 0) {
						one <- one[-one_na_idx,]
					}
					
				}

				target_idx_X2 <- intersect(which(mono[,2] == cell_name), which(mono[,4] == drug_name_X2))
				if(length(target_idx_X2) > 0 ){
					cell_drug_mat <- mono[target_idx_X2, c(5,6)]
					cell_drug_mat <- cell_drug_mat[order(cell_drug_mat$conc_first_agent),]
					two <- data.frame(concentration = cell_drug_mat$conc_first_agent, 
						viability = cell_drug_mat$first_agent_inhibition)
					two_na_idx <- which(is.na(two$viability))
					if(length(two_na_idx) > 0) {
						two <- two[-two_na_idx,]
					}
				}

				target_idx <- intersect(which(combo$cell_line == cell_name), 
							intersect(which(combo$first_agent == as.character(drug_name_X1)), 
								which(combo$second_agent == as.character(drug_name_X2))))
		
				if(length(target_idx) > 0) {
					combo_mat <- combo[target_idx,]
					combo_A_concentration <- unique(combo_mat$conc_first_agent)[order(unique(combo_mat$conc_first_agent))]
					combo_col_concentration <- combo_A_concentration
					combo_B_concentration <- unique(combo_mat$conc_second_agent)[order(unique(combo_mat$conc_second_agent))]
					combo_row_concentration <- combo_B_concentration
					combo_mat_z <- as.matrix(combo_mat[order(combo_mat$conc_first_agent, combo_mat$conc_second_agent), 12])
					Z_tmp <- array(combo_mat_z, dim=c(length(combo_B_concentration), length(combo_A_concentration)))
					colnames(Z_tmp) <- combo_A_concentration
					rownames(Z_tmp) <- combo_B_concentration

					Z <- cbind(unlist(two$viability), Z_tmp)
					Z <- rbind(c(0, one$viability), Z)

					colnames(Z) <- c(0, combo_A_concentration)
					rownames(Z) <- c(0, combo_B_concentration)

					meta <- data.frame(drug.col = drug_name_X1, drug.row = drug_name_X2, concUnit = "µM", blockIDs = 1)
					data <- list(dose.response.mats = list(Z), drug.pairs = meta)

					bliss_score <- NA; bliss_mat <- NA; bliss <- NA
					loewe_score <- NA; loewe_mat <- NA; loewe <- NA
					hsa_score <- NA; hsa_mat <- NA; hsa <- NA
					zip_score <- NA; zip_mat <- NA; zip <- NA

					bliss_exe <- tryCatch({ bliss <- CalculateSynergy(data=data, method="Bliss", correction=F)
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

					loewe_exe <- tryCatch({ loewe <- CalculateSynergy(data=data, method="Loewe", correction=F)
					}, error = function(err) { print("Loewe error") ; loewe <- NA
					})
					if(!is.na(loewe_exe)){
						loewe_score <- round(mean(loewe$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)
						loewe_mat <- loewe$scores[[1]]/100.0
					} else {
						loewe_mat <- matrix(rep(NA, length(c(Z))), ncol=ncol(Z), nrow=nrow(Z))
						rownames(loewe_mat) <- c(0,combo_row_concentration)
						colnames(loewe_mat) <- c(0,combo_col_concentration)
					}

					hsa_exe <- tryCatch({ 	hsa <- CalculateSynergy(data=data, method="HSA", correction=F)
					}, error = function(err) { print("HSA error"); hsa <- NA
					})
					if(!is.na(hsa_exe)){
						hsa_score <- round(mean(hsa$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)
						hsa_mat <- hsa$scores[[1]]/100.0
					} else {
						hsa_mat <- matrix(rep(NA, length(c(Z))), ncol=ncol(Z), nrow=nrow(Z))
						rownames(hsa_mat) <- c(0,combo_row_concentration)
						colnames(hsa_mat) <- c(0,combo_col_concentration)
					} 

					zip_exe <- tryCatch({ zip <- CalculateSynergy(data=data, method="ZIP", correction=F) 
					}, error = function(err) { print("ZIP error"); zip <- NA 
					})
					if(!is.na(zip_exe)){
						zip_score <- round(mean(zip$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)			
						zip_mat <- zip$scores[[1]]/100.0
					} else {
						zip_mat <- matrix(rep(NA, length(c(Z))), ncol=ncol(Z), nrow=nrow(Z))
						rownames(zip_mat) <- c(0,combo_row_concentration)
						colnames(zip_mat) <- c(0,combo_col_concentration)
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
			concA = as.numeric(rep(c(0,combo_row_concentration), each=nrow(bliss_mat))),
			concB = as.numeric(rep(c(0,combo_col_concentration), ncol(bliss_mat))),
			bliss = as.numeric(c(bliss_mat))
		)
	}
	bliss_matrix <- rbind(bliss_matrix, bliss_out)

	if(idDrugA < idDrugB) {
		loewe_out <- data.frame(
			name = rep(sample_name_id, length(c(loewe_mat))),
			drugA = rep(drug_name_X1_id, length(c(loewe_mat))),
			drugB = rep(drug_name_X2_id, length(c(loewe_mat))),
			concA = as.numeric(rep(c(0,combo_col_concentration), each=nrow(loewe_mat))),
			concB = as.numeric(rep(c(0,combo_row_concentration), ncol(loewe_mat))),
			loewe = as.numeric(c(loewe_mat))
		)
	} else {
		loewe_out <- data.frame(
			name = rep(sample_name_id, length(c(loewe_mat))),
			drugA = rep(drug_name_X2_id, length(c(loewe_mat))),
			drugB = rep(drug_name_X1_id, length(c(loewe_mat))),
			concA = as.numeric(rep(c(0,combo_row_concentration), each=nrow(loewe_mat))),
			concB = as.numeric(rep(c(0,combo_col_concentration), ncol(loewe_mat))),
			loewe = as.numeric(c(loewe_mat))
		)
	}
	loewe_matrix <- rbind(loewe_matrix, loewe_out)

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
			concA = as.numeric(rep(c(0,combo_row_concentration), each=nrow(hsa_mat))),
			concB = as.numeric(rep(c(0,combo_col_concentration), ncol(hsa_mat))),
			hsa = as.numeric(c(hsa_mat))
		)
	}
	hsa_matrix <- rbind(hsa_matrix, hsa_out)

	if(idDrugA < idDrugB) {
		zip_out <- data.frame(
			name = rep(sample_name_id, length(c(zip_mat))),
			drugA = rep(drug_name_X1_id, length(c(zip_mat))),
			drugB = rep(drug_name_X2_id, length(c(zip_mat))),
			concA = as.numeric(rep(c(0,combo_col_concentration), each=nrow(zip_mat))),
			concB = as.numeric(rep(c(0,combo_row_concentration), ncol(zip_mat))),
			zip = as.numeric(c(zip_mat))
		)
	} else {
		zip_out <- data.frame(
			name = rep(sample_name_id, length(c(zip_mat))),
			drugA = rep(drug_name_X2_id, length(c(zip_mat))),
			drugB = rep(drug_name_X1_id, length(c(zip_mat))),
			concA = as.numeric(rep(c(0,combo_row_concentration), each=nrow(zip_mat))),
			concB = as.numeric(rep(c(0,combo_col_concentration), ncol(zip_mat))),
			zip = as.numeric(c(zip_mat))
		)
	}
	zip_matrix <- rbind(zip_matrix, zip_out)
					

				} 
			}		
		}
	}
}

colnames(outfile_stats) <- c("idSample", "idDrugA", "idDrugB", "Bliss", "Loewe", "HSA", "ZIP")
write.table(outfile_stats, file=paste0("../results/tbl.yale-pdac.combo_summary.txt"), row.names=F, col.names=F, quote=F, sep="\t")

write.table(raw_matrix, file=paste("../results/tbl.yale-pdac.raw.txt",sep =""), row.names=F, col.names=T, quote=F, sep="\t")
write.table(bliss_matrix, file=paste("../results/tbl.yale-pdac.bliss.txt",sep =""), row.names=F, col.names=T, quote=F, sep="\t")
write.table(loewe_matrix, file=paste("../results/tbl.yale-pdac.loewe.txt",sep =""), row.names=F, col.names=T, quote=F, sep="\t")
write.table(hsa_matrix, file=paste("../results/tbl.yale-pdac.hsa.txt",sep =""), row.names=F, col.names=T, quote=F, sep="\t")
write.table(zip_matrix, file=paste("../results/tbl.yale-pdac.zip.txt",sep =""), row.names=F, col.names=T, quote=F, sep="\t")

quit("no")


