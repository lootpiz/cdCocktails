library(synergyfinder)

input <- read.table("../data/DECREASE/decrease_input.txt", header=F, sep="\t", stringsAsFactor=F)
colnames(input) <- c("Drug1", "Drug2", "Conc1", "Conc2", "Response", "ConcUnit", "Cell")
dim(input)

master_drugs <- read.table("../data/resource/master_mapping_drug.txt", header=T, sep="\t")
master_cells <- read.table("../data/resource/master_mapping_sample.txt", header=T, sep="\t")

cell_names <- unique(input[order(input[,7]),7])
drug_names <- unique(union(input[,1], input[,2]))

outfile_stats <- c()

raw_matrix <- c()
bliss_matrix <- c()
loewe_matrix <- c()
hsa_matrix <- c()
zip_matrix <- c()


cell_name <- "BT-549"
aDrug <- "Everolimus"
bDrug <- "Dactolisib"

drug_idx <- 1
drug_idx_x2 <- 32
drug_name_X1 <- drug_names[drug_idx]
drug_name_X2 <- drug_names[drug_idx_x2]

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

				target_idx <- intersect(which(input$Cell == cell_name), 
							intersect(which(input$Drug1 == as.character(drug_name_X1)), 
								which(input$Drug2 == as.character(drug_name_X2))))
				length(target_idx)

				if(length(target_idx) > 0) {
					combo_mat <- input[target_idx,]
					combo_col_concentration <- unique(combo_mat$Conc1)[order(unique(combo_mat$Conc1))]
					combo_row_concentration <- unique(combo_mat$Conc2)[order(unique(combo_mat$Conc2))]
					Z <- array(as.numeric(as.matrix(combo_mat[order(combo_mat$Conc1, combo_mat$Conc2), 5])),
						dim=c(length(combo_row_concentration), length(combo_col_concentration)))
					colnames(Z) <- combo_col_concentration
					rownames(Z) <- combo_row_concentration

					meta <- data.frame(drug.col = drug_name_X1, drug.row = drug_name_X2, concUnit = "µM", blockIDs = 1)
					data <- list(dose.response.mats = list(Z), drug.pairs = meta)

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

	if(idDrugA < idDrugB) {
		Z_out <- data.frame(
			name = rep(sample_name_id, length(c(Z))),
			drugA = rep(drug_name_X1_id, length(c(Z))),
			drugB = rep(drug_name_X2_id, length(c(Z))),
			concA = as.numeric(rep(colnames(Z), each=nrow(Z))),
			concB = as.numeric(rep(rownames(Z), ncol(Z))),
			experiments = as.numeric(c(Z))/100
		)
	} else {
		Z_out <- data.frame(
			name = rep(sample_name_id, length(c(Z))),
			drugA = rep(drug_name_X2_id, length(c(Z))),
			drugB = rep(drug_name_X1_id, length(c(Z))),
			concA = as.numeric(rep(rownames(Z), each=ncol(Z))),
			concB = as.numeric(rep(colnames(Z), nrow(Z))),
			experiments = as.numeric(c(Z))/100
		)
	}
	raw_matrix <- rbind(raw_matrix, Z_out)

	if(idDrugA < idDrugB) {
		bliss_out <- data.frame(
			name = rep(sample_name_id, length(c(bliss_mat))),
			drugA = rep(drug_name_X1_id, length(c(bliss_mat))),
			drugB = rep(drug_name_X2_id, length(c(bliss_mat))),
			concA = as.numeric(rep(c(combo_col_concentration), each=nrow(bliss_mat))),
			concB = as.numeric(rep(c(combo_row_concentration), ncol(bliss_mat))),
			bliss = as.numeric(c(bliss_mat))
		)
	} else {
		bliss_out <- data.frame(
			name = rep(sample_name_id, length(c(bliss_mat))),
			drugA = rep(drug_name_X2_id, length(c(bliss_mat))),
			drugB = rep(drug_name_X1_id, length(c(bliss_mat))),
			concA = as.numeric(rep(c(combo_row_concentration), each=nrow(bliss_mat))),
			concB = as.numeric(rep(c(combo_col_concentration), ncol(bliss_mat))),
			bliss = as.numeric(c(bliss_mat))
		)
	}
	bliss_matrix <- rbind(bliss_matrix, bliss_out)

	if(idDrugA < idDrugB) {
		loewe_out <- data.frame(
			name = rep(sample_name_id, length(c(loewe_mat))),
			drugA = rep(drug_name_X1_id, length(c(loewe_mat))),
			drugB = rep(drug_name_X2_id, length(c(loewe_mat))),
			concA = as.numeric(rep(c(combo_col_concentration), each=nrow(loewe_mat))),
			concB = as.numeric(rep(c(combo_row_concentration), ncol(loewe_mat))),
			loewe = as.numeric(c(loewe_mat))
		)
	} else {
		loewe_out <- data.frame(
			name = rep(sample_name_id, length(c(loewe_mat))),
			drugA = rep(drug_name_X2_id, length(c(loewe_mat))),
			drugB = rep(drug_name_X1_id, length(c(loewe_mat))),
			concA = as.numeric(rep(c(combo_row_concentration), each=nrow(loewe_mat))),
			concB = as.numeric(rep(c(combo_col_concentration), ncol(loewe_mat))),
			loewe = as.numeric(c(loewe_mat))
		)
	}
	loewe_matrix <- rbind(loewe_matrix, loewe_out)

	if(idDrugA < idDrugB) {
		hsa_out <- data.frame(
			name = rep(sample_name_id, length(c(hsa_mat))),
			drugA = rep(drug_name_X1_id, length(c(hsa_mat))),
			drugB = rep(drug_name_X2_id, length(c(hsa_mat))),
			concA = as.numeric(rep(c(combo_col_concentration), each=nrow(hsa_mat))),
			concB = as.numeric(rep(c(combo_row_concentration), ncol(hsa_mat))),
			hsa = as.numeric(c(hsa_mat))
		)
	} else {
		hsa_out <- data.frame(
			name = rep(sample_name_id, length(c(hsa_mat))),
			drugA = rep(drug_name_X2_id, length(c(hsa_mat))),
			drugB = rep(drug_name_X1_id, length(c(hsa_mat))),
			concA = as.numeric(rep(c(combo_row_concentration), each=nrow(hsa_mat))),
			concB = as.numeric(rep(c(combo_col_concentration), ncol(hsa_mat))),
			hsa = as.numeric(c(hsa_mat))
		)
	}
	hsa_matrix <- rbind(hsa_matrix, hsa_out)

	if(idDrugA < idDrugB) {
		zip_out <- data.frame(
			name = rep(sample_name_id, length(c(zip_mat))),
			drugA = rep(drug_name_X1_id, length(c(zip_mat))),
			drugB = rep(drug_name_X2_id, length(c(zip_mat))),
			concA = as.numeric(rep(c(combo_col_concentration), each=nrow(zip_mat))),
			concB = as.numeric(rep(c(combo_row_concentration), ncol(zip_mat))),
			zip = as.numeric(c(zip_mat))
		)
	} else {
		zip_out <- data.frame(
			name = rep(sample_name_id, length(c(zip_mat))),
			drugA = rep(drug_name_X2_id, length(c(zip_mat))),
			drugB = rep(drug_name_X1_id, length(c(zip_mat))),
			concA = as.numeric(rep(c(combo_row_concentration), each=nrow(zip_mat))),
			concB = as.numeric(rep(c(combo_col_concentration), ncol(zip_mat))),
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
write.table(outfile_stats, file="../output/tbl.decrease.combo_summary.txt", row.names=F, col.names=F, quote=F, sep="\t")

write.table(raw_matrix, file=paste("../output/tbl.decrease.raw.txt",sep =""), row.names=F, col.names=T, quote=F, sep="\t")
write.table(bliss_matrix, file=paste("../output/tbl.decrease.bliss.txt",sep =""), row.names=F, col.names=T, quote=F, sep="\t")
write.table(loewe_matrix, file=paste("../output/tbl.decrease.loewe.txt",sep =""), row.names=F, col.names=T, quote=F, sep="\t")
write.table(hsa_matrix, file=paste("../output/tbl.decrease.hsa.txt",sep =""), row.names=F, col.names=T, quote=F, sep="\t")
write.table(zip_matrix, file=paste("../output/tbl.decrease.zip.txt",sep =""), row.names=F, col.names=T, quote=F, sep="\t")

quit("no")
