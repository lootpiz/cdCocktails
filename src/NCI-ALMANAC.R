library(synergyfinder)
library(PharmacoGx)

getHill <- function (x, pars) 
{
    return(pars[[2]] + (1 - pars[[2]])/(1 + (10^x/10^pars[[3]])^pars[[1]]))
}


mono_ncol <- max(count.fields("../data/NCI-ALMANAC/nci-almanac_single.txt"))
mono <- read.table("../data/NCI-ALMANAC/nci-almanac_single.txt", header=F, row.names=1, sep="\t", fill=TRUE, stringsAsFactor=F, col.names=paste0('V', seq_len(mono_ncol)))

combo <- read.table("../data/NCI-ALMANAC/nci-almanac_combo.txt", header=F, row.names=1, sep="\t", fill=TRUE, stringsAsFactor=F)
colnames(combo) <- c("cellline", "drugA", "concentrationA", "drugB", "concentrationB", "viability") 

cell_names <-  unique(mono[order(mono[,1]),1])
drug_names <- unique(mono[order(mono[,2]),2])

outfile_stats <- c()

raw_matrix <- c()
bliss_matrix <- c()
loewe_matrix <- c()
hsa_matrix <- c()
zip_matrix <- c()

for(cell_name in cell_names){
for(drug_idx in c(1:length(drug_names))){
	drug_name_X1 <- drug_names[drug_idx]
	for(drug_idx_x2 in c(1:length(drug_names))) {
		drug_name_X2 <- drug_names[drug_idx_x2]

		if (drug_name_X1 != drug_name_X2) {
			label <- paste("Cellline: ", cell_name, ", Drug_X1: ", drug_name_X1, ", Drug_X2: ", drug_name_X2, sep="")
			print(label)

			target_idx_X1 <- intersect(which(mono[,1] == cell_name), which(mono[,2] == drug_name_X1))
			if(length(target_idx_X1) > 0 ){
				cell_drug_mat_X1 <- mono[target_idx_X1, c(3:ncol(mono))]
				cell_drug_mat_X1 <- cell_drug_mat_X1[order(cell_drug_mat_X1[,1]),]

				X1_tmp <- cell_drug_mat_X1[,1]
				Y1_tmp <- cell_drug_mat_X1[,2:ncol(cell_drug_mat_X1)]
				X1 <- rep(X1_tmp, ncol(Y1_tmp))
				Y1 <- as.vector(unlist(Y1_tmp))

				one_tmp <- data.frame(concentration = X1, viability = Y1)
				one_na_idx <- which(is.na(one_tmp$viability))
				if(length(one_na_idx) > 1) {
					one <- one_tmp[-one_na_idx,]
				} else {
					one <- one_tmp
				}
			}

			target_idx_X2 <- intersect(which(mono[,1] == cell_name), which(mono[,2] == drug_name_X2))
			if(length(target_idx_X2) > 0) {
				cell_drug_mat_X2 <- mono[target_idx_X2, c(3:ncol(mono)) ]
				cell_drug_mat_X2 <- cell_drug_mat_X2[order(cell_drug_mat_X2[,1]),]
		
				X2_tmp <- cell_drug_mat_X2[,1]
				Y2_tmp <- cell_drug_mat_X2[,2:ncol(cell_drug_mat_X2)]
				X2 <- rep(X2_tmp, ncol(Y2_tmp))
				Y2 <- as.vector(unlist(Y2_tmp))

				two_tmp <- data.frame(concentration = X2, viability = Y2)
				two_na_idx <- which(is.na(two_tmp$viability))
				if(length(two_na_idx) > 1) {
					two <- two_tmp[-two_na_idx,]
				} else {
					two <- two_tmp
				}
			}

			target_idx <- intersect(which(combo$cellline == cell_name), 
						intersect(which(combo$drugA == as.character(drug_name_X1)), 
							which(combo$drugB == as.character(drug_name_X2))))
			
			if(length(target_idx) > 0) {
				combo_mat <- combo[target_idx,]
				combo_A_concentration <- unique(combo_mat$concentrationA)[order(unique(combo_mat$concentrationA))]
				combo_B_concentration <- unique(combo_mat$concentrationB)[order(unique(combo_mat$concentrationB))]
				combo_mat_z <- as.matrix(combo_mat[order(combo_mat$concentrationB, combo_mat$concentrationA), 6])
				Z_tmp <- array(combo_mat_z, dim=c(length(combo_B_concentration), length(combo_A_concentration)))
				colnames(Z_tmp) <- combo_A_concentration
				rownames(Z_tmp) <- combo_B_concentration

				params_1 <- PharmacoGx::logLogisticRegression(conc = one$concentration, viability = one$viability, viability_as_pct=FALSE)
				params_1[[3]] <- log10(params_1[[3]])
				imputed_viability_A <- getHill(log10(combo_A_concentration), unlist(params_1))

				params_2 <- PharmacoGx::logLogisticRegression(conc = X2, viability = Y2, viability_as_pct=FALSE)
				params_2[[3]] <- log10(params_2[[3]])	
				imputed_viability_B <- getHill(log10(combo_B_concentration), unlist(params_2))

				Z <- rbind(imputed_viability_A, Z_tmp)
				Z <- cbind(c(1,imputed_viability_B), Z)

				colnames(Z) <- c(0, combo_A_concentration)
				rownames(Z) <- c(0, combo_B_concentration)

				meta <- data.frame(drug.col = drug_name_X1, drug.row = drug_name_X2, concUnit = "µM", blockIDs = 1)
				data <- list(dose.response.mats = list((1-Z)*100), drug.pairs = meta)

				bliss_score <- NA; bliss_mat <- NA; loewe_score <- NA; loewe_mat <- NA
				hsa_score <- NA; hsa_mat <- NA; zip_score <- NA; zip_mat <- NA

				bliss_exe <- tryCatch({ bliss <- CalculateSynergy(data=data, method="Bliss", correction=F)
					}, error = function(err) { print("Bliss error") ; bliss <- NA
				})
				if(!is.na(bliss_exe)){
					bliss_score <- round(mean(bliss$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)
					bliss_mat <- bliss$scores[[1]]/100.0
				} else {
					bliss_mat <- matrix(rep(NA, length(c(Z))), ncol=ncol(Z), nrow=nrow(Z))
					colnames(bliss_mat) <- c(0, combo_A_concentration)
					rownames(bliss_mat) <- c(0, combo_B_concentration)
				}

				loewe_exe <- tryCatch({ loewe <- CalculateSynergy(data=data, method="Loewe", correction=F)
				}, error = function(err) { print("Loewe error") ; loewe <- NA
				})
				if(!is.na(loewe_exe)){
					loewe_score <- round(mean(loewe$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)
					loewe_mat <- loewe$scores[[1]]/100.0
				} else {
					loewe_mat <- matrix(rep(NA, length(c(Z))), ncol=ncol(Z), nrow=nrow(Z))
					colnames(loewe_mat) <- c(0, combo_A_concentration)
					rownames(loewe_mat) <- c(0, combo_B_concentration)
				}

				hsa_exe <- tryCatch({ 	hsa <- CalculateSynergy(data=data, method="HSA", correction=F)
				}, error = function(err) { print("HSA error"); hsa <- NA
				})
				if(!is.na(hsa_exe)){
					hsa_score <- round(mean(hsa$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)
					hsa_mat <- hsa$scores[[1]]/100.0
				} else {
					hsa_mat <- matrix(rep(NA, length(c(Z))), ncol=ncol(Z), nrow=nrow(Z))
					colnames(hsa_mat) <- c(0, combo_A_concentration)
					rownames(hsa_mat) <- c(0, combo_B_concentration)
				}

				zip_exe <- tryCatch({ zip <- CalculateSynergy(data=data, method="ZIP", correction=F) 
				}, error = function(err) { print("ZIP error"); zip <- NA 
				})
				if(!is.na(zip_exe)){
					zip_score <- round(mean(zip$scores[[1]][-1,-1]/100.0, na.rm=TRUE), 4)
					zip_mat <- zip$scores[[1]]/100.0
				} else {
					zip_mat <- matrix(rep(NA, length(c(Z))), ncol=ncol(Z), nrow=nrow(Z))
					colnames(zip_mat) <- c(0, combo_A_concentration)
					rownames(zip_mat) <- c(0, combo_B_concentration)
				}

				synergy <- list(sample = cell_name, drug.row = drug_name_X1, drug.col = drug_name_X2, experiments = Z,
					bliss_score = bliss_score, bliss = bliss_mat,
					loewe_score = loewe_score, loewe = loewe_mat,
					hsa_score = hsa_score, hsa = hsa_mat,
					zip_score = zip_score, zip = zip_mat
				)

				result <- c(synergy$sample, synergy$drug.row, synergy$drug.col, 
					synergy$bliss_score, synergy$loewe_score, synergy$hsa_score, synergy$zip_score)

				outfile_stats <- rbind(outfile_stats, result)

				if(drug_name_X1 < drug_name_X2) {
					Z_out <- data.frame(
						name = rep(cell_name, length(c(Z))),
						drugA = rep(drug_name_X1, length(c(Z))),
						drugB = rep(drug_name_X2, length(c(Z))),
						concA = as.numeric(rep(colnames(Z), each=nrow(Z))),
						concB = as.numeric(rep(rownames(Z), ncol(Z))),
						experiments = as.numeric(c(Z))
					)
				} else {
					Z_out <- data.frame(
						name = rep(cell_name, length(c(Z))),
						drugA = rep(drug_name_X2, length(c(Z))),
						drugB = rep(drug_name_X1, length(c(Z))),
						concA = as.numeric(rep(rownames(Z), each=ncol(Z))),
						concB = as.numeric(rep(colnames(Z), nrow(Z))),
						experiments = as.numeric(c(Z))
					)
				}
				raw_matrix <- rbind(raw_matrix, Z_out)

				if(drug_name_X1 < drug_name_X2) {
					bliss_out <- data.frame(
						name = rep(cell_name, length(c(bliss_mat))),
						drugA = rep(drug_name_X1, length(c(bliss_mat))),
						drugB = rep(drug_name_X2, length(c(bliss_mat))),
						concA = as.numeric(rep(colnames(bliss_mat), each=nrow(bliss_mat))),
						concB = as.numeric(rep(rownames(bliss_mat), ncol(bliss_mat))),
						bliss = as.numeric(c(bliss_mat))
					)
				} else {
					bliss_out <- data.frame(
						name = rep(cell_name, length(c(bliss_mat))),
						drugA = rep(drug_name_X2, length(c(bliss_mat))),
						drugB = rep(drug_name_X1, length(c(bliss_mat))),
						concA = as.numeric(rep(rownames(bliss_mat), each=ncol(bliss_mat))),
						concB = as.numeric(rep(colnames(bliss_mat), nrow(bliss_mat))),
						bliss = as.numeric(c(bliss_mat))
					)
				}
				bliss_matrix <- rbind(bliss_matrix, bliss_out)

				if(drug_name_X1 < drug_name_X2) {
					loewe_out <- data.frame(
						name = rep(cell_name, length(c(loewe_mat))),
						drugA = rep(drug_name_X1, length(c(loewe_mat))),
						drugB = rep(drug_name_X2, length(c(loewe_mat))),
						concA = as.numeric(rep(colnames(loewe_mat), each=nrow(loewe_mat))),
						concB = as.numeric(rep(rownames(loewe_mat), ncol(loewe_mat))),
						loewe = as.numeric(c(loewe_mat))
					)
				} else {
					loewe_out <- data.frame(
						name = rep(cell_name, length(c(loewe_mat))),
						drugA = rep(drug_name_X2, length(c(loewe_mat))),
						drugB = rep(drug_name_X1, length(c(loewe_mat))),
						concA = as.numeric(rep(rownames(loewe_mat), each=ncol(loewe_mat))),
						concB = as.numeric(rep(colnames(loewe_mat), nrow(loewe_mat))),
						loewe = as.numeric(c(loewe_mat))
					)
				}
				loewe_matrix <- rbind(loewe_matrix, loewe_out)

				if(drug_name_X1 < drug_name_X2) {
					hsa_out <- data.frame(
						name = rep(cell_name, length(c(hsa_mat))),
						drugA = rep(drug_name_X1, length(c(hsa_mat))),
						drugB = rep(drug_name_X2, length(c(hsa_mat))),
						concA = as.numeric(rep(colnames(hsa_mat), each=nrow(hsa_mat))),
						concB = as.numeric(rep(rownames(hsa_mat), ncol(hsa_mat))),
						hsa = as.numeric(c(hsa_mat))
					)
				} else {
					hsa_out <- data.frame(
						name = rep(cell_name, length(c(hsa_mat))),
						drugA = rep(drug_name_X2, length(c(hsa_mat))),
						drugB = rep(drug_name_X1, length(c(hsa_mat))),
						concA = as.numeric(rep(rownames(hsa_mat), each=ncol(hsa_mat))),
						concB = as.numeric(rep(colnames(hsa_mat), nrow(hsa_mat))),
						hsa = as.numeric(c(hsa_mat))
					)
				}
				hsa_matrix <- rbind(hsa_matrix, hsa_out)

				if(drug_name_X1 < drug_name_X2) {
					zip_out <- data.frame(
						name = rep(cell_name, length(c(zip_mat))),
						drugA = rep(drug_name_X1, length(c(zip_mat))),
						drugB = rep(drug_name_X2, length(c(zip_mat))),
						concA = as.numeric(rep(colnames(zip_mat), each=nrow(zip_mat))),
						concB = as.numeric(rep(rownames(zip_mat), ncol(zip_mat))),
						zip = as.numeric(c(zip_mat))
					)
				} else {
					zip_out <- data.frame(
						name = rep(cell_name, length(c(zip_mat))),
						drugA = rep(drug_name_X2, length(c(zip_mat))),
						drugB = rep(drug_name_X1, length(c(zip_mat))),
						concA = as.numeric(rep(rownames(zip_mat), each=ncol(zip_mat))),
						concB = as.numeric(rep(colnames(zip_mat), nrow(zip_mat))),
						zip = as.numeric(c(zip_mat))
					)
				}
				zip_matrix <- rbind(zip_matrix, zip_out)
			} 
		}		
	}
}	
}

colnames(outfile_stats) <- c("Cellline", "DrugA", "DrugB", "Bliss", "Loewe", "HSA", "ZIP")
write.table(outfile_stats, file="../results/tbl.nci-almanac.combo_summary.txt", row.names=F, col.names=F, quote=F, sep="\t")

write.table(raw_matrix, file="../results/tbl.nci-almanac.raw.txt", row.names=F, col.names=T, quote=F, sep="\t")
write.table(bliss_matrix, file="../results/tbl.nci-almanac.bliss.txt", row.names=F, col.names=T, quote=F, sep="\t")
write.table(loewe_matrix, file="../results/tbl.nci-almanac.loewe.txt", row.names=F, col.names=T, quote=F, sep="\t")
write.table(hsa_matrix, file="../results/tbl.nci-almanac.hsa.txt", row.names=F, col.names=T, quote=F, sep="\t")
write.table(zip_matrix, file="../results/tbl.nci-almanac.zip.txt", row.names=F, col.names=T, quote=F, sep="\t")

quit("no")


