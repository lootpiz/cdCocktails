idDrugA <-
idDrugB <-
idSource <-
nameSource <-

library(wCI)

prefix <- "/Users/lootpiz/SYNERGxDB/biomarkers/"

data <- read.table(paste0(prefix, "profiles/", sourceName, "/p_", idDrugA, "_", idDrugB, ".txt"), header=F, stringsAsFactor=F, sep="\t")
colnames(data) <- c("idSample", "idDrugA", "idDrugB", "bliss", "loewe", "hsa", "zip", "gene", "fpkm")
genes <- unique(data$gene)

bliss_out <- c()
loewe_out <- c()
hsa_out <- c()
zip_out <- c()

for(gene in genes) {
	tuple_idx <- intersect(which(data$gene == gene),
		intersect(which(data$idDrugA == idDrugA), which(data$idDrugB == idDrugB)))
	subset <- data[tuple_idx,]

	for(method_idx in c(4:7)) {
		obj <- try(paired.concordance.index(subset$fpkm, subset[,method_idx], delta.obs=0, delta.pred=0), silent=T)
		if(is(obj, "try-error")) {
			concordanceIndex <- NA
			pValue <- NA
		} else {
			concordanceIndex <- obj$cindex
			pValue <- obj$p.value
		}
		result <- c(idSource, idDrugA, idDrugB, gene, concordanceIndex, pValue)

		if(method_idx == 4) {				
			bliss_out <- rbind(bliss_out, result)
		} else if(method_idx == 5) {
			loewe_out <- rbind(loewe_out, result)
		} else if(method_idx == 6) {
			hsa_out <- rbind(hsa_out, result)
		} else {
			zip_out <- rbind(zip_out, result)
		}
	}
}


write.table(bliss_out, file=paste0(prefix, sourceName, "/Biomarker__Bliss__", idSource, "__", idDrugA, "__", idDrugB, ".txt"), row.names=F, col.names=F, quote=F, sep="\t")
write.table(loewe_out, file=paste0(prefix, sourceName, "/Biomarker__Loewe__", idSource, "__", idDrugA, "__", idDrugB, ".txt"), row.names=F, col.names=F, quote=F, sep="\t")
write.table(hsa_out, file=paste0(prefix, sourceName, "/Biomarker__HSA__", idSource, "__", idDrugA, "__", idDrugB, ".txt"), row.names=F, col.names=F, quote=F, sep="\t")
write.table(zip_out, file=paste0(prefix, sourceName, "/Biomarker__ZIP__", idSource, "__", idDrugA, "__", idDrugB, ".txt"), row.names=F, col.names=F, quote=F, sep="\t")

q("no")
