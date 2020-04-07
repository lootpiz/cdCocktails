combo_file_name <- "combo"
combo_ncol <- max(count.fields(paste(combo_file_name, ".txt", sep="")))
combo <- read.table(paste0(combo_file_name, ".txt"), header=F, sep="\t", 
	fill=TRUE, stringsAsFactor=F, col.names=paste0('V', seq_len(combo_ncol)))

median <- apply(as.matrix(combo[,c(6:12)]), 1, median, na.rm=T)

write.table(cbind(combo[,c(1:5)], median), file=paste(combo_file_name, "_input.txt", sep=""), row.names=F, col.names=F, quote=F, sep="\t")

q("no")

