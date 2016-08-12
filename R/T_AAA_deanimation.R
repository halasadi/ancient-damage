

signature_counts <- get(load("../summary_data/signature-counts-clubbed-Lindo2016.rda"))

sig_set <- colnames(signature_counts);

sig_split <- do.call(rbind, lapply(colnames(signature_counts), function(x) strsplit(as.character(x), split="")[[1]]))

indices1 <- which(sig_split[,3]=="T" & sig_split[,6]=="C")
indices2 <- which(sig_split[,3]=="C" & sig_split[,6]=="T")

signature_counts_sub <- signature_counts[, -c(indices1,indices2)];

pr <- prcomp(t(limma::voom(t(signature_counts_sub))$E))

pc_data_frame <- data.frame("PC"=pr$x,
                            "labels"=c(rep("Ancient",25),
                                       rep("Modern",25)))

qplot(PC.PC1, PC.PC2,
      data=pc_data_frame,
      colour=labels)
