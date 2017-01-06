

library(aRchaic)
gossling_data <- get(load("../processed_data/annagosling2016-counts-table.rda"))
system.time(gossling_data_clubbed <- club_signature_counts(gossling_data))
gossling_data_clubbed <- gossling_data_clubbed[-28,];

modern_data <- get(load("../processed_data/1000Gmoderns-counts-table.rda"));
modern_data_clubbed <- club_signature_counts(modern_data)

pooled_names <- intersect(colnames(gossling_data_clubbed), colnames(modern_data_clubbed))
filtered_gossling <- gossling_data_clubbed[, match(pooled_names, colnames(gossling_data_clubbed))]
filtered_moderns <- modern_data_clubbed[, match(pooled_names, colnames(modern_data_clubbed))]

pooled_data <- rbind(filtered_gossling, filtered_moderns)

signature_set <- colnames(pooled_data)
sig_split <- t(sapply(1:length(signature_set), function(x) return(strsplit(signature_set[x], "")[[1]][1:8])))
new_sig_split <- matrix(0, dim(sig_split)[1], 5);
new_sig_split[,1] <- sig_split[,1]
new_sig_split[,2] <- sig_split[,2]
new_sig_split[,3] <- sapply(1:length(signature_set), function(x) return(paste(sig_split[x,3:6], collapse="")))
new_sig_split[,4] <- sig_split[,7]
new_sig_split[,5] <- sig_split[,8]

levels(new_sig_split[,1]) <- c("0", "1", "2", "3", "4")

pos <- t(sapply(1:length(signature_set), function(x)
{
  y = strsplit(signature_set[x], "")[[1]]
  return(paste(y[10:length(y)], collapse=""))
}))



mat <- matrix(0, dim(new_sig_split)[1], dim(new_sig_split)[2])
for(k in 1:dim(new_sig_split)[2]){
  temp <- as.factor(new_sig_split[,k])
  mat[,k] <- as.numeric(as.matrix(plyr::mapvalues(temp, from = levels(temp), to = 0:(length(levels(temp))-1))))
}

pos <- as.numeric(pos)
pos <- pos - min(pos)
pos <- factor(pos, levels = 0:22)

signatures <- mat;
signature_pos <- cbind.data.frame(signatures, pos)


out <- topics(pooled_data, K=5, tol=100, model="independent", signatures = signature_pos)
save(out, file="../processed_data/maptpx-runs/gosling2016-1000g-maptpx-independent-K-5.rda")

out <- get(load("../processed_data/maptpx-runs/gosling2016-1000g-maptpx-independent-K-5.rda"))


