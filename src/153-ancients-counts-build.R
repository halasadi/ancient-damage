

#########   153 ancient samples counts matrix building  #####################

ancient_files <- list.files("../summary_data/153-ancients-signature-counts/csv/")

length(ancient_files)

signature_ancient <- vector(mode="list")
signature_counts_ancient <- vector(mode="list")


for(num in 1:length(ancient_files)){
  data <- read.csv(paste0("../summary_data/153-ancients-signature-counts/csv/",ancient_files[num]), header=FALSE)
  signature_counts_ancient[[num]] <- data[,2];
  signature_ancient[[num]] <- as.character(data[,1]);
}

merged_signature_ancient <- signature_ancient[[1]]

for(num in 2:length(ancient_files)){
  merged_signature_ancient <- union(merged_signature_ancient, signature_ancient[[num]])
}

ancient_counts <- matrix(0, length(ancient_files), length(merged_signature_ancient))

for(num in 1:length(ancient_files)){
  ancient_counts[num, match(signature_ancient[[num]], merged_signature_ancient)] <- signature_counts_ancient[[num]]
}

ancient_files_filt <- as.character(sapply(ancient_files, function(x) strsplit(x, ".csv")[[1]][1]))

rownames(ancient_counts) <- ancient_files_filt

signature_split <- do.call(rbind, lapply(merged_signature_ancient, function(x) strsplit(as.character(x), split="")[[1]]))

indices1 <- which(signature_split[,3]==signature_split[,6])

indices2 <- numeric()
for(m in 1:8){
  indices2 <- c(indices2, which(signature_split[,m]=="N"));
}

indices <- union(indices1, indices2)

ancient_counts_filtered <- ancient_counts[,-indices]
rownames(ancient_counts_filtered) <- ancient_files_filt
colnames(ancient_counts_filtered) <- merged_signature_ancient[-indices]


save(ancient_counts_filtered,
            file="../summary_data/153-ancients-counts-table.rda")


