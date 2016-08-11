
###  Build signature counts matrix (Moderns + Ancients)

modern_files <- list.files("../summary_data/modern-signature-counts/")
ancient_files <- list.files("../summary_data/ancient-signature-counts/")

## check how many moderns and how many ancients we have 

length(modern_files)
length(ancient_files)

## Should be 25 for Lindo2016

signature_modern <- vector(mode="list")
signature_ancient <- vector(mode="list")
signature_counts_modern <- vector(mode="list")
signature_counts_ancient <- vector(mode="list")

for(num in 1:length(modern_files)){
  data <- read.csv(paste0("../summary_data/modern-signature-counts/",modern_files[num]), header=FALSE)
  signature_counts_modern[[num]] <- data[,2];
  signature_modern[[num]] <- as.character(data[,1]);
}

for(num in 1:length(ancient_files)){
  data <- read.csv(paste0("../summary_data/ancient-signature-counts/",ancient_files[num]), header=FALSE)
  signature_counts_ancient[[num]] <- data[,2];
  signature_ancient[[num]] <- as.character(data[,1]);
}

merged_signature_ancient <- signature_ancient[[1]]
merged_signature_modern <- signature_modern[[1]]
for(num in 2:length(ancient_files)){
  merged_signature_ancient <- union(merged_signature_ancient, signature_ancient[[num]])
}

for(num in 2:length(modern_files)){
  merged_signature_modern <- union(merged_signature_modern, signature_modern[[num]])
}

merged_signature <- union(merged_signature_modern, merged_signature_ancient)

modern_counts <- matrix(0, length(modern_files), length(merged_signature))
ancient_counts <- matrix(0, length(ancient_files), length(merged_signature))

for(num in 1:length(modern_files)){
  modern_counts[num, match(signature_modern[[num]], merged_signature)] <- signature_counts_modern[[num]]
}

for(num in 1:length(ancient_files)){
  ancient_counts[num, match(signature_ancient[[num]], merged_signature)] <- signature_counts_ancient[[num]]
}


signature_counts_data_pooled <- rbind(ancient_counts, modern_counts)

all_files <- c(ancient_files, modern_files)
all_files_filt <- as.character(sapply(all_files, function(x) strsplit(x, ".csv")[[1]][1]))
rownames(signature_counts_data_pooled) <- all_files_filt
colnames(signature_counts_data_pooled) <- merged_signature

merged_signature_split <- do.call(rbind, lapply(merged_signature, function(x) strsplit(as.character(x), split="")[[1]]))

indices1 <- which(merged_signature_split[,3]==merged_signature_split[,6])

indices2 <- numeric()
for(m in 1:8){
  indices2 <- c(indices2, which(merged_signature_split[,m]=="N"));
}

indices <- union(indices1, indices2)

signature_counts_data_pooled_filtered <- signature_counts_data_pooled[,-indices]

sig_split <- do.call(rbind, lapply(colnames(signature_counts_data_pooled_filtered), function(x) strsplit(as.character(x), split="")[[1]]))
# which(temp[,3]==temp[,6])

#save(signature_counts_data_pooled_filtered,
#            file="../summary_data/signature-counts-Lindo2016.rda")

