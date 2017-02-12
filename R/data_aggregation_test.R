
dir <- "../data/Sardinia_data_strand_flank/"
files = list.files(dir)
out <- damage_build_bin_counts(paste0(dir, files[1]), type=2)




file <- read.csv(paste0(dir, files[1]), header=FALSE)
breaks=NULL
min_dist_from_end <- apply(file[,2:3], 1, function(x) return(min(x)))


if(is.null(breaks)){
  bins <- c(-1, 5, 10, 15, 20, max(min_dist_from_end))
}else{
  message("breaks values provided : adjusting")
  bins <- c(intersect(breaks, (-1):max(min_dist_from_end)),max(min_dist_from_end))
}

bin_values <- .bincode(min_dist_from_end, bins, TRUE)

modified_file <- cbind.data.frame(file[,1], file[,4], file[,5], file[,6], file[,7], bin_values)
colnames(modified_file) <- c("pattern", "base.lsb", "base.rsb", "strand", "counts", "bin_values")

library(plyr)
df1 <- plyr::ddply(modified_file, .(pattern, base.lsb, base.rsb, strand, bin_values), summarise, newvar = sum(counts))
colnames(df1) <- c("pattern", "base.lsb", "base.rsb", "strand", "bin", "counts")

df2 <- cbind.data.frame(paste0(df1[,1], "_", df1[,4], "_", df1[,2], "_", df1[,3], "_", df1[,5]), df1[,6])
colnames(df2) <- c("pattern-strand-breaks", "counts")



ancient_files <- list.files(dir)
signature_ancient <- vector(mode="list")
signature_counts_ancient <- vector(mode="list")


for(num in 1:length(ancient_files)){
  tmp_dat <- damage_build_bin_counts(paste0(dir, ancient_files[num]),
                                     breaks=breaks,
                                     type=2)
  signature_counts_ancient[[num]] <- tmp_dat[,2];
  signature_ancient[[num]] <- as.character(tmp_dat[,1]);
  cat("Reading file ", num, "\n")
}

flanking_bases = 1
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

signature_split <- do.call(rbind, lapply(merged_signature_ancient, function(x) strsplit(as.character(x), split="")[[1]][1:(4+2*flanking_bases)]))

indices1 <- which(signature_split[,(flanking_bases+1)]==signature_split[,(flanking_bases+4)])

indices2 <- numeric()
for(m in 1:(4+2*flanking_bases)){
  indices2 <- c(indices2, which(signature_split[,m]=="N"));
}

indices <- union(indices1, indices2)

if(length(indices) > 0) {
  ancient_counts_filtered <- ancient_counts[,-indices]
}else{
  ancient_counts_filtered <- ancient_counts
}

rownames(ancient_counts_filtered) <- ancient_files_filt
if(length(indices) > 0){
  colnames(ancient_counts_filtered) <- merged_signature_ancient[-indices]
}else{
  colnames(ancient_counts_filtered) <- merged_signature_ancient
}


##########################################################################

##########################################################################


data <- get(load("../processed_data/sardinia2017-strand-flank.rda"))
signature_counts <- data

flanking_bases <- 1
signature_set <- colnames(signature_counts)

signature_set_split <- do.call(rbind, lapply(signature_set, function(x) strsplit(as.character(x), split="")[[1]][1:(4+2*flanking_bases)]))

indices_G <-  which(signature_set_split[,(flanking_bases+1)]=="G");
indices_A <-  which(signature_set_split[,(flanking_bases+1)]=="A");

indices <- c(indices_G, indices_A);

signature_set_2 <- signature_set
signature_set_2[indices] <- signatureclub(signature_set[indices])

signature_counts_pooled <- do.call(rbind, lapply(1:dim(signature_counts)[1], function(x) tapply(signature_counts[x,], signature_set_2, sum)))
rownames(signature_counts_pooled) <- rownames(signature_counts)
temp_split <- do.call(rbind, lapply(colnames(signature_counts_pooled), function(x) strsplit(as.character(x), split="")[[1]][1:(4+2*flanking_bases)]))


signatureclub <- function(signature_set){
  from <- c('A','T','G','C')
  to <- c('t','a','c','g');
  signature_set_mod <- array(0, length(signature_set));
  for(m in 1:length(signature_set)){
    temp <- toupper(gsub2(from, to, signature_set[m]))
    temp_split <- strsplit(as.character(temp), split="")[[1]]
    sign <- temp_split[6+2*flanking_bases]
    if(sign=="-"){sign = "+"}else if(sign=="+"){sign = "-"}
    leftb <- temp_split[10+2*flanking_bases]
    rightb <- temp_split[8+2*flanking_bases]
    temp_split[8+2*flanking_bases] <- leftb
    temp_split[10+2*flanking_bases] <- rightb
    temp_split[6+2*flanking_bases] <- sign
    temp_new <- paste0(temp_split, collapse = "")
    signature_set_mod[m] <- temp_new
  }
  return(signature_set_mod)
}


