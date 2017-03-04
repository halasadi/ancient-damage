

#############  Analysis of the Fu 2016 data  ##############################


#dir <- "../data/Fu_2016/";
#out <- aggregate_bin_counts(dir, breaks = c(-1, seq(1,20,1)))
#save(out,
#     file="../processed_data/Fu2016-strand-flank-2.rda")

data <- get(load("../processed_data/Fu2016-strand-flank-2.rda"))
clubbed_counts <- club_signature_counts_2(data, flanking_bases = 1)


##########  The read length  distribution of samples  ######################

files <- list.files("../data/Fu_2016/")
dir <- "../data/Fu_2016/"
tab <- list()
for(l in files){
  tab[[l]] <- read.csv(paste0(dir, l), header=FALSE)
  cat("We are at sample", l, "\n")
}

ancient_names <- as.character(sapply(files, function(x) return(strsplit(x, "[_]")[[1]][1])))

par(mfrow=c(3,3))
for(l in 1:length(tab)){
  tab1 <- tab[[l]]
  read_length <- tab1$V2 + tab1$V3
  indices <- grep("C->T", tab1$V1)
  indices <- c(indices, grep("G->A", tab1$V1))
  read_length_CtoT <- tab1[indices, ]$V2 + tab1[indices, ]$V3 ## read lengths for all C to T
  indices2 <- which(tab1$V2 < 5 | tab1$V3 < 5)
  indices_matched <- intersect(indices, indices2)
  read_length_3 <- (tab1[indices_matched, ]$V2 + tab1[indices_matched, ]$V3)
  plot(table(read_length)/sum(table(read_length)), type="o", col="red",
       main=paste0(ancient_names[l]), xlab="read pos", ylab="prop of occur", cex.main = 0.5)
  points(table(read_length_CtoT)/sum(table(read_length_CtoT)), type="o", col="green")
  points(table(read_length_3)/ sum(table(read_length_3)), type="o", col="blue")
  legend("topleft", fill=c("red", "green", "blue"), legend = c("all", "CtoT", "CtoT < 5"))
  cat("we are at sample", l, "\n")
}


##############   Topic model  fit   ######################################


