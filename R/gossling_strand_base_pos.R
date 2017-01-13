

############  Anna Gosling data with strand breaks + position + flanking bases  ######################

left_strand_breaks <- as.numeric()
right_strand_breaks <- as.numeric()

strand_breaks_gosling <- get(load("../processed_data/strand-breaks-gosling2016.rda"))
for(l in 1:length(strand_breaks_gosling)){
  left_strand_breaks <- rbind(left_strand_breaks, strand_breaks_gosling[[l]][[1]])
  right_strand_breaks <- rbind(right_strand_breaks, strand_breaks_gosling[[l]][[2]])
}

left_strand_breaks <- left_strand_breaks[-28,]
colnames(left_strand_breaks) <- paste0(colnames(left_strand_breaks), "_left")
right_strand_breaks <- right_strand_breaks[-28,]
colnames(right_strand_breaks) <- paste0(colnames(right_strand_breaks), "_right")


gosling_counts <- get(load("../processed_data/annagosling2016-counts-table.rda"))
temp <- club_signature_counts(gosling_counts)
temp <- temp[-28,];
gosling_filtered_counts <- filter_signatures_by_location(temp, max_pos=20, flanking_bases = 2)

gosling_pooled_data <- cbind.data.frame(gosling_filtered_counts, left_strand_breaks, right_strand_breaks)

names <- rownames(gosling_pooled_data);
control_indices <- c(grep("EXN", names), grep("Libneg", names), grep("PCRneg", names))
control_indices <- c(grep("EXN", names), grep("Libneg", names), grep("PCRneg", names))

labs <- character();
labs <- rep("ancient", dim(gosling_pooled_data)[1])
labs[control_indices] <- "controls"

signature_set <- colnames(gosling_filtered_counts)
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
pos <- factor(pos, levels = 0:20)

signatures <- mat;
signature_pos <- cbind.data.frame(signatures, pos)


library(CountClust)
topic_clus <- maptpx::topics(gosling_pooled_data, K=2, tol=100, model="full")
save(topic_clus, file="../processed_data/maptpx-runs/topic-clus-gosling-strand-2.rda")

library(gridBase)
damageLogo_pos_str(topic_clus$theta)


theta3 = rbind(c(0.25, 0.5), c(0.75, 0.5))
rownames(theta3) <- c("plus", "minus")
colnames(theta3) <- c("1", "2")

