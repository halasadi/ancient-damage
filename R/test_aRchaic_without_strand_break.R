

#########  test aRchaic without strand break  #######################

mat1 <- get(load("../data/Fu_2016/Fu_2016.rda"))
mat2 <- get(load("../data/moderns_lite/moderns_lite.rda"))

mat_pooled <- rbind(mat1, mat2)

mat_reduced <- filter_out_strand_break(mat_pooled)

signature_set <- colnames(mat_reduced)
sig_split <- t(sapply(1:length(signature_set), function(x) return(strsplit(signature_set[x], "")[[1]][1:(flanking_bases+5)])))
new_sig_split <- matrix(0, dim(sig_split)[1], 3);
new_sig_split[,1] <- sig_split[,1]
new_sig_split[,2] <- sapply(1:length(signature_set), function(x) return(paste(sig_split[x,(flanking_bases+1):(flanking_bases+4)], collapse="")))
new_sig_split[,3] <- sig_split[,(flanking_bases+5)]

pos <- sapply(signature_set, function(x) return(strsplit(x, "_")[[1]][2]))

pos <- as.numeric(pos)
pos <- pos - min(pos)
pos <- factor(pos, levels = 0:21)

