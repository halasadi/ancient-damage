
###  Signature Club :  ####################

signature_counts <- get(load("../summary_data/signature-counts-Lindo2016.rda"))

signature_counts <- get(load("../summary_data/153-ancients-counts-table.rda"))

signature_set <- colnames(signature_counts)

gsub2 <- function(pattern, replacement, x, ...) {
  for(i in 1:length(pattern))
    x <- chartr(pattern[i], replacement[i], x, ...)
  x
}

signatureclub <- function(signature_set){
  from <- c('A','T','G','C')
  to <- c('t','a','c','g');
  signature_set_mod <- array(0, length(signature_set));
  for(m in 1:length(signature_set)){
    signature_set_mod[m] <- toupper(gsub2(from, to, signature_set[m]))
  }
  return(signature_set_mod)
}

new_signature_set <- signatureclub(signature_set)

signature_set_split <- do.call(rbind, lapply(signature_set, function(x) strsplit(as.character(x), split="")[[1]]))

indices1 <-  which(signature_set_split[,3]=="C" & signature_set_split[,6]=="T")
indices2 <-  which(signature_set_split[,3]=="G" & signature_set_split[,6]=="A")

C_to_T_counts <- rowSums(signature_counts[,indices1]);
G_to_A_counts <- rowSums(signature_counts[,indices2]);

par(mar=c(5,4,4,4))
plot(C_to_T_counts, G_to_A_counts, xlab="C_to_T", ylab="G_to_A",
     pch=20, cex=1, col="red")
abline(0,1)


indices_G <-  which(signature_set_split[,3]=="G");
indices_A <-  which(signature_set_split[,3]=="A");

indices <- c(indices_G, indices_A);

signature_set_2 <- signature_set
signature_set_2[indices] <- signatureclub(signature_set[indices])

signature_counts_pooled <- do.call(rbind, lapply(1:dim(signature_counts)[1], function(x) tapply(signature_counts[x,], signature_set_2, sum)))
rownames(signature_counts_pooled) <- rownames(signature_counts)
temp_split <- do.call(rbind, lapply(colnames(signature_counts_pooled), function(x) strsplit(as.character(x), split="")[[1]]))

which(temp_split[,3]=="G") ## should be NA
which(temp_split[,3]=="A") ## should be NA

save(signature_counts_pooled, file="../summary_data/153-ancient-signature-counts-clubbed.rda")
