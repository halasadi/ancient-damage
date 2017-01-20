

###  clubbing the signature counts ################

lindo_ancients <- get(load("../processed_data/lindo2016ancients-counts-table-strand-flank.rda"))
lindo_moderns <- get(load("../processed_data/lindo2016moderns-counts-table-strand-flank.rda"))

lindo_moderns_club <- club_signature_counts(lindo_moderns, flanking_bases = 2)
lindo_ancients_club <- club_signature_counts(lindo_ancients, flanking_bases = 2)

save(lindo_moderns_club,
     file="../processed_data/lindo2016moderns-counts-table-strand-flank-clubbed.rda")
save(lindo_ancients_club,
     file="../processed_data/lindo2016ancients-counts-table-strand-flank-clubbed.rda")



pooled_names <- intersect(colnames(lindo_ancients_club),
                          colnames(lindo_moderns_club))

indices1 <- match(pooled_names, colnames(lindo_ancients_club))
indices2 <- match(pooled_names, colnames(lindo_moderns_club))

pooled_dat <- rbind(lindo_ancients_club[, indices1],
                    lindo_moderns_club[, indices2])

filtered_data <- filter_signatures_by_location(pooled_dat, max_pos = 20)
filtered_data_2 <- cbind(filtered_data,  pooled_dat[,(dim(pooled_dat)[2]-7):dim(pooled_dat)[2]])



leftflank <- grep("left", colnames(filtered_data_2))
rightflank <- grep("right", colnames(filtered_data_2))

signature_set <- colnames(filtered_data)
sig_split <- t(sapply(1:length(signature_set), function(x) return(strsplit(signature_set[x], "")[[1]][1:10])))
new_sig_split <- matrix(0, dim(sig_split)[1], 6);
new_sig_split[,1] <- sig_split[,1]
new_sig_split[,2] <- sig_split[,2]
new_sig_split[,3] <- sapply(1:length(signature_set), function(x) return(paste(sig_split[x,3:6], collapse="")))
new_sig_split[,4] <- sig_split[,7]
new_sig_split[,5] <- sig_split[,8]
new_sig_split[,6] <- sig_split[,10]
#indices_notCtoA <-  which(new_sig_split[,3] != "C->T")
#pooled_data <- pooled_data[, indices_notCtoA]

levels(new_sig_split[,1]) <- c("0", "1", "2", "3", "4")

pos <- t(sapply(1:length(signature_set), function(x)
{
  y = strsplit(signature_set[x], "")[[1]]
  return(paste(y[12:length(y)], collapse=""))
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


library(maptpx)
topic_clus <- topics(filtered_data_2, K=3, type = "independent",
                     signatures = signature_pos,
                     ind_model_indices = setdiff(1:dim(filtered_data_2)[2],c(leftflank, rightflank)),
                     tol=100)

topic_clus_2 <- topics(filtered_data_2, type="full", K=10, tol=100)

tpxThetaGroupInd(topic_clus_2$theta, ind_model_indices, signature_pos)

ind_model_indices = setdiff(1:dim(filtered_data_2)[2],c(leftflank, rightflank))
signatures = signature_pos

save(topic_clus_2, file="../processed_data/maptpx-runs/lindo2016_maptpx_full_K10_strand_break.rda")

labs <- c(rep("Ancient",25), rep("Modern",25))

omega <- topic_clus_2$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs)
)

cols1 <- c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
           "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
           "brown4","darkorchid","magenta","yellow", "azure1","azure4")

rownames(omega) <- annotation$sample_id;
plot.new()
CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = cols1,
                            yaxis_label = "Development Phase",
                            order_sample = FALSE,
                            figure_title = paste0("StructurePlot: K=", dim(omega)[2],", no C -> T / G -> A"),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))

damageLogo_pos_str(topic_clus_2$theta)


fac <- array(0, dim(topic_clus_2$theta)[1])
minus_indices <- grep("_-_", rownames(topic_clus_2$theta))
fac[minus_indices] <- "-"
fac[-minus_indices] <- "+"
out= tapply(topic_clus_2$theta[,3], factor(fac), mean)
out/sum(out)

