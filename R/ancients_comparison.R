

##########  Ancients comparison (Sardinia + Gossling Ancients + Sherpa) ###################


## Here we compare only the ancient sampled from three sources - Sardinia, Gossling Ancients and Sherpa
## We leave out the Lindo data because of it being clipped at the start
## We leave out the controls and the modern samples to remove any bias.

sardinia_counts <- get(load("../processed_data/sardinia2017.rda"))
temp <- club_signature_counts(sardinia_counts)
sardinia_filtered_counts <- filter_signatures_by_location(temp, max_pos=20, flanking_bases = 2)

gosling_counts <- get(load("../processed_data/annagosling2016-counts-table.rda"))
temp <- club_signature_counts(gosling_counts)
temp <- temp[-28,];
gosling_filtered_counts <- filter_signatures_by_location(temp, max_pos=20, flanking_bases = 2)

sherpa_counts <- get(load("../processed_data/sherpa2017.rda"))
temp <- club_signature_counts(sherpa_counts)
sherpa_filtered_counts <- filter_signatures_by_location(temp, max_pos=20, flanking_bases = 2)

names <- rownames(gosling_filtered_counts);
control_indices <- c(grep("EXN", names), grep("Libneg", names), grep("PCRneg", names))

gosling_controls_filtered_counts <- gosling_filtered_counts[control_indices,]
gosling_ancient_filtered_counts <- gosling_filtered_counts[-control_indices,]

matched_names <- intersect(colnames(gosling_ancient_filtered_counts),
                           intersect(colnames(sardinia_filtered_counts), colnames(sherpa_filtered_counts)))

gosling_filtered_counts_2 <- gosling_ancient_filtered_counts[,match(matched_names, colnames(gosling_ancient_filtered_counts))]
sherpa_filtered_counts_2 <- sherpa_filtered_counts[,match(matched_names, colnames(sherpa_filtered_counts))]
sardinia_filtered_counts_2 <- sardinia_filtered_counts[,match(matched_names, colnames(sardinia_filtered_counts))]

pooled_counts <- rbind(gosling_filtered_counts_2, sherpa_filtered_counts_2, sardinia_filtered_counts_2)

labs <- c(rep("Gosling-ancient", dim(gosling_filtered_counts_2)[1]),
          rep("sherpa", dim(sherpa_filtered_counts_2)[1]),
          rep("sardinia", dim(sardinia_filtered_counts_2)[1]))

gridPCA_signatures(pooled_counts, factor(labs))


signature_set <- colnames(pooled_counts)
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
topic_clus <- maptpx::topics(pooled_counts, K=3, tol=100, model="independent", signatures = signature_pos)
save(topic_clus, file="../processed_data/maptpx-runs/topic-clus-pooled-ancients-gosling-sherpa-sardinia-3.rda")


omega <- topic_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs)
)

CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = FALSE,
                            figure_title = paste0("StructurePlot: K=", dim(omega)[2],""),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))


damageLogo_pos(topic_clus$theta, max_pos=20)



