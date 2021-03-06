

##################   Sardinia + Sherpa + Gosling  ##############################

library(aRchaic)
gossling_data <- get(load("../processed_data/annagosling2016-counts-table.rda"))
system.time(gossling_data_clubbed <- club_signature_counts(gossling_data))
gossling_data_clubbed <- gossling_data_clubbed[-28,];

names <- rownames(gossling_data_clubbed);
control_indices <- c(grep("EXN", names), grep("Libneg", names), grep("PCRneg", names))

labs <- character();
labs <- rep("ancient", dim(gossling_data_clubbed)[1])
labs[control_indices] <- "controls"

indices <- which(labs == "ancient")

ancient_names = names[indices]
pop_names_1 <- as.factor(substring(ancient_names, 8, 8))
levels(pop_names_1)

levels(pop_names_1) = c("Chokhopani", "Kyang", "Rhirhi", "Mebrak", "Samdzong")
labs[indices] <- as.character(pop_names_1)

# names <- rownames(gossling_data_clubbed);
# control_indices <- c(grep("EXN", names), grep("Libneg", names), grep("PCRneg", names))
#
# labs <- character();
# labs <- rep("ancient", dim(gossling_data_clubbed)[1])
# labs[control_indices] <- "controls"
#
# indices <- which(labs == "ancient")

gossling_ancients <- gossling_data_clubbed[indices, ]
gossling_filtered_counts <- filter_signatures_by_location(gossling_ancients, max_pos=20, flanking_bases = 2)


sherpa_data <- get(load("../processed_data/sherpa2017.rda"))
sherpa_data <- club_signature_counts(sherpa_data)
sherpa_filtered_counts <- filter_signatures_by_location(sherpa_data, max_pos=20, flanking_bases = 2)

sardinia_data <- get(load("../processed_data/sardinia2017.rda"))
sardinia_data <- club_signature_counts(sardinia_data)
sardinia_filtered_counts <- filter_signatures_by_location(sardinia_data, max_pos=20, flanking_bases = 2)

hgdp_data <- get(load("../processed_data/HGDPmoderns-counts-table.rda"))
hgdp_data <- club_signature_counts(hgdp_data)
hgdp_filtered_counts <- filter_signatures_by_location(hgdp_data, max_pos=20, flanking_bases = 2)


pooled_names <- intersect(colnames(sherpa_filtered_counts),
                          intersect(colnames(sardinia_filtered_counts),
                                    intersect( colnames(hgdp_filtered_counts), colnames(gossling_filtered_counts))))


filtered_gossling <- gossling_filtered_counts[, match(pooled_names, colnames(gossling_filtered_counts))]
filtered_sherpa <- sherpa_filtered_counts[, match(pooled_names, colnames(sherpa_filtered_counts))]
filtered_sardinia <- sardinia_filtered_counts[, match(pooled_names, colnames(sardinia_filtered_counts))]
filtered_hgdp <- hgdp_filtered_counts[, match(pooled_names, colnames(hgdp_filtered_counts))]

pooled_data <- rbind(filtered_gossling, filtered_sherpa, filtered_sardinia, filtered_hgdp)



signature_set <- colnames(pooled_data)
sig_split <- t(sapply(1:length(signature_set), function(x) return(strsplit(signature_set[x], "")[[1]][1:8])))
new_sig_split <- matrix(0, dim(sig_split)[1], 5);
new_sig_split[,1] <- sig_split[,1]
new_sig_split[,2] <- sig_split[,2]
new_sig_split[,3] <- sapply(1:length(signature_set), function(x) return(paste(sig_split[x,3:6], collapse="")))
new_sig_split[,4] <- sig_split[,7]
new_sig_split[,5] <- sig_split[,8]

#indices_notCtoA <-  which(new_sig_split[,3] != "C->T")
#pooled_data <- pooled_data[, indices_notCtoA]

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

out <- maptpx::topics(pooled_data, K=4, tol=100, model="independent", signatures = signature_pos)
#out <- topics(pooled_data, K=3, tol=100, model="full")

save(out, file="../processed_data/maptpx-runs/sards-gosling-sherpa-hgdp-maptpx-independent-K-4.rda")

out <- get(load("../processed_data/maptpx-runs/sards-gosling-sherpa-hgdp-maptpx-independent-K-3.rda"))

labs1 <- c(labs[indices],
          rep("Sherpa", dim(filtered_sherpa)[1]),
          rep("Sardinia", dim(filtered_sardinia)[1]),
          rep("HGDP", dim(filtered_hgdp)[1]))

omega <- out$omega

cols1 <- c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
           "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
           "brown4","darkorchid","magenta","yellow", "azure1","azure4")

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs1, levels = c("Chokhopani", "Kyang", "Rhirhi", "Mebrak", "Samdzong",
                                          "Sherpa", "Sardinia", "HGDP"))
)

CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = cols1,
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = FALSE,
                            figure_title = paste0("StructurePlot: K=", dim(omega)[2],""),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))


damageLogo_pos(out$theta, renyi_alpha = 100)





