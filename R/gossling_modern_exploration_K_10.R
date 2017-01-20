
library(aRchaic)
gossling_data <- get(load("../processed_data/annagosling2016-counts-table.rda"))
system.time(gossling_data_clubbed <- club_signature_counts(gossling_data))
gossling_data_clubbed <- gossling_data_clubbed[-28,];

modern_data <- get(load("../processed_data/1000Gmoderns-counts-table.rda"));
modern_data_clubbed <- club_signature_counts(modern_data)

pooled_names <- intersect(colnames(gossling_data_clubbed), colnames(modern_data_clubbed))
filtered_gossling <- gossling_data_clubbed[, match(pooled_names, colnames(gossling_data_clubbed))]
filtered_moderns <- modern_data_clubbed[, match(pooled_names, colnames(modern_data_clubbed))]

filtered_gossling_2 <- filter_signatures_by_location(filtered_gossling, max_pos=20, flanking_bases = 2)
filtered_moderns_2 <- filter_signatures_by_location(filtered_moderns, max_pos=20, flanking_bases = 2)

names <- rownames(filtered_gossling_2);
control_indices <- c(grep("EXN", names), grep("Libneg", names), grep("PCRneg", names))

labs <- character();
labs <- rep("ancient", dim(filtered_gossling_2)[1])
labs[control_indices] <- "controls"
labs <- c(labs, rep("1000g", dim(filtered_moderns_2)[1]))

indices <- which(labs == "ancient")

gossling_ancients <- filtered_gossling_2[indices, ]

# ancient_names = names[indices]
# pop_names_1 <- as.factor(substring(ancient_names, 8, 8))
# levels(pop_names_1)
#
# levels(pop_names_1) = c("Chokhopani", "Kyang", "Rhirhi", "Mebrak", "Samdzong")
# labs[indices] <- as.character(pop_names_1)

pooled_data <- rbind(filtered_gossling_2, filtered_moderns_2)

#pooled_data <- filtered_gossling_2

signature_set <- colnames(pooled_data)
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


out <- topics(pooled_data, K=2, tol=100, model="independent", signatures = signature_pos)
save(out, file="../processed_data/maptpx-runs/gosling2016-1000g-maptpx-independent-K-2.rda")

out <- get(load("../processed_data/maptpx-runs/gosling2016-1000g-maptpx-independent-K-2.rda"))

# pop_names_1 <- as.factor(substring(, 1, 1))
# levels(pop_names_1)
#
# levels(pop_names_1) = c("Chokhopani", "Kyang", "Rhirhi", "Mebrak", "Samdzong")
#
#
# names <- rownames(filtered_gossling_2);
# control_indices <- c(grep("EXN", names), grep("Libneg", names), grep("PCRneg", names))

# labs <- character();
# labs <- rep("ancient", dim(filtered_gossling_2)[1])
# labs[control_indices] <- "controls"

#labs <- c(labs, rep("1000g", dim(filtered_moderns_2)[1]))

omega <- out$omega

cols1 <- c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
           "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
           "brown4","darkorchid","magenta","yellow", "azure1","azure4")

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs)
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

metadata <- get(load("../data/Gosling_metadata/metadata_ordered.rda"))
metadata$Sample[indices]
dat <- data.frame("omega" = out$omega[indices,1], "damage"= metadata$MapDamage[indices], "pop"= pop_names_1)
ggplot2::qplot(omega, damage, col=pop, data=dat, xlab="GoM cluster membership k=3 (clus1)",
               ylab="mapdamage C to T height ")

