

###########  Age and aRchaic  ############################

### In this script, we try to see how the age information of a sample (mainly
### aDNA sample) is reflected in aRchaic

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

gossling_ancients <- gossling_data_clubbed[indices, ]

metadata <- get(load("../data/Gosling_metadata/metadata_ordered.rda"))

signature_set <- colnames(gossling_ancients)
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
pos <- factor(pos, levels = 0:22)

signatures <- mat;
signature_pos <- cbind.data.frame(signatures, pos)


topic_clus <- topics(gossling_ancients, K=2, tol=10, type="independent", signatures = signature_pos)
save(topic_clus, file="../processed_data/maptpx-runs/topic-clus-pooled-ancients-gosling-2-independent.rda")


damageLogo_pos(topic_clus$theta)


topic_clus <- get(load("../processed_data/maptpx-runs/topic-clus-pooled-ancients-gosling-2-independent.rda"))

omega <- topic_clus$omega
labs <- rep("goss-ancient", dim(omega)[1])

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = 1:NROW(omega)
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


metadata$Sample[indices]
metadata$MapDamage[indices]
metadata$ContamMix.MAP.authentic[indices]

dim(omega)
plot(omega[,1], metadata$MapDamage[indices])

pop_names <- metadata$Sample[indices]
pop_names_1 <- as.factor(substring(pop_names, 1, 1))
levels(pop_names_1)

levels(pop_names_1) = c("Chokhopani", "Kyang", "Rhirhi", "Mebrak", "Samdzong")

dat <- data.frame("omega" = omega[,1], "damage"= metadata$Sex.determination[indices], "pop"= pop_names_1)
ggplot2::qplot(omega, damage, col=pop, data=dat, xlab="GoM cluster membership k=2 (clus1)",
               ylab="mapdamage C to T height ")

tmp <- as.character(metadata$endogenous.[indices])
endogenous <- as.numeric(as.character(sapply(tmp, function(x) return(substr(x, 1, nchar(x)-1)))))

dat <- data.frame("omega" = omega[,1], "endogenous"= endogenous, "pop"= pop_names_1)
ggplot2::qplot(omega, endogenous, col=pop, data=dat, xlab="GoM cluster membership k=2 (clus1)",
               ylab="percentage of endogenous DNA")
