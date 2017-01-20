

#####  clubbing the signature counts (Gossling et al data) #################

gossling_data <- get(load("../processed_data/annagosling2016-strand-flank.rda"))
gossling_club <- club_signature_counts(gossling_data, flanking_bases = 2)
gossling_club <- gossling_club[-28,];

sherpa_data <- get(load("../processed_data/sherpa2017-strand-flank.rda"))
sherpa_club <- club_signature_counts(sherpa_data, flanking_bases = 2)

pooled_names <- intersect(colnames(gossling_club),
                          colnames(sherpa_club))

indices1 <- match(pooled_names, colnames(gossling_club))
indices2 <- match(pooled_names, colnames(sherpa_club))

pooled_dat <- rbind(gossling_club[, indices1],
                    sherpa_club[, indices2])


names <- rownames(gossling_club);
control_indices <- c(grep("EXN", names), grep("Libneg", names), grep("PCRneg", names))

labs <- character();
labs <- rep("ancient", dim(gossling_club)[1])
labs[control_indices] <- "controls"

indices <- which(labs == "ancient")

ancient_names = names[indices]
pop_names_1 <- as.factor(substring(ancient_names, 8, 8))
levels(pop_names_1)

levels(pop_names_1) = c("Chokhopani", "Kyang", "Rhirhi", "Mebrak", "Samdzong")
labs[indices] <- as.character(pop_names_1)

labs <- c(labs, rep("sherpa", dim(sherpa_club)[1]))

library(maptpx)
topic_clus <- topics(pooled_dat, K=5, type = "full", tol=100)

save(topic_clus, file="../processed_data/maptpx-runs/topic_clus_maptpx_K_5_gossling_sherpa_strand_flank.rda")

omega <- topic_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs, levels=c("Chokhopani", "Kyang", "Mebrak", "Rhirhi",
                                       "Samdzong", "controls", "sherpa"))
)

rownames(omega) <- annotation$sample_id;
plot.new()
CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Development Phase",
                            order_sample = FALSE,
                            figure_title = paste0("StructurePlot: K=", dim(omega)[2],", no C -> T / G -> A"),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))

damageLogo_pos_str(topic_clus$theta)
