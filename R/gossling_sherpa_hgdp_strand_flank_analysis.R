

#########  Gossling + Sherpa + HGDP  #########################

gossling_data <- get(load("../processed_data/annagosling2016-strand-flank.rda"))
gossling_club <- club_signature_counts(gossling_data, flanking_bases = 2)
gossling_club <- gossling_club[-28,];
gossling_filtered <- filter_signatures_by_location(gossling_club, max_pos=20)
gossling_filtered_2 <- cbind(gossling_filtered, gossling_club[,(dim(gossling_club)[2]-7):dim(gossling_club)[2]])


sherpa_data <- get(load("../processed_data/sherpa2017-strand-flank.rda"))
sherpa_club <- club_signature_counts(sherpa_data, flanking_bases = 2)
sherpa_filtered <- filter_signatures_by_location(sherpa_club, max_pos = 20)
sherpa_filtered_2 <- cbind(sherpa_filtered, sherpa_club[,(dim(sherpa_club)[2]-7):dim(sherpa_club)[2]])


hgdp_data <- get(load("../processed_data/HGDPmoderns-counts-table-strand-flank.rda"))
hgdp_club <- club_signature_counts(hgdp_data, flanking_bases = 2)
hgdp_filtered <- filter_signatures_by_location(hgdp_club, max_pos = 20)
hgdp_filtered_2 <- cbind(hgdp_filtered, hgdp_club[,(dim(hgdp_club)[2]-7):dim(hgdp_club)[2]])


pooled_names <- intersect(colnames(gossling_filtered_2),
                          intersect(colnames(sherpa_filtered_2),
                          colnames(hgdp_filtered_2)))

indices1 <- match(pooled_names, colnames(gossling_filtered_2))
indices2 <- match(pooled_names, colnames(sherpa_filtered_2))
indices3 <- match(pooled_names, colnames(hgdp_filtered_2))

pooled_dat <- rbind(gossling_filtered_2[, indices1],
                    sherpa_filtered_2[, indices2],
                    hgdp_filtered_2[, indices3])

names <- rownames(gossling_filtered);
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

labs <- c(labs, rep("sherpa", dim(sherpa_filtered)[1]))
labs <- c(labs, rep("hgdp", dim(hgdp_filtered)[1]))

topic_clus <- topics(pooled_dat, K=10, type = "full", tol=100)
save(topic_clus, file="../processed_data/maptpx-runs/gossling_sherpa_hgdp_maptpx_independent_K10.rda")

omega <- topic_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs, levels=c("Chokhopani", "Kyang", "Mebrak", "Rhirhi",
                                       "Samdzong", "controls", "sherpa", "hgdp"))
)



rownames(omega) <- annotation$sample_id;

cols1 <- c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
           "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
           "brown4","darkorchid","magenta","yellow", "azure1","azure4")

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

damageLogo_pos_str(topic_clus$theta)

