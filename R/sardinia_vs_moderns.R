

##########  analysing Sardinia data : with moderns + sherpa + Gossling ##########################

library(aRchaic)

sardinia_data <- get(load("../processed_data/sardinia2017.rda"))
sherpa_data <- get(load("../processed_data/sherpa2017.rda"))
thousandg_data <- get(load("../processed_data/1000Gmoderns-counts-table.rda"));
hgdp_data <- get(load("../processed_data/HGDPmoderns-counts-table.rda"))
gossling_data <- get(load("../processed_data/annagosling2016-counts-table.rda"))

sardinia_data_clubbed <- club_signature_counts(sardinia_data)
sherpa_data_clubbed <- club_signature_counts(sherpa_data)
thousandg_data_clubbed <- club_signature_counts(thousandg_data)
hgdp_data_clubbed <- club_signature_counts(hgdp_data)
gossling_data_clubbed <- club_signature_counts(gossling_data)

matched_names <- Reduce(intersect, list(colnames(sardinia_data_clubbed), colnames(sherpa_data_clubbed),
                       colnames(thousandg_data_clubbed), colnames(hgdp_data_clubbed),
                       colnames(gossling_data_clubbed)))

sardinia_data_clubbed_2 <- sardinia_data_clubbed[, match(matched_names, colnames(sardinia_data_clubbed))]
sherpa_data_clubbed_2 <- sherpa_data_clubbed[, match(matched_names, colnames(sherpa_data_clubbed))]
thousandg_data_clubbed_2 <- thousandg_data_clubbed[, match(matched_names, colnames(thousandg_data_clubbed))]
hgdp_data_clubbed_2 <- hgdp_data_clubbed[, match(matched_names, colnames(hgdp_data_clubbed))]
gossling_data_clubbed_2 <- gossling_data_clubbed[, match(matched_names, colnames(gossling_data_clubbed))]

pooled_data <- rbind(sardinia_data_clubbed_2, sherpa_data_clubbed_2, thousandg_data_clubbed_2,
                     hgdp_data_clubbed_2, gossling_data_clubbed_2)

labs <- c(rep("sardinia", dim(sardinia_data_clubbed_2)[1]),
          rep("sherpa", dim(sherpa_data_clubbed_2)[1]),
          rep("1000-G", dim(thousandg_data_clubbed_2)[1]),
          rep("HGDP", dim(hgdp_data_clubbed_2)[1]))

names <- rownames(gossling_data_clubbed_2);
control_indices <- c(grep("EXN", names), grep("Libneg", names), grep("PCRneg", names))


labs1 <- character();
labs1 <- rep("Gossling-ancient", dim(gossling_data_filtered)[1])
labs1[control_indices] <- "Gossling-controls"

labs <- c(labs, labs1)

topic_clus <- maptpx::topics(pooled_data, K=3, tol=100)
save(topic_clus, file="../rda/maptpx-3-sherpa-sardinia-hgdp-1000g-gossling.rda")

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


topic_clus <- maptpx::topics(pooled_data, K=4, tol=500)
save(topic_clus, file="../rda/maptpx-4-sherpa-sardinia-hgdp-1000g-gossling.rda")

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
