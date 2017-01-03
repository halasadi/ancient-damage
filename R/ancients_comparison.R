

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

gosling_filtered_counts_2 <- gosling_filtered_counts[,match(matched_names, colnames(gosling_filtered_counts))]
sherpa_filtered_counts_2 <- sherpa_filtered_counts[,match(matched_names, colnames(sherpa_filtered_counts))]
sardinia_filtered_counts_2 <- sardinia_filtered_counts[,match(matched_names, colnames(sardinia_filtered_counts))]

pooled_counts <- rbind(gosling_filtered_counts_2, sherpa_filtered_counts_2, sardinia_filtered_counts_2)

labs <- c(rep("Gosling-ancient", dim(gosling_filtered_counts_2)[1]),
          rep("sherpa", dim(sherpa_filtered_counts_2)[1]),
          rep("sardinia", dim(sardinia_filtered_counts_2)[1]))

gridPCA_signatures(pooled_counts, factor(labs))

library(CountClust)
topic_clus <- maptpx::topics(pooled_counts, K=3, tol=100)
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



