

########  damage logo hist trial  #################

library(aRchaic)
par(mfrow=c(1,1))
signature_counts <- get(load("../processed_data/annagosling2016-counts-table.rda"))
validation_check <- club_signature_validation_plot(signature_counts, log=TRUE)
clubbed_counts <- club_signature_counts(signature_counts)
clubbed_counts <- clubbed_counts[-28,];

filtered_counts <- filter_signatures_by_location(clubbed_counts, max_pos = 15, flanking_bases = 2)

topic_clus <- maptpx::topics(filtered_counts, K=2, tol=10)


names <- rownames(clubbed_counts);
control_indices <- c(grep("EXN", names), grep("Libneg", names), grep("PCRneg", names))

labs <- character();
labs <- rep("ancient", dim(clubbed_counts)[1])
labs[control_indices] <- "controls"

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
                            figure_title = paste0("StructurePlot: K=", dim(omega)[2],": pmsignature: withput position information"),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))

damageLogo_pos(topic_clus$theta)

