

##########   Comparing the moderns data  ############################

library(aRchaic)

# dir <- "../data/1000Gmoderns/";
# out <- aggregate_bin_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))
# save(out,
#      file="../processed_data/1000Gmoderns-counts-table.rda")
#
#
# dir <- "../data/HGDPmoderns/";
# out <- aggregate_bin_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))
# save(out,
#      file="../processed_data/HGDPmoderns-counts-table.rda")
#
#
# dir <- "../data/Lindo2016moderns/";
# out <- aggregate_bin_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))
# save(out,
#      file="../processed_data/lindo2016moderns-counts-table.rda")


ThousandGmoderns <- get(load("../processed_data/1000Gmoderns-counts-table.rda"))
HGDPmoderns <- get(load("../processed_data/HGDPmoderns-counts-table.rda"))
Lindomoderns <- get(load("../processed_data/lindo2016moderns-counts-table.rda"))

clubbed_counts_thousandG <- club_signature_counts(ThousandGmoderns)
save(clubbed_counts_thousandG, file="../processed_data/1000Gmoderns-clubbed_counts-table.rda")
filtered_counts_thousandG <- filter_signatures_by_location(clubbed_counts_thousandG, max_pos=20, flanking_bases = 2)

clubbed_counts_HGDP <- club_signature_counts(HGDPmoderns)
filtered_counts_HGDP <- filter_signatures_by_location(clubbed_counts_HGDP, max_pos=20, flanking_bases = 2)

clubbed_counts_Lindomoderns <- club_signature_counts(Lindomoderns)
filtered_counts_Lindomoderns <- filter_signatures_by_location(clubbed_counts_Lindomoderns, max_pos=20, flanking_bases = 2)

matched_names <- intersect(colnames(filtered_counts_thousandG), intersect(colnames(filtered_counts_HGDP), colnames(filtered_counts_Lindomoderns)))

filtered_counts_thousandG_2 <- filtered_counts_thousandG[,match(matched_names, colnames(filtered_counts_thousandG))]
filtered_counts_HGDP_2 <- filtered_counts_HGDP[,match(matched_names, colnames(filtered_counts_HGDP))]
filtered_counts_Lindomoderns_2 <- filtered_counts_Lindomoderns[,match(matched_names, colnames(filtered_counts_Lindomoderns))]

labs <- c(rep("1000G", dim(filtered_counts_thousandG_2)[1]), rep("HGDP", dim(filtered_counts_HGDP_2)[1]), rep("Lindo-moderns", dim(filtered_counts_Lindomoderns_2)[1]))

pooled_counts <- rbind(filtered_counts_thousandG_2, filtered_counts_HGDP_2, filtered_counts_Lindomoderns_2)

gridPCA_signatures(pooled_counts, factor(labs))

library(CountClust)
topic_clus <- maptpx::topics(pooled_counts, K=2, tol=10)
save(topic_clus, file="../processed_data/maptpx-runs/topic-clus-pooled-moderns-2.rda")

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


damageLogo_pos(topic_clus$theta, max_pos=20)



