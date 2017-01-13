

############  Between moderns comparison through aRchaic  #########################


hgdp_data <- get(load("../processed_data/HGDPmoderns-counts-table.rda"))
hgdp_data <- club_signature_counts(hgdp_data)

thousandg_data <- get(load("../processed_data/1000Gmoderns-clubbed_counts-table.rda"))

pooled_names <- intersect(colnames(hgdp_data), colnames(thousandg_data))

filtered_hgdp <- hgdp_data[, match(pooled_names, colnames(hgdp_data))]
filtered_thousandg <- thousandg_data[, match(pooled_names, colnames(thousandg_data))]

pooled_data <- rbind(filtered_hgdp, filtered_thousandg)
pooled_data <- filter_signatures_wo_location(pooled_data)

labs <- c(rep("HGDP", dim(filtered_hgdp)[1]),
          as.character(sapply(rownames(thousandg_data), function(x) return(strsplit(x, "[.]")[[1]][5]))))

topic_clus <- maptpx::topics(pooled_data, K=5, tol=10)

omega <- topic_clus$omega

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


damageLogo(topic_clus$theta)
