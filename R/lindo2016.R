

####################   Lindo 2016 data analysis  #####################################

#################################   Ancients data  ################################


dir <- "../data/Lindo2016ancients/";
out <- aggregate_bin_counts(dir, breaks = c(-1,5,10,15))
save(out,
     file="../processed_data/lindo2016ancients-counts-table.rda")

signature_counts <- get(load("../processed_data/lindo2016ancients-counts-table.rda"))
validation_check <- club_signature_validation_plot(signature_counts, log=TRUE)
clubbed_counts <- club_signature_counts(signature_counts)
filtered_counts_ancient <- filter_signatures_wo_location(clubbed_counts);


################################   Moderns data  ######################################

dir <- "../data/Lindo2016moderns/";
out <- aggregate_bin_counts(dir, breaks = c(-1,5,10,15))
save(out,
     file="../processed_data/lindo2016moderns-counts-table.rda")

signature_counts <- get(load("../processed_data/lindo2016moderns-counts-table.rda"))
validation_check <- club_signature_validation_plot(signature_counts, log=TRUE)
clubbed_counts <- club_signature_counts(signature_counts)
filtered_counts_modern <- filter_signatures_wo_location(clubbed_counts);

filtered_counts_pooled <- rbind(filtered_counts_ancient, filtered_counts_modern)


signature_set <- colnames(filtered_counts_pooled)
sig_split <- t(sapply(1:length(signature_set), function(x) return(strsplit(signature_set[x], "")[[1]])))
new_sig_split <- matrix(0, dim(sig_split)[1], 5);
new_sig_split[,1] <- sig_split[,1]
new_sig_split[,2] <- sig_split[,2]
new_sig_split[,3] <- sapply(1:length(signature_set), function(x) return(paste(sig_split[x,3:6], collapse="")))
new_sig_split[,4] <- sig_split[,7]
new_sig_split[,5] <- sig_split[,8]

levels(new_sig_split[,1]) <- c("0", "1", "2", "3", "4")

mat <- matrix(0, dim(new_sig_split)[1], dim(new_sig_split)[2])
for(k in 1:dim(new_sig_split)[2]){
  temp <- as.factor(new_sig_split[,k])
  mat[,k] <- as.numeric(as.matrix(plyr::mapvalues(temp, from = levels(temp), to = 0:(length(levels(temp))-1))))
}

signatures <- mat;

topics_clus_1 <- maptpx::topics(filtered_counts_pooled, K=3, type="full", tol=10)
topics_clus_2 <- maptpx::topics(filtered_counts_pooled, K=3, type="independent", tol=10, signatures = signatures)

labs <- c(rep("Ancient",25), rep("Modern",25))

omega <- topics_clus_1$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs)
)

rownames(omega) <- annotation$sample_id;

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

signature_set <- colnames(filtered_counts_pooled)
apply(CountClust::ExtractTopFeatures(topics_clus_1$theta, top_features = 10, method="poisson", options="min"), c(1,2), function(x) signature_set[x])
