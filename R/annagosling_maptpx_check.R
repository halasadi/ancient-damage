

######################  test independent model on ancient DNA data  ###############################

library(aRchaic)
signature_counts <- get(load("../processed_data/annagosling2016-counts-table.rda"))
validation_check <- club_signature_validation_plot(signature_counts, log=TRUE)
clubbed_counts <- club_signature_counts(signature_counts)
clubbed_counts <- clubbed_counts[-28,];
filtered_counts <- filter_signatures_wo_location(clubbed_counts);

signature_set <- colnames(filtered_counts)
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

library(slam)
topics_clus <- maptpx::topics(filtered_counts, K=3, type="full", tol=10)
topics_clus <- maptpx::topics(filtered_counts, K=3, type="independent", tol=10, signatures = signatures)

#save(topics_clus, file="../rda/topicmodel_full_K_2.rda")
#save(topics_clus, file="../rda/topicmodel_independent_K_3.rda")

topics_clus <-  get(load("../rda/topicmodel_independent_K_2.rda"))

#topics_clus <- get(load("../rda/topicmodel_full_K_2.rda"))

omega <- topics_clus$omega

names <- rownames(filtered_counts);
control_indices <- c(grep("EXN", names), grep("Libneg", names), grep("PCRneg", names))

labs <- character();
labs <- rep("ancient", dim(filtered_counts)[1])
labs[control_indices] <- "controls"

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs)
)

par(mfrow=c(1,1))
CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = FALSE,
                            figure_title = paste0("StructurePlot: K=", dim(omega)[2],": pmsignature: independence"),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))

plot.new()
damageLogo(topics_clus$theta)


signature_set <- colnames(filtered_counts)
apply(CountClust::ExtractTopFeatures(topics_clus$theta, top_features = 10,
                                     method="poisson", options="min"), c(1,2), function(x) signature_set[x])



