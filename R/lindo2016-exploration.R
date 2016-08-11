
##  Exploration of the Lindo2016 data 

signature_counts <- get(load("../summary_data/signature-counts-clubbed-Lindo2016.rda"))

signature_set <- colnames(signature_counts)
###  number of mutations 

barplot(t(rowSums(signature_counts)), horiz=TRUE)



pr <- prcomp(t(limma::voom(t(signature_counts))$E))

pc_data_frame <- data.frame("PC"=pr$x,
                            "labels"=c(rep("Ancient",25),
                                       rep("Modern",25)))

qplot(PC.PC1, PC.PC2,
      data=pc_data_frame,
      colour=labels)

sig_split <- do.call(rbind, lapply(colnames(signature_counts), function(x) strsplit(as.character(x), split="")[[1]]))

indices1 <- which(sig_split[,3]=="C" & sig_split[,6]=="T")
indices2 <- which(sig_split[,3]=="G" & sig_split[,6]=="A")

indices <- union(indices1, indices2);

signature_counts_filter_C_to_T <- signature_counts[,-indices];

signature_set_filter <- colnames(signature_counts_filter_C_to_T)


pr <- prcomp(t(limma::voom(t(signature_counts_filter_C_to_T))$E))

pc_data_frame_filt <- data.frame("PC"=pr$x,
                                 "labels"=c(rep("Ancient",25),
                                            rep("Modern",25)))

qplot(PC.PC1, PC.PC2,
      data=pc_data_frame_filt,
      colour=labels)

qplot(PC.PC2, PC.PC3,
      data=pc_data_frame_filt,
      colour=labels)

plot(pr$center)

sort(pr$rotation[,1], decreasing=TRUE)[1:5]
sort(pr$rotation[,2], decreasing=TRUE)[1:5]

########   CountClust performance (with C-> T and G -> A) ########################

library(CountClust)
topics_clus <- FitGoM(signature_counts,
                      tol=0.1,
                      K=2:4)

 save(topics_clus, file="../rda/CountClust_output_Lindo2016_with_C_to_T.rda")

topics_clus <- get(load("../rda/CountClust_output_Lindo2016_with_C_to_T.rda"));

##############  clusters :  2#######################


omega <- topics_clus$clust_2$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(pc_data_frame$labels)
)

rownames(omega) <- annotation$sample_id;

StructureGGplot(omega = omega,
                annotation = annotation,
                palette = rev(RColorBrewer::brewer.pal(8, "Accent")),
                yaxis_label = "Development Phase",
                order_sample = FALSE,
                figure_title = paste0("StructurePlot: K=", dim(omega)[2],", no C -> T / G -> A"),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))

theta <- topics_clus$clust_2$theta;
sort(theta[,1])[1:10]
sort(theta[,2])[1:10]

apply(ExtractTopFeatures(theta, top_features = 10, method="poisson", options="min"), c(1,2), function(x) signature_set[x])


################## clusters: 3  ########################


omega <- topics_clus$clust_3$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(pc_data_frame$labels)
)

rownames(omega) <- annotation$sample_id;

StructureGGplot(omega = omega,
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

theta <- topics_clus$clust_3$theta;
sort(theta[,1])[1:10]
sort(theta[,2])[1:10]
sort(theta[,3])[1:10]

apply(ExtractTopFeatures(theta, top_features = 10, method="poisson", options="min"), c(1,2), function(x) signature_set[x])

###########  clusters: 4  ###########################
omega <- topics_clus$clust_4$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(pc_data_frame$labels)
)

rownames(omega) <- annotation$sample_id;

StructureGGplot(omega = omega,
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

theta <- topics_clus$clust_4$theta;
sort(theta[,1])[1:10]
sort(theta[,2])[1:10]
sort(theta[,3])[1:10]
sort(theta[,4])[1:10]

apply(ExtractTopFeatures(theta, top_features = 10, method="poisson", options="min"), c(1,2), function(x) signature_set[x])


########   CountClust performance (no C-> T and G -> A) ########################

library(CountClust)
topics_clus <- FitGoM(signature_counts_filter_C_to_T,
                      tol=0.1,
                      K=2:4)

save(topics_clus, file="../rda/CountClust_output_Lindo2016_without_C_to_T.rda")

topics_clus <- get(load("../rda/CountClust_output_Lindo2016_without_C_to_T.rda"))

##################  clusters:  2 ###################################

omega <- topics_clus$clust_2$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(pc_data_frame$labels)
)

rownames(omega) <- annotation$sample_id;

StructureGGplot(omega = omega,
                annotation = annotation,
                palette = rev(RColorBrewer::brewer.pal(8, "Accent")),
                yaxis_label = "Development Phase",
                order_sample = FALSE,
                figure_title = paste0("StructurePlot: K=", dim(omega)[2],", no C -> T / G -> A"),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))

theta <- topics_clus$clust_2$theta;
sort(theta[,1])[1:10]
sort(theta[,2])[1:10]

apply(ExtractTopFeatures(theta, top_features = 10, method="poisson", options="min"), c(1,2), function(x) signature_set_filter[x])


################## clusters: 3  ########################


omega <- topics_clus$clust_3$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(pc_data_frame$labels)
)

rownames(omega) <- annotation$sample_id;

StructureGGplot(omega = omega,
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

theta <- topics_clus$clust_3$theta;
sort(theta[,1])[1:10]
sort(theta[,2])[1:10]
sort(theta[,3])[1:10]

apply(ExtractTopFeatures(theta, top_features = 10, method="poisson", options="min"), c(1,2), function(x) signature_set_filter[x])

################## clusters: 4  ########################


omega <- topics_clus$clust_4$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(pc_data_frame$labels)
)

rownames(omega) <- annotation$sample_id;

StructureGGplot(omega = omega,
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


theta <- topics_clus$clust_4$theta;
sort(theta[,1])[1:10]
sort(theta[,2])[1:10]
sort(theta[,3])[1:10]
sort(theta[,4])[1:10]

apply(ExtractTopFeatures(theta, top_features = 10, method="poisson", options="min"), c(1,2), function(x) signature_set_filter[x])

