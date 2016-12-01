

##############    Anna Gosling 2016 DNA damage data study  ###########################

dir <- "../data/AnnaGosling2016data/";
out <- aggregate_bin_counts(dir, breaks = c(-1,5,10,15))
save(out,
     file="../processed_data/annagosling2016-counts-table.rda")

#########################   clubbing signatures  ######################################

signature_counts <- get(load("../processed_data/annagosling2016-counts-table.rda"))
validation_check <- club_signature_validation_plot(signature_counts, log=TRUE)
clubbed_counts <- club_signature_counts(signature_counts)

####################################################################

###################  Define labels for ancients and the controls ################

#####################################################################

names <- rownames(signature_counts);
control_indices <- c(grep("EXN", names), grep("Libneg", names), grep("PCRneg", names))

labs <- character();
labs <- rep("ancient", dim(out)[1])
labs[control_indices] <- "controls"

###################################################################

########################   PCA plot  ###############################

###################################################################


gridPCA_signatures(clubbed_counts, labs)

topics_clus <- maptpx::topics(clubbed_counts, K=2, tol=0.1);
save(topics_clus, file="../rda/annagosling2016/topics-annagosling2016-k-2.rda")

topics_clus <- maptpx::topics(clubbed_counts, K=3, tol=0.1);
save(topics_clus, file="../rda/annagosling2016/topics-annagosling2016-k-3.rda")

topics_clus <- maptpx::topics(clubbed_counts, K=4, tol=1);
save(topics_clus, file="../rda/annagosling2016/topics-annagosling2016-k-4.rda")

topics_clus <- maptpx::topics(clubbed_counts, K=5, tol=1);
save(topics_clus, file="../rda/annagosling2016/topics-annagosling2016-k-5.rda")

topics_clus <- get(load("../rda/annagosling2016/topics-annagosling2016-k-3.rda"))

omega <- topics_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs)
)

rownames(omega) <- annotation$sample_id;

CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = FALSE,
                            figure_title = paste0("StructurePlot: K=", dim(omega)[2],": pmsignature: with C->T/G->A"),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))

signature_set <- colnames(clubbed_counts)
apply(CountClust::ExtractTopFeatures(topics_clus$theta, top_features = 50, method="poisson", options="min"), c(1,2), function(x) signature_set[x])


####  it seemed the location of the mutation drove two clusters ######

### As a result, we could not distinguish the controls from ancients ####

### So we now get rid of the location of the reads and just focus on ##

### mutation signatures

filtered_counts <- filter_signatures_wo_location(clubbed_counts);

gridPCA_signatures(filtered_counts, labs)

topics_clus <- maptpx::topics(filtered_counts, K=2, tol=0.1);
save(topics_clus, file="../rda/annagosling2016/topics-annagosling2016-k-2-filtered.rda")

topics_clus <- maptpx::topics(filtered_counts, K=3, tol=0.1);
save(topics_clus, file="../rda/annagosling2016/topics-annagosling2016-k-3-filtered.rda")

topics_clus <- maptpx::topics(filtered_counts, K=4, tol=1);
save(topics_clus, file="../rda/annagosling2016/topics-annagosling2016-k-4-filtered.rda")

topics_clus <- maptpx::topics(filtered_counts, K=5, tol=1);
save(topics_clus, file="../rda/annagosling2016/topics-annagosling2016-k-5-filtered.rda")



topics_clus <- get(load("../rda/annagosling2016/topics-annagosling2016-k-3-filtered.rda"))

omega <- topics_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs)
)

rownames(omega) <- annotation$sample_id;

CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = FALSE,
                            figure_title = paste0("StructurePlot: K=", dim(omega)[2],": pmsignature: with C->T/G->A"),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))

signature_set <- colnames(filtered_counts)
apply(CountClust::ExtractTopFeatures(topics_clus$theta, top_features = 50, method="poisson", options="min"), c(1,2), function(x) signature_set[x])


##################   filter mutation profiles  ###########################

filtered_counts_2 <- filter_signatures_w_mutation(clubbed_counts);

gridPCA_signatures(filtered_counts_2, labs)

topics_clus <- maptpx::topics(filtered_counts_2, K=2, tol=0.1);
save(topics_clus, file="../rda/annagosling2016/topics-annagosling2016-k-2-filtered-mutation.rda")

topics_clus <- maptpx::topics(filtered_counts_2, K=3, tol=0.1);
save(topics_clus, file="../rda/annagosling2016/topics-annagosling2016-k-3-filtered-mutation.rda")

topics_clus <- maptpx::topics(filtered_counts_2, K=4, tol=1);
save(topics_clus, file="../rda/annagosling2016/topics-annagosling2016-k-4-filtered-mutation.rda")

topics_clus <- maptpx::topics(filtered_counts_2, K=5, tol=1);
save(topics_clus, file="../rda/annagosling2016/topics-annagosling2016-k-5-filtered-mutation.rda")

topics_clus <- get(load("../rda/annagosling2016/topics-annagosling2016-k-3-filtered-mutation.rda"))

omega <- topics_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs)
)

rownames(omega) <- annotation$sample_id;

CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = FALSE,
                            figure_title = paste0("StructurePlot: K=", dim(omega)[2],": pmsignature: with C->T/G->A"),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))

signature_set <- colnames(filtered_counts_2)
apply(CountClust::ExtractTopFeatures(topics_clus$theta, top_features = 1, method="poisson", options="min"), c(1,2), function(x) signature_set[x])

################   Focusing on ADR-T2-C2.dup.q30 ####################

###  this is the sample that showed very different behavior from the
###  other ancient samples

topics_clus <- get(load("../rda/annagosling2016/topics-annagosling2016-k-3-filtered.rda"))

which.max(topics_clus$omega[,3])
plot(clubbed_counts[28, ], col="red")

#########  this sample has just one mutation detected  ####################

apply(clubbed_counts, 1, function(x) return(sum(x)))

###########   this sample was eliminated   ###########################

index <- which.max(topics_clus$omega[,3])

clubbed_counts_reduced <- clubbed_counts[-index,];
labs_reduced <- labs[-index]

###########  extracting mutation signatures wo location #################

filtered_counts <- filter_signatures_wo_location(clubbed_counts_reduced);

gridPCA_signatures(filtered_counts, labs_reduced)

topics_clus <- maptpx::topics(filtered_counts, K=2, tol=0.1);
save(topics_clus, file="../rda/annagosling2016/topics-annagosling2016-k-2-filtered_reduced.rda")

topics_clus <- maptpx::topics(filtered_counts, K=3, tol=0.1);
save(topics_clus, file="../rda/annagosling2016/topics-annagosling2016-k-3-filtered_reduced.rda")

topics_clus <- maptpx::topics(filtered_counts, K=4, tol=1);
save(topics_clus, file="../rda/annagosling2016/topics-annagosling2016-k-4-filtered_reduced.rda")

topics_clus <- maptpx::topics(filtered_counts, K=5, tol=1);
save(topics_clus, file="../rda/annagosling2016/topics-annagosling2016-k-5-filtered_reduced.rda")

topics_clus <- get(load("../rda/annagosling2016/topics-annagosling2016-k-3-filtered_reduced.rda"))

omega <- topics_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs_reduced)
)

rownames(omega) <- annotation$sample_id;

CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = FALSE,
                            figure_title = paste0("StructurePlot: K=", dim(omega)[2],": pmsignature: with C->T/G->A"),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))

signature_set <- colnames(filtered_counts)
apply(CountClust::ExtractTopFeatures(topics_clus$theta, top_features = 50, method="poisson", options="min"), c(1,2), function(x) signature_set[x])


###########  extracting mutation signatures w mutation #################

filtered_counts_2 <- filter_signatures_w_mutation(clubbed_counts_reduced);

gridPCA_signatures(filtered_counts_2, labs_reduced)

topics_clus <- maptpx::topics(filtered_counts_2, K=2, tol=0.1);
save(topics_clus, file="../rda/annagosling2016/topics-annagosling2016-k-2-filtered_reduced_w_mutation.rda")

topics_clus <- maptpx::topics(filtered_counts_2, K=3, tol=0.1);
save(topics_clus, file="../rda/annagosling2016/topics-annagosling2016-k-3-filtered_reduced_w_mutation.rda")

topics_clus <- maptpx::topics(filtered_counts_2, K=4, tol=1);
save(topics_clus, file="../rda/annagosling2016/topics-annagosling2016-k-4-filtered_reduced_w_mutation.rda")

topics_clus <- maptpx::topics(filtered_counts_2, K=5, tol=1);
save(topics_clus, file="../rda/annagosling2016/topics-annagosling2016-k-5-filtered_reduced_w_mutation.rda")

topics_clus <- get(load("../rda/annagosling2016/topics-annagosling2016-k-3-filtered_reduced_w_mutation.rda"))

omega <- topics_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs_reduced)
)

rownames(omega) <- annotation$sample_id;

CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = FALSE,
                            figure_title = paste0("StructurePlot: K=", dim(omega)[2],": pmsignature: with C->T/G->A"),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))

signature_set <- colnames(filtered_counts_2)
apply(CountClust::ExtractTopFeatures(topics_clus$theta, top_features = 2, method="poisson", options="min"), c(1,2), function(x) signature_set[x])

names <- rownames(topics_clus$omega);
control_indices <- c(grep("EXN", names), grep("Libneg", names), grep("PCRneg", names))
topics_clus$omega[control_indices,]


#######   proportional mutation signature content ####################

filtered_counts_2 <- filter_signatures_w_mutation(clubbed_counts_reduced);

prop_substitute_sig <- t(apply(filtered_counts_2, 1, function(x) return(x/sum(x))))

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(prop_substitute_sig))),
  tissue_label = factor(labs_reduced)
)

rownames(prop_substitute_sig) <- annotation$sample_id;

StructureGGplot(omega = prop_substitute_sig,
                annotation = annotation,
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Development Phase",
                order_sample = FALSE,
                figure_title = paste0("StructurePlot: 6 substituion patterns (with C -> T)"),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))

############   remove the X from the counts data  #######################


colnames(filtered_counts)[grep("X", colnames(filtered_counts))]
filtered_counts_2_no_X <- filtered_counts_2[,-grep("X", colnames(filtered_counts))]

prop_substitute_sig_no_X <- t(apply(filtered_counts_2_no_X, 1, function(x) return(x/sum(x))))

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(prop_substitute_sig_no_X))),
  tissue_label = factor(labs_reduced)
)

rownames(prop_substitute_sig_no_X) <- annotation$sample_id;

StructureGGplot(omega = prop_substitute_sig_no_X,
                annotation = annotation,
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Development Phase",
                order_sample = FALSE,
                figure_title = paste0("StructurePlot: 6 substituion patterns (with C -> T)"),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))


###  remove the C->T not at the ends   ##########################

filtered_counts <- filter_signatures_wo_location(clubbed_counts_reduced);

CtoTnoXindices <- setdiff( grep("C->T", colnames(filtered_counts)), grep("X", colnames(filtered_counts)))
colnames(filtered_counts)[CtoTnoXindices]

filtered_counts_CtoTjustX <- filtered_counts[,CtoTnoXindices]
filtered_counts_2_CtoTjustX <- filter_signatures_w_mutation(filtered_counts_CtoTjustX);

prop_substitute_sig_CtoTjustX <- t(apply(filtered_counts_2_CtoTjustX, 1, function(x) return(x/sum(x))))

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(prop_substitute_sig_CtoTjustX))),
  tissue_label = factor(labs_reduced)
)

rownames(prop_substitute_sig_CtoTjustX) <- annotation$sample_id;

StructureGGplot(omega = prop_substitute_sig_CtoTjustX,
                annotation = annotation,
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Development Phase",
                order_sample = FALSE,
                figure_title = paste0("StructurePlot: 6 substituion patterns (with just extremal C -> T counted)"),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))



########  Comparative study between extremal C->T and non extremal C->T  #########

filtered_counts <- filter_signatures_wo_location(clubbed_counts_reduced);

names <- colnames(filtered_counts);

names_mod <- as.vector(sapply(names, function(x) {

  if(!is.element("X", strsplit(x, "")[[1]])){
    return(paste(strsplit(x, "")[[1]][3:6], collapse=""))
   } else{
  return(paste(c(strsplit(x, "")[[1]][3:6], " (X)"), collapse=""))
   }
}))

filtered_counts_w_mutation <- numeric();
for(num in 1:dim(filtered_counts)[1]){
  filtered_counts_w_mutation <- rbind(filtered_counts_w_mutation,
                             tapply(filtered_counts[num,],
                                    factor(names_mod), sum))
}

rownames(filtered_counts_w_mutation) <- rownames(filtered_counts)

prop_substitute_sig <- t(apply(filtered_counts_w_mutation, 1, function(x) return(x/sum(x))))

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(prop_substitute_sig))),
  tissue_label = factor(labs_reduced)
)

rownames(prop_substitute_sig) <- annotation$sample_id;

color=c("red","blue","cornflowerblue","black","cyan","darkblue",
          "brown4","burlywood","darkgoldenrod1","darkgray","deepskyblue","darkkhaki",
          "firebrick","darkorchid","hotpink","green","magenta","yellow", "azure1","azure4");

color1 = sample(color, 12)

StructureGGplot(omega = prop_substitute_sig,
                annotation = annotation,
                palette = color1,
                yaxis_label = "Development Phase",
                order_sample = FALSE,
                figure_title = paste0("StructurePlot: 6 substituion patterns (X stands extremal occurrences)"),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))

###  we now run Structure plot on this data  #####################

topics_clus <- maptpx::topics(filtered_counts_w_mutation, K=2, tol=0.1)

omega <- topics_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs_reduced)
)

rownames(omega) <- annotation$sample_id;

CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = FALSE,
                            figure_title = paste0("StructurePlot: K=", dim(omega)[2],": pmsignature: with C->T/G->A"),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))

signature_set <- colnames(filtered_counts_w_mutation)
apply(CountClust::ExtractTopFeatures(topics_clus$theta, top_features = 2, method="poisson", options="min"), c(1,2), function(x) signature_set[x])

##############  maptpx on the data without X   ##################

filtered_counts <- filter_signatures_wo_location(clubbed_counts_reduced);

CtoTnoXindices <- setdiff( grep("C->T", colnames(filtered_counts)), grep("X", colnames(filtered_counts)))
colnames(filtered_counts)[CtoTnoXindices]

filtered_counts_noX <- filtered_counts[,-grep("X", colnames(filtered_counts))]
colnames(filtered_counts_noX)

topics_clus <- maptpx::topics(filtered_counts_noX, K=2, tol=0.1)
save(topics_clus, file="../rda/annagosling2016/topics-annagosling2016-k-2-filtered_reduced_wo_position_noX.rda")

topics_clus <- maptpx::topics(filtered_counts_noX, K=3, tol=0.1)
save(topics_clus, file="../rda/annagosling2016/topics-annagosling2016-k-3-filtered_reduced_wo_position_noX.rda")

topic_clus <- get(load("../rda/annagosling2016/topics-annagosling2016-k-3-filtered_reduced_wo_position_noX.rda"))

omega <- topics_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs_reduced)
)

rownames(omega) <- annotation$sample_id;

CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = FALSE,
                            figure_title = paste0("StructurePlot: K=", dim(omega)[2],": pmsignature: no extremal signatures"),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))

signature_set <- colnames(filtered_counts_noX)
apply(CountClust::ExtractTopFeatures(topics_clus$theta, top_features = 50, method="poisson", options="min"), c(1,2), function(x) signature_set[x])

################  just extremal signatures  #########################


