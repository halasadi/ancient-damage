

############   Lindo 2016 classtpx application  ###################

library(slam)
library(maptpx)
source("../../classtpx/R/class_count.R")
source("../../classtpx/R/class_tpx.R")
source("../../classtpx/R/class_topics.R")
source("../../classtpx/R/class_varselect.R")
source("../../classtpx/R/tpx.R")
source("../../classtpx/R/count.R")



signature_counts <- get(load("../summary_data/signature-counts-clubbed-Lindo2016.rda"))
signature_set <- colnames(signature_counts)


class_labs <- rep(1,25)
known_samples <- 26:50

Topic_clus <- class_topics(
  signature_counts, 
  K=2, 
  known_samples = known_samples,
  class_labs = class_labs,
  method="theta.fix",
  shrink=FALSE,
  tol=0.01,
  prior_omega = c(0.8,0.2),
  ord=FALSE)

omega <- Topic_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(c(rep("Ancient",25), rep("Modern",25)))
)

rownames(omega) <- annotation$sample_id;

CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = FALSE,
                            figure_title = "StructurePlot",
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))

indices <- CountClust::ExtractTopFeatures(Topic_clus$theta, top_features=10, method="poisson", options="min")
imp_features <- apply(indices, c(1,2), function(x) signature_set[x])

imp_features

Topic_clus <- class_topics(
  signature_counts, 
  K=3, 
  known_samples = known_samples,
  class_labs = class_labs,
  method="theta.fix",
  shrink=FALSE,
  tol=0.01,
  prior_omega = c(0.95,0.025, 0.025),
  ord=FALSE)

omega <- Topic_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(c(rep("Ancient",25), rep("Modern",25)))
)

rownames(omega) <- annotation$sample_id;

CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = FALSE,
                            figure_title = "StructurePlot",
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))

indices <- CountClust::ExtractTopFeatures(Topic_clus$theta, top_features=10, method="poisson", options="min")
imp_features <- apply(indices, c(1,2), function(x) signature_set[x])

imp_features


sig_split <- do.call(rbind, lapply(colnames(signature_counts), function(x) strsplit(as.character(x), split="")[[1]]))

indices <- which(sig_split[,3]=="C" & sig_split[,6]=="T")
signature_counts_noCtoT <- signature_counts[,-indices];
signature_set_noCtoT <- signature_set[-indices];


class_labs <- rep(1,25)
known_samples <- 26:50

Topic_clus <- class_topics(
  signature_counts_noCtoT, 
  K=2, 
  known_samples = known_samples,
  class_labs = class_labs,
  method="theta.fix",
  shrink=FALSE,
  tol=0.01,
  prior_omega = c(0.5,0.5),
  ord=FALSE)

omega <- Topic_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(c(rep("Ancient",25), rep("Modern",25)))
)

rownames(omega) <- annotation$sample_id;

CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = FALSE,
                            figure_title = "StructurePlot",
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))

Topic_clus <- class_topics(
  signature_counts_noCtoT, 
  K=3, 
  known_samples = known_samples,
  class_labs = class_labs,
  method="theta.fix",
  shrink=FALSE,
  tol=0.01,
  prior_omega = c(0.9999,0.001,0.001),
  ord=FALSE)

omega <- Topic_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(c(rep("Ancient",25), rep("Modern",25)))
)

rownames(omega) <- annotation$sample_id;

CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = FALSE,
                            figure_title = "StructurePlot",
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))

indices <- CountClust::ExtractTopFeatures(Topic_clus$theta, top_features=10, method="poisson", options="min")
imp_features <- apply(indices, c(1,2), function(x) signature_set_noCtoT[x])

imp_features

Topic_clus <- class_topics(
  signature_counts_noCtoT, 
  K=3, 
  known_samples = known_samples,
  class_labs = class_labs,
  method="theta.fix",
  optional_theta = rep(1/length(signature_set_noCtoT), length(signature_set_noCtoT)),
  shrink=FALSE,
  tol=0.01,
  prior_omega = c(0.9999,0.001,0.001),
  ord=FALSE)

omega <- Topic_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(c(rep("Ancient",25), rep("Modern",25)))
)

rownames(omega) <- annotation$sample_id;

CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = FALSE,
                            figure_title = "StructurePlot",
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))


Topic_clus <- classtpx::class_topics(
  signature_counts_noCtoT, 
  K=4, 
  known_samples = known_samples,
  class_labs = class_labs,
  method="theta.fix",
  optional_theta = rep(1/length(signature_set_noCtoT), length(signature_set_noCtoT)),
  shrink=FALSE,
  tol=0.01,
  prior_omega = c(0.9999,0.001,0.001,0.001),
  ord=FALSE)

omega <- Topic_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(c(rep("Ancient",25), rep("Modern",25)))
)

rownames(omega) <- annotation$sample_id;

CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = FALSE,
                            figure_title = "StructurePlot",
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))

indices <- CountClust::ExtractTopFeatures(Topic_clus$theta, top_features=10, method="poisson", options="min")
imp_features <- apply(indices, c(1,2), function(x) signature_set_noCtoT[x])

imp_features

