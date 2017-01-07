

##############   PCA plot of C->T patterns along reads  ##########################

library(aRchaic)
sardinia_counts <- get(load("../processed_data/sardinia2017.rda"))
temp <- club_signature_counts(sardinia_counts)
sardinia_counts_CtoT <- filter_signatures_per_substitution(temp, pattern="C->T", use_prop = FALSE)
sardinia_counts_CtoT <- filter_signatures_only_location(sardinia_counts_CtoT, max_pos = 10)

lindo_ancients_counts <- get(load("../processed_data/lindo2016ancients-counts-table.rda"))
temp <- club_signature_counts(lindo_ancients_counts)
lindo_ancients_counts_CtoT <- filter_signatures_per_substitution(temp, pattern="C->T", use_prop=FALSE)
lindo_ancients_counts_CtoT <- filter_signatures_only_location(lindo_ancients_counts_CtoT, max_pos = 10)


lindo_moderns_counts <- get(load("../processed_data/lindo2016moderns-counts-table.rda"))
temp <- club_signature_counts(lindo_moderns_counts)
lindo_moderns_counts_CtoT <- filter_signatures_per_substitution(temp, pattern="C->T", use_prop = FALSE)
lindo_moderns_counts_CtoT <- filter_signatures_only_location(lindo_moderns_counts_CtoT, max_pos = 10)

gosling_counts <- get(load("../processed_data/annagosling2016-counts-table.rda"))
temp <- club_signature_counts(gosling_counts)
temp <- temp[-28,];
gosling_counts_CtoT <- filter_signatures_per_substitution(temp, pattern="C->T", use_prop = FALSE)
gosling_counts_CtoT <- filter_signatures_only_location(gosling_counts_CtoT, max_pos = 10)

sherpa_counts <- get(load("../processed_data/sherpa2017.rda"))
temp <- club_signature_counts(sherpa_counts)
sherpa_counts_CtoT <- filter_signatures_per_substitution(temp, pattern="C->T", use_prop = FALSE)
sherpa_counts_CtoT <- filter_signatures_only_location(sherpa_counts_CtoT, max_pos = 10)

thousandg_counts_clubbed <- get(load("../processed_data/1000Gmoderns-clubbed_counts-table.rda"))
temp <- thousandg_counts_clubbed
thousandg_counts_CtoT <- filter_signatures_per_substitution(temp, pattern="C->T", use_prop = FALSE)
thousandg_counts_CtoT <- filter_signatures_only_location(thousandg_counts_CtoT, max_pos = 10)

hgdp_counts <- get(load("../processed_data/HGDPmoderns-counts-table.rda"))
temp <- club_signature_counts(hgdp_counts)
hgdp_counts_CtoT <- filter_signatures_per_substitution(temp, pattern="C->T", use_prop = FALSE)
hgdp_counts_CtoT <- filter_signatures_only_location(hgdp_counts_CtoT, max_pos = 10)

names <- rownames(gosling_counts);
control_indices <- c(grep("EXN", names), grep("Libneg", names), grep("PCRneg", names))

gosling_controls_counts_CtoT <- gosling_counts_CtoT[control_indices,]
gosling_ancient_counts_CtoT <- gosling_counts_CtoT[-control_indices,]

I_counts <- get(load("../processed_data/Idata-counts-table.rda"))
temp <- club_signature_counts(I_counts)
I_counts_CtoT <- filter_signatures_per_substitution(temp, pattern="C->T", use_prop = FALSE)
I_counts_CtoT <- filter_signatures_only_location(I_counts_CtoT, max_pos = 10)

RISE_counts <- get(load("../processed_data/RISE-counts-table.rda"))
temp <- club_signature_counts(RISE_counts)
RISE_counts_CtoT <- filter_signatures_per_substitution(temp, pattern="C->T", use_prop = FALSE)
RISE_counts_CtoT <- filter_signatures_only_location(RISE_counts_CtoT, max_pos = 10)



ll <- list()
#ll[["lindo-moderns"]] <- lindo_moderns_counts_CtoT
#ll[["lindo-old"]] <- lindo_ancients_counts_CtoT
ll[["sardinia"]] <- sardinia_counts_CtoT
ll[["sherpa"]] <- sherpa_counts_CtoT
ll[["gosling-control"]] <- gosling_controls_counts_CtoT
ll[["gosling-ancient"]] <- gosling_ancient_counts_CtoT
ll[["1000G"]] <- thousandg_counts_CtoT[1:100,]
ll[["hgdp"]] <- hgdp_counts_CtoT
#ll[["I"]] <- I_counts_CtoT
#ll[["RISE"]] <- RISE_counts_CtoT

library(aRchaic)
out <- aRchaic::gridPCA_position_signatures_combo(ll,
                                  input_pos=1:20,
                                  normalize=FALSE,
                                  cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
                                           "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
                                           "brown4","darkorchid","magenta","yellow", "azure1","azure4"))

mat <- do.call(rbind, ll);
labs <- unlist(lapply(1:length(ll), function(x) return (rep(names(ll)[x], dim(ll[[x]])[1]))))
tmp <- gridPCA_signatures(mat, factor(labs), normalize = TRUE, cols=c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
                                                               "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
                                                               "brown4","darkorchid","magenta","yellow", "azure1","azure4"))


topic_clus <- maptpx::topics(mat, K=4, tol=0.01)
omega <- topic_clus$omega;
annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs))


CountClust::StructureGGplot(omega = as.matrix(omega),
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Tissue type",
                            order_sample = FALSE,
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))

par(mfrow=c(2,2))
plot_graph(topic_clus$theta[,1], max_pos=10, max_prob=max(topic_clus$theta[,1]), main="cluster 1")
plot_graph(topic_clus$theta[,2], max_pos=10, max_prob=max(topic_clus$theta[,2]), main="cluster 2")
plot_graph(topic_clus$theta[,3], max_pos=10, max_prob=max(topic_clus$theta[,3]), main="cluster 3")
plot_graph(topic_clus$theta[,4], max_pos=10, max_prob=max(topic_clus$theta[,4]), main="cluster 4")

#################  Apply svm classifier to get contamination rate  #####################################


thousandg_names <- colnames(ll[["1000G"]])
hgdp_names <- colnames(ll[["hgdp"]])
gosling_ancient_names <-  colnames(ll[["gosling-ancient"]])
gosling_control_names <- colnames(ll[["gosling-control"]])
sherpa_names <- colnames(ll[["sherpa"]])
sardinia_names <- colnames(ll[["sardinia"]])
#lindo_moderns_names <- colnames(ll[["lindo-moderns"]])
#lindo_ancient_names <- colnames(ll[["lindo-old"]])

matched_names <- Reduce(intersect, list(thousandg_names, hgdp_names, gosling_ancient_names,
                                        gosling_control_names, sherpa_names, sardinia_names))
                                        #lindo_moderns_names, lindo_ancient_names))

trainX1 <- rbind(ll[["1000G"]][, match(matched_names, thousandg_names)],
                ll[["hgdp"]][, match(matched_names, hgdp_names)])

trainX2 <- rbind(ll[["sherpa"]][, match(matched_names, sherpa_names)],
                 ll[["gosling-ancient"]][, match(matched_names, gosling_ancient_names)])



trainX <- rbind(trainX1, trainX2)

testX <- rbind(ll[["sardinia"]][, match(matched_names, sardinia_names)],
               ll[["gosling-control"]][, match(matched_names, gosling_control_names)])

y <- c(rep("modern", dim(trainX1)[1]), rep("ancient", dim(trainX2)[1]))
y <- factor(y)

data <- cbind.data.frame(y, trainX);
library(e1071)
out <- svm(y ~ ., data=data, probability=TRUE)
predict(out, testX, probability=TRUE)

########  Classification Grade of Memebrship Models  ##############################

library(classtpx)

pooled_data <- rbind(trainX, testX);
class_labs <- c(rep(1, dim(trainX1)[1]), rep(2, dim(trainX2)[1]))
known_samples <- 1:(dim(trainX1)[1] + dim(trainX2)[1])

pooled_data <- floor(10^12*pooled_data);
Topic_clus <- class_topics(
  as.matrix(pooled_data),
  K=2,
  known_samples = known_samples,
  class_labs = class_labs,
  method="theta.fix",
  shrink=FALSE,
  shrink.method = 1,
  tol=0.0001,
  ord=FALSE)

labs <- c(rep("1000g", dim(ll[["1000G"]])[1]), rep("hgdp", dim(ll[["hgdp"]])[1]),
          rep("sherpa", dim(ll[["sherpa"]])[1]), rep("gosling-ancient", dim(ll[["gosling-ancient"]])[1]),
          rep("sardinia", dim(ll[["sardinia"]])[1]), rep("gosling-control", dim(ll[["gosling-control"]])[1]))

omega <- Topic_clus$omega;
annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs))


CountClust::StructureGGplot(omega = as.matrix(omega),
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Tissue type",
                            order_sample = FALSE,
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))

par(mfrow=c(1,2))
plot_graph(Topic_clus$theta[,1], max_pos=10, max_prob=max(Topic_clus$theta[,1]), main="cluster 1")
plot_graph(Topic_clus$theta[,2], max_pos=10, max_prob=max(Topic_clus$theta[,2]), main="cluster 2")


pooled_data <- rbind(trainX, testX);
class_labs <- c(rep(1, dim(trainX1)[1]))
known_samples <- 1:(dim(trainX1)[1])
#pooled_data <- floor(10^12*pooled_data);
Topic_clus <- class_topics(
  as.matrix(pooled_data),
  K=2,
  known_samples = known_samples,
  class_labs = class_labs,
  method="theta.fix",
  shrink=FALSE,
  shrink.method = 1,
  tol=0.001,
  ord=FALSE)

omega <- Topic_clus$omega;
annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs))


CountClust::StructureGGplot(omega = as.matrix(omega),
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Tissue type",
                            order_sample = TRUE,
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))


temp_omega <- tail(omega, 12)
