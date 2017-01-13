library(aRchaic)
sardinia_counts <- get(load("../processed_data/sardinia2017.rda"))
temp <- club_signature_counts(sardinia_counts)
sardinia_counts_CtoT <- filter_signatures_per_substitution(temp, pattern="C->T", use_prop = FALSE)
sardinia_counts_CtoT <- filter_signatures_only_location(sardinia_counts_CtoT, max_pos = 10)

left_strand_breaks <- as.numeric()
right_strand_breaks <- as.numeric()

strand_breaks <- get(load("../processed_data/strand-breaks-sardinia2017.rda"))
for(l in 1:length(strand_breaks)){
  left_strand_breaks <- rbind(left_strand_breaks, strand_breaks[[l]][[1]])
  right_strand_breaks <- rbind(right_strand_breaks, strand_breaks[[l]][[2]])
}
colnames(left_strand_breaks) <- paste0(colnames(left_strand_breaks), "_left")
colnames(right_strand_breaks) <- paste0(colnames(right_strand_breaks), "_right")


sardinia_filtered <- cbind(sardinia_counts_CtoT, left_strand_breaks, right_strand_breaks)


gosling_counts <- get(load("../processed_data/annagosling2016-counts-table.rda"))
temp <- club_signature_counts(gosling_counts)
temp <- temp[-28,];
gosling_counts_CtoT <- filter_signatures_per_substitution(temp, pattern="C->T", use_prop = FALSE)
gosling_counts_CtoT <- filter_signatures_only_location(gosling_counts_CtoT, max_pos = 10)

left_strand_breaks <- as.numeric()
right_strand_breaks <- as.numeric()

strand_breaks <- get(load("../processed_data/strand-breaks-gosling2016.rda"))
for(l in 1:length(strand_breaks)){
  left_strand_breaks <- rbind(left_strand_breaks, strand_breaks[[l]][[1]])
  right_strand_breaks <- rbind(right_strand_breaks, strand_breaks[[l]][[2]])
}
colnames(left_strand_breaks) <- paste0(colnames(left_strand_breaks), "_left")
colnames(right_strand_breaks) <- paste0(colnames(right_strand_breaks), "_right")


gosling_filtered <- cbind(gosling_counts_CtoT, left_strand_breaks[-28,], right_strand_breaks[-28,])

sherpa_counts <- get(load("../processed_data/sherpa2017.rda"))
temp <- club_signature_counts(sherpa_counts)
sherpa_counts_CtoT <- filter_signatures_per_substitution(temp, pattern="C->T", use_prop = FALSE)
sherpa_counts_CtoT <- filter_signatures_only_location(sherpa_counts_CtoT, max_pos = 10)

left_strand_breaks <- as.numeric()
right_strand_breaks <- as.numeric()

strand_breaks <- get(load("../processed_data/strand-breaks-sherpa.rda"))
for(l in 1:length(strand_breaks)){
  left_strand_breaks <- rbind(left_strand_breaks, strand_breaks[[l]][[1]])
  right_strand_breaks <- rbind(right_strand_breaks, strand_breaks[[l]][[2]])
}
colnames(left_strand_breaks) <- paste0(colnames(left_strand_breaks), "_left")
colnames(right_strand_breaks) <- paste0(colnames(right_strand_breaks), "_right")

sherpa_filtered <- cbind(sherpa_counts_CtoT, left_strand_breaks, right_strand_breaks)

thousandg_counts_clubbed <- get(load("../processed_data/1000Gmoderns-clubbed_counts-table.rda"))
temp <- thousandg_counts_clubbed
thousandg_counts_CtoT <- filter_signatures_per_substitution(temp, pattern="C->T", use_prop = FALSE)
thousandg_counts_CtoT <- filter_signatures_only_location(thousandg_counts_CtoT, max_pos = 10)

left_strand_breaks <- as.numeric()
right_strand_breaks <- as.numeric()

strand_breaks <- get(load("../processed_data/strand-breaks-1000g.rda"))
indices <- match(as.character(paste0(rownames(thousandg_counts_clubbed), ".csv")), as.character(names(strand_breaks)))
for(l in 1:length(strand_breaks)){
  left_strand_breaks <- rbind(left_strand_breaks, strand_breaks[[l]][[1]])
  right_strand_breaks <- rbind(right_strand_breaks, strand_breaks[[l]][[2]])
}
colnames(left_strand_breaks) <- paste0(colnames(left_strand_breaks), "_left")
colnames(right_strand_breaks) <- paste0(colnames(right_strand_breaks), "_right")

left_strand_breaks <- left_strand_breaks[indices, ]
right_strand_breaks <- right_strand_breaks[indices, ]

thousandg_filtered <- cbind(thousandg_counts_CtoT, left_strand_breaks, right_strand_breaks)

hgdp_counts <- get(load("../processed_data/HGDPmoderns-counts-table.rda"))
temp <- club_signature_counts(hgdp_counts)
hgdp_counts_CtoT <- filter_signatures_per_substitution(temp, pattern="C->T", use_prop = FALSE)
hgdp_counts_CtoT <- filter_signatures_only_location(hgdp_counts_CtoT, max_pos = 10)

left_strand_breaks <- as.numeric()
right_strand_breaks <- as.numeric()

strand_breaks <- get(load("../processed_data/strand-breaks-hgdp.rda"))
for(l in 1:length(strand_breaks)){
  left_strand_breaks <- rbind(left_strand_breaks, strand_breaks[[l]][[1]])
  right_strand_breaks <- rbind(right_strand_breaks, strand_breaks[[l]][[2]])
}
colnames(left_strand_breaks) <- paste0(colnames(left_strand_breaks), "_left")
colnames(right_strand_breaks) <- paste0(colnames(right_strand_breaks), "_right")

hgdp_filtered <- cbind(hgdp_counts_CtoT, left_strand_breaks, right_strand_breaks)


names <- rownames(gosling_counts);
control_indices <- c(grep("EXN", names), grep("Libneg", names), grep("PCRneg", names))

gosling_controls_filtered <- gosling_filtered[control_indices,]
gosling_ancient_filtered <- gosling_filtered[-control_indices,]



ll <- list()

ll[["sardinia"]] <- sardinia_filtered
ll[["sherpa"]] <- sherpa_filtered
ll[["gosling-control"]] <- gosling_controls_filtered
ll[["gosling-ancient"]] <- gosling_ancient_filtered
ll[["1000G"]] <- thousandg_filtered[1:100,]
ll[["hgdp"]] <- hgdp_filtered


mat2 <- do.call(rbind, ll);
labs <- unlist(lapply(1:length(ll), function(x) return (rep(names(ll)[x], dim(ll[[x]])[1]))))
tmp <- gridPCA_signatures(mat2[,1:10], factor(labs), normalize = TRUE, cols=c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
                                                                      "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
                                                                      "brown4","darkorchid","magenta","yellow", "azure1","azure4"))



topic_clus <- maptpx::topics(mat2, K=4, tol=0.01)
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
plot_graph(topic_clus$theta[1:10,1], max_pos=10, max_prob=max(topic_clus$theta), main="cluster 1")
plot_graph(topic_clus$theta[1:10,2], max_pos=10, max_prob=max(topic_clus$theta), main="cluster 2")
plot_graph(topic_clus$theta[1:10,3], max_pos=10, max_prob=max(topic_clus$theta), main="cluster 3")
plot_graph(topic_clus$theta[1:10,4], max_pos=10, max_prob=max(topic_clus$theta), main="cluster 4")

topic_clus$theta[11:14,]

thousandg_names <- colnames(ll[["1000G"]])
hgdp_names <- colnames(ll[["hgdp"]])
gosling_ancient_names <-  colnames(ll[["gosling-ancient"]])
gosling_control_names <- colnames(ll[["gosling-control"]])
sherpa_names <- colnames(ll[["sherpa"]])
sardinia_names <- colnames(ll[["sardinia"]])

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
temp2 <- attr(predict(out, testX, probability=TRUE), "probabilities")

library(classtpx)

pooled_data <- rbind(trainX, testX);
class_labs <- c(rep(1, dim(trainX1)[1]), rep(2, dim(trainX2)[1]))
known_samples <- 1:(dim(trainX1)[1] + dim(trainX2)[1])

#pooled_data <- floor(10^12*pooled_data);
Topic_clus <- class_topics(
  as.matrix(pooled_data),
  K=2,
  known_samples = known_samples,
  class_labs = class_labs,
  method="omega.fix",
  shrink=FALSE,
  shrink.method = 1,
  tol=0.001,
  ord=FALSE)


labs <- c(rep("1000g", dim(ll[["1000G"]])[1]), rep("hgdp", dim(ll[["hgdp"]])[1]),
          rep("sherpa", dim(ll[["sherpa"]])[1]), rep("gosling-ancient", dim(ll[["gosling-ancient"]])[1]),
          rep("sardinia", dim(ll[["sardinia"]])[1]), rep("gosling-control", dim(ll[["gosling-control"]])[1]))

labs <- labs[-known_samples]
omega <- Topic_clus$omega[-known_samples, ];
omega <- temp2
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
plot_graph(Topic_clus$theta[1:10,1], max_pos=10, max_prob=max(Topic_clus$theta), main="modern cluster")
plot_graph(Topic_clus$theta[1:10,2], max_pos=10, max_prob=max(Topic_clus$theta), main="ancient cluster")


cols1 <- RColorBrewer::brewer.pal(8, "Accent")
par(mfrow=c(1,2))
plot(Topic_clus$theta[c(11,13,12,14), 2]/ sum(Topic_clus$theta[11:14, 2]), col=cols1[2],
     xaxt="n", pch=20, cex=3, main="left strand break")
axis(1, at=1:4, c("A", "G", "C", "T"))
points(Topic_clus$theta[c(11,13,12,14), 1]/ sum(Topic_clus$theta[11:14, 1]),
     col=cols1[1], xaxt="n", pch=20, cex=3)
axis(1, at=1:4, c("A", "G", "C", "T"))




plot(Topic_clus$theta[c(15,17,16,18), 2]/ sum(Topic_clus$theta[c(15,17,16,18), 2]), col=cols1[2],
     xaxt="n", pch=20, cex=3, main="right strand break")
axis(1, at=1:4, c("A", "G", "C", "T"))
points(Topic_clus$theta[c(15,17,16,18), 1]/ sum(Topic_clus$theta[c(15,17,16,18), 1]),
     col=cols1[1], xaxt="n", pch=20, cex=3)
axis(1, at=1:4, c("A", "G", "C", "T"))


lines(Topic_clus$theta[11:14, 2], col="red", type="b", xaxt="n")

omega[-known_samples, ]
