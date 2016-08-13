

#######   Class preparation for pm signature  ####################

####  This script is aimed at creating the MutationFeatureData class
####  from the Lindo2016 data that we shall use for pmsignature model.


signature_counts <- get(load("../summary_data/signature-counts-clubbed-Lindo2016.rda"))

signature_set <- colnames(signature_counts)
sig_split <- do.call(rbind, lapply(colnames(signature_counts), function(x) strsplit(as.character(x), split="")[[1]]))

sig_split[,1]

gsub2 <- function(pattern, replacement, x, ...) {
  for(i in 1:length(pattern))
    x <- gsub(pattern[i], replacement[i], x, ...)
  x
}

site_left_2 <- gsub2(c("A","C","G","T"), c(1,2,3,4), x=sig_split[,1])
site_left_1 <- gsub2(c("A","C","G","T"), c(1,2,3,4), x=sig_split[,2])
site_right_1 <- gsub2(c("A","C","G","T"), c(1,2,3,4), x=sig_split[,7])
site_right_2 <- gsub2(c("A","C","G","T"), c(1,2,3,4), x=sig_split[,8])

sub_pattern <- sapply(1:dim(sig_split)[1], function(x) paste(sig_split[x,3:6], collapse=""))

from = c("C->A", "C->G", "C->T", "T->A", "T->C", "T->G")
to = c(1,2,3,4,5,6)

substitute <- gsub2(from, to, sub_pattern)

mutation_features <- sapply(1:dim(sig_split)[1], 
                           function(m) paste0(substitute[m], ",", site_left_2[m], ",", site_left_1[m], ",", site_right_1[m], ",", site_right_2[m]))


temp <- get(load("../summary_data/signature-counts-Lindo2016.rda"))
sample_names <- rownames(temp);

lookupSampleInd <- 1:length(sample_names)
names(lookupSampleInd) <- sample_names

lookupFeatInd <- 1:length(mutation_features)
names(lookupFeatInd) <- mutation_features


w <- which(tableCount > 0, arr.ind=TRUE)

row_col_indices <- which(signature_counts > 0, arr.ind=TRUE)
out <- slam::simple_triplet_matrix(signature_counts)
out$

library(slam)
CheckCounts <- function(counts){
    if(class(counts)[1] == "TermDocumentMatrix"){ counts <- t(counts) }
    if(is.null(dimnames(counts)[[1]])){ dimnames(counts)[[1]] <- paste("doc",1:nrow(counts)) }
    if(is.null(dimnames(counts)[[2]])){ dimnames(counts)[[2]] <- paste("wrd",1:ncol(counts)) }
    empty <- row_sums(counts) == 0
    if(sum(empty) != 0){
      counts <- counts[!empty,]
      cat(paste("Removed", sum(empty), "blank documents.\n")) }
    return(as.simple_triplet_matrix(counts))
}

stm <- CheckCounts(signature_counts)

proc_count <- cbind(stm$j, stm$i, stm$v)
rownames(proc_count) <- NULL

mut_features_mat <- data.frame(cbind(substitute, site_left_2, site_left_1, site_right_1, site_right_2))
rownames(mut_features_mat) <- mutation_features
colnames(mut_features_mat)=NULL
mut_features_mat <- as.matrix(apply(mut_features_mat, c(1,2), function(x) as.numeric(x)))

type <- "independent"
flankingBasesNum = as.integer(numBases)
trDir <- FALSE
fdim <- c(6, rep(4, numBases - 1), rep(2, as.integer(trDir)))

library(pmsignature)
new(Class = "MutationFeatureData", 
    type = type,
    flankingBasesNum = as.integer(numBases),
    transcriptionDirection = trDir,
    possibleFeatures = as.integer(fdim),
    featureVectorList = t(mut_features_mat),
    sampleList = sample_names,
    countData = t(proc_count),
    mutationPosition = data.frame()
)

G <- new(Class = "MutationFeatureData", 
    type = type,
    flankingBasesNum = as.integer(numBases),
    transcriptionDirection = trDir,
    possibleFeatures = as.integer(fdim),
    featureVectorList = t(mut_features_mat),
    sampleList = sample_names,
    countData = t(proc_count),
    mutationPosition = data.frame()
)

out <- slot(G, "countData")

Param <- getPMSignature(G, K = 4)
save(Param, file="../rda/pmsignature_fit_K_4.rda")
visPMSignature(Param, 1)
visPMSignature(Param, 2)
visPMSignature(Param, 3)
visPMSignature(Param, 4)

omega <- slot(Param, "sampleSignatureDistribution")
library(CountClust)
annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(c(rep("Ancient",25), rep("Modern",25)))
)

rownames(omega) <- annotation$sample_id;

StructureGGplot(omega = omega,
                annotation = annotation,
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Development Phase",
                order_sample = FALSE,
                figure_title = paste0("StructurePlot: K=", dim(omega)[2],": pmsignature: with C->T/G->A"),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))


#########   Reperform topic model removing C -> T  ###########


indices_noCtoT <- which(mut_features_mat[,1]!=3);

dim(signature_counts)
signature_counts_noCtoT <- signature_counts[,indices_noCtoT];

mut_features_mat_noCtoT <-  mut_features_mat[indices_noCtoT,];

stm <- CheckCounts(signature_counts_noCtoT)

proc_count_noCtoT <- cbind(stm$j, stm$i, stm$v)
rownames(proc_count_noCtoT) <- NULL

G <- new(Class = "MutationFeatureData", 
         type = "custom",
         flankingBasesNum = as.integer(numBases),
         transcriptionDirection = trDir,
         possibleFeatures = as.integer(fdim),
         featureVectorList = t(mut_features_mat_noCtoT),
         sampleList = sample_names,
         countData = t(proc_count_noCtoT),
         mutationPosition = data.frame()
)

Param <- getPMSignature(G, K = 4, numInit = 5)
visPMSignature(Param, 1)
visPMSignature(Param, 2)
visPMSignature(Param, 3)
visPMSignature(Param, 4)

omega <- slot(Param, "sampleSignatureDistribution")
library(CountClust)
annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(c(rep("Ancient",25), rep("Modern",25)))
)

rownames(omega) <- annotation$sample_id;

StructureGGplot(omega = omega,
                annotation = annotation,
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Development Phase",
                order_sample = FALSE,
                figure_title = paste0("StructurePlot: K=", dim(omega)[2],": pmsignature: with C->T/G->A"),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))

###########################################################

signature_counts_noCtoT 

pr <- prcomp(t(limma::voom(t(signature_counts_noCtoT))$E))

pc_data_frame <- data.frame("PC"=pr$x,
                            "labels"=c(rep("Ancient",25),
                                       rep("Modern",25)))

qplot(PC.PC1, PC.PC2,
      data=pc_data_frame,
      colour=labels)

library(CountClust)
topics_clus <- FitGoM(signature_counts_noCtoT,
                      tol=0.1,
                      K=2:4)

save(topics_clus, file="../rda/CountClust_output_Lindo2016_without_C_to_T.rda")

topics_clus <- get(load("../rda/CountClust_output_Lindo2016_without_C_to_T.rda"));

omega <- topics_clus$clust_4$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = c(rep("Ancient",25),rep("Modern",25))
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

theta <- topics_clus$clust_4$theta;
sort(theta[,1])[1:10]
sort(theta[,2])[1:10]

signature_set_noCtoT <- signature_set[indices_noCtoT]
apply(ExtractTopFeatures(theta, top_features = 10, method="poisson", options="min"), c(1,2), function(x) signature_set_noCtoT [x])

