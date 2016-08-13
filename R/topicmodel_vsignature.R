
###  Visulaization of topic model on damage data

signature_counts <- get(load("../summary_data/signature-counts-clubbed-Lindo2016.rda"))

signature_set <- colnames(signature_counts)
sig_split <- do.call(rbind, lapply(colnames(signature_counts), function(x) strsplit(as.character(x), split="")[[1]]))

signature_counts_noCtoT <- signature_counts[,-indices];
signature_set_noCtoT <- signature_set[-indices];

sig_split_noCtoT <- sig_split[-indices,];


site_left_2 <- gsub2(c("A","C","G","T"), c(1,2,3,4), x=sig_split_noCtoT[,1])
site_left_1 <- gsub2(c("A","C","G","T"), c(1,2,3,4), x=sig_split_noCtoT[,2])
site_right_1 <- gsub2(c("A","C","G","T"), c(1,2,3,4), x=sig_split_noCtoT[,7])
site_right_2 <- gsub2(c("A","C","G","T"), c(1,2,3,4), x=sig_split_noCtoT[,8])

sub_pattern <- sapply(1:dim(sig_split_noCtoT)[1], function(x) paste(sig_split_noCtoT[x,3:6], collapse=""))

from = c("C->A", "C->G", "C->T", "T->A", "T->C", "T->G")
to = c(1,2,3,4,5,6)

substitute <- gsub2(from, to, sub_pattern)

mutation_features_noCtoT <- sapply(1:dim(sig_split_noCtoT)[1], 
                                   function(m) paste0(substitute[m], ",", site_left_2[m], ",", site_left_1[m], ",", site_right_1[m], ",", site_right_2[m]))

mut_features_mat_noCtoT <- data.frame(cbind(substitute, site_left_2, site_left_1, site_right_1, site_right_2))
rownames(mut_features_mat_noCtoT) <- mutation_features_noCtoT
colnames(mut_features_mat_noCtoT)=NULL
mut_features_mat_noCtoT <- as.matrix(apply(mut_features_mat_noCtoT, c(1,2), function(x) as.numeric(x)))

topics_clus <- FitGoM(signature_counts_noCtoT,
                      tol=0.1,
                      K=2:4)

save(topics_clus, file="pmsignature_fit_noCtoT_2_to_4.rda")

topics_clus <- get(load("pmsignature_fit_noCtoT_2_to_4.rda"))

theta <- topics_clus$clust_4$theta;
omega <- topics_clus$clust_4$omega

K <- dim(omega)[2]
F <- array(0, c(K, length(fdim), max(fdim)))
for(k in 1:K){
  for(kk in 1:length(fdim)){
    temp <- tapply(theta[,k], mut_features_mat_noCtoT[,kk], sum)
    F[k,kk,as.numeric(names(temp))] <- as.numeric(temp)
  }
}

isBG <- TRUE
Param <- new(Class = "EstimatedParameters", 
           type = slot(mutationFeatureData, "type"),
           flankingBasesNum = slot(mutationFeatureData, "flankingBasesNum"),
           transcriptionDirection = slot(mutationFeatureData, "transcriptionDirection"),
           possibleFeatures = slot(mutationFeatureData, "possibleFeatures"),
           sampleList = slot(mutationFeatureData, "sampleList"),
           signatureNum = as.integer(K),
           isBackGround = isBG,
           signatureFeatureDistribution = F,
           sampleSignatureDistribution = omega,
           loglikelihood = -100000)

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
                figure_title = paste0("StructurePlot: K=", dim(omega)[2],": pmsignature: without C->T/G->A"),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))

apply(ExtractTopFeatures(theta, top_features = 10, method="poisson", options="min"), c(1,2), function(x) signature_set_noCtoT [x])

