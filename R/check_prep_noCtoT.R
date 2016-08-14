

signature_counts <- get(load("../summary_data/signature-counts-clubbed-Lindo2016.rda"))

signature_set <- colnames(signature_counts)
sig_split <- do.call(rbind, lapply(colnames(signature_counts), function(x) strsplit(as.character(x), split="")[[1]]))

indices <- which(sig_split[,3]=="C" & sig_split[,6]=="T")

signature_counts_noCtoT <- signature_counts[,-indices];
signature_set_noCtoT <- signature_set[-indices];

sig_split_noCtoT <- sig_split[-indices,];

gsub2 <- function(pattern, replacement, x, ...) {
  for(i in 1:length(pattern))
    x <- gsub(pattern[i], replacement[i], x, ...)
  x
}

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

temp <- get(load("../summary_data/signature-counts-Lindo2016.rda"))
sample_names <- rownames(temp);

lookupSampleInd <- 1:length(sample_names)
names(lookupSampleInd) <- sample_names

lookupFeatInd <- 1:length(mutation_features_noCtoT)
names(lookupFeatInd) <- mutation_features_noCtoT

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

stm <- CheckCounts(signature_counts_noCtoT)

proc_count_noCtoT <- cbind(stm$j, stm$i, stm$v)
rownames(proc_count_noCtoT) <- NULL

mut_features_mat_noCtoT <- data.frame(cbind(substitute, site_left_2, site_left_1, site_right_1, site_right_2))
rownames(mut_features_mat_noCtoT) <- mutation_features_noCtoT
colnames(mut_features_mat_noCtoT)=NULL
mut_features_mat_noCtoT <- as.matrix(apply(mut_features_mat_noCtoT, c(1,2), function(x) as.numeric(x)))

type <- "independent"
flankingBasesNum = as.integer(numBases)
trDir <- FALSE
fdim <- c(6, rep(4, numBases - 1), rep(2, as.integer(trDir)))

library(pmsignature)
G <- new(Class = "MutationFeatureData", 
    type = type,
    flankingBasesNum = as.integer(numBases),
    transcriptionDirection = trDir,
    possibleFeatures = as.integer(fdim),
    featureVectorList = t(mut_features_mat_noCtoT),
    sampleList = sample_names,
    countData = t(proc_count_noCtoT),
    mutationPosition = data.frame()
 )


out <- slot(G, "countData")

Param <- getPMSignature(G, K = 4)

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



