


##################   Sardinia 2017 data  ################################

library(aRchaic)
dir <- "../data/Sardinia2017/";
out <- aggregate_bin_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))

save(out,
     file="../processed_data/sardinia2017.rda")

signature_counts <- get(load("../processed_data/sardinia2017.rda"))
validation_check <- club_signature_validation_plot(signature_counts, log=TRUE)
clubbed_counts_sardinia_ancient <- club_signature_counts(signature_counts)
filtered_counts_sardinia_ancient <- filter_signatures_wo_location(clubbed_counts_sardinia_ancient);


signature_counts <- get(load("../processed_data/lindo2016ancients-counts-table.rda"))
validation_check <- club_signature_validation_plot(signature_counts, log=TRUE)
clubbed_counts_lindo_ancient <- club_signature_counts(signature_counts)
filtered_counts_lindo_ancient <- filter_signatures_wo_location(clubbed_counts_lindo_ancient);

signature_counts <- get(load("../processed_data/lindo2016moderns-counts-table.rda"))
validation_check <- club_signature_validation_plot(signature_counts, log=TRUE)
clubbed_counts_lindo_modern <- club_signature_counts(signature_counts)
filtered_counts_lindo_modern <- filter_signatures_wo_location(clubbed_counts_lindo_modern);


common_signatures <- intersect(colnames(clubbed_counts_lindo_modern),
                               intersect(colnames(clubbed_counts_lindo_ancient),
                               colnames(clubbed_counts_sardinia_ancient)))

clubbed_counts_lindo_modern_1 <- clubbed_counts_lindo_modern[,match(common_signatures, colnames(clubbed_counts_lindo_modern))]
clubbed_counts_lindo_ancient_1 <- clubbed_counts_lindo_ancient[,match(common_signatures, colnames(clubbed_counts_lindo_ancient))]
clubbed_counts_sardinia_ancient_1 <- clubbed_counts_sardinia_ancient[,match(common_signatures, colnames(clubbed_counts_sardinia_ancient))]


clubbed_counts_pooled <- rbind(clubbed_counts_lindo_ancient_1, clubbed_counts_lindo_modern_1,
                                clubbed_counts_sardinia_ancient_1)



signature_set <- colnames(clubbed_counts_pooled)
sig_split <- t(sapply(1:length(signature_set), function(x) return(strsplit(signature_set[x], "")[[1]][1:8])))
new_sig_split <- matrix(0, dim(sig_split)[1], 6);
new_sig_split[,1] <- sig_split[, 1]
new_sig_split[,2] <- sig_split[, 2]
new_sig_split[,3] <- sapply(1:length(signature_set), function(x) return(paste(sig_split[x,3:6], collapse="")))
new_sig_split[,4] <- sig_split[, 7]
new_sig_split[,5] <- sig_split[, 8]

position_sig <- sapply(1:length(signature_set), function(x)
{
  tmp <- strsplit(signature_set[x], "")[[1]];
  return(paste(tmp[10:length(tmp)], collapse=""))
})

new_sig_split[,6]  <- factor(position_sig, levels=0:max(as.numeric(position_sig)))

mat <- matrix(0, dim(new_sig_split)[1], dim(new_sig_split)[2])
for(k in 1:dim(new_sig_split)[2]){
  if(k < dim(new_sig_split)[2]){
      temp <- as.factor(new_sig_split[,k])
      mat[,k] <- as.numeric(as.matrix(plyr::mapvalues(temp, from = levels(temp), to = 0:(length(levels(temp))-1))))
  }else{
      mat[,k] <- as.numeric(as.matrix(new_sig_split[,k]))
  }
}

signatures <- mat;

topics_clus_1 <- maptpx::topics(clubbed_counts_pooled, K=5, type="full", tol=10)
topics_clus_2 <- maptpx::topics(clubbed_counts_pooled, K=5, type="independent", tol=10, signatures = signatures)

labs <- c(rep("Lindo Ancient",25), rep("Lindo Modern",25), rep("Sardinia Ancient",5))

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

signature_set <- colnames(clubbed_counts_pooled)
apply(CountClust::ExtractTopFeatures(topics_clus_1$theta, top_features = 10, method="poisson", options="min"), c(1,2), function(x) signature_set[x])

signature_counts <- get(load("../processed_data/sardinia2017.rda"))

file <- "../data/Sardinia2017/ISB001.A0201_S0_L008_R1_001.fastq.fq.qF.sorted.cleaned_rmdup.q30.csv"


###############################   C ->  T  ###################################################


pattern_plot("../data/Sardinia2017/ISB001.A0201_S0_L008_R1_001.fastq.fq.qF.sorted.cleaned_rmdup.q30.csv", pattern="C->T", plot_type = "left", use_prop = TRUE)
pattern_plot("../data/Sardinia2017/LON001.A0201_S0_L008_R1_001.fastq.fq.qF.sorted.cleaned_rmdup.q30.csv", pattern="C->T", plot_type = "left", use_prop = TRUE)
pattern_plot("../data/Sardinia2017/SUA001.A0201_S0_L008_R1_001.fastq.fq.qF.sorted.cleaned_rmdup.q30.csv", pattern="C->T", plot_type = "left", use_prop = TRUE)
pattern_plot("../data/Sardinia2017/SUA002.A0201_S0_L008_R1_001.fastq.fq.qF.sorted.cleaned_rmdup.q30.csv", pattern="C->T", plot_type = "left", use_prop = TRUE)
pattern_plot("../data/Sardinia2017/SUA003.A0201_S0_L008_R1_001.fastq.fq.qF.sorted.cleaned_rmdup.q30.csv", pattern="C->T", plot_type = "left", use_prop = TRUE)

pattern_plot("../data/Lindo2016ancients/125_all_chr.q30.csv", pattern="C->T", plot_type = "left", use_prop = TRUE)
pattern_plot("../data/Lindo2016ancients/158_all_chr.q30.csv", pattern="C->T", plot_type = "left", use_prop = TRUE)
pattern_plot("../data/Lindo2016ancients/300_all_chr.q30.csv", pattern="C->T", plot_type = "left", use_prop = TRUE)
pattern_plot("../data/Lindo2016ancients/939_all_chr.q30.csv", pattern="C->T", plot_type = "left", use_prop = TRUE)


pattern_plot("../data/Lindo2016moderns/S001_all_chr.q30.csv", pattern="C->T", plot_type = "left", use_prop = TRUE)
pattern_plot("../data/Lindo2016moderns/S002_all_chr.q30.csv", pattern="C->T", plot_type = "left", use_prop = TRUE)
pattern_plot("../data/Lindo2016moderns/T008_all_chr.q30.csv", pattern="C->T", plot_type = "left", use_prop = TRUE)

##########################################   T -> A    #####################################################

pattern_plot("../data/Sardinia2017/ISB001.A0201_S0_L008_R1_001.fastq.fq.qF.sorted.cleaned_rmdup.q30.csv", pattern="T->A", plot_type = "right", use_prop = TRUE)
pattern_plot("../data/Sardinia2017/LON001.A0201_S0_L008_R1_001.fastq.fq.qF.sorted.cleaned_rmdup.q30.csv", pattern="T->A", plot_type = "right", use_prop = TRUE)
pattern_plot("../data/Sardinia2017/SUA001.A0201_S0_L008_R1_001.fastq.fq.qF.sorted.cleaned_rmdup.q30.csv", pattern="T->A", plot_type = "right", use_prop = TRUE)
pattern_plot("../data/Sardinia2017/SUA002.A0201_S0_L008_R1_001.fastq.fq.qF.sorted.cleaned_rmdup.q30.csv", pattern="T->A", plot_type = "right", use_prop = TRUE)
pattern_plot("../data/Sardinia2017/SUA003.A0201_S0_L008_R1_001.fastq.fq.qF.sorted.cleaned_rmdup.q30.csv", pattern="T->A", plot_type = "right", use_prop = TRUE)

pattern_plot("../data/Lindo2016ancients/125_all_chr.q30.csv", pattern="T->A", plot_type = "left", use_prop = TRUE)
pattern_plot("../data/Lindo2016ancients/158_all_chr.q30.csv", pattern="T->A", plot_type = "left", use_prop = TRUE)
pattern_plot("../data/Lindo2016ancients/300_all_chr.q30.csv", pattern="T->A", plot_type = "left", use_prop = TRUE)
pattern_plot("../data/Lindo2016ancients/939_all_chr.q30.csv", pattern="T->A", plot_type = "right", use_prop = TRUE)


pattern_plot("../data/Lindo2016moderns/S001_all_chr.q30.csv", pattern="T->A", plot_type = "left", use_prop = TRUE)
pattern_plot("../data/Lindo2016moderns/S002_all_chr.q30.csv", pattern="T->A", plot_type = "left", use_prop = TRUE)
pattern_plot("../data/Lindo2016moderns/T008_all_chr.q30.csv", pattern="T->A", plot_type = "left", use_prop = TRUE)


########  PCA plot of the Sardinian mutation profile data #########################


gridPCA_signatures(clubbed_counts_pooled, labs)

###########   Full pattern plot analysis  ##################################################

color=c("red","blue","cornflowerblue","black","cyan","darkblue",
        "brown4","burlywood","darkgoldenrod1","darkgray","deepskyblue","darkkhaki",
        "firebrick","darkorchid","hotpink","green","magenta","yellow", "azure1","azure4");


pattern_plot_full(file="../data/Sardinia2017/ISB001.A0201_S0_L008_R1_001.fastq.fq.qF.sorted.cleaned_rmdup.q30.csv",
                  pattern = c("C->T", "C->G", "C->A", "T->G", "T->A", "T->C",
                              "G->A", "G->C", "G->T", "A_>G", "A->T", "A->C"),
                  plot_type="left",
                  sample_name = "ISB001",
                  cols = color)

pattern_plot_full(file="../data/Sardinia2017/LON001.A0201_S0_L008_R1_001.fastq.fq.qF.sorted.cleaned_rmdup.q30.csv",
                  pattern = c("C->T", "C->G", "C->A", "T->G", "T->A", "T->C",
                              "G->A", "G->C", "G->T", "A_>G", "A->T", "A->C"),
                  plot_type="left",
                  sample_name = "LON001",
                  cols = color)

pattern_plot_full(file="../data/Sardinia2017/SUA001.A0201_S0_L008_R1_001.fastq.fq.qF.sorted.cleaned_rmdup.q30.csv",
                  pattern = c("C->T", "C->G", "C->A", "T->G", "T->A", "T->C",
                              "G->A", "G->C", "G->T", "A_>G", "A->T", "A->C"),
                  plot_type="left",
                  sample_name = "SUA001",
                  cols = color)

pattern_plot_full(file="../data/Sardinia2017/SUA002.A0201_S0_L008_R1_001.fastq.fq.qF.sorted.cleaned_rmdup.q30.csv",
                  pattern = c("C->T", "C->G", "C->A", "T->G", "T->A", "T->C",
                              "G->A", "G->C", "G->T", "A_>G", "A->T", "A->C"),
                  plot_type="left",
                  sample_name = "SUA002",
                  cols = color)

pattern_plot_full(file="../data/Sardinia2017/SUA003.A0201_S0_L008_R1_001.fastq.fq.qF.sorted.cleaned_rmdup.q30.csv",
                  pattern = c("C->T", "C->G", "C->A", "T->G", "T->A", "T->C",
                              "G->A", "G->C", "G->T", "A_>G", "A->T", "A->C"),
                  plot_type="left",
                  sample_name = "SUA002",
                  cols = color)

############################   C to T comaparisons  ################################################

signature_set <- colnames(clubbed_counts_pooled)
sig_split <- t(sapply(1:length(signature_set), function(x) return(strsplit(signature_set[x], "")[[1]][1:8])))
substitution_pattern <- sapply(1:length(signature_set), function(x) return(paste(sig_split[x,3:6], collapse="")))

clubbed_counts_C_to_T <- clubbed_counts_pooled[, grep("C->T", substitution_pattern)];

labs <- c(rep("Lindo Ancient",25), rep("Lindo Modern",25), rep("Sardinia Ancient",5))
indices1 <- grep("Modern", labs)
labs1 <- labs[-indices1];
clubbed_counts_C_to_T_1 <- clubbed_counts_C_to_T[-indices1,]


topics_clus_3 <- maptpx::topics(clubbed_counts_C_to_T_1, K=2, type="full", tol=10)

omega <- topics_clus_3$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs1)
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

signature_set <- colnames(clubbed_counts_C_to_T_1)
apply(CountClust::ExtractTopFeatures(topics_clus_3$theta, top_features = 10, method="poisson", options="min"), c(1,2), function(x) signature_set[x])


####################   More signature patterns exploration  ##########################################

topics_clus_1 <- maptpx::topics(clubbed_counts_pooled, K=5, type="full", tol=10)

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

signature_set <- colnames(clubbed_counts_pooled)
apply(CountClust::ExtractTopFeatures(topics_clus_1$theta, top_features = 10, method="poisson", options="min"), c(1,2), function(x) signature_set[x])

theta <- topics_clus_1$theta
mutation_pattern <- substring(as.character(rownames(theta)), 1, 8)

library(dplyr)
new_theta <- tbl_df(theta) %>% mutate(mutation_pattern) %>% group_by(mutation_pattern) %>% summarise_each(funs(sum)) %>% as.data.frame()
rownames(new_theta) <- new_theta[,1]
new_theta <- new_theta[,-1]

par(mfrow=c(2,1))
damageLogo(new_theta)

position <- as.numeric(as.matrix(sapply(rownames(theta), function(x) return (substring(as.character(x), 10, nchar(as.character(x)))))))

new_theta_pos <- tbl_df(theta) %>% mutate(position) %>% group_by(position) %>% summarise_each(funs(sum)) %>% as.data.frame()

plot_graph(new_theta_pos[,1])

par(mfrow=c(1,2))
l1 <- pattern_plot("../data/Sardinia2017/ISB001.A0201_S0_L008_R1_001.fastq.fq.qF.sorted.cleaned_rmdup.q30.csv", pattern="T->A", plot_type = "left", use_prop = TRUE)
l2 <- pattern_plot("../data/AnnaGosling2016data/ADR-T1-M241.dup.q30.csv", pattern="T->A", plot_type = "left", use_prop = TRUE)
l3 <- pattern_plot("../data/Lindo2016ancients/300_all_chr.q30.csv", pattern="T->A", plot_type = "left", use_prop = TRUE)
l4 <- pattern_plot("../data/Lindo2016moderns/S002_all_chr.q30.csv", pattern="T->A", plot_type = "left", use_prop = TRUE)
l5 <- pattern_plot("../data/AnnaGosling2016data/ADR-T2-PCRneg.dup.q30.csv", pattern="T->A", plot_type = "left", use_prop = TRUE)

plot(l1, col="red", type="l")
lines(l2, col="green")
lines(l3, col="blue")
lines(l4, col="black")
lines(l5, col="orange")

legend("topright", legend=c("sardinia", "gosling ancient", "lindo ancient", "lindo modern", "control"),
       fill=c("red", "green", "blue", "black", "orange"))


l2 <- pattern_plot("../data/Lindo2016ancients/507_all_chr.q30.csv", pattern="C->T", plot_type = "left", use_prop = TRUE)
l3 <- pattern_plot("../data/Lindo2016moderns/S002_all_chr.q30.csv", pattern="C->T", plot_type = "left", use_prop = TRUE)
