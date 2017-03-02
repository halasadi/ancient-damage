

dir <- "../data/Skoglund/";

###  aggregate the counts of signatures
###  per read position bin (defined by the breaks)
###


out <- aggregate_signature_counts(dir,
                                  breaks = c(-1, seq(1,20,1)),
                                  flanking_bases = 1,
                                  pattern = "Ajv59.hs37d5.fa.merged.q30.csv")


###  Read length distribution for a single file



par(mfrow=c(1,1))
par(mar=c(2,2,2,2))
read_length_distribution(dir,
                         pattern = "Ajv59.hs37d5.fa.merged.q30.csv",
                         end_break = 5,
                         plot_layout = c(1,1),
                         cex.main = 0.8,
                         cex_legend = 0.7)


file <- "../data/Skoglund/Ajv59.hs37d5.fa.merged.q30.csv"



############  Mutation trends from the ends  ######################

color=c("red","blue","cornflowerblue","black","cyan","darkblue",
        "brown4","burlywood","darkgoldenrod1","darkgray","deepskyblue","darkkhaki",
        "firebrick","darkorchid","hotpink","green","magenta","yellow", "azure1","azure4");

par(mfrow=c(1,2))
pattern_plot_full(file="../data/Skoglund/Ajv59.hs37d5.fa.merged.q30.csv",
                  pattern = c("C->T", "C->G", "C->A", "T->G", "T->A", "T->C",
                              "G->A", "G->C", "G->T", "A_>G", "A->T", "A->C"),
                  plot_type="left",
                  sample_name = "Skoglund",
                  cols = color)

pattern_plot_full(file="../data/Skoglund/Ajv59.hs37d5.fa.merged.q30.csv",
                  pattern = c("C->T", "C->G", "C->A", "T->G", "T->A", "T->C",
                              "G->A", "G->C", "G->T", "A_>G", "A->T", "A->C"),
                  plot_type="right",
                  sample_name = "Skoglund",
                  cols = color)


par(mfrow=c(1,2))
pattern_plot_full(file="../data/Skoglund/Ajv59.hs37d5.fa.merged.q30.csv",
                  pattern = c("C->A", "T->G", "T->A", "T->C"),
                  plot_type="right",
                  sample_name = "Skoglund",
                  cols = color)

pattern_plot_full(file="../data/Skoglund/Ajv59.hs37d5.fa.merged.q30.csv",
                  pattern = c("C->A", "T->G", "T->A", "T->C"),
                  plot_type="left",
                  sample_name = "Skoglund",
                  cols = color)


##########  aRchaic logo plot  ##############################


clubbed_data <- club_signature_counts(out, flanking_bases = 1)
clubbed_data_normed <- clubbed_data/sum(clubbed_data)
damageLogo_five(t(clubbed_data_normed), renyi_alpha=100, xaxis_fontsize=12,
                y_fontsize=12, pop_names = "Skoglund aDNA")




dir <- "../data/Sardinia_data_strand_flank/";
out <- aggregate_signature_counts(dir, breaks = c(-1, seq(1,20,1)), flanking_bases = 1)
save(out, file = "../processed_data/Sardinian_data_strand_flank.rda")
out <- get(load("../processed_data/Sardinian_data_strand_flank.rda"))

input_file <- "../data/Sardinia_data_strand_flank/ISB001snp_sorted_deduped_chrALL.q30.csv"


read_length_distribution(dir, end_break = 5, plot_layout = c(3,3), cex.main = 0.8)
clubbed_data <- club_signature_counts(out, flanking_bases = 1)
clubbed_data_normed <- apply(clubbed_data, 1, function(x) return(x/sum(x)))

library(grid)
library(gridBase)
library(Logolas)
library(ggplot2)

damageLogo_five(clubbed_data_normed[,c(1,2, 5, 25)], renyi_alpha=100,
                xaxis_fontsize=12,
                y_fontsize=12, pop_names = c("ISB001", "LON001", "MA108", "SUA002"))


dir <- "../data/1000Gmoderns_data_strand_flank_2/";
out <- aggregate_signature_counts(dir, breaks = c(-1, seq(1,20,1)), flanking_bases = 1)

out <- get(load("../processed_data/1000G_data_strand_flank.rda"))
clubbed_data <- club_signature_counts(out, flanking_bases = 1)
clubbed_data_normed <- apply(clubbed_data, 1, function(x) return(x/sum(x)))

save(out, file = "../processed_data/1000G_data_strand_flank.rda")


damageLogo_five(clubbed_data_normed[,c(20, 30)], renyi_alpha=100, xaxis_fontsize=12,
                y_fontsize=12, pop_names = c("Modern 1", "Modern 2"))



###########  combining data from Sardinia and 1000G and topic model  #######


data1 <- get(load("../processed_data/Sardinian_data_strand_flank.rda"))
data2 <- get(load("../processed_data/1000G_data_strand_flank.rda"))

clubbed_data_1 <- club_signature_counts(data1, flanking_bases = 1)
clubbed_data_2 <- club_signature_counts(data2, flanking_bases = 1)

pooled_names <- intersect(colnames(clubbed_data_1), colnames(clubbed_data_2))

filtered_data_1 <- clubbed_data_1[,match(pooled_names, colnames(clubbed_data_1))]
filtered_data_2 <- clubbed_data_2[,match(pooled_names, colnames(clubbed_data_2))]

pooled_data <- rbind(filtered_data_1, filtered_data_2)



signature_set <- colnames(pooled_data)
sig_split <- t(sapply(1:length(signature_set), function(x) return(strsplit(signature_set[x], "")[[1]][1:8])))
new_sig_split <- matrix(0, dim(sig_split)[1], 3);
new_sig_split[,1] <- sig_split[,1]
new_sig_split[,2] <- sapply(1:length(signature_set), function(x) return(paste(sig_split[x,2:5], collapse="")))
new_sig_split[,3] <- sig_split[,6]

levels(new_sig_split[,1]) <- c("0", "1", "2", "3", "4")

pos <- t(sapply(1:length(signature_set), function(x)
{
  y = strsplit(signature_set[x], "")[[1]]
  return(paste(y[12:length(y)], collapse=""))
}))



mat <- matrix(0, dim(new_sig_split)[1], dim(new_sig_split)[2])
for(k in 1:dim(new_sig_split)[2]){
  temp <- as.factor(new_sig_split[,k])
  mat[,k] <- as.numeric(as.matrix(plyr::mapvalues(temp, from = levels(temp), to = 0:(length(levels(temp))-1))))
}

pos <- as.numeric(pos)
pos <- pos - min(pos)
pos <- factor(pos, levels = 0:22)

signatures <- mat;
signature_pos <- cbind.data.frame(signatures, pos)

topic_clus <- maptpx::topics(pooled_data, K=3, tol=100, model="independent", signatures = signature_pos)

library(grid)
library(gridBase)
library(ggplot2)
library(Logolas)
library(CountClust)

damageLogo_five(topic_clus$theta, renyi_alpha=100, xaxis_fontsize=15,
                y_fontsize=15, pop_names=paste0("Cluster ",1:dim(topic_clus$theta)[2]))

omega <- topic_clus$omega
labs <- c("ISB", "LON", rep("MA", 21), rep("SUA", 3), rep("moderns", 50))
annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs, levels = rev(c("moderns", "MA", "SUA", "ISB", "LON")))
)

cols1 <- c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
           "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
           "brown4","darkorchid","magenta","yellow", "azure1","azure4")


par(mfrow=c(1,1))
StructureGGplot(omega = omega,
                annotation = annotation,
                palette = cols1,
                yaxis_label = "Moderns vs Ancients",
                order_sample = FALSE,
                figure_title = paste0("  StructurePlot: K=", dim(omega)[2],""),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"),
                legend_title_size = 12,
                legend_text_size = 9)
