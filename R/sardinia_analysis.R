

#############   applying aRchaic on Sardinian data  ####################

data1 <- get(load("../processed_data/sardinia2017-strand-flank.rda"))
clubbed_counts_1 <- club_signature_counts_2(data1, flanking_bases = 1)
data2 <- get(load("../processed_data/1000g-strand-flank-2.rda"))
clubbed_counts_2 <- club_signature_counts_2(data2, flanking_bases = 1)

pooled_names <- intersect(colnames(clubbed_counts_1), colnames(clubbed_counts_2))
filtered_sardinia <- clubbed_counts_1[, match(pooled_names, colnames(clubbed_counts_1))]
filtered_moderns <- clubbed_counts_2[, match(pooled_names, colnames(clubbed_counts_2))]

pooled_data <- rbind(filtered_sardinia, filtered_moderns)



labs <- c(substring(rownames(clubbed_counts_1), 1, 2), rep("moderns", 50))
library(maptpx)

topic_clus <- maptpx::topics(clubbed_counts, K=3, tol=100, type = "full")
save(topic_clus, file="../processed_data/maptpx-runs/sardinia-26-full-maptpx.rda")

omega <- topic_clus$omega

cols1 <- c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
           "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
           "brown4","darkorchid","magenta","yellow", "azure1","azure4")

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs)
)

StructureGGplot(omega = omega,
                annotation = annotation,
                palette = cols1,
                yaxis_label = "Sardinia",
                order_sample = FALSE,
                figure_title = paste0("StructurePlot: K=", dim(omega)[2],""),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))


damageLogo_three(topic_clus$theta, flanking_bases = 1)

indices <- CountClust::ExtractTopFeatures(topic_clus$theta, top_features = 10, method="poisson", options="min")

apply(indices, c(1,2), function(x) rownames(topic_clus$theta)[x])




signature_set <- colnames(pooled_data)
sig_split <- t(sapply(1:length(signature_set), function(x) return(strsplit(signature_set[x], "")[[1]][1:12])))
new_sig_split <- matrix(0, dim(sig_split)[1], 3);
new_sig_split[,1] <- sig_split[,1]
new_sig_split[,2] <- sapply(1:length(signature_set), function(x) return(paste(sig_split[x,2:5], collapse="")))
new_sig_split[,3] <- sig_split[,6]

breakbase <- sig_split[,10]
strand <- sig_split[,8]

new_sig_split <- cbind.data.frame(new_sig_split, breakbase, strand)

pos <- t(sapply(1:length(signature_set), function(x)
{
  y = strsplit(signature_set[x], "")[[1]]
  return(paste(y[12:length(y)], collapse=""))
}))

pos <- as.numeric(pos)
pos <- pos - min(pos)
pos <- factor(pos, levels = 0:20)

mat <- matrix(0, dim(new_sig_split)[1], dim(new_sig_split)[2])
for(k in 1:dim(new_sig_split)[2]){
  temp <- as.factor(new_sig_split[,k])
  mat[,k] <- as.numeric(as.matrix(plyr::mapvalues(temp, from = levels(temp), to = 0:(length(levels(temp))-1))))
}


signatures <- mat;
signature_pos <- cbind.data.frame(signatures, t(pos))

topic_clus <- maptpx::topics(pooled_data, K=2, tol=100, model="independent", signatures = signature_pos)

save(topic_clus, file="../processed_data/maptpx-runs/sardinia-26-moderns-50-ind-maptpx-K2.rda")

omega <- topic_clus$omega
annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs)
)

cols1 <- c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
           "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
           "brown4","darkorchid","magenta","yellow", "azure1","azure4")


StructureGGplot(omega = omega,
                annotation = annotation,
                palette = cols1,
                yaxis_label = "Sardinia",
                order_sample = FALSE,
                figure_title = paste0("StructurePlot: K=", dim(omega)[2],""),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))



damageLogo_five(topic_clus$theta, flanking_bases = 1)



plot_logo(breaks_theta_vec = breaks_theta[,1, drop=FALSE])




theta_pool <- topic_clus$theta

sig_names = NULL
ic.scale=TRUE
max_pos = 20
flanking_bases=1
yscale_change = TRUE
xaxis=TRUE
yaxis=TRUE
xaxis_fontsize=5
xlab_fontsize=10
y_fontsize=10
mut_width=2
start=0.0001
renyi_alpha = 1
pop_names=paste0("Cluster ",1:dim(theta_pool)[2])
logoport_x = 0.24
logoport_y= 0.50
logoport_width= 0.30
logoport_height= 0.40
lineport_x = 0.70
lineport_y=0.50
lineport_width=0.20
lineport_height=0.25
stackbarport_x = 1.00
stackbarport_y = 5
stackbarport_width=0.30
stackbarport_height=0.7
barport_x = 0.58
barport_y=0.65
barport_width=0.25
barport_height=0.25







test  <- data.frame(person=c("A", "B", "C", "D", "E"),
                    value1=c(100,150,120,80,150),
                    value2=c(25,30,45,30,30) ,
                    value3=c(100,120,150,150,200))

library(reshape2) # for melt

melted <- melt(test, "person")

melted$cat <- ''
melted[melted$variable == 'value1',]$cat <- "first"
melted[melted$variable != 'value1',]$cat <- "second"

ggplot(melted, aes(x = cat, y = value, fill = variable)) +
  geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ person)

A = c(0.8, 0.3)
G = c(0.2, 0.1)
C = c(0.4, 0.4)
T = c(0.2, 0.6)
type=c(rep("purine",4), rep("pyrimidine",4))

test  <- data.frame(person=c("5' strand break", "3' strand break"),
                    A=A,  G=G,  C=C, T=T)

melted <- melt(test, "person")
melted$type <- type

ggplot(melted, aes(x = type, y = value, fill = variable)) +
  geom_bar(stat = 'identity', position = 'stack') +
  facet_wrap(~ person, nrow=2) +
  labs(x = "") + labs(y="") + labs(fill="") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        panel.border = element_rect(colour = "white")) +
  theme(plot.margin = unit(c(0.01,0.01,0.01,0.01), "cm"))




breaks_theta_vec <- breaks_theta[,2, drop=FALSE]
mat <- matrix(0, nrow=4, ncol=2)
mat[c(1,3),1] <- breaks_theta_vec[c(1,3),1]
mat[c(2,4),2] <- breaks_theta_vec[c(2,4),1]
rownames(mat) <- c("A", "C", "G", "T")
colnames(mat) <- c("purine", "pyrimidine")

Logolas::logomaker(mat,
                   hist=TRUE,
                   frame_width = 1,
                   cols= c("blue", "red"),
          ic.scale = TRUE,
          yscale_change = TRUE,
          pop_name = "5' strand break",
          xlab = "Papers",
          ylab = "",
          yaxis=FALSE,
          cols_per_column = TRUE,
          col_line_split="black")

mat <- matrix(0, nrow=4, ncol=2)
mat[c(1,3),1] <- breaks_theta_vec[c(5,7),1]
mat[c(2,4),2] <- breaks_theta_vec[c(6,8),1]
rownames(mat) <- c("A", "C", "G", "T")
colnames(mat) <- c("purine", "pyrimidine")

logomaker(mat,
          hist=TRUE,
          frame_width = 1,
          cols= c("blue", "red"),
          ic.scale = TRUE,
          yscale_change = TRUE,
          pop_name = "3' strand break",
          xlab = "Papers",
          ylab = "",
          yaxis=FALSE,
          cols_per_column = TRUE,
          col_line_split="black")


ic_computer(mat, hist=TRUE)

