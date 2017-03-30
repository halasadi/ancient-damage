

#############  clipped first two bases (moderns_lite, nepal)  ###########################

end_breaks <- 2

data1 <- get(load("../data/moderns_lite/moderns_lite.rda"))
data2 <- get(load("../data/Nepal/Nepal.rda"))

pos <- sapply(colnames(data1), function(x) return(strsplit(x, "_")[[1]][4]))

data1_filtered <- data1[, which(as.numeric(pos) > end_breaks)]
data2_filtered <- data2[, which(as.numeric(pos) > end_breaks)]

data_pool_filtered <- rbind(data1_filtered, data2_filtered)

split_profile <- do.call(rbind, lapply(colnames(data_pool_filtered), function(x) return(strsplit(x, "_")[[1]])))
split_profile[,3] = "A"
paste_profile <- apply(split_profile, 1, function(x) return(paste0(x, collapse = "_")))
colnames(data_pool_filtered) <- paste_profile


split_profile_T <- do.call(rbind, lapply(unique(colnames(data_pool_filtered)), function(x) return(strsplit(x, "_")[[1]])))
split_profile_T[,3] = "T"

split_profile_A <- do.call(rbind, lapply(unique(colnames(data_pool_filtered)), function(x) return(strsplit(x, "_")[[1]])))
split_profile_A[,3] = "A"

split_profile_C <- do.call(rbind, lapply(unique(colnames(data_pool_filtered)), function(x) return(strsplit(x, "_")[[1]])))
split_profile_C[,3] = "C"

split_profile_G <- do.call(rbind, lapply(unique(colnames(data_pool_filtered)), function(x) return(strsplit(x, "_")[[1]])))
split_profile_G[,3] = "G"

paste_profile_T <- apply(split_profile_T, 1, function(x) return(paste0(x, collapse = "_")))
paste_profile_A <- apply(split_profile_A, 1, function(x) return(paste0(x, collapse = "_")))
paste_profile_G <- apply(split_profile_G, 1, function(x) return(paste0(x, collapse = "_")))
paste_profile_C <- apply(split_profile_C, 1, function(x) return(paste0(x, collapse = "_")))

paste_profile_new <- c(rbind(paste_profile_T, paste_profile_A, paste_profile_G, paste_profile_C))

new_data <- matrix(0, dim(data_pool_filtered)[1], dim(data_pool_filtered)[2])

for(m in 1:dim(data_pool_filtered)[1]){
  group_data <- tapply(data_pool_filtered[m,], colnames(data_pool_filtered), sum)
  group_data_extended <- rep(round(group_data/4), each=4)
  unique_names <- unique(colnames(data_pool_filtered))
  new_data[m, ] <- group_data_extended
}

colnames(new_data) <- paste_profile_new

signature_set <- colnames(new_data)
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
pos <- factor(pos, levels = 0:21)

signatures <- mat;
signature_pos <- cbind.data.frame(signatures, pos)

topic_clus <- maptpx::topics(new_data, K=3, type="independent", tol=10,
                             signatures = signature_pos)


files <- list.files("../data/Nepal/", pattern = ".csv")
labs <- c(rep("moderns", 50), paste0(substring(paste0(files), 1, 3), "_", substring(paste0(files), 8, 12)))
levels <- unique(labs)

omega <- topic_clus$omega
annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs, levels = levels)
)


output_dir <- "../utilities/clipped_moderns_Nepal_clip_2/"

topic_cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
               "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
               "brown4","darkorchid","magenta","yellow", "azure1","azure4")

structure.control = list()
logo.control = list()
topics.control = list()

structure_width = 5
structure_height = 8

plot.new()
grid.newpage()
do.call(StructureGGplot, append(list(omega= omega,
                                     annotation = annotation,
                                     palette = topic_cols,
                                     order_sample=FALSE),
                                structure.control))

ggplot2::ggsave(paste0(output_dir, "structure.png"), width = structure_width,
                height = structure_height)

plot.new()
damageLogo_five(topic_clus$theta, save_plot=TRUE, output_dir = output_dir)
graphics.off()

###################  Lets take a look at a second example #############



##########   clipped first two bases (moderns_lite, sardinia)  ##########################


end_breaks <- 10

data1 <- get(load("../data/moderns_lite/moderns_lite.rda"))
data2 <- get(load("../data/Sardinia2017/Sardinia2017.rda"))

pos <- sapply(colnames(data1), function(x) return(strsplit(x, "_")[[1]][4]))

data1_filtered <- data1[, which(as.numeric(pos) > end_breaks)]
data2_filtered <- data2[, which(as.numeric(pos) > end_breaks)]

data_pool_filtered <- rbind(data1_filtered, data2_filtered)


signature_set <- colnames(data_pool_filtered)
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
pos <- factor(pos, levels = 0:21)

signatures <- mat;
signature_pos <- cbind.data.frame(signatures, pos)

topic_clus <- maptpx::topics(data_pool_filtered, K=3, type="independent", tol=100,
                             signatures = signature_pos)


labs_sard <- c("ISB", "LON", rep("MA", 21), rep("SUA", 3))
labs <- c(rep("moderns", 50), labs_sard)
levels <- unique(labs)

omega <- topic_clus$omega
annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs, levels = levels)
)


output_dir <- "../utilities/clipped_moderns_Sardinia_clip_2/"


plot.new()
grid.newpage()
do.call(StructureGGplot, append(list(omega= omega,
                                     annotation = annotation,
                                     palette = topic_cols),
                                structure.control))

ggplot2::ggsave(paste0(output_dir, "structure.png"), width = structure_width,
                height = structure_height)

plot.new()
damageLogo_five(topic_clus$theta, save_plot=TRUE, output_dir = output_dir)
graphics.off()



#######################################################################

theta_pool <- topic_clus$theta
sig_names = NULL
ic.scale=TRUE
max_pos = 20
flanking_bases=1
yscale_change = TRUE
xaxis=TRUE
yaxis=TRUE
xlab = " "
xaxis_fontsize=5
xlab_fontsize=10
title_aligner = 18
y_fontsize=10
title_fontsize = 20
mut_width=2
start=0.0001
renyi_alpha = 1
pop_names=paste0("Cluster ",1:dim(theta_pool)[2])
logoport_x = 0.25
logoport_y= 0.50
logoport_width= 0.28
logoport_height= 0.40
lineport_x = 0.9
lineport_y=0.40
lineport_width=0.32
lineport_height=0.28
breaklogoport_x = 1.00
breaklogoport_y = 0.40
breaklogoport_width=0.40
breaklogoport_height=0.50
barport_x = 0.58
barport_y=0.60
barport_width=0.25
barport_height=0.25
output_dir = NULL
output_width = 1200
output_height = 700
save_plot=TRUE
