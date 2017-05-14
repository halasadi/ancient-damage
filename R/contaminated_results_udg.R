


####################    contaminated moderns + ancients  ###########################


ancient <- get(load("../data/Pinhasi/Pinhasi.rda"))
moderns <- get(load("../data/moderns_lite/moderns_lite.rda"))

rowSums(ancient)
rowSums(moderns)

samp_1 <- rmultinom(1, 10000, ancient[1,]) + rmultinom(1, 20000, moderns[1,])
samp_2 <- rmultinom(1, 10000, ancient[2,]) + rmultinom(1, 10000, moderns[20,])
samp_3 <- rmultinom(1, 10000, ancient[3,]) + rmultinom(1, 5000, moderns[30,])
samp_4 <- rmultinom(1, 10000, ancient[4,]) + rmultinom(1, 50000, moderns[40,])
samp_5 <- rmultinom(1, 5000, ancient[5,]) + rmultinom(1, 500, moderns[30,])
samp_6 <- rmultinom(1, 500, ancient[2,]) + rmultinom(1, 700, moderns[20,])
samp_7 <- rmultinom(1, 1000, ancient[3,]) + rmultinom(1, 5000, moderns[40,])
samp_8 <- rmultinom(1, 10000, ancient[5,]) + rmultinom(1, 100, moderns[20,])
samp_9 <- rmultinom(1, 200, ancient[3,]) + rmultinom(1, 1700, moderns[50,])
samp_10 <- rmultinom(1, 20, ancient[3,]) + rmultinom(1, 10000, moderns[20,])

true_prop <- c(1/3, 1/2, 10/15, 1/6, 5000/5500, 5/12, 1/6, 100/101, 2/19,  2/1000)


contam_samp <- t(cbind(samp_1, samp_2, samp_3, samp_4, samp_5, samp_6, samp_7, samp_8, samp_9, samp_10))


pooled_data <- rbind(contam_samp, ancient, moderns)

labs <- c(rep("contaminated", 10), rep("ancient", dim(ancient)[1]),
          rep("moderns", 50))
levels <- unique(labs)

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
pos <- factor(pos, levels = 0:21)

signatures <- mat;
signature_pos <- cbind.data.frame(signatures, pos)

topic_clus <- maptpx::topics(pooled_data, K=2, type="independent", tol=10,
                             signatures = signature_pos)



labs <- c(rep("contaminated", 5), rep("lazaridis", dim(ancient)[1]),
          rep("moderns", 50))
levels <- unique(labs)

omega <- topic_clus$omega
annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs, levels = levels)
)

rownames(omega) <- annotation$sample_id
output_dir <- "../utilities/contaminated_modern_lazaridis/"

topic_cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
               "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
               "brown4","darkorchid","magenta","yellow", "azure1","azure4")
structure_width = 5
structure_height = 8
structure.control <- list()
structure.control.default <- list(yaxis_label = "aRchaic pops",
                                  order_sample = FALSE,
                                  figure_title = paste0("  StructurePlot: K=", K,""),
                                  axis_tick = list(axis_ticks_length = .1,
                                                   axis_ticks_lwd_y = .1,
                                                   axis_ticks_lwd_x = .1,
                                                   axis_label_size = 7,
                                                   axis_label_face = "bold"),
                                  legend_title_size = 8,
                                  legend_key_size = 0.4,
                                  legend_text_size = 5)
structure.control <- modifyList(structure.control, structure.control.default)


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

