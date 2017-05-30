

############################  moderns + UDG + non UDG  ################################

left_base_composition <- get(load("../utilities/left_right_base_composition_ancient.rda"))


mathieson_data <- get(load("../data/Mathieson-2015/Mathieson-2015.rda"))
mathieson_data_sampled <- get(load("../data/Mathieson-2015-subsampled/Mathieson-2015-subsampled.rda"))
names <- as.character(sapply(rownames(mathieson_data), function(x) return(strsplit(x, "[.]")[[1]][1])))

meta_table <- read.csv("../utilities/Mathieson-sub/meta_table.csv")
names_meta <- as.character(meta_table[,1])
meta_table_filtered <- meta_table[match(names, names_meta),]

topic_clus <- get(load("../utilities/Mathieson-sub/clus_2/model.rda"))
omega <- topic_clus$omega

metadata <- meta_table_filtered$Max.Date[-which(meta_table_filtered$Max.Date == "..")]
omega1 <- omega[1:163,]
omega2 <- omega1[-which(meta_table_filtered$Max.Date == ".."),]

plot( omega2[,1], as.numeric(as.character(metadata)), pch=20, cex = 1.2, xlab= "1st cluster GoM",
      ylab = "Min Date", col="red")

plot( omega[1:163,1], meta_table_filtered$Y.derived.SNPs.supporting.haplogroup.determination, pch=20, cex = 1, xlab= "1st cluster GoM",
      ylab = "Min Date")


##########################   Lindo 2016   ####################################

folders <- c("../data/Lindo2016/")
labs <- c(rep("ancient", 25),
          rep("moderns", 25))
levels <- unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 3,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            output_dir = "../utilities/Lindo2016/clus_3/")


lindo_metadata <- read.table("../utilities/Lindo_metadata.txt")
folders <- c("../data/Lindo2016moderns/")
labs <- lindo_metadata[,2]
levels <- unique(labs)



clus_out <- aRchaic_cluster(folders = folders,
                            K = 3,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            output_dir = "../utilities/Lindo2016moderns/clus_3/")

clus_out <- aRchaic_cluster_beta(folders = folders,
                                 type = "mutation-pos",
                                 pattern = "C->T",
                                 K = 2,
                                 labs = labs,
                                 levels = levels,
                                 tol = 0.1,
                                 run_from = "plot",
                                 output_dir = "../utilities/Lindo2016moderns/mutation-pos/clus_2/")



base_composition <- get(load("../utilities/left_right_base_composition_ancient.rda"))

names <- rownames(base_composition$right)
labs <- as.character(sapply(names, function(x) return(strsplit(x, "[/]")[[1]][3])))

indices1 <- which(labs == "Lazaridis")
indices2 <- which(labs == "Reich")
indices3 <- which(labs == "neanderthal")
indices4 <- which(labs == "Pinhasi")
indices5 <- which(labs == "Allentoft")
indices6 <- which(labs == "Sherpa_data")

pooled_indices <- c(indices1, indices2, indices3, indices4, indices5, indices6)
labs <- c(rep("Lazaridis", length(indices1)), rep("Mathieson", length(indices2)),
          rep("neanderthal", length(indices3)), rep("Pinhasi", length(indices4)),
          rep("Allentoft", length(indices5)), rep("Sherpa", length(indices6)))
levels <- c("Lazaridis", "Mathieson", "neanderthal", "Pinhasi", "Allentoft", "Sherpa")

omega <- base_composition$right[pooled_indices,]
omega <- t(apply(omega,1, function(x) return(x/sum(x))))

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs, levels = levels)
)

topic_cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
               "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
               "brown4","darkorchid","magenta","yellow", "azure1","azure4")

rownames(omega) <- annotation$sample_id
CountClust::StructureGGplot(omega = omega,
                annotation = annotation,
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Ancients",
                order_sample = FALSE,
                figure_title = paste0("right flanking base composition"),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 10,
                                 axis_label_face = "bold"))

