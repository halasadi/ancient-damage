

###############  Test aRchaic on Nepal samples  ############################

source("../../aRchaic.site/R/aggregate_signature_counts.R")
source("../../aRchaic.site/R/club_signature_counts.R")
source("../../aRchaic.site/R/aRchaic_view.R")
source("../../aRchaic.site/R/filter_signatures_only_location.R")
source("../../aRchaic.site/R/aRchaic_cluster.R")
source("../../aRchaic.site/R/aRchaic_pool.R")
source("../../aRchaic.site/R/damageLogo_5.R")
source("../../CountClust/R/StructureGGplot.R")


folders <- c("../data/Nepal/", "../data/moderns_lite/", "../data/Sardinia2017/")
files <- list.files(folders[1], pattern = ".csv")
labs_sard <- c("ISB", "LON", rep("MA", 21), rep("SUA", 3))
labs <- c(paste0(substring(paste0(files), 1, 3), "_", substring(paste0(files), 8, 12)),
          rep("moderns", 50), labs_sard)
levels <- unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 3,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            positive_logo = FALSE,
                            output_dir = "../utilities/Nepal_sardinia_moderns/clus_3/")


#############   Lazaridis + Moderns  ####################################

folders <- c("../data/moderns_lite/", "../data/Lazaridis/")

labs <- c(rep("moderns", 50), c("LBK", "Loschbour", rep("Motala", 7)))
levels <- unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 3,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            output_dir = "../utilities/Lazaridis_moderns/clus_3/",
                            save_plot = TRUE)


##########   Lazaridis + Nepal + Sardinia + Moderns #######################

folders <- c("../data/moderns_lite/", "../data/Lazaridis/", "../data/Nepal/", "../data/Sardinia2017/")
labs_moderns <- rep("moderns", 50)
labs_laz <- c("LBK", "Loschbour", rep("Motala", 7))
labs_nepal <- paste0(substring(paste0(files), 1, 3), "_", substring(paste0(files), 8, 12))
labs_sard <- c("ISB", "LON", rep("MA", 21), rep("SUA", 3))
labs <- c(labs_moderns, labs_laz, labs_nepal, labs_sard)
levels <- unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 5,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            output_dir = "../utilities/Lazaridis_nepal_sardinia_moderns/clus_5/",
                            save_plot = TRUE)


############   contaminated  + moderns + Allentoft  ###########################

folders <-c("../data/contaminated/", "../data/moderns_lite/", "../data/Allentoft/")
labs <- c(rep("contaminated", 6), rep("moderns", 50), rep("RISE", 42))
levels <- unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 4,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            output_dir = "../utilities/contaminated_moderns_Allentoft/clus_4/",
                            save_plot = TRUE)


############   Lazaridis  + moderns + Allentoft  ###########################

folders <-c("../data/Lazaridis/", "../data/moderns_lite/", "../data/Allentoft/")
labs_laz <- c("LBK", "Loschbour", rep("Motala", 7))
labs <- c(labs_laz, rep("moderns", 50), rep("RISE", 42))
levels <- unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 3,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "gom",
                            output_dir = "../utilities/Lazaridis_moderns_Allentoft/clus_3/",
                            save_plot = TRUE)

#############   Lazaridis + Iceman + Skoglund + moderns  #######################

folders <-c("../data/Iceman/", "../data/Skoglund/","../data/moderns_lite/", "../data/Lazaridis/")
labs_laz <- c("LBK", "Loschbour", rep("Motala", 7))
labs <- c(rep("Iceman", 3), "Skoglund", rep("moderns", 50), labs_laz)
levels <- unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 2,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "gom",
                            output_dir = "../utilities/Iceman_Skoglund_moderns_Lazaridis/clus_2/",
                            save_plot = TRUE)


##############   moderns + Sherpa  ##########################

folders <- c("../data/moderns_lite/", "../data/Sherpa_data/")
labs <- c(rep("moderns", 50), rep("sherpa", 5))
levels <- unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 2,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            output_dir = "../utilities/moderns_Sherpa/clus_2/",
                            save_plot = TRUE)



############   Fu  +  moderns   ################################

folders <- c("../data/moderns_lite/", "../data/Fu_2016/")
labs <- c(rep("moderns", 50), rep("Fu", 44))
levels <- unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 5,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            output_dir = "../utilities/moderns_Fu/clus_5/",
                            save_plot = TRUE)

#############  allentoft and moderns  #######################


folders <- c("../data/Allentoft/", "../data/moderns_lite/")
labs <- c(rep("Allentoft", length(list.files("../data/Allentoft/", pattern = ".csv"))),
          rep("moderns", length(list.files("../data/moderns_lite/", pattern = ".csv"))))

levels = unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 3,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            output_dir = "../utilities/moderns_Allentoft/clus_3/",
                            save_plot = TRUE)

topic_clus <- get(load("../utilities/moderns_Allentoft/clus_2/model.rda"))

library(readxl)
metadata_allentoft <- read.csv("../utilities/moderns_Allentoft/nature_paper.csv")
mt.contamination <- metadata_allentoft$MT.contamination.....

contam_estimates <- as.numeric(sapply(mt.contamination, function(x) return(strsplit(as.character(x), "[(]")[[1]][1])))


names <- as.character(rownames(topic_clus$omega))
names_mod <- sapply(names, function(x) return(strsplit(as.character(x), "[.]")[[1]][1]))
names_mod <- as.vector(names_mod)[1:42]

matched_indices <- match(names_mod, as.character(metadata_allentoft[,1]))

metadata_allentoft_mod <- metadata_allentoft[matched_indices,]
mt.contamination.mod <- contam_estimates[matched_indices]

plot(topic_clus$omega[1:42,1], log(mt.contamination.mod), cex=1, pch=20)



#############  Jones and moderns  #######################


folders <- c("../data/Jones_2015/", "../data/moderns_lite/")
labs <- c(rep("Jones", length(list.files("../data/Jones_2015/", pattern = ".csv"))),
          rep("moderns", length(list.files("../data/moderns_lite/", pattern = ".csv"))))

levels = unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 2,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            output_dir = "../utilities/modern_Jones/clus_2/",
                            save_plot = TRUE)


#############  Skoglund and moderns  #######################

folders <- c("../data/Skoglund/", "../data/moderns_lite/")
files <- list.files(folders[1], pattern = ".csv")

labs_skoglund <- substring(paste0(files),1,3)
labs <- c(labs_skoglund, rep("moderns", 50))
levels <- unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 4,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            output_dir = "../utilities/modern_Skoglund/clus_4/",
                            save_plot = TRUE)


##################  Skoglund, Pinhasi and moderns  #######################

folders <- c("../data/Skoglund/", "../data/Pinhasi/", "../data/moderns_lite/")
files <- list.files(folders[1], pattern = ".csv")
labs_skoglund <- substring(paste0(files),1,3)
labs <- c(labs_skoglund, rep("Pinhasi", 12), rep("moderns", 50))
levels <- unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 3,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            output_dir = "../utilities/modern_Skoglund_Pinhasi/clus_3/",
                            save_plot = TRUE)


##################   moderns + Fu (X chromosome)   ##########################

folders <- c("../data/X_moderns_lite/", "../data/X_Fu_2016/")
files <- list.files(folders[2], pattern = ".csv")
labs <- c(rep("X_moderns", 50), rep("X_Fu", length(files)))
levels <- unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 4,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            output_dir = "../utilities/X_modern_Fu/clus_4/",
                            save_plot = TRUE)

##################   moderns + Fu (MT-DNA chromosome)   ##########################

folders <- c("../data/MT_moderns_lite/", "../data/MT_Fu_2016/")
files <- list.files(folders[2], pattern = ".csv")
labs <- c(rep("MT_moderns", 50), rep("MT_Fu", length(files)))
levels <- unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 4,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            output_dir = "../utilities/MT_modern_Fu/clus_4/",
                            save_plot = TRUE)

##################   moderns + Fu (nuclear DNA)   ##########################

folders <- c("../data/moderns_lite/", "../data/Fu_2016/")
files <- list.files(folders[2], pattern = ".csv")
labs <- c(rep("moderns", 50), rep("Fu", length(files)))
levels <- unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 3,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            output_dir = "../utilities/moderns_Fu/clus_3/")

##################  moderns + Raghavan  ############################

folders <- c("../data/moderns_lite/", "../data/Raghavan/")
files <- list.files(folders[2], pattern = ".csv")
labs_raghavan <- substring(files, 1, 2)
labs <- c(rep("moderns", 50), labs_raghavan)
levels <- unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 3,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            output_dir = "../utilities/moderns_Raghavan/clus_3/",
                            save_plot = TRUE)


##################  Skoglund, Lazaridis and moderns  #######################

folders <- c("../data/Pinhasi/", "../data/Lazaridis/", "../data/moderns_lite/")
files <- list.files(folders[1], pattern = ".csv")
labs_skoglund <- substring(paste0(files),1,3)
labs <- c(rep("Pinhasi (non UDG)", 12),
          rep("Lazaridis (UDG)", 9), rep("moderns", 50))
levels <- unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 3,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "gom",
                            output_dir = "../utilities/modern_Pinhasi_Lazaridis/clus_3/")


files <- list.files()


######################   Lindo 2016  + moderns  ############################

folders <- c("../data/Lindo2016/")
labs <- c(rep("ancient", 25),
          rep("moderns", 25))
levels <- unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 4,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            output_dir = "../utilities/Lindo2016/clus_4/",
                            save_plot = TRUE)



dat <- get(load("../data/Lindo2016/Lindo2016.rda"))
dat <- dat[-25, ]
labs <- labs[-25]

signature_set <- colnames(dat)
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

topic_clus <- maptpx::topics(dat, K=2, type="independent", tol=50,
                             signatures = signature_pos)


omega <- topic_clus$omega
#
levels <- unique(labs)
annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs, levels = levels)
)

rownames(omega) <- annotation$sample_id
output_dir <- "../utilities/Lindo2016_wo_939/clus_2/"

topic_cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
               "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
               "brown4","darkorchid","magenta","yellow", "azure1","azure4")
structure_width = 5
structure_height = 8
structure.control <- list()
K=2
structure.control.default <- list(yaxis_label = "aRchaic pops",
                                  order_sample = FALSE,
                                  figure_title = paste0("  StructurePlot: K=", K,""),
                                  axis_tick = list(axis_ticks_length = .1,
                                                   axis_ticks_lwd_y = .1,
                                                   axis_ticks_lwd_x = .1,
                                                   axis_label_size = 10,
                                                   axis_label_face = "bold"),
                                  legend_title_size = 10,
                                  legend_key_size = 0.7,
                                  legend_text_size = 8)
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

damageLogo_five(topic_clus$theta, save_plot=TRUE,
                inflation_factor = c(2,1,2), output_dir = output_dir)
graphics.off()

###################   Lindo moderns analysis  ##############################

lindo_metadata <- read.table("../utilities/Lindo_metadata.txt")
folders <- c("../data/Lindo2016moderns/")
labs <- lindo_metadata[,2]
levels <- unique(labs)



clus_out <- aRchaic_cluster(folders = folders,
                            K = 2,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "gom",
                            output_dir = "../utilities/Lindo2016moderns/clus_2/")

###############   contaminated  samples  ##################################

