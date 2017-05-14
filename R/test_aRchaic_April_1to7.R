

###############  test aRchaic on new Sardinian data  #################################


out <- aRchaic_pool(folders = "../data/Sardinia2017/")



##################   Sardinia vs moderns : aRchaic_cluster()   ##########################

library(grid)
library(gridBase)
library(gridExtra)
library(ggplot2)

source("../../aRchaic.site/R/aggregate_signature_counts.R")
source("../../aRchaic.site/R/aRchaic_cluster.R")
source("../../aRchaic.site/R/filter_signatures_only_location.R")

folders <- c("../data/moderns_lite/", "../data/Sardinia2017/")
labs_moderns <- rep("moderns", 50)
labs_sard <- c("ISB", "LON", rep("MA", 21), rep("SUA", 3))
labs <- c(labs_moderns, labs_sard)
levels <- c("moderns", "MA", "ISB", "LON", "SUA")

clus_out <- aRchaic_cluster(folders = folders,
                            K = 3,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            output_dir = "../utilities/moderns_Sardinia/clus_3/")


###################  Mathieson (subsampled + non subsampled) ###########################

folders <- c("../data/Mathieson-2015/", "../data/Mathieson-2015-subsampled/")
labs <- c(rep("Mathieson", length(list.files(folders[1], pattern = ".csv"))),
          rep("Mathieson-subsampled", length(list.files(folders[1], pattern = ".csv"))))
levels <- unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 3,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            output_dir = "../utilities/Mathieson-sub/clus_3/")



###################  Lazaridis + moderns + Pinhasi    ###########################

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

#################  Lazaridis + Skoglund + moderns  ##########################

folders <-c("../data/Skoglund/","../data/moderns_lite/", "../data/Lazaridis/")
labs_laz <- c("LBK", "Loschbour", rep("Motala", 7))
labs <- c(rep("Skoglund",9), rep("moderns", 50), labs_laz)
levels <- unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 3,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            output_dir = "../utilities/moderns_Lazaridis_Skoglund/clus_3/")


##################   Mathieson (original + subsampled) ###########################

folders <-c("../data/Mathieson-2015/","../data/Mathieson-2015-subsampled/")
labs <- c(rep("Mathieson", length(list.files(folders[1], pattern = ".csv"))),
        rep("Mathieson-sub", length(list.files(folders[2], pattern = ".csv"))))
levels <- unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 3,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            output_dir = "../utilities/Mathieson-sub/clus_3/")


#################   moderns + Allentoft   #############################

folders <-c("../data/Allentoft/","../data/moderns_lite/")
labs <- c(rep("Allentoft", length(list.files(folders[1], pattern = ".csv"))),
          rep("moderns", length(list.files(folders[2], pattern = ".csv"))))
levels <- unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 3,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            output_dir = "../utilities/moderns_Allentoft/clus_3/")


