

###############   Generating Fig 1  ######################################

folders <- c("../data/moderns_lite/", "../data/Pinhasi/", "../data/Lazaridis/")
pinhasi_files <- list.files("../data/Pinhasi/", pattern = ".csv")
lazaridis_files <- list.files("../data/Lazaridis/", pattern = ".csv")

labs <- c(rep("moderns", 50),
          rep("Pinhasi (non UDG)", length(pinhasi_files)),
          rep("Lazaridis (UDG)", length(lazaridis_files)))

levels <- unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 5,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "gom",
                            output_dir = "../utilities/figure_1/clus_5/",
                            save_plot = TRUE)


folders <- c("../data/moderns_lite/", "../data/Lazaridis/")
pinhasi_files <- list.files("../data/Pinhasi/", pattern = ".csv")
lazaridis_files <- list.files("../data/Lazaridis/", pattern = ".csv")

labs <- c(rep("moderns", 50),
          rep("Lazaridis (UDG)", length(lazaridis_files)))

levels <- unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 2,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "gom",
                            output_dir = "../utilities/figure_1_udg/clus_2/",
                            save_plot = TRUE)
