folders <- c("../data/Lazaridis/", "../data/moderns_lite/")
files <- list.files(folders[1], pattern = ".csv")
labs_skoglund <- substring(paste0(files),1,3)
labs <- c(rep("Lazaridis (UDG)", 9), rep("moderns", 50))
levels <- unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 3,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            logo.control = list(use_log = TRUE),
                            output_dir = "../utilities/modern_Pinhasi_Lazaridis/clus_3/")

