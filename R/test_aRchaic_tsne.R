

#############  test aRchaic_tsne  ################################


data1 <- get(load("../data/Fu_2016/Fu_2016.rda"))
data2 <- get(load("../data/moderns_lite/moderns_lite.rda"))

pooled_data <- rbind(data1, data2)
pooled_data <- pooled_data[,1:1000]
labs <- c(rep("fu", 44), rep("moderns", 50))
dims = 3

################  Multiple folders (aRchaic pca) ##########################

################### ###################################################

folders <- c("../data/Jones_2015/", "../data/moderns_lite/")
labs <- c(rep("Jones", length(list.files("../data/Jones_2015/", pattern = ".csv"))),
          rep("moderns", length(list.files("../data/moderns_lite/", pattern = ".csv"))))

clus_out <- aRchaic_tsne(folders = folders,
                        labs = labs,
                        run_from = "plot",
                        normalize = TRUE,
                        dims = 5,
                        filter_indices = 1:1000,
                        output_dir = "../utilities/modern_Jones/tsne/")
