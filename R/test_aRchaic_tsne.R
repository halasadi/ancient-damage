

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


mat1 <- get(load("../data/Fu_2016/Fu_2016.rda"))
mat2 <- get(load("../data/Pinhasi/Pinhasi.rda"))
mat3 <- get(load("../data/Lazaridis/Lazaridis.rda"))
mat4 <- get(load("../data/moderns_lite/moderns_lite.rda"))

pooled_mat <- rbind(mat1, mat2, mat3, mat4)

labs <- c(rep("Fu", dim(mat1)[1]), rep("Pinhasi", dim(mat2)[1]),
          rep("Lazaridis", dim(mat3)[1]), rep("moderns", dim(mat4)[1]))

clus_out <- aRchaic_tsne_beta(mat = pooled_mat,
                              labs = labs,
                              type = "wo-strand",
                              normalize = TRUE,
                              dims = 3,
                              output_dir = "../utilities/pca_filter/",
                              output_name = "tsne-wo-strand")


