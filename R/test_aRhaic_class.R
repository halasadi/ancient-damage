


########################   test  aRchaic_clus()   ###################################


out <- aRchaic_pool(folders = "../data/Pinhasi/")
aRchaic_pool(folders = "../data/moderns_lite/")
aRchaic_pool(folders = "../data/Sardinia2017/")

folders <- c("../data/Pinhasi/", "../data/moderns_lite/", "../data/Gosling-controls/")

data1 <- get(load("../data/Pinhasi/Pinhasi.rda"))
data2 <- get(load("../data/moderns_lite/moderns_lite.rda"))

class_labs <- c(rep(1, dim(data1)[1]),
                rep(2, dim(data2)[1]),
                rep(NA, length(list.files(folders[3], pattern = ".csv"))))

tab <- aRchaic_class(folders,
                     class_labs = class_labs,
                     class_method = "SVM",
                     run_from = "class",
                     normalize = FALSE,
                     classtpx.control = list())

labs <- c(rep("Pinhasi", dim(data1)[1]),
          rep("moderns", dim(data2)[1]),
          rep("controls", length(list.files(folders[3], pattern = ".csv"))))
levels <- unique(labs)

clus_out <- aRchaic_cluster(folders = folders,
                            K = 2,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            output_dir = "../utilities/moderns_Pinhasi_controls/clus_2/")


data1 <- get(load("../data/Pinhasi/Pinhasi.rda"))
data2 <- get(load("../data/moderns_lite/moderns_lite.rda"))
data3 <- get(load("../data/Gosling-controls/Gosling-controls.rda"))

datalist <- list()
datalist[[1]] <- data1
datalist[[2]] <- data2
datalist[[3]] <- data3

sig_names <- colnames(datalist[[1]])
row_names_pool <- rownames(datalist[[1]])
if(length(datalist) >= 2){
  for(num in 2:length(datalist)){
    sig_names <- union(sig_names, colnames(datalist[[num]]))
    row_names_pool <- c(row_names_pool, rownames(datalist[[num]]))
  }
}

pooled_data <- matrix(0, length(row_names_pool), length(sig_names))
rownames(pooled_data) <- row_names_pool
colnames(pooled_data) <- sig_names

for(num in 1:length(datalist)){
  pooled_data[match(rownames(datalist[[num]]), rownames(pooled_data)), match(colnames(datalist[[num]]), sig_names)] <- datalist[[num]]
}

tab <- aRchaic_class_beta(pooled_data,
                          type = "mutation-flank",
                          class_labs = class_labs,
                          class_method = "SVM",
                          normalize = FALSE,
                          classtpx.control = list())


library(gridBase)
library(gridExtra)
library(ggplot2)
labs <- c(rep("Pinhasi", dim(data1)[1]),
          rep("moderns", dim(data2)[1]),
          rep("controls", length(list.files(folders[3], pattern = ".csv"))))
aRchaic_pca(folders = folders,
            labs = labs,
            run_from = "plot",
            normalize = TRUE,
            pcs_to_plot = c("PC1", "PC2", "PC3", "PC4", "PC5"),
            output_dir = "../utilities/moderns_Pinhasi_controls/")

aRchaic_tsne_beta(pooled_data,
                  labs = labs,
                  type = "mutation-flank",
                  normalize = TRUE,
                  dims = 3,
                  output_dir = "../utilities/moderns_Pinhasi_controls/")



