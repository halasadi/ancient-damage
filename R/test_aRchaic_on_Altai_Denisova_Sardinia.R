

##############   test aRchaic on Altai + Denisova + Sardinia ###################


source("../../aRchaic.site/R/aggregate_signature_counts.R")
source("../../aRchaic.site/R/club_signature_counts.R")
source("../../aRchaic.site/R/aRchaic_view.R")
source("../../aRchaic.site/R/filter_signatures_only_location.R")
source("../../aRchaic.site/R/aRchaic_cluster.R")
source("../../aRchaic.site/R/aRchaic_pool.R")
source("../../aRchaic.site/R/damageLogo_5.R")


folders <- c("../data/Denisova/", "../data/Altai/", "../data/Sardinia2017/")

labs <- c("Denisova", "Altai", c("ISB", "LON", rep("MA", 21), rep("SUA", 3)))

clus_out <- aRchaic_cluster(folders = folders,
                            K = 4,
                            labs = labs,
                            tol = 0.5,
                            run_from = "gom",
                            output_dir = "../utilities/denisova_altai_sardinia/clus_4/",
                            save_plot = TRUE)



folders <- "../data/Sardinia2017/"
labs <- c("ISB", "LON", rep("MA", 21), rep("SUA", 3))

clus_out <- aRchaic_cluster(folders = folders,
                            K = 2,
                            labs = labs,
                            tol = 0.5,
                            run_from = "gom",
                            output_dir = "../utilities/sardinia/clus_2/",
                            save_plot = TRUE)

folders <- c("../data/Denisova/", "../data/Altai/", "../data/Sardinia2017/", "../data/moderns_lite/")

labs <- c("Denisova", "Altai", c("ISB", "LON", rep("MA", 21), rep("SUA", 3)), rep("moderns", 50))

clus_out <- aRchaic_cluster(folders = folders,
                            K = 2,
                            labs = labs,
                            tol = 0.5,
                            run_from = "gom",
                            output_dir = "../utilities/denisova_altai_sardinia_moderns/clus_2/",
                            save_plot = TRUE)



folders <- c("../data/wolves/", "../data/moderns_lite/")

labs <- c(rep("RKW", 2), rep("RSRW", 5), rep("moderns", 50))
levels <- c("RKW", "RSRW", "moderns")

clus_out <- aRchaic_cluster(folders = folders,
                            K = 3,
                            labs = labs,
                            tol = 0.5,
                            run_from = "gom",
                            output_dir = "../utilities/wolves_moderns/clus_3/",
                            save_plot = TRUE)



#######################    Gosling  data   ################################

names <- list.files("../data/AnnaGosling2016data/", pattern = ".csv")
control_indices <- c(grep("EXN", names), grep("Libneg", names), grep("PCRneg", names))
labs <- character();
labs <- rep("ancient", length(names))
labs[control_indices] <- "controls"

indices <- which(labs == "ancient")

ancient_names = names[indices]
pop_names_1 <- as.factor(substring(ancient_names, 8, 8))
levels(pop_names_1)

levels(pop_names_1) = c("Chokhopani", "Kyang", "Rhirhi", "Mebrak", "Samdzong")
labs[indices] <- as.character(pop_names_1)

levels <- levels(pop_names_1)

folders <- c("../data/AnnaGosling2016data/")

clus_out <- aRchaic_cluster(folders = folders,
                            K = 2,
                            labs = labs,
                            levels = levels,
                            tol = 1,
                            run_from = "gom",
                            output_dir = "../utilities/Gosling2016/clus_2/",
                            save_plot = TRUE)


#############   Gosling + moderns data  ############################

names <- list.files("../data/AnnaGosling2016data/", pattern = ".csv")
control_indices <- c(grep("EXN", names), grep("Libneg", names), grep("PCRneg", names))
labs <- character();
labs <- rep("ancient", length(names))
labs[control_indices] <- "controls"

indices <- which(labs == "ancient")

ancient_names = names[indices]
pop_names_1 <- as.factor(substring(ancient_names, 8, 8))
levels(pop_names_1)

levels(pop_names_1) = c("Chokhopani", "Kyang", "Rhirhi", "Mebrak", "Samdzong")
labs[indices] <- as.character(pop_names_1)

labs2 <- c(labs, rep("moderns", 50))

levels <- c(levels(pop_names_1), "moderns")

folders <- c("../data/AnnaGosling2016data/", "../data/moderns_lite/")

clus_out <- aRchaic_cluster(folders = folders,
                            K = 3,
                            labs = labs2,
                            levels = levels,
                            tol = 1,
                            run_from = "plot",
                            output_dir = "../utilities/Gosling_moderns/clus_3/",
                            save_plot = TRUE)


#################    Reich  + Allentoft  + moderns  ########################

folders <- c("../data/Reich/", "../data/Allentoft/", "../data/moderns_lite/")
labs <- c(rep("Reich", length(list.files("../data/Reich/", pattern = ".csv"))),
          rep("Allentoft", length(list.files("../data/Allentoft/", pattern = ".csv"))),
          rep("moderns", length(list.files("../data/moderns_lite/", pattern = ".csv"))))

levels = NULL

clus_out <- aRchaic_cluster(folders = folders,
                            K = 2,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "gom",
                            output_dir = "../utilities/Reich_Allentoft_moderns/clus_2/",
                            save_plot = TRUE)


topic_clus <- get(load(paste0("../utilities/Reich_Allentoft_moderns/clus_4/model.rda")))

I_gom <- topic_clus$omega[1:110,]
meta_table <- read.csv("../utilities/Reich_Allentoft_moderns/meta_table.csv")
coverage <- meta_table$Coverage

I_id <- as.vector(sapply(rownames(I_gom), function(x) return(strsplit(x, "[.]")[[1]][1])))
I_full_id <- as.character(meta_table$Unique.ID)
meta_table_reorder <- meta_table[match(I_id, I_full_id),]

plot(I_gom[,2], meta_table_reorder$Min.Date, col="red",
     xlab= "Blue cluster gom", ylab = "Min Date",
     main = "Cluster memberships wrt metadata", pch=20,cex=1)


