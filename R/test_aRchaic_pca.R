


source("../../aRchaic.site/R/aggregate_signature_counts.R")
source("../../aRchaic.site/R/club_signature_counts.R")
source("../../aRchaic.site/R/aRchaic_view.R")
source("../../aRchaic.site/R/filter_signatures_only_location.R")
source("../../aRchaic.site/R/aRchaic_cluster.R")
source("../../aRchaic.site/R/aRchaic_pool.R")
source("../../aRchaic.site/R/damageLogo_5.R")
source("../../CountClust/R/StructureGGplot.R")


###################   Lindo moderns analysis  ##############################

lindo_metadata <- read.table("../utilities/Lindo_metadata.txt")
folders <- c("../data/Lindo2016moderns/")
labs <- lindo_metadata[,2]

out <- aRchaic_pca(folders = folders,
                   labs = labs,
                   run_from = "plot",
                   normalize = TRUE,
                   pcs_to_plot = c("PC1", "PC2", "PC3", "PC4", "PC5"),
                   output_dir = "../utilities/Lindo2016moderns/pca/")

################  Multiple folders (aRchaic pca) ##########################

################### ###################################################

folders <- c("../data/Jones_2015/", "../data/moderns_lite/")
labs <- c(rep("Jones", length(list.files("../data/Jones_2015/", pattern = ".csv"))),
          rep("moderns", length(list.files("../data/moderns_lite/", pattern = ".csv"))))

clus_out <- aRchaic_pca(folders = folders,
                        labs = labs,
                        run_from = "plot",
                        normalize = TRUE,
                        pcs_to_plot = c("PC3", "PC4", "PC5"),
                        output_dir = "../utilities/modern_Jones/pca/")

pca_out <- get(load("../utilities/modern_Jones/pca/pca.rda"))
screeplot(pca_out)

##############   Multiple folders (aRchaic pca)  ########################
##########################################################################



