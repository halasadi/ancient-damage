


###########  aRchaic_cluster_beta   :  different types  ###################


#############  test aRchaic_cluster_beta_mutation()  #############################


mat1 <- get(load("../data/Fu_2016/Fu_2016.rda"))
mat2 <- get(load("../data/moderns_lite/moderns_lite.rda"))

mat_pooled <- rbind(mat1, mat2)
labs <- c(rep("fu", 44), rep("moderns", 50))

aRchaic_cluster_beta_mutation(mat_pooled,
                              K=3,
                              tol=0.01,
                              labs = labs,
                              output_dir = "../utilities/structure-mutation/moderns_Fu/")


mat1 <- get(load("../data/Pinhasi/Pinhasi.rda"))
mat2 <- get(load("../data/moderns_lite/moderns_lite.rda"))
labs <- c(rep("Pinhasi", 12), rep("moderns", 50))

mat_pooled <- rbind(mat1, mat2)

aRchaic_cluster_beta_mutation(mat_pooled,
                              K=2,
                              tol=0.01,
                              labs = labs,
                              output_dir = "../utilities/structure-mutation/moderns_Pinhasi/")



################  mutation + flank   ###########################

mat1 <- get(load("../data/Fu_2016/Fu_2016.rda"))
mat2 <- get(load("../data/moderns_lite/moderns_lite.rda"))

mat_pooled <- rbind(mat1, mat2)
labs <- c(rep("fu", 44), rep("moderns", 50))

aRchaic_cluster_beta_mutation_flank(mat_pooled,
                                    K=3,
                                    tol=0.01,
                                    labs = labs,
                                    output_dir = "../utilities/structure_mutation_flank/moderns_Fu/")


mat1 <- get(load("../data/Pinhasi/Pinhasi.rda"))
mat2 <- get(load("../data/moderns_lite/moderns_lite.rda"))
labs <- c(rep("Pinhasi", 12), rep("moderns", 50))

mat_pooled <- rbind(mat1, mat2)

aRchaic_cluster_beta_mutation_flank(mat_pooled,
                                    K=2,
                                    tol=0.01,
                                    labs = labs,
                                    output_dir = "../utilities/structure_mutation_flank/moderns_Pinhasi/")


#############  test aRchaic_cluster_beta_pos()  #############################


source('~/Documents/CountClust/R/StructureGGplot.R')
source('~/Documents/aRchaic.site/R/filter_by_pos_pattern.R')

mat1 <- get(load("../data/Fu_2016/Fu_2016.rda"))
mat2 <- get(load("../data/moderns_lite/moderns_lite.rda"))

mat_pooled <- rbind(mat1, mat2)
labs <- c(rep("fu", 44), rep("moderns", 50))

aRchaic_cluster_beta_pos(mat_pooled,
                         pattern = "C->T",
                         K=3,
                         max_pos = 20,
                         tol=0.01,
                         labs = labs,
                         output_dir = "../utilities/structure_pos/moderns_Fu/")


mat1 <- get(load("../data/Pinhasi/Pinhasi.rda"))
mat2 <- get(load("../data/moderns_lite/moderns_lite.rda"))
labs <- c(rep("Pinhasi", 12), rep("moderns", 50))

mat_pooled <- rbind(mat1, mat2)

aRchaic_cluster_beta_pos(mat_pooled,
                         pattern = "C->T",
                         K=3,
                         max_pos = 20,
                         tol=0.01,
                         labs = labs,
                         output_dir = "../utilities/structure_pos/moderns_Pinhasi/")


##########  aRchaic_cluster_beta: mutation + flank + pos ################

source('~/Documents/aRchaic.site/R/filter_by_pos.R')
source('~/Documents/aRchaic.site/R/damageLogo_5.R')
source('~/Documents/aRchaic.site/R/damageLogo_1.R')
source('~/Documents/aRchaic.site/R/aRchaic_cluster_beta_mutation_flank_pos.R')

mat1 <- get(load("../data/Fu_2016/Fu_2016.rda"))
mat2 <- get(load("../data/moderns_lite/moderns_lite.rda"))

mat_pooled <- rbind(mat1, mat2)
labs <- c(rep("fu", 44), rep("moderns", 50))

aRchaic_cluster_beta_mutation_flank_pos(mat_pooled,
                                        K=3,
                                        max_pos = 20,
                                        tol=0.01,
                                        labs = labs,
                                        output_dir = "../utilities/structure_mutation_flank_pos/moderns_Fu/")


mat1 <- get(load("../data/Pinhasi/Pinhasi.rda"))
mat2 <- get(load("../data/moderns_lite/moderns_lite.rda"))
labs <- c(rep("Pinhasi", 12), rep("moderns", 50))

mat_pooled <- rbind(mat1, mat2)

aRchaic_cluster_beta_mutation_flank_pos(mat_pooled,
                                        K=2,
                                        max_pos = 20,
                                        tol=0.1,
                                        labs = labs,
                                        output_dir = "../utilities/structure_mutation_flank_pos/moderns_Pinhasi/")

##########  aRchaic_cluster_beta: mutation + flank + pos + strand  ################


