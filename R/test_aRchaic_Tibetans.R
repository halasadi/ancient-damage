

###########  test aRchaic on Tibetans  #######################

folders <- c("../data/Tibetans/", "../data/moderns_lite/")
files <- list.files(folders[1], pattern = ".csv")
labs1 <- as.character(sapply(files, function(x) return (substring(as.character(x),1,1))))
labs <- factor(c(labs1, rep("moderns", 50)), levels = c("K", "R", "S", "M", "moderns"))
levels <- unique(labs)

modfiles <- list.files("../data/moderns_lite_lite/", pattern = ".csv")
bases_freq_mixcomp_left_mat <- c()
bases_freq_mixcomp_right_mat <- c()
mutations_freq_mixcomp_mat <- c()

for(f in 1:length(modfiles)){
  file <- read.csv(paste0("../data/moderns_lite_lite/", modfiles[f]),header=FALSE)
  muts <- file[,1]
  leftflankcomp <- substring(muts, 1, 1)
  bases_freq <- tapply(file[,7], leftflankcomp, sum)
  bases_freq <- bases_freq[match(c("A", "C", "G", "T"), names(bases_freq))]
  bases_freq_mixcomp_left <- bases_freq/sum(bases_freq)
  bases_freq_mixcomp_left_mat <- rbind(bases_freq_mixcomp_left_mat, bases_freq_mixcomp_left)

  rightflankcomp <- substring(muts, 6, 6)
  bases_freq <- tapply(file[,7], rightflankcomp, sum)
  bases_freq <- bases_freq[match(c("A", "C", "G", "T"), names(bases_freq))]
  bases_freq_mixcomp_right <- bases_freq/sum(bases_freq)
  bases_freq_mixcomp_right_mat <- rbind(bases_freq_mixcomp_right_mat, bases_freq_mixcomp_right)


  mutations <- substring(muts, 2, 5)
  mutations_freq <- tapply(file[,7], mutations, sum)
  mutations_freq <- mutations_freq[match(c("C->T", "C->G", "C->A", "T->C", "T->A", "T->G"), names(mutations_freq))]
  mutations_freq_mixcomp <- mutations_freq/sum(mutations_freq)
  mutations_freq_mixcomp_mat <- rbind(mutations_freq_mixcomp_mat, mutations_freq_mixcomp)
}

bases_freq_mixcomp_left_mean <- colMeans(bases_freq_mixcomp_left_mat)
bases_freq_mixcomp_right_mean <- colMeans(bases_freq_mixcomp_right_mat)
mutations_freq_mixcomp_mean <- colMeans(mutations_freq_mixcomp_mat)
base_probs_temp <- list()
base_probs_temp[[1]] <- bases_freq_mixcomp_left_mean
base_probs_temp[[3]] <- bases_freq_mixcomp_right_mean
base_probs_temp[[2]] <- mutations_freq_mixcomp_mean
#names(base_probs_temp[[1]]) <- c("A", "C", "G", "T")
#names(base_probs_temp[[3]]) <- c("A", "C", "G", "T")
#names(base_probs_temp[[2]]) <- c("C->T", "C->A", "C->G", "T->A", "T->C", "T->G")

save(base_probs_temp, file = "../rda/base_probs_temp.rda")

base_probs_temp <- get(load("../rda/base_probs_temp.rda"))
levels <- c("R", "M", "K", "S", "moderns")
clus_out <- aRchaic_cluster(folders = folders,
                            K = 2,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            logo.control = list(renyi_alpha = 2,
                                                use_log = FALSE,
                                                base_probs_list = base_probs_temp),
                            output_dir = "../utilities/modern_Tibetans/clus_2/")


