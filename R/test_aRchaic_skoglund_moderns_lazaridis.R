

###########  test aRchaic negative on Skoglund + moderns + Lazaridis  ##############

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
                            run_from = "plot",
                            logo.control = list(use_log = TRUE),
                            output_dir = "../utilities/modern_Pinhasi_Lazaridis/clus_3/")


topic_clus <- get(load("../utilities/modern_Pinhasi_Lazaridis/clus_3/model.rda"))
theta <- topic_clus$theta

theta_pool <- theta
sig_names = NULL
ic.scale=TRUE
use_log = FALSE
max_pos = 20
flanking_bases=1
yscale_change = TRUE
xaxis=TRUE
yaxis=TRUE
xlab = " "
xaxis_fontsize=10
xlab_fontsize=20
title_aligner = 15
y_fontsize=20
title_fontsize = 20
mut_width=2
start=0.0001
renyi_alpha = 1
inflation_factor = c(2,1,2)
pop_names=paste0("Cluster ",1:dim(theta_pool)[2])
logoport_x = 0.25
logoport_y= 0.50
logoport_width= 0.28
logoport_height= 0.40
lineport_x = 0.9
lineport_y=0.40
lineport_width=0.32
lineport_height=0.28
breaklogoport_x = 1.00
breaklogoport_y = 0.40
breaklogoport_width=0.40
breaklogoport_height=0.50
barport_x = 0.58
barport_y=0.60
barport_width=0.25
barport_height=0.40
output_dir = NULL
output_width = 1200
output_height = 700



pwm <- prop_patterns_list

alpha=1
inflation_factor = c(1,1,1)
base_probs = NULL



folders <- c("../data/Pinhasi/", "../data/Lazaridis/", "../data/moderns_lite/")
files <- list.files(folders[1], pattern = ".csv")
labs_skoglund <- substring(paste0(files),1,3)
labs <- c(rep("Pinhasi (non UDG)", 12),
          rep("Lazaridis (UDG)", 9), rep("moderns", 50))
levels <- unique(labs)


files <- list.files("../data/moderns_lite_lite/", pattern = ".csv")
bases_freq_mixcomp_left_mat <- c()
bases_freq_mixcomp_right_mat <- c()
mutations_freq_mixcomp_mat <- c()

for(f in 1:length(files)){
  file <- read.csv(paste0("../data/moderns_lite_lite/", files[f]),header=FALSE)
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



clus_out <- aRchaic_cluster(folders = folders,
                            K = 3,
                            labs = labs,
                            levels = levels,
                            tol = 10,
                            run_from = "plot",
                            logo.control = list(renyi_alpha = 2,
                                                use_log = FALSE,
                                                base_probs_list = base_probs_temp),
                            output_dir = "../utilities/modern_Pinhasi_Lazaridis/clus_3/")




out <- damage.ic(pwm, base_probs_list = base_probs_temp)
