

##############  Strand breaks  ########################

files <- list.files("../data/AnnaGosling2016data/")
str <- list()
for(m in files){
  str[[m]] <- strand_breaks_composition(paste0("../data/AnnaGosling2016data/", m))
  cat("at iteration", m, "\n")
}

save(str, file="../processed_data/strand-breaks-gosling2016.rda")

files <- list.files("../data/Sardinia2017/")
str <- list()
for(m in 1:length(files)){
  str[[m]] <- strand_breaks_composition(paste0("../data/Sardinia2017/", files[m]))
  cat("at iteration", m, "\n")
}

save(str, file="../processed_data/strand-breaks-sardinia2017.rda")

files <- list.files("../data/HGDPmoderns/")
str <- list()
for(m in 1:length(files)){
  str[[m]] <- strand_breaks_composition(paste0("../data/HGDPmoderns/", files[m]))
  cat("at iteration", m, "\n")
}

save(str, file="../processed_data/strand-breaks-hgdp.rda")


files <- list.files("../data/Lindo2016moderns/")
str <- list()
for(m in 1:length(files)){
  str[[m]] <- strand_breaks_composition(paste0("../data/Lindo2016moderns/", files[m]))
  cat("at iteration", m, "\n")
}

save(str, file="../processed_data/strand-breaks-lindomoderns.rda")


files <- list.files("../data/Lindo2016ancients/")
str <- list()
for(m in 1:length(files)){
  str[[m]] <- strand_breaks_composition(paste0("../data/Lindo2016ancients/", files[m]))
  cat("at iteration", m, "\n")
}

save(str, file="../processed_data/strand-breaks-lindoancients.rda")


files <- list.files("../data/1000Gmoderns/")
str <- list()
for(m in 1:length(files)){
  str[[m]] <- strand_breaks_composition(paste0("../data/1000Gmoderns/", files[m]))
  cat("at iteration", m, "\n")
}

save(str, file="../processed_data/strand-breaks-1000g.rda")



files <- list.files("../csv/Sherpa_data/")
str <- list()
for(m in 1:length(files)){
  str[[m]] <- strand_breaks_composition(paste0("../csv/Sherpa_data/", files[m]))
  cat("at iteration", m, "\n")
}

save(str, file="../processed_data/strand-breaks-sherpa.rda")



#############################   Strand breaks analysis  ##############################################


left_strand_breaks <- as.numeric()
right_strand_breaks <- as.numeric()
labs <- character();

strand_breaks_thousandG <- get(load("../processed_data/strand-breaks-1000g.rda"))
for(l in 1:length(strand_breaks_thousandG)){
  left_strand_breaks <- rbind(left_strand_breaks, strand_breaks_thousandG[[l]][[1]])
  right_strand_breaks <- rbind(right_strand_breaks, strand_breaks_thousandG[[l]][[2]])
}
labs <- c(labs, rep("1000G", length(strand_breaks_thousandG)))

strand_breaks_lindomoderns <- get(load("../processed_data/strand-breaks-lindomoderns.rda"))
for(l in 1:length(strand_breaks_lindomoderns)){
  left_strand_breaks <- rbind(left_strand_breaks, strand_breaks_lindomoderns[[l]][[1]])
  right_strand_breaks <- rbind(right_strand_breaks, strand_breaks_lindomoderns[[l]][[2]])
}

labs <- c(labs, rep("Lindo-moderns", length(strand_breaks_lindomoderns)))


strand_breaks_hgdp <- get(load("../processed_data/strand-breaks-hgdp.rda"))
for(l in 1:length(strand_breaks_hgdp)){
  left_strand_breaks <- rbind(left_strand_breaks, strand_breaks_hgdp[[l]][[1]])
  right_strand_breaks <- rbind(right_strand_breaks, strand_breaks_hgdp[[l]][[2]])
}

labs <- c(labs, rep("HGDP", length(strand_breaks_hgdp)))


strand_breaks_lindoancients <- get(load("../processed_data/strand-breaks-lindoancients.rda"))
for(l in 1:length(strand_breaks_lindoancients)){
  left_strand_breaks <- rbind(left_strand_breaks, strand_breaks_lindoancients[[l]][[1]])
  right_strand_breaks <- rbind(right_strand_breaks, strand_breaks_lindoancients[[l]][[2]])
}

labs <- c(labs, rep("Lindo-ancient", length(strand_breaks_lindoancients)))


strand_breaks_sardinia <- get(load("../processed_data/strand-breaks-sardinia2017.rda"))
for(l in 1:length(strand_breaks_sardinia)){
  left_strand_breaks <- rbind(left_strand_breaks, strand_breaks_sardinia[[l]][[1]])
  right_strand_breaks <- rbind(right_strand_breaks, strand_breaks_sardinia[[l]][[2]])
}

labs <- c(labs, rep("Sardinia", length(strand_breaks_sardinia)))



strand_breaks_sherpa <- get(load("../processed_data/strand-breaks-sherpa.rda"))
for(l in 1:length(strand_breaks_sherpa)){
  left_strand_breaks <- rbind(left_strand_breaks, strand_breaks_sherpa[[l]][[1]])
  right_strand_breaks <- rbind(right_strand_breaks, strand_breaks_sherpa[[l]][[2]])
}

labs <- c(labs, rep("Sherpa", length(strand_breaks_sherpa)))


strand_breaks_gosling <- get(load("../processed_data/strand-breaks-gosling2016.rda"))
for(l in 1:length(strand_breaks_gosling)){
  left_strand_breaks <- rbind(left_strand_breaks, strand_breaks_gosling[[l]][[1]])
  right_strand_breaks <- rbind(right_strand_breaks, strand_breaks_gosling[[l]][[2]])
}


signature_counts <- get(load("../processed_data/annagosling2016-counts-table.rda"))
names <- rownames(signature_counts);
control_indices <- c(grep("EXN", names), grep("Libneg", names), grep("PCRneg", names))

labs1 <- character();
labs1 <- rep("Gossling-ancient", dim(signature_counts)[1])
labs1[control_indices] <- "Gossling-controls"

labs <- c(labs, labs1)

library(aRchaic)

aRchaic::gridPCA_signatures(left_strand_breaks, factor(labs), normalize = TRUE)
aRchaic::gridPCA_signatures(right_strand_breaks, factor(labs), normalize = TRUE)

all_strand_breaks <- cbind(left_strand_breaks, right_strand_breaks)
aRchaic::gridPCA_signatures(all_strand_breaks, factor(labs), normalize = TRUE)

omega_left <- t(apply(left_strand_breaks, 1, function(x) return(x/sum(x))));
#omega_left <- left_strand_breaks

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega_left))),
  tissue_label = factor(labs)
)

rownames(omega_left) <- annotation$sample_id;

CountClust::StructureGGplot(omega = omega_left,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = FALSE,
                            figure_title = paste0("Composition of nucleotides at left strand breaks"),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))


omega_right <- t(apply(right_strand_breaks, 1, function(x) return(x/sum(x))));
#omega_left <- left_strand_breaks

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega_right))),
  tissue_label = factor(labs)
)

rownames(omega_right) <- annotation$sample_id;

CountClust::StructureGGplot(omega = omega_right,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = FALSE,
                            figure_title = paste0("Composition of nucleotides at right strand breaks"),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))

library(dplyr)
omega_left_filt <- tbl_df(as.matrix(left_strand_breaks)) %>% mutate(labs) %>% group_by(labs) %>% summarise_each(funs(sum)) %>% as.data.frame()

mat <- as.matrix(omega_left_filt[,-1])
mat <- t(mat)
colnames(mat) <- omega_left_filt[,1]

library(Logolas)
logomaker(mat, cols= RColorBrewer::brewer.pal(dim(mat)[1],name = "Spectral"), ic.scale = FALSE)


omega_right_filt <- tbl_df(as.matrix(right_strand_breaks)) %>% mutate(labs) %>% group_by(labs) %>% summarise_each(funs(sum)) %>% as.data.frame()

mat <- as.matrix(omega_right_filt[,-1])
mat <- t(mat)
colnames(mat) <- omega_right_filt[,1]

library(Logolas)
logomaker(mat, cols= RColorBrewer::brewer.pal(dim(mat)[1],name = "Spectral"), ic.scale = FALSE)
