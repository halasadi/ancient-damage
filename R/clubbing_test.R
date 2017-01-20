

###  clubbing the signature counts ################

lindo_ancients <- get(load("../processed_data/lindo2016ancients-counts-table-strand-flank.rda"))
lindo_moderns <- get(load("../processed_data/lindo2016moderns-counts-table-strand-flank.rda"))

club_signature_validation_plot(lindo_moderns)
lindo_moderns_club <- club_signature_counts(lindo_moderns, flanking_bases = 2)
lindo_ancients_club <- club_signature_counts(lindo_ancients, flanking_bases = 2)

save(lindo_moderns_club,
     file="../processed_data/lindo2016moderns-counts-table-strand-flank-clubbed.rda")
save(lindo_ancients_club,
     file="../processed_data/lindo2016ancients-counts-table-strand-flank-clubbed.rda")



pooled_names <- intersect(colnames(lindo_ancients_club),
                          colnames(lindo_moderns_club))

indices1 <- match(pooled_names, colnames(lindo_ancients_club))
indices2 <- match(pooled_names, colnames(lindo_moderns_club))

pooled_dat <- rbind(lindo_ancients_club[, indices1],
                    lindo_moderns_club[, indices2])

library(maptpx)
topic_clus <- topics(pooled_dat, K=2, type = "full", tol=100)

labs <- c(rep("Ancient",25), rep("Modern",25))

omega <- topic_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs)
)

rownames(omega) <- annotation$sample_id;
plot.new()
CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Development Phase",
                            order_sample = FALSE,
                            figure_title = paste0("StructurePlot: K=", dim(omega)[2],", no C -> T / G -> A"),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))


signature_set <- rownames(topic_clus$theta)

indices_left <- grep("left", signature_set)
indices_right <- grep("right", signature_set)

breaks_theta <- topic_clus$theta[c(indices_left, indices_right),]

theta_new <- topic_clus$theta[-c(indices_left, indices_right),]
theta_new<- apply(theta_new, 2, function(x) return(x/sum(x)))

indices_minus <- grep("_-_", signature_set)
strand_theta <- data.frame("minus" = colSums(theta_new[indices_minus,]),
                           "plus" = colSums(theta_new[-indices_minus,]))
strand_theta_vec <- data.frame(strand_theta[1,])

signature_set_new <- rownames(theta_new)
signature_set_split <- do.call(rbind, lapply(signature_set_new, function(x) {
                                y = strsplit(as.character(x), split="")[[1]][-c((4+2*flanking_bases + 1):(4+2*flanking_bases + 2))]
                                return(paste0(y, collapse=""))
                                }))
tapply(theta_new[,1], signature_set_split, sum)

library(dplyr)
theta <- tbl_df(theta_new) %>% mutate(sig = signature_set_split) %>% group_by(sig) %>% summarise_each(funs(sum)) %>% as.data.frame()
rownames(theta) <-  theta[,1]
theta <- theta[,-1]


damageLogo_pos_str(topic_clus$theta)

sig_names = NULL
ic.scale=TRUE
max_pos = 20
flanking_bases=2
yscale_change = TRUE
xaxis=TRUE
yaxis=TRUE
xaxis_fontsize=5
xlab_fontsize=10
y_fontsize=10
mut_width=2
start=0.0001
renyi_alpha = 1
pop_names=paste0("Cluster ",1:dim(theta)[2])
logoport_x = 0.24
logoport_y= 0.50
logoport_width= 0.30
logoport_height= 0.40
lineport_x = 0.70
lineport_y=0.50
lineport_width=0.20
lineport_height=0.25
pieport_x = 1.18
pieport_y = 5
pieport_width=0.70
pieport_height=0.6
barport_x = 0.58
barport_y=0.65
barport_width=0.25
barport_height=0.25
