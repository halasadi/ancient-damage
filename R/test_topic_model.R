


################   Checking Admixture Performance  #######################

omega_modern <- c(0.8, 0.1, 0.1);
omega_ancient_1 <- c(0.6, 0.3, 0.1);
omega_ancient_2 <- c(0.7, 0.1, 0.2);

#omega_modern <- c(0.8,0.2)
#omega_ancient_1 <- c(0.6,0.4);

theta1 <- rep(1/6,6); ## equal prob to 6 types of mutations
theta2 <- c(1,rep(0,5))
theta3 <- c(rep(0,5),1)

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

omega_modern_set <- rep.row(omega_modern,25);
omega_ancient_1_set <- rep.row(omega_ancient_1, 20);
omega_ancient_2_set <- rep.row(omega_ancient_2, 5);

omega_set <- rbind(omega_modern_set,
                   omega_ancient_1_set,
                   omega_ancient_2_set)

theta_set <- cbind(theta1, theta2, theta3)

p_set <- omega_set%*%t(theta_set)

counts <- matrix(0, dim(p_set)[1], dim(p_set)[2])

for(n in 1:dim(counts)[1]){
  counts[n,] <- rmultinom(1, 100, prob=p_set[n,])
}


topics_clus <- maptpx::topics(counts, K=3, tol=0.001)

w <- topics_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(w))),
  tissue_label = factor(c(rep("Modern",25), rep("Ancient",25)))
)

rownames(w) <- annotation$sample_id;

CountClust::StructureGGplot(omega = w,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = TRUE,
                            figure_title = "StructurePlot",
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))

library(maptpx)
library(slam)

source("../../classtpx/R/class_topics.R")
source("../../classtpx/R/class_tpx.R")
source("../../classtpx/R/class_count.R")
source("../../classtpx/R/class_varselect.R")


known_samples <- 1:25;
class_labs <- rep(1,25);

K <- 3
known_samples=known_samples
class_labs=class_labs
method="theta.fix"
shrink=FALSE
mash_user=NULL
shape=NULL
initopics=NULL
tol=0.1
bf=FALSE
kill=2
ord=TRUE 
verb=1
tmax=10000 
wtol=10^(-4)
qn=100
grp=NULL 
admix=TRUE 
prior_omega = c(0.8,0.1,0.1)
nonzero=FALSE 
dcut=-10

Topic_clus <- class_topics(
  counts, 
  K=3, 
  known_samples = known_samples,
  class_labs = class_labs,
  method="theta.fix",
  shrink=FALSE,
  tol=0.01,
  prior_omega = c(0.8,0.1,0.1),
  ord=FALSE)

omega <- Topic_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(c(rep("Modern",25), rep("Ancient",25)))
)

rownames(omega) <- annotation$sample_id;

CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = TRUE,
                            figure_title = "StructurePlot",
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))

Topic_clus <- class_topics(
  counts, 
  K=3, 
  known_samples = known_samples,
  class_labs = class_labs,
  method="theta.fix",
  optional_theta = rep(1/6,6),
  shrink=FALSE,
  tol=0.01,
  prior_omega = c(0.8,0.1,0.1),
  ord=FALSE)

omega <- Topic_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(c(rep("Modern",25), rep("Ancient",25)))
)

rownames(omega) <- annotation$sample_id;

CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = FALSE,
                            figure_title = "StructurePlot",
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))
