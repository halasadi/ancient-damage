


#############  test independent model in maptpx  ###########################

library(aRchaic)
HGDPmoderns <- get(load("../processed_data/HGDPmoderns-counts-table.rda"))
HGDPmoderns_clubbed <- club_signature_counts(HGDPmoderns)

signature_set <- colnames(HGDPmoderns_clubbed)
sig_split <- t(sapply(1:length(signature_set), function(x) return(strsplit(signature_set[x], "")[[1]][1:8])))
new_sig_split <- matrix(0, dim(sig_split)[1], 5);
new_sig_split[,1] <- sig_split[,1]
new_sig_split[,2] <- sig_split[,2]
new_sig_split[,3] <- sapply(1:length(signature_set), function(x) return(paste(sig_split[x,3:6], collapse="")))
new_sig_split[,4] <- sig_split[,7]
new_sig_split[,5] <- sig_split[,8]

levels(new_sig_split[,1]) <- c("0", "1", "2", "3", "4")

pos <- t(sapply(1:length(signature_set), function(x)
                                          {
                                              y = strsplit(signature_set[x], "")[[1]]
                                              return(paste(y[10:length(y)], collapse=""))
}))



mat <- matrix(0, dim(new_sig_split)[1], dim(new_sig_split)[2])
for(k in 1:dim(new_sig_split)[2]){
  temp <- as.factor(new_sig_split[,k])
  mat[,k] <- as.numeric(as.matrix(plyr::mapvalues(temp, from = levels(temp), to = 0:(length(levels(temp))-1))))
}

pos <- as.numeric(pos)
pos <- pos - min(pos)
pos <- factor(pos, levels = 0:22)

signatures <- mat;
signature_pos <- cbind.data.frame(signatures, pos)

out <- topics(HGDPmoderns_clubbed, K=2, type="full", signatures = signature_pos)

omega <- out$omega
labs <- rep("hgdp", dim(omega)[1])

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs)
)

CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = FALSE,
                            figure_title = paste0("StructurePlot: K=", dim(omega)[2],""),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))






counts <- HGDPmoderns_clubbed
K <- 2
shape=NULL
initopics=NULL
tol=10
bf=FALSE
kill=2
ord=TRUE
verb=1
admix=TRUE
nbundles=1
use_squarem=FALSE
init.adapt=FALSE
type="independent"
signatures=NULL
light=1
method_admix=1
sample_init=TRUE
tmax=10000
wtol=10^(-4)
qn=100
grp=NULL
nonzero=FALSE
dcut=-10
top_genes=100
burn_in=5
signatures=signature_pos



X <- CheckCounts(counts)
p <- ncol(X)

#############in tpx run these #############

alpha <- shape
theta <- initopics


