

###############  Non-negative Matrix Factorization  ############################


signature_counts <- get(load("../summary_data/signature-counts-clubbed-Lindo2016.rda"))

## Apply NMF

install.packages("NMF")

meth <- nmfAlgorithm(version = "R")
meth <- c(names(meth), meth)
meth

library(NMF)
out <- nmf(t(signature_counts), rank=2)

w <- basis(out)
dim(w)

h <- t(coef(out))
dim(h)

hnorm <- t(apply(h, 1, function(x) return(x/sum(x))))

library(CountClust)

omega <- hnorm

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(c(rep("Ancient",25), rep("Modern",25)))
)

StructureGGplot(omega = omega,
                annotation = annotation,
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Moderns/Ancients",
                order_sample = FALSE,
                figure_title = paste0("NMF StructurePlot: K=", dim(omega)[2],": pmsignature: with C->T/G->A"),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))


library(NMF)
out <- nmf(t(signature_counts), rank=3)

w <- basis(out)
dim(w)

h <- t(coef(out))
dim(h)

hnorm <- t(apply(h, 1, function(x) return(x/sum(x))))

library(CountClust)

omega <- hnorm

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(c(rep("Ancient",25), rep("Modern",25)))
)

StructureGGplot(omega = omega,
                annotation = annotation,
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Moderns/Ancients",
                order_sample = FALSE,
                figure_title = paste0("NMF StructurePlot: K=", dim(omega)[2],": pmsignature: with C->T/G->A"),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))



