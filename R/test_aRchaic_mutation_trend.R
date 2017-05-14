

################  test aRchaic mutation trend  #########################

file <- "../data/Jones_2015/Bichon.sort.rmdup.IR.q30.mapDamage.subsampled.q30.csv"


par(mfrow=c(1,2))
aRchaic_mutation_trend(file,
                       plot_type = "left",
                       strand = "both",
                       layout_left = "topright",
                       layout_right = "topleft",
                       legend_cex = 0.3)

aRchaic_mutation_trend(file,
                       plot_type = "right",
                       strand = "both",
                       layout_left = "topright",
                       layout_right = "topleft",
                       legend_cex = 0.3)



file <- "../data/Pinhasi/KO1.q30.csv"


par(mfrow=c(1,2))
aRchaic_mutation_trend(file,
                       plot_type = "left",
                       strand = "both",
                       layout_left = "topright",
                       layout_right = "topleft",
                       legend_cex = 0.3,
                       lab_cex = 0.7)

aRchaic_mutation_trend(file,
                       plot_type = "right",
                       strand = "both",
                       layout_left = "topright",
                       layout_right = "topleft",
                       legend_cex = 0.3,
                       lab_cex = 0.7)



file <- "../data/AnnaGosling2016data/ADR-T1-CNE1.dup.q30.csv"


par(mfrow=c(1,2))
aRchaic_mutation_trend(file,
                       pattern = c("C->T", "G->A", "T->C", "C->A"),
                       plot_type = "left",
                       strand = "both",
                       layout_left = "topright",
                       layout_right = "topleft",
                       legend_cex = 0.3,
                       lab_cex = 0.7)

aRchaic_mutation_trend(file,
                       pattern = c("C->T", "G->A", "T->C", "C->A"),
                       plot_type = "right",
                       strand = "both",
                       layout_left = "topright",
                       layout_right = "topleft",
                       legend_cex = 0.3,
                       lab_cex = 0.7)



file <- "../data/AnnaGosling2016data/ADR-T1-KS20.dup.q30.csv"


par(mfrow=c(1,2))
aRchaic_mutation_trend(file,
                       plot_type = "left",
                       strand = "both",
                       layout_left = "topright",
                       layout_right = "topleft",
                       legend_cex = 0.3,
                       lab_cex = 0.7)

aRchaic_mutation_trend(file,
                       plot_type = "right",
                       strand = "both",
                       layout_left = "topright",
                       layout_right = "topleft",
                       legend_cex = 0.3,
                       lab_cex = 0.7)

