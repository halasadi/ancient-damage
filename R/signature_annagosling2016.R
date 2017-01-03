

############   signature plot Anna Gosling 2016  ########################


topics_clus <- get(load(file="../rda/annagosling2016/topics-annagosling2016-k-3-filtered.rda"))
theta <- topics_clus$theta

damageLogo(theta, ic.scale = FALSE, pop_names = paste0("Cluster ",1:dim(theta)[2]))

omega <- topics_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs)
)


CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = FALSE,
                            figure_title = paste0("StructurePlot: K=", dim(omega)[2],": pmsignature: with position information"),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))

signature_set <- rownames(theta)
apply(CountClust::ExtractTopFeatures(theta, top_features = 10, method="poisson", options="min"), c(1,2), function(x) signature_set[x])


##############   just mutation   #################################

filtered_counts_2 <- filter_signatures_w_mutation(clubbed_counts);
topics_clus <- maptpx::topics(filtered_counts_2, K=2, tol=0.1);

theta <- topics_clus$theta
omega <- topics_clus$omega
CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = FALSE,
                            figure_title = paste0("StructurePlot: K=", dim(omega)[2],": pmsignature: withput position information"),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))
signature_set <- colnames(filtered_counts_2)
apply(CountClust::ExtractTopFeatures(topics_clus$theta, top_features = 1, method="poisson", options="min"), c(1,2), function(x) signature_set[x])


##############  no X  data   ######################################

topics_clus <- get(load(file="../rda/annagosling2016/topics-annagosling2016-k-3-filtered_reduced_wo_position_noX.rda"))

theta <- topics_clus$theta

damageLogo(theta, ic.scale = FALSE, pop_names = paste0("Cluster ",1:dim(theta)[2]))

omega <- topics_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs)
)

rownames(omega) <- annotation$sample_id;

CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = FALSE,
                            figure_title = paste0("StructurePlot: K=", dim(omega)[2],": pmsignature: wo position, wo X"),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))

signature_set <-rownames(theta)
apply(CountClust::ExtractTopFeatures(topics_clus$theta, top_features = 10, method="poisson", options="min"), c(1,2), function(x) signature_set[x])


################   pattern plot damage   ##########################


pattern_plot_full(file="../data/AnnaGosling2016data/ADR-T2-EXN3.dup.q30.csv",
             pattern = c("C->T", "C->G", "C->A", "T->G", "T->A", "T->C"),
             plot_type="left",
             sample_name="EXN3")

color=c("red","blue","cornflowerblue","black","cyan","darkblue",
        "brown4","burlywood","darkgoldenrod1","darkgray","deepskyblue","darkkhaki",
        "firebrick","darkorchid","hotpink","green","magenta","yellow", "azure1","azure4");


pattern_plot_full(file="../data/AnnaGosling2016data/ADR-T1-S34.dup.q30.csv",
             pattern = c("C->T", "C->G", "C->A", "T->G", "T->A", "T->C",
                         "G->A", "G->C", "G->T", "A_>G", "A->T", "A->C"),
             plot_type="left",
             sample_name = "S34",
             cols = color)


pattern_plot_full(file="../data/AnnaGosling2016data/ADR-T2-Libneg3.dup.q30.csv",
                  pattern = c("C->T", "C->G", "C->A", "T->G", "T->A", "T->C",
                              "G->A", "G->C", "G->T", "A_>G", "A->T", "A->C"),
                  plot_type="left",
                  sample_name = "Libneg3",
                  cols = color)

pattern_plot_full(file="../data/AnnaGosling2016data/ADR-T2-Libneg2.dup.q30.csv",
                  pattern = c("C->T", "C->G", "C->A", "T->G", "T->A", "T->C",
                              "G->A", "G->C", "G->T", "A_>G", "A->T", "A->C"),
                  plot_type="left",
                  sample_name = "Libneg2",
                  cols=color)

pattern_plot_full(file="../data/AnnaGosling2016data/ADR-T2-Libneg1.dup.q30.csv",
                  pattern = c("C->T", "C->G", "C->A", "T->G", "T->A", "T->C",
                              "G->A", "G->C", "G->T", "A_>G", "A->T", "A->C"),
                  plot_type="left",
                  sample_name = "Libneg1",
                  cols=color)

pattern_plot_full(file="../data/AnnaGosling2016data/ADR-T2-PCRneg.dup.q30.csv",
                  pattern = c("C->T", "C->G", "C->A", "T->G", "T->A", "T->C",
                              "G->A", "G->C", "G->T", "A_>G", "A->T", "A->C"),
                  plot_type="left",
                  sample_name = "PCRneg",
                  cols = color)

color=c("red","blue","cornflowerblue","black","cyan","darkblue",
        "brown4","burlywood","darkgoldenrod1","darkgray","deepskyblue","darkkhaki",
        "firebrick","darkorchid","hotpink","green","magenta","yellow", "azure1","azure4");

color1 <- color;
color1[7:9] <- color[1:3]
color1[1:3] <- color[7:9]

par(mfrow=c(1,2))

pattern_plot_full(file="../data/AnnaGosling2016data/ADR-T1-S20.dup.q30.csv",
                  pattern = c("C->T", "C->G", "C->A", "T->G", "T->A", "T->C",
                              "G->A", "G->C", "G->T", "A_>G", "A->T", "A->C"),
                  plot_type="left",
                  sample_name = "S20",
                  cols = color)

pattern_plot_full(file="../data/AnnaGosling2016data/ADR-T1-S20.dup.q30.csv",
             pattern = c("G->A", "G->C", "G->T",  "T->G", "T->A", "T->C",
                         "C->T", "C->G", "C->A",  "A->G", "A->T", "A->C"),
             plot_type="right",
             sample_name = "S20",
             cols = color1)



par(mfrow=c(1,2))

pattern_plot_full(file="../data/AnnaGosling2016data/ADR-T1-S34.dup.q30.csv",
                  pattern = c("C->T", "C->G", "C->A", "T->G", "T->A", "T->C",
                              "G->A", "G->C", "G->T", "A_>G", "A->T", "A->C"),
                  plot_type="left",
                  sample_name = "S34",
                  cols = color)

pattern_plot_full(file="../data/AnnaGosling2016data/ADR-T1-S34.dup.q30.csv",
                  pattern = c("G->A", "G->C", "G->T",  "T->G", "T->A", "T->C",
                              "C->T", "C->G", "C->A",  "A->G", "A->T", "A->C"),
                  plot_type="right",
                  sample_name = "S34",
                  cols = color1)


par(mfrow=c(1,2))

pattern_plot_full(file="../data/AnnaGosling2016data/ADR-T2-EXN3.dup.q30.csv",
                  pattern = c("C->T", "C->G", "C->A", "T->G", "T->A", "T->C",
                              "G->A", "G->C", "G->T", "A_>G", "A->T", "A->C"),
                  plot_type="left",
                  sample_name = "EXN3",
                  cols = color)

pattern_plot_full(file="../data/AnnaGosling2016data/ADR-T2-EXN3.dup.q30.csv",
                  pattern = c("G->A", "G->C", "G->T",  "T->G", "T->A", "T->C",
                              "C->T", "C->G", "C->A",  "A->G", "A->T", "A->C"),
                  plot_type="right",
                  sample_name = "EXN3",
                  cols = color1)


par(mfrow=c(1,2))

pattern_plot_full(file="../data/AnnaGosling2016data/ADR-T2-Libneg1.dup.q30.csv",
                  pattern = c("C->T", "C->G", "C->A", "T->G", "T->A", "T->C",
                              "G->A", "G->C", "G->T", "A_>G", "A->T", "A->C"),
                  plot_type="left",
                  sample_name = "Libneg1",
                  cols = color)

pattern_plot_full(file="../data/AnnaGosling2016data/ADR-T2-Libneg1.dup.q30.csv",
                  pattern = c("G->A", "G->C", "G->T",  "T->G", "T->A", "T->C",
                              "C->T", "C->G", "C->A",  "A->G", "A->T", "A->C"),
                  plot_type="right",
                  sample_name = "Libneg1",
                  cols = color1)


par(mfrow=c(1,2))

pattern_plot_full(file="../data/AnnaGosling2016data/ADR-T2-Libneg2.dup.q30.csv",
                  pattern = c("C->T", "C->G", "C->A", "T->G", "T->A", "T->C",
                              "G->A", "G->C", "G->T", "A_>G", "A->T", "A->C"),
                  plot_type="left",
                  sample_name = "Libneg2",
                  cols = color)

pattern_plot_full(file="../data/AnnaGosling2016data/ADR-T2-Libneg2.dup.q30.csv",
                  pattern = c("G->A", "G->C", "G->T",  "T->G", "T->A", "T->C",
                              "C->T", "C->G", "C->A",  "A->G", "A->T", "A->C"),
                  plot_type="right",
                  sample_name = "Libneg2",
                  cols = color1)


par(mfrow=c(1,2))

pattern_plot_full(file="../data/AnnaGosling2016data/ADR-T2-EXN1.dup.q30.csv",
                  pattern = c("C->T", "C->G", "C->A", "T->G", "T->A", "T->C",
                              "G->A", "G->C", "G->T", "A_>G", "A->T", "A->C"),
                  plot_type="left",
                  sample_name = "EXN1",
                  cols = color)

pattern_plot_full(file="../data/AnnaGosling2016data/ADR-T2-EXN1.dup.q30.csv",
                  pattern = c("G->A", "G->C", "G->T",  "T->G", "T->A", "T->C",
                              "C->T", "C->G", "C->A",  "A->G", "A->T", "A->C"),
                  plot_type="right",
                  sample_name = "EXN1",
                  cols = color1)


par(mfrow=c(1,2))

pattern_plot_full(file="../data/AnnaGosling2016data/ADR-T2-EXN2.dup.q30.csv",
                  pattern = c("C->T", "C->G", "C->A", "T->G", "T->A", "T->C",
                              "G->A", "G->C", "G->T", "A_>G", "A->T", "A->C"),
                  plot_type="left",
                  sample_name = "EXN2",
                  cols = color)

pattern_plot_full(file="../data/AnnaGosling2016data/ADR-T2-EXN2.dup.q30.csv",
                  pattern = c("G->A", "G->C", "G->T",  "T->G", "T->A", "T->C",
                              "C->T", "C->G", "C->A",  "A->G", "A->T", "A->C"),
                  plot_type="right",
                  sample_name = "EXN2",
                  cols = color1)


par(mfrow=c(1,2))

pattern_plot_full(file="../data/AnnaGosling2016data/ADR-T1-S32.dup.q30.csv",
                  pattern = c("C->T", "C->G", "C->A", "T->G", "T->A", "T->C",
                              "G->A", "G->C", "G->T", "A_>G", "A->T", "A->C"),
                  plot_type="left",
                  sample_name = "S32",
                  cols = color)

pattern_plot_full(file="../data/AnnaGosling2016data/ADR-T1-S32.dup.q30.csv",
                  pattern = c("G->A", "G->C", "G->T",  "T->G", "T->A", "T->C",
                              "C->T", "C->G", "C->A",  "A->G", "A->T", "A->C"),
                  plot_type="right",
                  sample_name = "S32",
                  cols = color1)
