

###############    aRchaic on wolves   ##############################


##########   aRchiac.view  :  View contents of a file  #####################


source("../../aRchaic.site/R/aggregate_signature_counts.R")
source("../../aRchaic.site/R/club_signature_counts.R")
source("../../aRchaic.site/R/aRchaic_view.R")
source("../../aRchaic.site/R/filter_signatures_only_location.R")
source("../../aRchaic.site/R/aRchaic_cluster.R")
source("../../aRchaic.site/R/aRchaic_pool.R")
source("../../aRchaic.site/R/damageLogo_5.R")


getwd()

dir = "../data/wolves/"
pattern = "chr01_RKW1548.subsampled.q30.csv"
breaks = c(-1, seq(1, 20, 1))
flanking_bases = 1
out <- aggregate_signature_counts(dir = dir,
                                  pattern = pattern,
                                  breaks = breaks,
                                  flanking_bases = flanking_bases)

clubbed_counts <- club_signature_counts(out, flanking_bases = 1)
clubbed_counts_norm <- clubbed_counts/ sum(clubbed_counts)


out <- aRchaic_view(dir = "../data/wolves/",
                    file = "chr01_RKW1548.subsampled.q30",
                    breaks = c(-1, seq(1, 20, 1)),
                    flanking_bases = 1,
                    logo.control = list(),
                    output_dir = "../utilities/wolves/",
                    save_plot = TRUE)


plot.new()
clus_out <- aRchaic_cluster(folders = "../data/wolves/",
                            K = 4,
                            labs = c(rep("RKW", 2), rep("RSRW", 5)),
                            tol = 0.5,
                            run_from = "plot",
                            output_dir = "../utilities/wolves/clus_4/",
                            save_plot = TRUE)



folders = "../data/wolves/"
K = 3
tol=0.5
labs = c(rep("RKW", 2), rep("RSRW", 5))
levels = NULL
run_from = "plot"
run_index = 1:length(folders),
breaks = c(-1, seq(1,20,1))
flanking_bases = 1
gom_method = "independent"
topic_cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
               "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
               "brown4","darkorchid","magenta","yellow", "azure1","azure4")
structure.control = list()
logo.control = list()
topics.control = list()
output_dir = "../utilities/wolves/clus_3/"
save_plot = TRUE
structure_width = 5
structure_height = 8





data <- get(load(file = "../data/wolves/wolves.rda"))
topic_clus  <- maptpx::topics(data, K=2, model = "full")

topic_clus <- get(load(file = "../utilities/wolves/clus_4/model.rda"))


K <- 4
levels <- unique(labs)
omega <- topic_clus$omega
annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs, levels = levels)
)

png(paste0("../utilities/wolves/clus_4/", "structure.png"), width = 4, height = 10)
do.call(CountClust::StructureGGplot, append(list(omega= omega,
                                                 annotation = annotation,
                                                 palette = topic_cols),
        structure.control))
ggsave(paste0("../utilities/wolves/clus_4/", "structure.png"), width=4, height=10)
graphics.off()

structure.control = list()
logo.control = list()
topics.control = list()







