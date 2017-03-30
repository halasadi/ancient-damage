

########   test aRchaic view  ###############

aRchaic_view (file = "../data/Skoglund/Ajv52.hs37d5.fa.merged.q30.csv",
              breaks = c(-1, seq(1,20,1)),
              flanking_bases =1,
              logo.control = list(),
              title = "Skoglund aDNA",
              output_dir = NULL)

grid.newpage()
par(new=TRUE)
par(mar=c(1,1,1,1))
aRchaic_view (file = "../data/Skoglund/Ajv59.hs37d5.fa.merged.q30.csv",
              breaks = c(-1, seq(1,20, 2)),
              flanking_bases =1,
              logo.control = list(),
              title = "Skoglund aDNA",
              output_dir = NULL,
              save_plot = FALSE)


aRchaic_view (file = "../data/Siberia/MA1_1stextraction.hg19.q30.csv",
              breaks = c(-1, seq(1,20,1)),
              flanking_bases =1,
              logo.control = list(renyi_alpha = 200),
              title = "Siberia aDNA",
              output_dir = NULL,
              save_plot = TRUE)

file <- "../data/Skoglund/Ajv59.hs37d5.fa.merged.q30.csv"
title <- head(strsplit(rev((as.vector(strsplit(file, "/" )[[1]])))[1], ".csv")[[1]],1)
file <- get(load("../data/Skoglund/Skoglund.rda"))

read_length_distribution(dir = "../data/Skoglund/",
                         pattern = "Ajv59.hs37d5.fa.merged.q30.csv",
                         end_break = 5,
                         plot_layout = c(1,1),
                         cols = c("red", "green", "blue"),
                         cex_legend = 0.5,
                         cex.main = 0.5)


