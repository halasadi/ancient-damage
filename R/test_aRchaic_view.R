
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
              output_dir = NULL)


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

###########################   Test aRchaic view   #############################
dir <- "../data/Tyrolean_Iceman/"

file = "teresa.srt.subsampled.q30"
breaks = c(-1, seq(1,20,1))
flanking_bases =1
logo.control = list()
output_dir <- NULL

sig_names = NULL
ic.scale=TRUE
max_pos = 20
flanking_bases=1
yscale_change = TRUE
xaxis=TRUE
yaxis=TRUE
xlab = " "
xaxis_fontsize=5
xlab_fontsize=10
title_aligner = 8
y_fontsize=10
title_fontsize = 20
mut_width=2
start=0.0001
renyi_alpha = 1
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
barport_height=0.25
output_dir = NULL
output_width = 1200
output_height = 700
save_plot=TRUE

dat <- get(load("../data/Mathieson-2015-subsampled/Mathieson-2015-subsampled.rda"))
club <- dat[1,]/sum(dat[1, ])
temp <- t(club)
temp <- t(temp)

theta_pool <- temp




library(grid)
library(gridBase)
library(Logolas)
damageLogo_five(temp)

pattern <- "I0011.1240k.subsampled.q30.csv"
out <- aggregate_signature_counts(dir = "../data/Mathieson-2015-subsampled/",
                                  pattern = pattern,
                                  breaks = c(-1, seq(1,20,1)),
                                  flanking_bases = 1)
clubbed_counts <- club_signature_counts(out, flanking_bases = 1)
clubbed_counts_norm <- clubbed_counts/ sum(clubbed_counts)

damageLogo_five(theta_pool = temp, output_dir = "../utilities/")

theta_pool = temp
output_dir = output_dir
save_plot = save_plot


aRchaic_view(dir = "../data/Tyrolean_Iceman/",
             file = "teresa.srt.subsampled.q30",
             breaks = c(-1, seq(1,20,1)),
             flanking_bases =1,
             logo.control = list(),
             output_dir = NULL,
             save_plot = FALSE)


file <- "../data/Mathieson/I0011.1240k.subsampled.q30.csv"



