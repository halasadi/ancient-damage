

###########################   Test aRchaic view   ###############################################

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



