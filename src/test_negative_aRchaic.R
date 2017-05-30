

######################## test aRchaic negative plots  #####################################


topic_clus <- get(load("../utilities/Nepal_sardinia_moderns/clus_5/model.rda"))
theta <- topic_clus$theta

K <- 5
logo.control <- list()

logo.control.default <- list(sig_names = NULL, ic.scale=TRUE,
                             max_pos = 20, flanking_bases=1,
                             yscale_change = TRUE, xaxis=TRUE,
                             yaxis=TRUE, xlab = " ", xaxis_fontsize=30,
                             xlab_fontsize=10, title_aligner = 11,
                             y_fontsize=27, title_fontsize = 35,
                             mut_width=2, start=0.0001,
                             renyi_alpha = 5, inflation_factor = c(3,1,3),
                             pop_names = paste0("Cluster : ", 1:K),
                             logoport_x = 0.25, logoport_y= 0.50, logoport_width= 0.25, logoport_height= 0.50,
                             lineport_x = 0.9, lineport_y=0.40, lineport_width=0.32, lineport_height=0.28,
                             breaklogoport_x = 0.94, breaklogoport_y = 0.40, breaklogoport_width=0.30, breaklogoport_height=0.45,
                             barport_x = 0.60, barport_y=0.60, barport_width=0.25, barport_height=0.35,
                             output_width = 1200, output_height = 700)

logo.control <- modifyList(logo.control.default, logo.control)


do.call(damageLogo_five_neg, append(list(theta_pool = theta,
                                     output_dir = "../utilities/Nepal_sardinia_moderns/clus_5/"),
                                logo.control))

theta_pool <- theta
sig_names = NULL
ic.scale=TRUE
max_pos = 20
flanking_bases=1
yscale_change = TRUE
xaxis=TRUE
yaxis=TRUE
xlab = " "
xaxis_fontsize=10
xlab_fontsize=20
title_aligner = 15
y_fontsize=20
title_fontsize = 20
mut_width=2
start=0.0001
renyi_alpha = 1
inflation_factor = c(2,1,2)
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
barport_height=0.40
output_dir = NULL
output_width = 1200
output_height = 700

pwm <- prop_patterns_list[[2]]
ic <- ic[,2]
probs = prob_mutation[2,]
breaks_theta_vec = breaks_theta[,2, drop=FALSE]
strand_theta_vec = strand_theta[2,]
