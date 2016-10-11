

#######  testing logo plots on Lindo data #############


topics_clus <- get(load("../rda/CountClust_output_Lindo2016_with_C_to_T.rda"));


#######  Try the default set up  #################


theta <- topics_clus$clust_2$theta
damageLogo(theta)

theta <- topics_clus$clust_3$theta
damageLogo(theta)

theta <- topics_clus$clust_4$theta
damageLogo(theta)


#######  suppose you wanna play with mutation width ####


damageLogo(theta, mut_width = 3)
damageLogo(theta, mut_width = 1.5)


#######  suppose you want the scales to be same for all ###


damageLogo(theta, mut_width = 2, yscale_change = FALSE)

######  suppose you want the height to be 1 for all ###

damageLogo(theta, ic.scale = FALSE)

## change the axis label size

damageLogo(theta, xaxis_fontsize = 11)
