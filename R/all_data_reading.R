

##############   All data processing  ###########################

library(aRchaic)

dir <- "../data/Lindo2016ancients_strand_flank/";
out <- aggregate_signature_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))
save(out,
     file="../processed_data/lindo2016ancients-counts-table-strand-flank.rda")


dir <- "../data/Lindo2016moderns_strand_flank/";
out <- aggregate_signature_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))
save(out,
     file="../processed_data/lindo2016moderns-counts-table-strand-flank.rda")

dir <- "../data/Sardinia_data_strand_flank/";
out <- aggregate_signature_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))

save(out,
     file="../processed_data/sardinia2017-strand-flank.rda")

dir <- "../data/Sherpa_data_strand_flank/";
out <- aggregate_signature_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))

save(out,
     file="../processed_data/sherpa2017-strand-flank.rda")

dir <- "../data/AnnaGosling2016_strand_flank//";
out <- aggregate_signature_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))

save(out,
     file="../processed_data/annagosling2016-strand-flank.rda")


dir <- "../data/1000Gmoderns_data_strand_flank/";
out <- aggregate_signature_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))
save(out,
     file="../processed_data/1000Gmoderns-counts-table.rda")

##################   I  data   ##################################

dir <- "../data/I_data/"
out <- aggregate_signature_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))
save(out,
     file="../processed_data/Idata-counts-table.rda")

##############   RISE data  #######################################
dir <- "../data/RISE_data/"
out <- aggregate_signature_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))
save(out,
     file="../processed_data/RISE-counts-table.rda")

###############   HGDP moderns data  ##############################

dir <- "../data/HGDPmoderns_data_strand_flank/";
out <- aggregate_signature_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))
save(out,
     file="../processed_data/HGDPmoderns-counts-table-strand-flank.rda")

################  Sardinia ancients data  ##########################

dir <- "../data/Sardinia_data_strand_flank/";
out <- aggregate_signature_counts(dir, breaks = c(-1, seq(1,20,1)))
save(out,
     file="../processed_data/sardinia2017-strand-flank.rda")


##################  1000 Genomes moderns data   ##########################

dir <- "../data/1000Gmoderns_data_strand_flank_2/";
out <- aggregate_signature_counts(dir, breaks = c(-1, seq(1,20,1)))
save(out,
     file="../processed_data/1000g-strand-flank-2.rda")

####################   Fu data reading and loading  #######################

dir <- "../data/Fu_2016/";
out <- aggregate_signature_counts(dir, breaks = c(-1, seq(1,20,1)))
save(out,
     file="../processed_data/Fu2016-strand-flank-2.rda")

clubbed_counts <- club_signature_counts_2(out, flanking_bases = 1)


####################   Allentoft   ################################


dir <- "../data/Allentoft/";
out <- aggregate_signature_counts(dir, breaks = c(-1, seq(1,20,1)), flanking_bases = 1)
save(out,
     file="../processed_data/allentoft-counts-table-strand-flank.rda")

dir <- "../data/Skoglund/";
out <- aggregate_signature_counts(dir, breaks = c(-1, seq(1,20,1)), flanking_bases = 1)

