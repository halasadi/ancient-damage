

##############   All data processing  ###########################

library(aRchaic)

dir <- "../data/Lindo2016ancients_strand_flank/";
out <- aggregate_bin_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))
save(out,
     file="../processed_data/lindo2016ancients-counts-table-strand-flank.rda")


dir <- "../data/Lindo2016moderns_strand_flank/";
out <- aggregate_bin_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))
save(out,
     file="../processed_data/lindo2016moderns-counts-table-strand-flank.rda")

dir <- "../data/Sardinia_data_strand_flank/";
out <- aggregate_bin_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))

save(out,
     file="../processed_data/sardinia2017-strand-flank.rda")

dir <- "../data/Sherpa_data_strand_flank/";
out <- aggregate_bin_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))

save(out,
     file="../processed_data/sherpa2017-strand-flank.rda")

dir <- "../data/AnnaGosling2016_strand_flank//";
out <- aggregate_bin_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))

save(out,
     file="../processed_data/annagosling2016-strand-flank.rda")


dir <- "../data/1000Gmoderns_data_strand_flank/";
out <- aggregate_bin_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))
save(out,
     file="../processed_data/1000Gmoderns-counts-table.rda")

dir <- "../data/I_data/"
out <- aggregate_bin_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))
save(out,
     file="../processed_data/Idata-counts-table.rda")

dir <- "../data/RISE_data/"
out <- aggregate_bin_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))
save(out,
     file="../processed_data/RISE-counts-table.rda")

dir <- "../data/HGDPmoderns_data_strand_flank/";
out <- aggregate_bin_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))
save(out,
     file="../processed_data/HGDPmoderns-counts-table-strand-flank.rda")

