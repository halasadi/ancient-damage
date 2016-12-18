

##############   All data processing  ###########################

library(aRchaic)

dir <- "../data/Lindo2016ancients/";
out <- aggregate_bin_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))
save(out,
     file="../processed_data/lindo2016ancients-counts-table.rda")


dir <- "../data/Lindo2016moderns/";
out <- aggregate_bin_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))
save(out,
     file="../processed_data/lindo2016moderns-counts-table.rda")

dir <- "../data/Sardinia2017/";
out <- aggregate_bin_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))

save(out,
     file="../processed_data/sardinia2017.rda")

dir <- "../data/Sherpa_data/";
out <- aggregate_bin_counts(dir, breaks = c(-1, seq(1,20,1), 25, 35))

save(out,
     file="../processed_data/sherpa2017.rda")
