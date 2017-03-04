

R CMD BATCH --vanilla '--args dir = "../data/Skoglund/" breaks = c(-1, seq(1,20,1)) flanking_bases = 1 output_rda = "../processed_data/skoglund-counts-table-strand-flank.rda" ' ../../aRchaic/R/aggregate_signature_counts_cmd.R &


R CMD BATCH --vanilla '--args dir="../data/Skoglund/" breaks=c(-1,seq(1,20,1)) flanking_bases=1 output_rda="../processed_data/skoglund-counts-table-strand-flank.rda"' ../../aRchaic/R/aggregate_signature_counts_cmd.R &
