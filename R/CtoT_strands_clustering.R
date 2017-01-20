library(aRchaic)
sardinia_counts <- get(load("../processed_data/sardinia2017.rda"))
temp <- club_signature_counts(sardinia_counts)
sardinia_counts_CtoT <- filter_signatures_per_substitution(temp, pattern="C->T", use_prop = FALSE)
sardinia_counts_CtoT <- filter_signatures_only_location(sardinia_counts_CtoT, max_pos = 10)

left_strand_breaks <- as.numeric()
right_strand_breaks <- as.numeric()

strand_breaks_RISE <- get(load("../processed_data/strand-breaks-sardinia2017.rda"))
for(l in 1:length(strand_breaks_RISE)){
  left_strand_breaks <- rbind(left_strand_breaks, strand_breaks_RISE[[l]][[1]])
  right_strand_breaks <- rbind(right_strand_breaks, strand_breaks_RISE[[l]][[2]])
}
colnames(left_strand_breaks) <- paste0(colnames(left_strand_breaks), "_left")
colnames(right_strand_breaks) <- paste0(colnames(right_strand_breaks), "_right")


sardinia_filtered <- cbind(sardinia_counts_CtoT, left_strand_breaks, right_strand_breaks)
