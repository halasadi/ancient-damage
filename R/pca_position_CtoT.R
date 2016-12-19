

##############   PCA plot of C->T patterns along reads  ##########################

library(aRchaic)
sardinia_counts <- get(load("../processed_data/sardinia2017.rda"))
temp <- club_signature_counts(sardinia_counts)
sardinia_counts_CtoT <- filter_signatures_per_substitution(temp, pattern="C->T", use_prop = TRUE)

lindo_ancients_counts <- get(load("../processed_data/lindo2016ancients-counts-table.rda"))
temp <- club_signature_counts(lindo_ancients_counts)
lindo_ancients_counts_CtoT <- filter_signatures_per_substitution(temp, pattern="C->T", use_prop=TRUE)

lindo_moderns_counts <- get(load("../processed_data/lindo2016moderns-counts-table.rda"))
temp <- club_signature_counts(lindo_moderns_counts)
lindo_moderns_counts_CtoT <- filter_signatures_per_substitution(temp, pattern="C->T", use_prop = TRUE)

gosling_counts <- get(load("../processed_data/annagosling2016-counts-table.rda"))
temp <- club_signature_counts(gosling_counts)
temp <- temp[-28,];
gosling_counts_CtoT <- filter_signatures_per_substitution(temp, pattern="C->T", use_prop = TRUE)

sherpa_counts <- get(load("../processed_data/sherpa2017.rda"))
temp <- club_signature_counts(sherpa_counts)
sherpa_counts_CtoT <- filter_signatures_per_substitution(temp, pattern="C->T", use_prop = TRUE)

names <- rownames(gosling_counts_CtoT);
control_indices <- c(grep("EXN", names), grep("Libneg", names), grep("PCRneg", names))

gosling_controls_counts_CtoT <- gosling_counts_CtoT[control_indices,]
gosling_ancient_counts_CtoT <- gosling_counts_CtoT[-control_indices,]

ll <- list()
ll[["lindo-moderns"]] <- lindo_moderns_counts_CtoT
ll[["lindo-old"]] <- lindo_ancients_counts_CtoT
ll[["sardinia"]] <- sardinia_counts_CtoT
ll[["sherpa"]] <- sherpa_counts_CtoT
ll[["gosling-control"]] <- gosling_controls_counts_CtoT
ll[["gosling-ancient"]] <- gosling_ancient_counts_CtoT

out <- gridPCA_signatures_combo(ll,
                         input_pos=1:10,
                         normalize=FALSE)

