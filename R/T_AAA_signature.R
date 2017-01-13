

##########   T->AAA  signatures  ##########################

library(aRchaic)

lindoancient_data <- get(load("../processed_data/lindo2016ancients-counts-table.rda"))
lindomoderns_data <- get(load("../processed_data/lindo2016moderns-counts-table.rda"))
sardinia_data <- get(load("../processed_data/sardinia2017.rda"))
sherpa_data <- get(load("../processed_data/sherpa2017.rda"))
thousandg_data <- get(load("../processed_data/1000Gmoderns-counts-table.rda"));
hgdp_data <- get(load("../processed_data/HGDPmoderns-counts-table.rda"))
gossling_data <- get(load("../processed_data/annagosling2016-counts-table.rda"))


lindoancient_data_clubbed <- club_signature_counts(lindoancient_data)
lindomoderns_data_clubbed <- club_signature_counts(lindomoderns_data)
sardinia_data_clubbed <- club_signature_counts(sardinia_data)
sherpa_data_clubbed <- club_signature_counts(sherpa_data)
thousandg_data_clubbed <- club_signature_counts(thousandg_data)
hgdp_data_clubbed <- club_signature_counts(hgdp_data)
gossling_data_clubbed <- club_signature_counts(gossling_data)

lindoancient_data_filtered <- filter_signatures_wo_location(lindoancient_data_clubbed)
lindomoderns_data_filtered <- filter_signatures_wo_location(lindomoderns_data_clubbed)
sardinia_data_filtered <- filter_signatures_wo_location(sardinia_data_clubbed)
sherpa_data_filtered <- filter_signatures_wo_location(sherpa_data_clubbed)
thousandg_data_filtered <- filter_signatures_wo_location(thousandg_data_clubbed)
hgdp_data_filtered <- filter_signatures_wo_location(hgdp_data_clubbed)
gossling_data_filtered <- filter_signatures_wo_location(gossling_data_clubbed)

sig_T_AAA_counts <- c();
sig_T_AAA_counts <- c(sig_T_AAA_counts, rowSums(lindoancient_data_filtered[,grep("T->AAA", colnames(lindoancient_data_filtered))])/ rowSums(lindoancient_data_filtered))
sig_T_AAA_counts <- c(sig_T_AAA_counts, rowSums(lindomoderns_data_filtered[,grep("T->AAA", colnames(lindomoderns_data_filtered))])/ rowSums(lindomoderns_data_filtered))
sig_T_AAA_counts <- c(sig_T_AAA_counts, rowSums(sardinia_data_filtered[,grep("T->AAA", colnames(sardinia_data_filtered))])/ rowSums(sardinia_data_filtered))
sig_T_AAA_counts <- c(sig_T_AAA_counts, rowSums(sherpa_data_filtered[,grep("T->AAA", colnames(sherpa_data_filtered))])/ rowSums(sherpa_data_filtered))
sig_T_AAA_counts <- c(sig_T_AAA_counts, rowSums(thousandg_data_filtered[,grep("T->AAA", colnames(thousandg_data_filtered))])/ rowSums(thousandg_data_filtered))
sig_T_AAA_counts <- c(sig_T_AAA_counts, rowSums(hgdp_data_filtered[,grep("T->AAA", colnames(hgdp_data_filtered))])/ rowSums(hgdp_data_filtered))
sig_T_AAA_counts <- c(sig_T_AAA_counts, rowSums(gossling_data_filtered[,grep("T->AAA", colnames(gossling_data_filtered))])/ rowSums(gossling_data_filtered))

sig_T_AAA_counts <- c();
sig_T_AAA_counts <- c(sig_T_AAA_counts, lindoancient_data_filtered[,grep("AAT->AAA", colnames(lindoancient_data_filtered))]/rowSums(lindoancient_data_filtered))
sig_T_AAA_counts <- c(sig_T_AAA_counts, lindomoderns_data_filtered[,grep("AAT->AAA", colnames(lindomoderns_data_filtered))]/rowSums(lindomoderns_data_filtered))
sig_T_AAA_counts <- c(sig_T_AAA_counts, sardinia_data_filtered[,grep("AAT->AAA", colnames(sardinia_data_filtered))]/rowSums(sardinia_data_filtered))
sig_T_AAA_counts <- c(sig_T_AAA_counts, sherpa_data_filtered[,grep("AAT->AAA", colnames(sherpa_data_filtered))]/rowSums(sherpa_data_filtered))
sig_T_AAA_counts <- c(sig_T_AAA_counts, thousandg_data_filtered[,grep("AAT->AAA", colnames(thousandg_data_filtered))]/rowSums(thousandg_data_filtered))
sig_T_AAA_counts <- c(sig_T_AAA_counts, hgdp_data_filtered[,grep("AAT->AAA", colnames(hgdp_data_filtered))]/rowSums(hgdp_data_filtered))
sig_T_AAA_counts <- c(sig_T_AAA_counts, gossling_data_filtered[,grep("AAT->AAA", colnames(gossling_data_filtered))]/rowSums(gossling_data_filtered))

labs <- c(rep("Lindo-ancient", dim(lindoancient_data_filtered)[1]),
          rep("Lindo-modern", dim(lindomoderns_data_filtered)[1]),
          rep("Sardinia", dim(sardinia_data_filtered)[1]),
          rep("Sherpa", dim(sherpa_data_filtered)[1]),
          rep("1000G", dim(thousandg_data_filtered)[1]),
          rep("HGDP", dim(hgdp_data_filtered)[1]))

names <- rownames(gossling_data_filtered);
control_indices <- c(grep("EXN", names), grep("Libneg", names), grep("PCRneg", names))


labs1 <- character();
labs1 <- rep("Gossling-ancient", dim(gossling_data_filtered)[1])
labs1[control_indices] <- "Gossling-controls"

labs <- c(labs, labs1)

cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
         "hotpink","burlywood","darkkhaki","yellow","darkgray","deepskyblue",
         "brown4","darkorchid","magenta", "azure1","azure4")

par(mfrow=c(1,1))
par(mar=c(2,2,2,2))
plot(sig_T_AAA_counts, col=cols[as.numeric(factor(labs))], pch=20, main="prop. of AAT->AAA across aDNA sources")
legend("topright", legend=levels(factor(labs)),fill=cols[1:length(levels(factor(labs)))], cex=0.5)

which(sig_T_AAA_counts > 0.01)
#########  total damage distribution ##################


totdamage <- c()
totdamage <- c(totdamage, rowSums(lindoancient_data_filtered))
totdamage <- c(totdamage, rowSums(lindomoderns_data_filtered))
totdamage <- c(totdamage, rowSums(sardinia_data_filtered))
totdamage <- c(totdamage, rowSums(sherpa_data_filtered))
totdamage <- c(totdamage, rowSums(thousandg_data_filtered))
totdamage <- c(totdamage, rowSums(hgdp_data_filtered))
totdamage <- c(totdamage, rowSums(gossling_data_filtered))


par(mfrow=c(1,1))
par(mar=c(2,2,2,2))
plot(log(totdamage), col=cols[as.numeric(factor(labs))], pch=20, main="Log no. of damages across aDNA sources")
legend("bottomleft", legend=levels(factor(labs)),fill=cols[1:length(levels(factor(labs)))], cex=0.5)
