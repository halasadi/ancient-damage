

############## Modern mutation proportion and nucleotide comp  ##########################


file <- read.csv("../data/moderns_lite/HG00100.mapped.ILLUMINA.bwa.GBR.low_coverage.20130415.q30.csv",
                 header=FALSE)

muts <- file[,1]
leftflankcomp <- substring(muts, 1, 1)
bases_freq <- tapply(file[,7], leftflankcomp, sum)
bases_freq <- bases_freq[match(c("A", "C", "G", "T"), names(bases_freq))]
bases_freq_mixcomp_left <- bases_freq/sum(bases_freq)

rightflankcomp <- substring(muts, 6, 6)
bases_freq <- tapply(file[,7], rightflankcomp, sum)
bases_freq <- bases_freq[match(c("A", "C", "G", "T"), names(bases_freq))]
bases_freq_mixcomp_right <- bases_freq/sum(bases_freq)


mutations <- substring(muts, 2, 5)
mutations_freq <- tapply(file[,7], mutations, sum)
mutations_freq <- mutations_freq[match(c("C->T", "C->G", "C->A", "T->C", "T->A", "T->G"), names(mutations_freq))]
mutations_freq_mixcomp <- mutations_freq/sum(mutations_freq)

topic_clus <- get(load("../utilities/Nepal_sardinia_moderns/clus_5/model.rda"))
theta <- topic_clus$theta

out <- theta_breakdown(theta)
out[[2]]$`mismatch-flank`[-(1:4),2]
mutations_freq_mixcomp
