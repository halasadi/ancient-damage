

################ Test read memberships in aRchaic  ###################

library(aRchaic)

reads_Pinhasi <- read.csv(file = "../data/Pinhasi/NE2.realigned.sorted.q30.csv",
                          header=FALSE)


reads_moderns <- read.csv(file = "../data/moderns_lite/NA19795.mapped.ILLUMINA.bwa.MXL.low_coverage.20130415.realigned.sorted.q30.csv",
                          header=FALSE)

model <- get(load("../utilities/moderns_Pinhasi/clus_2/model.rda"))

system.time(output <- read_memberships(model,
                           reads_moderns,
                           subset = 1000,
                           nostrand = TRUE))

fit <- model
reads_file <- reads_moderns
subset <- 10
nostrand <- TRUE
verbose <- TRUE

reads_Lazaridis <- read.csv(file = "../utilities/reads_data/LBK.hg19_1000g.subsampled.q30.csv",
                            header=FALSE)

model <- get(load("../utilities/Lazaridis_moderns/clus_2/model.rda"))

system.time(output <- read_memberships(model,
                                       reads_Lazaridis,
                                       subset = 1000,
                                       nostrand = TRUE))


median(as.numeric(output$`cl- 1`))
median(as.numeric(output$`cl- 2`))
plot(density((as.numeric(output$`cl- 1`))))
plot(density((as.numeric(output$`cl- 2`))))


system.time(output <- read_memberships(model,
                           reads_Pinhasi,
                           subset = 5000,
                           nostrand = TRUE))

median(as.numeric(output$`cl- 1`))
median(as.numeric(output$`cl- 2`))
plot(density((as.numeric(output$`cl- 1`))))
plot(density((as.numeric(output$`cl- 2`))))

signatures <- output$sig

sig_split <- sapply(signatures, function(l) return(strsplit(l, ";")))
counts_sig_split <- as.vector(unlist(lapply(sig_split, function(x) return(length(x)))))

indices_1 <- which(counts_sig_split == 1)
indices_2 <- which(counts_sig_split == 2)




########     check for the independence of the double damage patterns   #######

count_CtoT <- as.vector(unlist(lapply(sig_split_2, function(l){
  return(length(grep("C->T", l)))
})))

prop_double_CtoT <- length(which(count_CtoT == 2))/length(count_CtoT)
prop_single_CtoT <- length(which(count_CtoT == 1))/length(count_CtoT)
marginal_CtoT <- sum(count_CtoT)/(2*length(count_CtoT))
marginal_CtoT*marginal_CtoT
prop_double_CtoT

###########     check for the independence of single base damages     ###########

count_CtoT <- as.vector(unlist(lapply(sig_split_1, function(l){
  return(length(grep("C->T", l)))
})))

prop_single_CtoT <- length(which(count_CtoT == 1))/length(count_CtoT)


