
# tab_pooled <- c();
# library(dplyr)
# for(chr in 1:22){
#   tab <- read.csv(paste0("patt_chr_", chr, ".csv"), header=FALSE);
#   tab_new = tbl_df(tab) %>% group_by(V1) %>% mutate( percentage = V3/sum(V3)) %>%  as.data.frame()
#   tab1 = tab_new[grep("AAT->AAA", tab_new[,2]), ];
#   tab2 = tab_new[grep("TTA->TTT", tab_new[,2]), ];
#   tab_pool = cbind.data.frame(rbind(tab1, tab2), paste0('chr-', chr))
#   tab_pooled = rbind(tab_pooled, tab_pool)
#   cat("we are at chromosome", chr, "\n")
# }
#
# colnames(tab_pooled) <- c("pop", "sig", "count", "percentage", "chr")


######## We first read the data from of the AAT->AAA or TTA->TTT patterns  ###################


vcf_patterns_data <- get(load("../data/vcf_moderns/vcf_patterns_prop_across_chr.rda"))

library(dplyr)

tbl1 <- tbl_df(vcf_patterns_data) %>% filter(sig == "AAT->AAA") %>% select(pop, percentage) %>% group_by(pop) %>% summarise(sumprop=sum(percentage))
tbl2 <- tbl_df(vcf_patterns_data) %>% filter(sig == "TTA->TTT") %>% select(pop, percentage) %>% group_by(pop) %>% summarise(sumprop=sum(percentage))
par(mfrow=c(1,2))

library(readxl)
pops = read.table(file = '../data/vcf_moderns/igsr_samples.tsv', sep = '\t', header = TRUE)
pop_ids <- pops$Sample.name
pop_names <- pops$Population.code

populations <- pop_names[match(tbl1$pop, pop_ids)]

dat1 <- data.frame("pop"=populations, "samp"=tbl1$pop, "count"= tbl1$sumprop, "id"=1:dim(tbl1)[1])
dat2 <- data.frame("pop"=populations, "samp"=tbl2$pop, "count"= tbl2$sumprop, "id"=1:dim(tbl2)[1])

par(mfrow=c(1,2))
ggplot2::qplot(id, count, main="AAT->AAA", data=dat1, colour = populations)
ggplot2::qplot(id, count, main="TTA->TTT", data=dat2, colour = populations)

thousandg_data <- get(load("../processed_data/1000Gmoderns-counts-table.rda"));
signature_pops_thousandg <- as.character(sapply(rownames(thousandg_data), function(x) strsplit(x, "[.]")[[1]][1]))
pops_sig = as.character(sapply(rownames(thousandg_data), function(x) strsplit(x, "[.]")[[1]][5]))

pops4 = pops[match(signature_pops_thousandg, pops$Sample.name), ];
pop_ids2 <- pops4$Sample.name
pop_names2 <- pops_sig

indices <- which(is.na(match(pops4$Sample.name, tbl1$pop)))
tbl12 = tbl1[which(!is.na(match(pops4$Sample.name, tbl1$pop))),]
tbl22 = tbl2[which(!is.na(match(pops4$Sample.name, tbl2$pop))),]
populations2 <- pop_names2[-indices]

dat12 = data.frame("pop"=populations2, "samp"=tbl12$pop, "count"= tbl12$sumprop, "id"=1:dim(tbl12)[1])
dat22 = data.frame("pop"=populations2, "samp"=tbl22$pop, "count"= tbl22$sumprop, "id"=1:dim(tbl22)[1])

par(mfrow=c(1,2))
ggplot2::qplot(id, count, main="AAT->AAA", data=dat12, colour = pop)
ggplot2::qplot(id, count, main="TTA->TTT", data=dat22, colour = pop)



plot(tbl12$sumcount, col="red", pch=20)
plot(tbl22$sumcount, col="red", pch=20)




tbl <-  tbl_df(vcf_patterns_data) %>% filter(sig == "AAT->AAA")  %>% mutate(pop=pops3$Population.code)
plot(tbl$count, col="red", pch=20)
tbl$pop[order(tbl$count, decreasing=TRUE)[1:5]]



pops2 = pops[which(!is.na(match(pop_ids,  as.character(tbl$pop)))), ]

pops3 = pops2[match(tbl$pop, pops2$Sample.name), ]



tbl1 <- tbl_df(vcf_patterns_data) %>% filter(sig == "AAT->AAA") %>% select(pop, count) %>% group_by(pop) %>% summarise(sumcount=sum(count))
plot(tbl1$sumcount, col="red", pch=20)
tbl1$pop[order(tbl1$sumcount, decreasing=TRUE)[1:5]]


pops3[match(tbl1$pop[order(tbl1$sumcount, decreasing=TRUE)[1:100]], pops3$Sample.name), 4]



thousandg_data <- get(load("../processed_data/1000Gmoderns-counts-table.rda"));
signature_pops_thousandg <- as.character(sapply(rownames(thousandg_data), function(x) strsplit(x, "[.]")[[1]][1]))

pops = read.table(file = '../data/vcf_moderns/igsr_samples.tsv', sep = '\t', header = TRUE)

pops4 = pops[match(signature_pops_thousandg, pops$Sample.name), ];
pops4$Sample.name

tbl2 = tbl1[match(pops4$Sample.name, tbl1$pop),]
plot(tbl2$sumcount, col="red", pch=20)
pops4[match(tbl2$pop[order(tbl2$sumcount, decreasing=TRUE)[1:10]], pops4$Sample.name), 4]

tbl2$pop[order(tbl2$sumcount, decreasing=TRUE)[1:5]]
