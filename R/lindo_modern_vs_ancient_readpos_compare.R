

########  Lindo et al data analysis ####################

tab <- read.csv("../data/Lindo2016ancients/125_all_chr.q30.csv", header=FALSE);
all_tab <- tab[,2] + tab[,3]
mean(all_tab)
min(all_tab)
table(tab[,3])


tab <- read.csv("../data/Lindo2016moderns/S001_all_chr.q30.csv", header=FALSE);
tbl_df(tab) %>% filter(V1 %in% c("AAT->AAA", "TTA->TTT")) %>% arrange(desc(V4))

tab <- read.csv("../data/Lindo2016ancients/163_all_chr.q30.csv", header=FALSE)
tbl_df(tab) %>% filter(V1 %in% c("AAT->AAA", "TTA->TTT")) %>% arrange(desc(V4))

tab <- read.csv("../data/Lindo2016ancients/939_all_chr.q30.csv", header=FALSE)
tbl_df(tab) %>% filter(V1 %in% c("AAT->AAA", "TTA->TTT")) %>% arrange(desc(V4))

tab <- read.csv("../data/Lindo2016ancients/468_all_chr.q30.csv", header=FALSE)
tbl_df(tab) %>% filter(V1 %in% c("AAT->AAA", "TTA->TTT")) %>% arrange(desc(V4))


tab<- read.csv("../data/1000Gmoderns/NA19675.mapped.ILLUMINA.bwa.MXL.low_coverage.20120522.q30.csv", header=FALSE)
tbl_df(tab) %>% filter(V1 %in% c("AAT->AAA", "TTA->TTT")) %>% arrange(desc(V4))

tab<- read.csv("../data/1000Gmoderns/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.q30.csv", header=FALSE)
tbl_df(tab) %>% filter(V1 %in% c("AAT->AAA", "TTA->TTT")) %>% arrange(desc(V4))

tab<- read.csv("../data/1000Gmoderns/HG00599.mapped.ILLUMINA.bwa.CHS.low_coverage.20121211.q30.csv", header=FALSE)
tbl_df(tab) %>% filter(V1 %in% c("AAT->AAA", "TTA->TTT")) %>% arrange(desc(V4))

tab<- read.csv("../data/1000Gmoderns/NA19717.mapped.ILLUMINA.bwa.MXL.low_coverage.20120522.q30.csv", header=FALSE)
tbl_df(tab) %>% filter(V1 %in% c("AAT->AAA", "TTA->TTT")) %>% arrange(desc(V4))

tab<- read.csv("../data/1000Gmoderns/HG01047.mapped.ILLUMINA.bwa.PUR.low_coverage.20120522.q30.csv", header=FALSE)
tbl_df(tab) %>% filter(V1 %in% c("AAT->AAA", "TTA->TTT")) %>% arrange(desc(V4))
