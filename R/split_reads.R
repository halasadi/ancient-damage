

#########  split the reads in a mix ancient - modern samples  ##########


reads_KO1 <- read.csv(file = "../utilities/reads_data/KO1.realigned.sorted.q30.csv",
                          header=FALSE)
d <- sample(1:dim(reads_KO1)[1], dim(reads_KO1)[1], replace=FALSE)
chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))

out <- chunk(d, 10)

for(num in 1:10){
  reads_1 <- reads_KO1[out[[num]],]
  write.csv(reads_1,
            file = paste0("../utilities/reads_data/contam_data/KO1_split_", num, ".csv"),
            row.names = FALSE)
}

library(aRchaic)

mff <- aRchaic::aRchaic_pool(folders = "../utilities/reads_data/contam_data/")
file <- read.csv(file = paste0("../utilities/reads_data/contam_data/KO1_split_", 1, ".csv"), header=TRUE)
file <- file[,-7]
colnames(file) <- c("mut", "leftpos", "rightpos", "lsb", "rsb", "strand")

breaks <- c(-1, seq(1,20,1))
tab <- tbl_df(file) %>% group_by(mut, leftpos, rightpos, lsb, rsb, strand) %>% summarize(n= n())
tab2 <- data.frame(tab)

mff <- aRchaic_pool(folders = "../utilities/reads_data/contam_data/")


file <- read.csv(file = paste0("../utilities/reads_data/contam_data/KO1_split_", 1, ".csv"))
file <- file[,-7]


tmp <- aggregate_signature_counts("../utilities/reads_data/contam_data/")

