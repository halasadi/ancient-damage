
dir <- "../data/Lindo2016ancients/";

tab <- read.csv("../data/Lindo2016ancients/181_all_chr.q30.csv", header=FALSE)

dim(tab)

head(tab)

read_length <- tab$V2 + tab$V3

indices <- grep("C->T", tab$V1)
pos_CtoT <- pmin(tab[indices, ]$V2, tab[indices, ]$V3)
read_length_CtoT <- tab[indices, ]$V2 + tab[indices, ]$V3

mean(read_length)
mean(pos_CtoT)

plot(table(read_length), type="o", col="red", main="read length distributions (all vs CtoT)")
points(table(read_length_CtoT), type="o", col="blue")
legend("topleft", fill=c("red", "blue"), legend = c("all", "CtoT"))

indices <- grep("C->T", tab$V1)
indices <- c(indices, grep("G->A", tab$V1))

table(tab$V2)

table(tab$V3)

indices2 <- which(tab$V2 < 5 | tab$V3 < 5)
length(indices2)

length(grep("G->A", tab$V1))

length(grep("C->T", tab$V1))

indices_matched <- intersect(indices, indices2)

length(indices_matched)
read_length_3 <- (tab[indices_matched, ]$V2 + tab[indices_matched, ]$V3)
table(read_length_3)

plot(table(read_length)/sum(table(read_length)), type="o", col="red", 
     main="read length distributions (all vs CtoT)")
points(table(read_length_CtoT)/sum(table(read_length_CtoT)), type="o", col="green")
points(table(read_length_3)/ sum(table(read_length_3)), type="o", col="blue")
legend("topleft", fill=c("red", "green", "blue"), legend = c("all", "CtoT", "CtoT < 5"))

files <- list.files("../data/Lindo2016ancients/")
dir <- "../data/Lindo2016ancients/"
tab <- list()
for(l in files){
    tab[[l]] <- read.csv(paste0(dir, l), header=FALSE)
    cat("We are at sample", l, "\n")
}

ancient_names <- as.character(sapply(files, function(x) return(strsplit(x, "[_]")[[1]][1])))

par(mfrow=c(3,3))
for(l in 1:length(tab)){
    tab1 <- tab[[l]]
    read_length <- tab1$V2 + tab1$V3
    indices <- grep("C->T", tab1$V1)
    indices <- c(indices, grep("G->A", tab1$V1))
    read_length_CtoT <- tab1[indices, ]$V2 + tab1[indices, ]$V3 ## read lengths for all C to T
    indices2 <- which(tab1$V2 < 5 | tab1$V3 < 5)
    indices_matched <- intersect(indices, indices2)
    read_length_3 <- (tab1[indices_matched, ]$V2 + tab1[indices_matched, ]$V3)
    plot(table(read_length)/sum(table(read_length)), type="o", col="red", 
     main=paste0(ancient_names[l]), xlab="read pos", ylab="prop of occur")
    points(table(read_length_CtoT)/sum(table(read_length_CtoT)), type="o", col="green")
    points(table(read_length_3)/ sum(table(read_length_3)), type="o", col="blue")
    legend("topleft", fill=c("red", "green", "blue"), legend = c("all", "CtoT", "CtoT < 5"))
    cat("we are at sample", l, "\n")
}

files <- list.files("../data/Lindo2016moderns/")
dir <- "../data/Lindo2016moderns/"
tab <- list()
for(l in files){
    tab[[l]] <- read.csv(paste0(dir, l), header=FALSE)
    cat("We are at sample", l, "\n")
}

modern_names <- as.character(sapply(files, function(x) return(strsplit(x, "[_]")[[1]][1])))

par(mfrow=c(3,3))
for(l in 1:length(tab)){
    tab1 <- tab[[l]]
    read_length <- tab1$V2 + tab1$V3
    indices <- grep("C->T", tab1$V1)
    indices <- c(indices, grep("G->A", tab1$V1))
    read_length_CtoT <- tab1[indices, ]$V2 + tab1[indices, ]$V3 ## read lengths for all C to T
    indices2 <- which(tab1$V2 < 5 | tab1$V3 < 5)
    indices_matched <- intersect(indices, indices2)
    read_length_3 <- (tab1[indices_matched, ]$V2 + tab1[indices_matched, ]$V3)
    plot(table(read_length)/sum(table(read_length)), type="o", col="red", 
     main=paste0(modern_names[l]), xlab="read pos", ylab="prop of occur")
    points(table(read_length_CtoT)/sum(table(read_length_CtoT)), type="o", col="green")
    points(table(read_length_3)/ sum(table(read_length_3)), type="o", col="blue")
    legend("topleft", fill=c("red", "green", "blue"), legend = c("all", "CtoT", "CtoT < 5"))
    cat("we are at sample", l, "\n")
}

files <- list.files("../data/Sardinia2017/")
dir <- "../data/Sardinia2017/"
tab <- list()
for(l in files){
    tab[[l]] <- read.csv(paste0(dir, l), header=FALSE)
    cat("We are at sample", l, "\n")
}

par(mfrow=c(2,3))
for(l in 1:length(tab)){
    tab1 <- tab[[l]]
    read_length <- tab1$V2 + tab1$V3
    indices <- grep("C->T", tab1$V1)
    indices <- c(indices, grep("G->A", tab1$V1))
    read_length_CtoT <- tab1[indices, ]$V2 + tab1[indices, ]$V3 ## read lengths for all C to T
    indices2 <- which(tab1$V2 < 5 | tab1$V3 < 5)
    indices_matched <- intersect(indices, indices2)
    read_length_3 <- (tab1[indices_matched, ]$V2 + tab1[indices_matched, ]$V3)
    plot(table(read_length)/sum(table(read_length)), type="o", col="red", 
     main="read length dist", xlab="read pos", ylab="prop of occur")
    points(table(read_length_CtoT)/sum(table(read_length_CtoT)), type="o", col="green")
    points(table(read_length_3)/ sum(table(read_length_3)), type="o", col="blue")
    legend("topleft", fill=c("red", "green", "blue"), legend = c("all", "CtoT", "CtoT < 5"))
    cat("we are at sample", l, "\n")
}

files <- list.files("../data/AnnaGosling2016data/")
dir <- "../data/AnnaGosling2016data/"
tab <- list()
for(l in files){
    tab[[l]] <- read.csv(paste0(dir, l), header=FALSE)
    cat("We are at sample", l, "\n")
}

names <- as.character(sapply(files, function(x) return(strsplit(x, c("[-]"))[[1]][3])))

names

par(mfrow=c(3,3))
for(l in 1:length(tab)){
    tab1 <- tab[[l]]
    read_length <- tab1$V2 + tab1$V3
    indices <- grep("C->T", tab1$V1)
    indices <- c(indices, grep("G->A", tab1$V1))
    read_length_CtoT <- tab1[indices, ]$V2 + tab1[indices, ]$V3 ## read lengths for all C to T
    indices2 <- which(tab1$V2 < 5 | tab1$V3 < 5)
    indices_matched <- intersect(indices, indices2)
    read_length_3 <- (tab1[indices_matched, ]$V2 + tab1[indices_matched, ]$V3)
    plot(table(read_length)/sum(table(read_length)), type="o", col="red", 
     main=paste0(names[l]), xlab="read pos", ylab="prop of occur")
    points(table(read_length_CtoT)/sum(table(read_length_CtoT)), type="o", col="green")
    points(table(read_length_3)/ sum(table(read_length_3)), type="o", col="blue")
    legend("topleft", fill=c("red", "green", "blue"), legend = c("all", "CtoT", "CtoT < 5"), cex=0.5)
    cat("we are at sample", l, "\n")
}

files <- list.files("../data/HGDPmoderns/")
dir <- "../data/HGDPmoderns/"
tab <- list()
for(l in files){
    tab[[l]] <- read.csv(paste0(dir, l), header=FALSE)
    cat("We are at sample", l, "\n")
}

par(mfrow=c(2,3))
for(l in 1:length(tab)){
    tab1 <- tab[[l]]
    read_length <- tab1$V2 + tab1$V3
    indices <- grep("C->T", tab1$V1)
    indices <- c(indices, grep("G->A", tab1$V1))
    read_length_CtoT <- tab1[indices, ]$V2 + tab1[indices, ]$V3 ## read lengths for all C to T
    indices2 <- which(tab1$V2 < 5 | tab1$V3 < 5)
    indices_matched <- intersect(indices, indices2)
    read_length_3 <- (tab1[indices_matched, ]$V2 + tab1[indices_matched, ]$V3)
    plot(table(read_length)/sum(table(read_length)), type="o", col="red", 
     main="read length dist", xlab="read pos", ylab="prop of occur")
    points(table(read_length_CtoT)/sum(table(read_length_CtoT)), type="o", col="green")
    points(table(read_length_3)/ sum(table(read_length_3)), type="o", col="blue")
    legend("topleft", fill=c("red", "green", "blue"), legend = c("all", "CtoT", "CtoT < 5"), cex=0.5)
    cat("we are at sample", l, "\n")
}

files <- sample(list.files("../data/1000Gmoderns/"))[1:50]
dir <- "../data/1000Gmoderns/"
tab <- list()
for(l in files){
    tab[[l]] <- read.csv(paste0(dir, l), header=FALSE)
    cat("We are at sample", l, "\n")
}

modern_names <- as.character(sapply(files, function(x) return(strsplit(x, c("[.]"))[[1]][5])))

par(mfrow=c(2,3))
for(l in 1:length(tab)){
    tab1 <- tab[[l]]
    read_length <- tab1$V2 + tab1$V3
    indices <- grep("C->T", tab1$V1)
    indices <- c(indices, grep("G->A", tab1$V1))
    read_length_CtoT <- tab1[indices, ]$V2 + tab1[indices, ]$V3 ## read lengths for all C to T
    indices2 <- which(tab1$V2 < 5 | tab1$V3 < 5)
    indices_matched <- intersect(indices, indices2)
    read_length_3 <- (tab1[indices_matched, ]$V2 + tab1[indices_matched, ]$V3)
    plot(table(read_length)/sum(table(read_length)), type="o", col="red", 
     main=paste0(modern_names[l]), xlab="read pos", ylab="prop of occur")
    points(table(read_length_CtoT)/sum(table(read_length_CtoT)), type="o", col="green")
    points(table(read_length_3)/ sum(table(read_length_3)), type="o", col="blue")
    legend("topleft", fill=c("red", "green", "blue"), legend = c("all", "CtoT", "CtoT < 5"), cex=0.5)
    cat("we are at sample", l, "\n")
}


files <- list.files("../data/RISE_data/")
dir <- "../data/RISE_data/"
tab <- list()
for(l in files){
    tab[[l]] <- read.csv(paste0(dir, l), header=FALSE)
    cat("We are at sample", l, "\n")
}

files <- list.files("../data/RISE_data/")
length(files)

par(mfrow=c(2,3))
for(l in 1:length(tab)){
    tab1 <- tab[[l]]
    read_length <- tab1$V2 + tab1$V3
    indices <- grep("C->T", tab1$V1)
    indices <- c(indices, grep("G->A", tab1$V1))
    read_length_CtoT <- tab1[indices, ]$V2 + tab1[indices, ]$V3 ## read lengths for all C to T
    indices2 <- which(tab1$V2 < 5 | tab1$V3 < 5)
    indices_matched <- intersect(indices, indices2)
    read_length_3 <- (tab1[indices_matched, ]$V2 + tab1[indices_matched, ]$V3)
    plot(table(read_length)/sum(table(read_length)), type="o", col="red", 
     main="read length distribution", xlab="read pos", ylab="prop of occur")
    points(table(read_length_CtoT)/sum(table(read_length_CtoT)), type="o", col="green")
    points(table(read_length_3)/ sum(table(read_length_3)), type="o", col="blue")
    legend("topleft", fill=c("red", "green", "blue"), legend = c("all", "CtoT", "CtoT < 5"), cex=0.5)
    cat("we are at sample", l, "\n")
}
