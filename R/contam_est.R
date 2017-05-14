

###########  contamination studies across board   #######################

model <- get(load("../utilities/Gosling_moderns/clus_2/model.rda"))
omega <- model$omega

contamination_metadata <- read.csv("../utilities/moderns_Pinhasi/contamination_metadata.csv", header=TRUE)

omega1 <- omega[51:62,]
rownames(omega1) <- c("IR1", "KO1", "KO2", "NE1", "NE2", "NE3", "NE4", "NE5", "NE6", "NE7", "BR1", "BR2")


contamination_metadata_filtered <- contamination_metadata[match(rownames(omega1),contamination_metadata[,1]),]

plot(omega1[,1], contamination_metadata_filtered[,6], pch=20, col="red")




###################  moderns + Gosling  ##############################


model <- get(load("../utilities/Gosling_moderns/clus_2/model.rda"))
omega <- model$omega
omega1 <- omega[1:54,]
omega1_names <- rownames(omega1)

names2 <- as.character(sapply(omega1_names, function(x) return(strsplit(strsplit(x, "-")[[1]][3], "[.]")[[1]][1])))

metadata <- read.csv("../utilities/Gosling_moderns/gosling_metadata.csv", header=TRUE)

metadata1 <- metadata[match(names2, as.character(metadata[,1])), ]

omega1_filtered <- omega1[-unique(which(is.na(metadata1), arr.ind=TRUE)[,1]),]
metadata1_filtered  <- metadata1[-unique(which(is.na(metadata1), arr.ind=TRUE)[,1]),]

omega2_filtered <- omega1_filtered[-c(22,23,24,30,31,35),]
metadata2_filtered <- metadata1_filtered[-c(22,23,24,30,31,35),]


plot(omega2_filtered[,1], metadata2_filtered$Cont., pch=20)
