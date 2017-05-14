

modern_fu <- get(load("../utilities/moderns_Fu/clus_4/model.rda"))
modern_fu$omega[grep("Pavlov1", rownames(modern_fu$omega)),]
dim(modern_fu$omega)
