
############   Create  a reads matrix for prediction    ###################


reads_Pinhasi <- read.csv(file = "../data/Pinhasi/NE2.realigned.sorted.q30.csv",
                          header=FALSE)


reads_moderns <- read.csv(file = "../data/moderns_lite/NA19795.mapped.ILLUMINA.bwa.MXL.low_coverage.20130415.realigned.sorted.q30.csv",
                          header=FALSE)

model <- get(load("../utilities/moderns_Pinhasi/clus_2/model.rda"))

output <- read_memberships(model,
                             reads_moderns,
                             subset = 10000,
                             nostrand = TRUE)

median(as.numeric(output$`cl- 1`))
median(as.numeric(output$`cl- 2`))
plot(density((as.numeric(output$`cl- 1`))))
plot(density((as.numeric(output$`cl- 2`))))

###################      Exploratory analysis       ###########################

theta <- model$theta
rownames <- rownames(theta)
topic_break <- theta_breakdown(theta)

reads_collection <- reads_moderns

read_id <- as.list(unique(reads_collection[,7])[1:10000])
table <- list()

cl <- parallel::makeCluster(parallel::detectCores(),type=ifelse(.Platform$OS.type=="unix","FORK","PSOCK"))
print(cl)

read_member_mat <- as.data.frame(do.call(rbind, parLapply(cl, read_id, read_probability)))
colnames(read_member_mat) <- c(paste("cl-", 1:dim(theta)[2]), "sig", "ID")


system.time(
for(l in 1:length(read_id)){
  idx <- which(reads_collection[, 7] == read_id[l])
  reads <- reads_collection[idx,]
  dist <-  as.numeric(apply(as.matrix(reads), 1, function(x) return(min(x[2], x[3]))))
  modf <- cbind.data.frame(paste0(reads[,1], "_", reads[,6], "_", reads[,4], "_", reads[,5], "_", dist))
  colnames(modf) <- "signature"
  modf2 <- as.character(modf[,1])
  symbol_vec <- signatureclub2(modf2, flanking_bases =  1)

  sym_prob <- array(1, length(topic_break))
  for(v in 1:length(symbol_vec)){
    sym_prob <- tryCatch(sym_prob * get_prob_symbol(topic_break, as.character(symbol_vec[v])),
                         error = function(e){sym_prob <- sym_prob})
  }
  sym_prob <- sym_prob/sum(sym_prob)
  sym_prob <- round(sym_prob, round_off)
  sym_prob <- sym_prob/sum(sym_prob)

  sym_prob_2 <- c(sym_prob, paste0(symbol_vec, collapse = " ; "))
  table[[l]] <- sym_prob_2
})


read_probability <- function(read_x){
  idx <- which(reads_collection[, 7] == read_x)
  reads <- reads_collection[idx,]
  dist <-  as.numeric(apply(as.matrix(reads), 1, function(x) return(min(x[2], x[3]))))
  modf <- cbind.data.frame(paste0(reads[,1], "_", reads[,6], "_", reads[,4], "_", reads[,5], "_", dist))
  colnames(modf) <- "signature"
  modf2 <- as.character(modf[,1])
  symbol_vec <- signatureclub2(modf2, flanking_bases =  1)

  sym_prob <- array(1, length(topic_break))
  for(v in 1:length(symbol_vec)){
    sym_prob <- tryCatch(sym_prob * get_prob_symbol(topic_break, as.character(symbol_vec[v])),
                         error = function(e){sym_prob <- sym_prob})
  }
  sym_prob <- sym_prob/sum(sym_prob)
  sym_prob <- round(sym_prob, round_off)
  sym_prob <- sym_prob/sum(sym_prob)

  sym_prob_2 <- c(sym_prob, paste0(symbol_vec, collapse = " ; "), read_x)
  return(sym_prob_2)
}


read <- reads_Pinhasi[1,]
modf <- cbind.data.frame(paste0(read[,1], "_", read[,6], "_", read[,4], "_", read[,5], "_", min(read[,2], read[,3])))
colnames(modf) <- "signature"
modf2 <- as.character(modf[,1])
symbol_vec <- signatureclub2(modf2, flanking_bases =  1)

theta <- model$theta
rownames <- rownames(theta)
topic_break <- theta_breakdown(theta)
sym_prob <- array(1, length(topic_break))
for(v in 1:length(symbol_vec)){
  sym_prob <- sym_prob * get_prob_symbol(topic_break, as.character(symbol_vec[v]))
}
sym_prob <- sym_prob/sum(sym_prob)
sym_prob <- round(sym_prob, round_off)
sym_prob <- sym_prob/sum(sym_prob)

sym_prob_2 <- c(sym_prob, paste0(symbol_vec, collapse = " ; "))



if(read[,6] == "+"){
  pos <- read[2]
  string <- read[1]

}



file = "../data/moderns_lite/NA19795.mapped.ILLUMINA.bwa.MXL.low_coverage.20130415.realigned.sorted.q30.csv"
