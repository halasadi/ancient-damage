

##########  R function to read damage files  ############################

damage.build.counts =  function(file,
                                breaks=NULL,
                                type=1)
{
  file <-  read.csv(file=file);
  file[which(file[,3]==-1), 3] <- 0
  
  min_dist_from_end <- apply(file[,2:3], 1, function(x) return(min(x)))
  
  
  if(is.null(breaks)){
    bins <- c(-1, 5, 10, 15, 20, 30, 40, max(min_dist_from_end))
  }else{
    message("one of the breaks values is beyond the range of read lengths: adjusting")
    bins <- c(intersect(breaks, (-1):max(min_dist_from_end)),max(min_dist_from_end)) 
  }
  
  bin_values <- .bincode(min_dist_from_end, bins, TRUE)
  
  modified_file <- cbind.data.frame(file[,1], file[,4], bin_values)
  colnames(modified_file) <- c("pattern", "counts", "bin_values")
  
  library(plyr)
  df1 <- ddply(modified_file, .(pattern, bin_values), summarise, newvar = sum(counts))
  colnames(df1) <- c("pattern", "bin", "counts")
  
  if(type==2){
    df2 <- cbind.data.frame(paste0(df1[,1], "_", df1[,2]), df1[,3])
    colnames(df2) <- c("pattern-bins", "counts")
    out <- df2
  }else{
    out <- df1
  }
  
  return(out)
}


# test example 
## out <- damage.build.counts(file="../summary_data/I0026.q30.csv")
