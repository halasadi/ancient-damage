

#################  pattern plot across the read length ########################
#######  tracking how a damage signature/pattern changes across a read ############

pattern_plot <- function(file,
                         pattern,
                         plot_type=c("left", "right", "both"))
{
  data <- read.csv(file, header=FALSE)
  pattern_data <- data[grep(pattern, data[,1]),]
  pattern_data[which(pattern_data[,3]==-1), 3] <- 0
  
  pattern_data <- data[grep(pattern, data[,1]),]
  
  pattern_data[which(pattern_data[,3]==-1), 3] <- 0
  

  tab_pattern_left <- tapply(pattern_data[,4], pattern_data[,2], sum);
  
  lo_left <- loess(as.numeric(tab_pattern_left)~as.numeric(names(tab_pattern_left)))
  
  tab_pattern_right <- tapply(pattern_data[,4], pattern_data[,3], sum);
  
  lo_right <- loess(as.numeric(tab_pattern_right)~as.numeric(names(tab_pattern_right)))
  
  if(plot_type=="left"){
  par(mfrow=c(1,1))
  plot(predict(lo_left), col='red', lwd=2, type="l", 
       ylab=paste0("predicted no. of damages: ", pattern),
       xlab="read positions (from left)",
       main=paste0("Damage plot across reads: ", pattern) 
  )}
  
  if(plot_type=="right"){
    par(mfrow=c(1,1))
    plot(predict(lo_right), col='red', lwd=2, type="l", 
         ylab=paste0("predicted no. of damages: ", pattern),
         xlab="read positions (from right)",
         main=paste0("Damage plot across reads: ", pattern), 
         xlim=rev(range(as.numeric(names(tab_pattern_right)))))
  }
  
  if(plot_type=="both"){
    par(mfrow=c(1,2))
    plot(predict(lo_left), col='red', lwd=2, type="l", 
         ylab=paste0("predicted no. of damages: ", pattern),
         xlab="read positions (from left)",
         main=paste0("Damage plot across reads: ", pattern) 
    )
    plot(predict(lo_right), col='red', lwd=2, type="l", 
         ylab=paste0("predicted no. of damages: ", pattern),
         xlab="read positions (from right)",
         main=paste0("Damage plot across reads: ", pattern), 
         xlim=rev(range(as.numeric(names(tab_pattern_right)))))
  }
}

pattern_plot(file="../summary_data/I0026.q30.csv",
             pattern="C->T",
             plot_type="left")

