

##########  pooled Damage signature plots  #######################
#  this is the interface connecting damage logo plots
#  to the theta matrix of proportions of different signnatures 
#  as would be obtained from a topic model output.

# theta_mat: matrix obtained from a GoM model output 

source("damagelogo.skeleton.R")

damageLogo <- function(theta,
                       sig_names = NULL,
                       ic.scale=TRUE,
                       yscale_change = TRUE,
                       xaxis=TRUE, 
                       yaxis=TRUE, 
                       xaxis_fontsize=10,
                       xlab_fontsize=15,
                       y_fontsize=15,
                       mut_width=2,
                       start=0.0001,
                       pop_names=paste0("Cluster ",1:dim(theta)[2])){
  if(is.null(sig_names))
    sig_names <- rownames(theta)

  sig_split <- do.call(rbind, 
                     lapply(sig_names, 
                            function(x) strsplit(as.character(x), split="")[[1]]))

  ncol_sig <- ncol(sig_split)
  flanks <- (ncol_sig - 4)/2;

   if(flanks%%1 != 0){
     stop("flanking bases not evenly distributed")
   }
  
  
   sub_pattern <- sapply(1:dim(sig_split)[1], 
                        function(x) paste(sig_split[x,(flanks+1):(flanks+4)], collapse=""))

   new_sig_split <- cbind(sig_split[,1:flanks], sub_pattern, sig_split[,((ncol_sig - flanks +1):ncol_sig)])
   colnames(new_sig_split) = NULL
   
   prop_patterns_list <- list()
   
   for(l in 1:dim(theta)[2]){
     prop_patterns_list[[l]] <- numeric();
      for(j in 1:ncol(new_sig_split)){
           temp <- tapply(theta[,l], factor(new_sig_split[,j], levels=c("A", "C", "G", "T", 
                                                        "C->T", "C->A", "C->G", 
                                                        "T->A", "T->C", "T->G")), sum)
   
           temp[is.na(temp)]=0
           prop_patterns_list[[l]] <- cbind(prop_patterns_list[[l]], temp)
      }
   }
   
   ic <- damage.ic(prop_patterns_list)
  
   grob_list <- list()
   for(l in 1:length(prop_patterns_list)){
    damageLogo.skeleton(prop_patterns_list[[l]], 
                        ic = ic[,l], 
                        ic.scale = ic.scale,
                        yscale_change = yscale_change,
                        xaxis=xaxis, 
                        yaxis=yaxis, 
                        xaxis_fontsize=xaxis_fontsize,
                        xlab_fontsize=xlab_fontsize,
                        y_fontsize=y_fontsize,
                        mut_width=mut_width,
                        start=start,
                        pop_name = pop_names[l])
     
   }
}
   
pwm2ic<-function(pwm) {
  npos<-ncol(pwm)
  ic<-numeric(length=npos)
  for (i in 1:npos) {
    ic[i]<- log(nrow(pwm), base=2) + sum(sapply(pwm[, i], function(x) { 
      if (x > 0) { x*log2(x) } else { 0 }
    }))
  }    
  ic
}


damage.ic<-function(pwm) {
  npos<-ncol(pwm[[1]])
  ic<- matrix(0, npos, length(pwm))
  
  for(i in 1:npos){
    mat <- numeric()
    for(j in 1:length(pwm)){
      mat <- cbind(mat, pwm[[j]][,i])
    }
    mat_clean <- mat[rowSums(mat) != 0,]
    ic[i,] <- pwm2ic(mat_clean)
  }
  
  return(ic)
}


   
   