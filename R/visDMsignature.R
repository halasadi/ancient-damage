

visDMSignature <- function(theta,
                           pattern,
                           sample_names,
                           numBases=5,
                           flankingBasesNum,
                           trDir,
                           fdim,
                           title_size=10,
                           panel_title_size=10,
                           panel_title_font=4,
                           layout=c(6,1)){
  
  graphList <- list(mode="vector")
  
  index <- c("C->A", "C->G", "C->T", "T->A", "T->C", "T->G")
  
  prop_clus_sig <- do.call(cbind, lapply(1:dim(theta)[2], function(l) return(tapply(theta[,l], sub_pattern, sum))))
  
  for(num in 1:length(index)){
    theta_filt <- matrix(0, dim(theta)[1], dim(theta)[2])
    theta_filt[which(sub_pattern==index[num]),] <- theta[which(sub_pattern==index[num]),];
    theta_filt_norm <- apply(theta_filt, 2, function(x) return(normalize(x)))
    
    isBG <- TRUE
    K <- dim(theta)[2]
    F <- array(0, c(K, length(fdim), max(fdim)))
    for(k in 1:K){
      for(kk in 1:length(fdim)){
        temp <- tapply(theta_filt_norm[,k], mut_features_mat[,kk], sum)
        F[k,kk,as.numeric(names(temp))] <- as.numeric(temp)
      }
    }
    
    Param <- new(Class = "EstimatedParameters", 
                 type = type,
                 flankingBasesNum = flankingBasesNum,
                 transcriptionDirection = trDir,
                 possibleFeatures = as.integer(fdim),
                 sampleList = sample_names,
                 signatureNum = as.integer(K),
                 isBackGround = isBG,
                 signatureFeatureDistribution = F,
                 sampleSignatureDistribution = omega,
                 loglikelihood = -100000)
    
    graphList[[num]] <- list(mode="vector")
    for(k in 1:dim(theta)[2]){
      b <- pmsignature::visPMSignature(Param,k);
      b <- b +ggplot2::ggtitle(paste("prop frequency: ", round(prop_clus_sig[num,k],2))) + theme(plot.title = element_text(size = 8, face = "bold"))
      graphList[[num]][[k]] <- b
    }
  }
      
  for(k in 1:dim(theta)[2]){
    graphviz_list <- lapply(1:length(index), function(n) graphList[[n]][[k]])
    
    
    library(grid)
    library(gridExtra)
    do.call("grid.arrange", 
            args = list(grobs=graphviz_list,
                        ncol = layout[2],
                        nrow = layout[1],
                        top=textGrob(paste0("Cluster:", k),                                gp=gpar(fontsize=15,                                         font=3))))
  }
}