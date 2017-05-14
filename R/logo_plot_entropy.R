

################   Fixing the information content   #########################


topic_clus <- get(load("../utilities/Nepal_moderns/clus_2/model.rda"))
theta_pool <- topic_clus$theta

sig_names = NULL
ic.scale=TRUE
max_pos = 20
flanking_bases=1
yscale_change = TRUE
xaxis=TRUE
yaxis=TRUE
xlab = " "
xaxis_fontsize=5
xlab_fontsize=10
title_aligner = 18
y_fontsize=10
title_fontsize = 20
mut_width=2
start=0.0001
renyi_alpha = 1
pop_names=paste0("Cluster ",1:dim(theta_pool)[2])



signature_set <- rownames(theta_pool)
signature_patterns <- substring(signature_set, 1, 4+2*flanking_bases)
library(dplyr)
theta <- dplyr::tbl_df(data.frame(theta_pool)) %>% dplyr::mutate(sig = signature_patterns) %>% dplyr::group_by(sig) %>% dplyr::summarise_each(funs(sum)) %>% as.data.frame()
rownames(theta) <-  theta[,1]
theta <- theta[,-1]

indices_minus <- grep("_-_", signature_set)
strand_theta <- data.frame("minus" = colSums(theta_pool[indices_minus,]),
                           "plus" = colSums(theta_pool[-indices_minus,]))

if(flag == 1){
  strand_theta <- data.frame("minus" = colSums(matrix(theta_pool[indices_minus,])),
                             "plus" = colSums(matrix(theta_pool[-indices_minus,])))
  strand_theta <- strand_theta/2;
}
breakbase <- substring(signature_set, 8+2*flanking_bases,  8+2*flanking_bases)

theta_break <- dplyr::tbl_df(data.frame(theta_pool)) %>% dplyr::mutate(sig = breakbase) %>% dplyr::group_by(sig) %>% dplyr::summarise_each(funs(sum)) %>% as.data.frame()
rownames(theta_break) <- theta_break[,1]
theta_break <- theta_break[,-1]

theta_break <- theta_break[match(c("A", "C", "G", "T"), rownames(theta_break)),]
breaks_theta <- theta_break


if(is.null(sig_names))
  sig_names <- rownames(theta)

prob_mutation <- filter_signatures_only_location(t(theta_pool), max_pos = max_pos, flanking_bases = flanking_bases)
prob_mutation <- t(apply(prob_mutation, 1, function(x) {
  y <- x[!is.na(x)];
  return(y/sum(y))
}))
max_prob <- max(prob_mutation);

sig_split <- do.call(rbind,
                     lapply(sig_names,
                            function(x) strsplit(as.character(x), split="")[[1]][1:(4+2*flanking_bases)]))

ncol_sig <- (4+2*flanking_bases)

if(flanking_bases%%1 != 0){
  stop("flanking bases not evenly distributed")
}


sub_pattern <- sapply(1:dim(sig_split)[1],
                      function(x) paste(sig_split[x,(flanking_bases+1):(flanking_bases+4)], collapse=""))

new_sig_split <- cbind(sig_split[,1:flanking_bases], sub_pattern, sig_split[,((ncol_sig - flanking_bases +1):ncol_sig)])
colnames(new_sig_split) = NULL

prop_patterns_list <- list()


for(l in 1:dim(theta)[2]){
  prop_patterns_list[[l]] <- numeric();
  for(j in 1:ncol(new_sig_split)){
    temp <- tapply(theta[,l], factor(new_sig_split[,j], levels=c("A", "C", "G", "T", "X",
                                                                 "C->T", "C->A", "C->G",
                                                                 "T->A", "T->C", "T->G")), sum)

    temp[is.na(temp)]=0
    prop_patterns_list[[l]] <- cbind(prop_patterns_list[[l]], temp)
  }
}

ic <- damage.ic(prop_patterns_list, alpha=renyi_alpha)

mat <- prop_patterns_list[[1]]
i <- 1
alpha = 3000
log(length(which(mat[,i]!=0.00)), base=2) - (1/(1-alpha))* log (sum(mat[,i]^{alpha}), base=2)

i <- 1
ll <- 6
(1/(1-alpha))* log (sum((1/ll)^{alpha}), base=2) - (1/(1-alpha))* log (sum(mat[,i]^{alpha}), base=2)
