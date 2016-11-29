

#source("https://bioconductor.org/biocLite.R")
#biocLite("seqLogo")

source("damagelogo.R")
source("pwm.R")

library(seqLogo)
mFile <- system.file("Exfiles/pwm1", package="seqLogo")
m <- read.table(mFile)
p <- makePWM(m)
p
mat1 <- cbind(p@pwm[,c(1,2)], rep(0,4), p@pwm[,c(7,8)]);
colnames(mat1) <- c("-2", "-1", "0", "1", "2")
mat2 <- cbind(rep(0,6), rep(0,6),
              c(0.5, 0.2, 0.2, 0.05, 0.05, 0),
              rep(0,6), rep(0,6))
rownames(mat2) <- c("C->T", "C->A", "C->G",
                    "T->A", "T->C", "T->G")

pwm1 <- rbind(mat1, mat2)
colSums(pwm1)

mat1 <- cbind(p@pwm[,c(3,4)], rep(0,4), p@pwm[,c(5,6)]);
colnames(mat1) <- c("-2", "-1", "0", "1", "2")
mat2 <- cbind(rep(0,6), rep(0,6),
              c(0.3, 0.2, 0.2, 0.1, 0.1, 0.1),
              rep(0,6), rep(0,6))
rownames(mat2) <- c("C->T", "C->A", "C->G",
                    "T->A", "T->C", "T->G")

pwm2 <- rbind(mat1, mat2)
colSums(pwm2)

mat1 <- cbind(p@pwm[,c(3,4)], rep(0,4), p@pwm[,c(7,8)]);
colnames(mat1) <- c("-2", "-1", "0", "1", "2")
mat2 <- cbind(rep(0,6), rep(0,6),
              c(0.4, 0.4, 0.05, 0.05, 0.05, 0.05),
              rep(0,6), rep(0,6))
rownames(mat2) <- c("C->T", "C->A", "C->G",
                    "T->A", "T->C", "T->G")

pwm3 <- rbind(mat1, mat2)
colSums(pwm3)

mat1 <- cbind(p@pwm[,c(1,3)], rep(0,4), p@pwm[,c(5,6)]);
colnames(mat1) <- c("-2", "-1", "0", "1", "2")
mat2 <- cbind(rep(0,6), rep(0,6),
              c(0.4, 0.1, 0, 0, 0.25, 0.25),
              rep(0,6), rep(0,6))
rownames(mat2) <- c("C->T", "C->A", "C->G",
                    "T->A", "T->C", "T->G")

pwm4 <- rbind(mat1, mat2)
colSums(pwm4)




pwm <- list()
pwm[[1]] <- pwm1
pwm[[2]] <- pwm2
pwm[[3]] <- pwm3
pwm[[4]] <- pwm4


ic.scale=FALSE
xaxis=TRUE
yaxis=TRUE
xaxis_fontsize=5
xlab_fontsize=15
y_fontsize=15

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

damage_ic <- damage.ic(pwm)

damageLogo.skeleton(pwm[[1]], ic = damage_ic[,1], ic.scale = TRUE)
damageLogo.skeleton(pwm[[2]], ic = damage_ic[,2], ic.scale = TRUE)
damageLogo.skeleton(pwm[[3]], ic = damage_ic[,3], ic.scale = TRUE)
damageLogo.skeleton(pwm[[4]], ic = damage_ic[,4], ic.scale = TRUE)


topics_clus <- get(load("../rda/CountClust_output_Lindo2016_with_C_to_T.rda"));
theta <- topics_clus$clust_4$theta

damageLogo(theta)





