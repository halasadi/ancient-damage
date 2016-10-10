
###################################################################
#####################  Damage Logos  ##############################
###################################################################



###################  letter  A   #################################


letterA <- function(x.pos,y.pos,ht,wt,id=NULL){
  
  x <- c(0,4,6,10,8,6.8,3.2,2,0,3.6,5,6.4,3.6)
  y <- c(0,10,10,0,0,3,3,0,0,4,7.5,4,4)
  x <- 0.1*x
  y <- 0.1*y
  
  x <- x.pos + wt*x
  y <- y.pos + ht*y
  
  if (is.null(id)){
    id <- c(rep(1,9),rep(2,4))
  }else{
    id <- c(rep(id,9),rep(id+1,4))
  }
  
  fill <- c("green","white")
  
  list(x=x,y=y,id=id,fill=fill)
}

#grid.newpage()

#grid.polygon(x=a_let$x, y=a_let$y, gp=gpar(fill=a_let$fill,
#             col="transparent"))

################   letter  T   ##################################

letterT <- function(x.pos,y.pos,ht,wt,id=NULL){
  
  x <- c(0,10,10,6,6,4,4,0)
  y <- c(10,10,9,9,0,0,9,9)
  x <- 0.1*x
  y <- 0.1*y
  
  x <- x.pos + wt*x
  y <- y.pos + ht*y
  
  if (is.null(id)){
    id <- rep(1,8)
  }else{
    id <- rep(id,8)
  }
  
  fill <- "red"
  
  list(x=x,y=y,id=id,fill=fill)
}

#################### letter  C #####################################

letterC <- function(x.pos,y.pos,ht,wt,id=NULL){
  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)
  
  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)
  
  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))
  
  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)
  
  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)
  
  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)
  
  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))
  
  x <- c(x,rev(x1))
  y <- c(y,rev(y1))
  
  x <- x.pos + wt*x
  y <- y.pos + ht*y
  
  if (is.null(id)){
    id <- rep(1,length(x))
  }else{
    id <- rep(id,length(x))
  }
  
  fill <- "blue"
  
  list(x=x,y=y,id=id,fill=fill)
}


##################  letter G  ######################################

letterG <- function(x.pos,y.pos,ht,wt,id=NULL){
  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)
  
  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)
  
  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))
  
  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)
  
  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)
  
  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)
  
  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))
  
  x <- c(x,rev(x1))
  y <- c(y,rev(y1))
  
  h1 <- max(y.l1)
  r1 <- max(x.l1)
  
  h1 <- 0.4
  x.add <- c(r1,0.5,0.5,r1-0.2,r1-0.2,r1,r1)
  y.add <- c(h1,h1,h1-0.1,h1-0.1,0,0,h1)
  
  
  
  if (is.null(id)){
    id <- c(rep(1,length(x)),rep(2,length(x.add)))
  }else{
    id <- c(rep(id,length(x)),rep(id+1,length(x.add)))
  }
  
  x <- c(rev(x),x.add)
  y <- c(rev(y),y.add)
  
  x <- x.pos + wt*x
  y <- y.pos + ht*y
  
  
  fill <- c("orange","orange")
  
  list(x=x,y=y,id=id,fill=fill)
  
}



################  letter C to T  ###############################

letter_C_to_T <- function(x.pos,y.pos,ht,wt,id=NULL){
  
  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)
  
  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)
  
  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))
  
  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)
  
  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)
  
  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)
  
  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))
  
  x <- c(x,rev(x1))
  y <- c(y,rev(y1))
  x1 <- 0.4*x
  y1 <- 1*y
  
  
  x <- c( 0.4, 0.4, 0, 0, 1, 1, 0.6, 0.6)
  y <- c( 0, 0.85, 0.85, 1, 1, 0.85, 0.85, 0)
  x2 <- 0.6 + 0.4*x
  y2 <- 1*y
  
  
  x3 <- c(0.42, 0.42, 0.55, 0.55, 0.60, 0.55, 0.55)
  y3 <- c(0.45, 0.55, 0.55, 0.60, 0.50, 0.40, 0.45)
  
  xpool <- c(x1, x2, x3)
  ypool <- c(y1, y2, y3)
  
  if(!is.null(id)){
    id_pool <- id + c(rep(1,length(x1)), rep(3, length(x2)), 
                      rep(4, length(x3)))
  }else{
    id_pool <- c(rep(1,length(x1)), rep(3, length(x2)), 
                 rep(4, length(x3)))
  }
  
  x <- x.pos + wt*xpool
  y <- y.pos + ht*ypool
  
  fill=c("blue","red","grey80")
  
  list(x=x,y=y,id=id_pool,fill=fill)
  
}

###############  letter C to G  ###################################

letter_C_to_G <- function(x.pos,y.pos,ht,wt,id=NULL){
  
    angle1 <- seq(0.3+pi/2,pi,length=100)
    angle2 <- seq(pi,1.5*pi,length=100)
    x.l1 <- 0.5 + 0.5*sin(angle1)
    y.l1 <- 0.5 + 0.5*cos(angle1)
    x.l2 <- 0.5 + 0.5*sin(angle2)
    y.l2 <- 0.5 + 0.5*cos(angle2)

    x.l <- c(x.l1,x.l2)
    y.l <- c(y.l1,y.l2)

    x <- c(x.l,rev(x.l))
    y <- c(y.l,1-rev(y.l))

    x.i1 <- 0.5 +0.35*sin(angle1)
    y.i1 <- 0.5 +0.35*cos(angle1)
    x.i1 <- x.i1[y.i1<=max(y.l1)]
    y.i1 <- y.i1[y.i1<=max(y.l1)]
    y.i1[1] <- max(y.l1)

    x.i2 <- 0.5 +0.35*sin(angle2)
    y.i2 <- 0.5 +0.35*cos(angle2)

    x.i <- c(x.i1,x.i2)
    y.i <- c(y.i1,y.i2)

    x1 <- c(x.i,rev(x.i))
    y1 <- c(y.i,1-rev(y.i))

    x <- c(x,rev(x1))
    y <- c(y,rev(y1))
    x1_pool <- 0.4*x
    y1_pool <- 1*y

    id1_pool <- rep(1,length(x1_pool))



    angle1 <- seq(0.3+pi/2,pi,length=100)
    angle2 <- seq(pi,1.5*pi,length=100)
    x.l1 <- 0.5 + 0.5*sin(angle1)
    y.l1 <- 0.5 + 0.5*cos(angle1)
    x.l2 <- 0.5 + 0.5*sin(angle2)
    y.l2 <- 0.5 + 0.5*cos(angle2)

    x.l <- c(x.l1,x.l2)
    y.l <- c(y.l1,y.l2)

    x <- c(x.l,rev(x.l))
    y <- c(y.l,1-rev(y.l))

    x.i1 <- 0.5 +0.35*sin(angle1)
    y.i1 <- 0.5 +0.35*cos(angle1)
    x.i1 <- x.i1[y.i1<=max(y.l1)]
    y.i1 <- y.i1[y.i1<=max(y.l1)]
    y.i1[1] <- max(y.l1)

    x.i2 <- 0.5 +0.35*sin(angle2)
    y.i2 <- 0.5 +0.35*cos(angle2)

    x.i <- c(x.i1,x.i2)
    y.i <- c(y.i1,y.i2)

    x1 <- c(x.i,rev(x.i))
    y1 <- c(y.i,1-rev(y.i))

    x <- c(x,rev(x1))
    y <- c(y,rev(y1))

    h1 <- max(y.l1)
    r1 <- max(x.l1)

    h1 <- 0.4
    x.add <- c(r1,0.5,0.5,r1-0.2,r1-0.2,r1,r1)
    y.add <- c(h1,h1,h1-0.1,h1-0.1,0,0,h1)



   id2_pool <- c(rep(2,length(x)),rep(3,length(x.add)))


    x2_pool <- 0.6 + 0.4*c(rev(x),x.add)
    y2_pool <- c(rev(y),y.add)


    x3_pool <- c(0.42, 0.42, 0.55, 0.55, 0.60, 0.55, 0.55)
    y3_pool <- c(0.45, 0.55, 0.55, 0.60, 0.50, 0.40, 0.45)

    xpool <- c(x1_pool, x2_pool, x3_pool)
    ypool <- c(y1_pool, y2_pool, y3_pool)

    if(!is.null(id)){
        id_pool <- id + c(id1_pool, id2_pool, rep(4, length(x3_pool)))
    }else{
        id_pool <- c(id1_pool, id2_pool, rep(4, length(x3_pool)))
    }

    x <- x.pos + wt*xpool
    y <- y.pos + ht*ypool
    
    fill <- c("blue","orange", "orange", "grey80")
    
    list(x=x,y=y,id=id_pool,fill=fill)
}

##############  letter C to A  ####################################

letter_C_to_A <- function(x.pos,y.pos,ht,wt,id=NULL){
  
  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)
  
  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)
  
  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))
  
  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)
  
  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)
  
  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)
  
  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))
  
  x <- c(x,rev(x1))
  y <- c(y,rev(y1))
  x1 <- 0.4*x
  y1 <- 1*y
  
  
  x <- 0.1* (c(0,4,6,10,8,6.8,3.2,2,0,3.6,5,6.4,3.6))
  y <- 0.1*(c(0,10,10,0,0,3,3,0,0,4,7.5,4,4))
  x2 <- 0.6 + 0.4*x
  y2 <- 1*y
  
  x3 <- c(0.42, 0.42, 0.55, 0.55, 0.60, 0.55, 0.55)
  y3 <- c(0.45, 0.55, 0.55, 0.60, 0.50, 0.40, 0.45)
  
  xpool <- c(x1, x2, x3)
  ypool <- c(y1, y2, y3)
  
  if(!is.null(id)){
    id_pool <- id +  c(rep(1,length(x1)), c(rep(2,9),rep(3,4)), 
                       rep(4, length(x3)))
  }else{
    id_pool <- c(rep(1,length(x1)), c(rep(2,9),rep(3,4)), 
                 rep(4, length(x3)))
  }
  
  x <- x.pos + wt*xpool
  y <- y.pos + ht*ypool
  
  fill=c("blue","green", "white", "grey80")
  
  list(x=x,y=y,id=id_pool,fill=fill)
}

################  letter T to G  ################################


letter_T_to_G <- function(x.pos,y.pos,ht,wt,id=NULL){
  
  x <- c( 0.4, 0.4, 0, 0, 1, 1, 0.6, 0.6)
  y <- c( 0, 0.85, 0.85, 1, 1, 0.85, 0.85, 0)
  x1_pool <- 0.4*x
  y1_pool <- 1*y
  id1_pool <- rep(1,length(x1_pool))
  
  
  x3_pool <- c(0.42, 0.42, 0.55, 0.55, 0.60, 0.55, 0.55)
  y3_pool <- c(0.45, 0.55, 0.55, 0.60, 0.50, 0.40, 0.45)
  
  
  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)
  
  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)
  
  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))
  
  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)
  
  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)
  
  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)
  
  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))
  
  x <- c(x,rev(x1))
  y <- c(y,rev(y1))
  
  h1 <- max(y.l1)
  r1 <- max(x.l1)
  
  h1 <- 0.4
  x.add <- c(r1,0.5,0.5,r1-0.2,r1-0.2,r1,r1)
  y.add <- c(h1,h1,h1-0.1,h1-0.1,0,0,h1)
  
  
  
  id2_pool <- c(rep(2,length(x)),rep(3,length(x.add)))
  
  
  x2_pool <- 0.6 + 0.4*c(rev(x),x.add)
  y2_pool <- c(rev(y),y.add)
  
  if(!is.null(id)){
    id_pool <- id + c(id1_pool, id2_pool, rep(4, length(x3_pool)))
  }else{
    id_pool <- c(id1_pool, id2_pool, rep(4, length(x3_pool)))
  }
  
  xpool <- c(x1_pool, x2_pool, x3_pool)
  ypool <- c(y1_pool, y2_pool, y3_pool)
  
  x <- x.pos + wt*xpool
  y <- y.pos + ht*ypool
  
  fill=c("red","orange","orange", "grey80")
  
  list(x=x,y=y,id=id_pool,fill=fill)
}

################ letter  T to C  ################################

letter_T_to_C <- function(x.pos,y.pos,ht,wt,id=NULL){
  
  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)
  
  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)
  
  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))
  
  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)
  
  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)
  
  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)
  
  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))
  
  x <- c(x,rev(x1))
  y <- c(y,rev(y1))
  x2 <- 0.6 + 0.4*x
  y2 <- 1*y
  
  
  x <- c( 0.4, 0.4, 0, 0, 1, 1, 0.6, 0.6)
  y <- c( 0, 0.85, 0.85, 1, 1, 0.85, 0.85, 0)
  x1 <- 0.4*x
  y1 <- 1*y
  
  
  x3 <- c(0.42, 0.42, 0.55, 0.55, 0.60, 0.55, 0.55)
  y3 <- c(0.45, 0.55, 0.55, 0.60, 0.50, 0.40, 0.45)
  
  xpool <- c(x1, x2, x3)
  ypool <- c(y1, y2, y3)
  
  if(!is.null(id)){
    id_pool <- id +c(rep(1,length(x1)), rep(3, length(x2)), 
                     rep(4, length(x3)))
  }else{
    id_pool <- c(rep(1,length(x1)), rep(3, length(x2)),
                 rep(4, length(x3)))
  }
  
  x <- x.pos + wt*xpool
  y <- y.pos + ht*ypool
  
  fill=c("blue","red","grey80")
  
  list(x=x,y=y,id=id_pool,fill=fill)
  
}

#################  letter T to A  ####################################

letter_T_to_A <- function(x.pos,y.pos,ht,wt,id=NULL){
  x <- c( 0.4, 0.4, 0, 0, 1, 1, 0.6, 0.6)
  y <- c( 0, 0.85, 0.85, 1, 1, 0.85, 0.85, 0)
  x1 <- 0.4*x
  y1 <- 1*y
  
  
  x <- 0.1* (c(0,4,6,10,8,6.8,3.2,2,0,3.6,5,6.4,3.6))
  y <- 0.1*(c(0,10,10,0,0,3,3,0,0,4,7.5,4,4))
  x2 <- 0.6 + 0.4*x
  y2 <- 1*y
  
  
  x3 <- c(0.42, 0.42, 0.55, 0.55, 0.60, 0.55, 0.55)
  y3 <- c(0.45, 0.55, 0.55, 0.60, 0.50, 0.40, 0.45)
  
  xpool <- c(x1, x2, x3)
  ypool <- c(y1, y2, y3)
  
  if(!is.null(id)){
    id_pool <- id + c( rep(3, length(x1)), c(rep(1,9),rep(2,4)), 
                       rep(4, length(x3)))
  }else{
    id_pool <- c( rep(3, length(x1)), c(rep(1,9),rep(2,4)),
                  rep(4, length(x3)))
  }
  
  fill=c("green","white", "red", "grey80")
  
  x <- x.pos + wt*xpool
  y <- y.pos + ht*ypool
  
  list(x=x,y=y,id=id_pool,fill=fill)
}

################# add letters to logo plot #########################


addLetter <- function(letters,which,x.pos,y.pos,ht,wt){
  
  if (which == "A"){
    letter <- letterA(x.pos,y.pos,ht,wt)
  }else if (which == "C"){
    letter <- letterC(x.pos,y.pos,ht,wt)    
  }else if (which == "G"){
    letter <- letterG(x.pos,y.pos,ht,wt)    
  }else if (which == "T"){
    letter <- letterT(x.pos,y.pos,ht,wt)    
  }else if (which == "C->T"){
    letter <- letter_C_to_T(x.pos,y.pos,ht,wt)    
  }else if (which == "C->A"){
    letter <- letter_C_to_A(x.pos,y.pos,ht,wt)    
  }else if (which == "C->G"){
    letter <- letter_C_to_G(x.pos,y.pos,ht,wt)    
  }else if (which == "T->A"){
    letter <- letter_T_to_A(x.pos,y.pos,ht,wt)    
  }else if (which == "T->G"){
    letter <- letter_T_to_G(x.pos,y.pos,ht,wt)    
  }else if (which == "T->C"){
    letter <- letter_T_to_C(x.pos,y.pos,ht,wt)    
  }else{
    stop("which must be one of A,C,G,T, C->T, 
         C->G, C->A, T->A, T->C, T->G")
  }
  
  letters$x <- c(letters$x,letter$x)
  letters$y <- c(letters$y,letter$y)
  
  lastID <- ifelse(is.null(letters$id),0,max(letters$id))
  letters$id <- c(letters$id,lastID+letter$id)
  letters$fill <- c(letters$fill,letter$fill)
  letters
}

damageLogo <- function(pwm, 
                       ic,
                       ic.scale=TRUE, 
                       xaxis=TRUE, 
                       yaxis=TRUE, 
                       xaxis_fontsize=10,
                       xlab_fontsize=15,
                       y_fontsize=15){
  
  if (class(pwm) == "data.frame"){
    pwm <- as.matrix(pwm)
  }else if (class(pwm) != "matrix"){
    stop("pwm must be of class matrix or data.frame")
  }
  
  if (any(abs(1 - apply(pwm,2,sum)) > 0.01))
    stop("Columns of PWM must add up to 1.0")
  
  chars <- c("A", "C", "G", "T", 
             "C->T", "C->A", "C->G", 
             "T->A", "T->C", "T->G")
  
  letters <- list(x=NULL,y=NULL,id=NULL,fill=NULL)
  npos <- ncol(pwm)
  
  
  if (ic.scale){
    ylim <- 2
    ylab <- "Information content"
    facs <- ic;
  }else{
    ylim <- 1
    ylab <- "Probability"
    facs <- rep(1, npos)
  }
  
  wt <- c(rep(1, floor(npos/2)),2,rep(1, floor(npos/2)))
  
  x.pos <- 0  
  for (j in 1:npos){
    
    column <- pwm[,j]
    hts <- 0.99*column*facs[j]
    letterOrder <- order(hts)
    
        y.pos <- 0    
        for (i in 1:10){
        letter <- chars[letterOrder[i]]
        ht <- hts[letterOrder[i]]
        if (ht>0) letters <- addLetter(letters,letter,x.pos,y.pos,ht,wt[j])
            y.pos <- y.pos + ht + 0.0001
        }
        x.pos <- x.pos + wt[j]
    
  }
  
  grid.newpage()
  bottomMargin = ifelse(xaxis, 2 + xfontsize/3.5, 3)
  leftMargin = ifelse(yaxis, 2 + yfontsize/3.5, 3)
  pushViewport(plotViewport(c(bottomMargin,leftMargin,3,3)))
  pushViewport(dataViewport(0:ncol(pwm),0:ylim,name="vp1"))
  grid.polygon(x=unit(letters$x,"native"), y=unit(letters$y,"native"),
               id=letters$id, gp=gpar(fill=letters$fill,col="transparent"))
  grid.polygon(x=unit(letters$x,"native"), y=unit(letters$y,"native"),
               id=letters$id,
               gp=gpar(fill=letters$fill,col="transparent"))
  
  xlim <- c(0.5, 1.5, 3, 4.5, 5.5)
  low_xlim <- c(xlim - 0.5*wt, xlim[length(xlim)]+0.5*wt[length(xlim)])
  
  for(n in 1:length(xlim)){
  grid.lines(x = unit(low_xlim[n], "native"),
             y = unit(c(0, 2), "native"),
             gp=gpar(col="grey80"))
  }
             
  if (xaxis){
    grid.xaxis(at=xlim,
               label=(c("\n left \n flank \n -2", 
                      "\n left \n flank \n -1", 
                      "\n mutation", 
                      "\n right \n flank \n +1",
                      "\n right \n flank \n +2")), 
               gp=gpar(fontsize=xaxis_fontsize))
    grid.text("Position",y=unit(-3,"lines"), 
              gp=gpar(fontsize=xlab_fontsize))
  }
  if (yaxis){
    grid.yaxis(gp=gpar(fontsize=yfontsize))
    grid.text(ylab,x=unit(-3,"lines"),rot=90, 
              gp=gpar(fontsize=y_fontsize))
  }
  popViewport()
  popViewport()
  par(ask=FALSE)
}


  
