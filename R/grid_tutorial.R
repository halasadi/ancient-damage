

#########   grid construct   ##################

# Grid provides tools to draw and arrange basic shapes. 
# It is a very low level graphics package. 
# Grid provides functions that allow the plot canvas to be 
# accessed programmatically. 
#Viewports create a context for drawing. 
# The basic shapes available are polygons, curves, 
# raster images, and data symbols.



library(grid)
grid.newpage()
vp <- viewport(x=0.5,y=0.5,width=0.9, height=0.9)

pushViewport(vp)
grid.rect(gp=gpar(lty="dashed"))
grid.circle(x=0.6, y=0.4, r=0.3)

stickperson <- function() {
  grid.circle(x=.5, y=.8, r=.1, gp=gpar(fill="yellow"))
  grid.lines(c(.5,.5), c(.7,.2)) # vertical line for body
  grid.lines(c(.5,.7), c(.6,.7)) # right arm
  grid.lines(c(.5,.3), c(.6,.7)) # left arm
  grid.lines(c(.5,.65), c(.2,0)) # right leg
  grid.lines(c(.5,.35), c(.2,0)) # left leg
}

stickperson()

vp1 <- viewport(x=0.5, y=0.75, width=0.6, height=0.3)
pushViewport(vp1)

grid.circle(gp=gpar(col="blue"))
# plot the outline of vp1:
grid.rect()

upViewport()

pushViewport(viewport())
grid.lines(c(.05, .95), c(.95, .05))
grid.lines(c(.05, .95), c(.05, .95))
for (i in 1:100) {
  vp <- viewport(h=.9, w=.9)
  pushViewport(vp)
  grid.rect()
}

viewport(x = 0.5, y = 0.5, width = 0.5, height = 0.25, angle=45)
popViewport()

grid.rect(gp = gpar(lty = "dashed"))
vp1 <- viewport(x = 0, y = 0.5, w = 0.5, h = 0.5,
                just = c("left", "bottom"), name = "vp1")
vp2 <- viewport(x = 0.5, y = 0, w = 0.5, h = 0.5,
                just = c("left", "bottom"))

pushViewport(vp1)
grid.rect(gp = gpar(col = "grey"))
grid.text("Some drawing in graphics region 1", y = 0.8)
upViewport()
pushViewport(vp2)
grid.rect(gp = gpar(col = "grey"))
grid.text("Some drawing in graphics region 2", y = 0.8)
upViewport()
downViewport("vp1")
grid.text("MORE drawing in graphics region 1", y = 0.2)
popViewport()

grid.rect(gp = gpar(lty = "dashed"))
vp <- viewport(width = 0.5, height = 0.5)
pushViewport(vp)
grid.rect(gp = gpar(col = "grey", fill="blue"))
grid.text("quarter of the page", y = 0.85)
pushViewport(vp)
grid.rect()
popViewport()

# npc: Normalised parent coordinates
# native: relative to current x-, y-scale
# 

##  Paul Murrell R graphics codes

data("pressure")
plot(pressure)
text(150, 600, 
     "Pressure (mm Hg)\nversus\nTemperature (Celsius)")


###############   Draw a picture   in  R   #########################

grid.newpage()
pushViewport(viewport(xscale=c(0, 1), yscale=c(0.5, 1),
                      clip=TRUE))

res <- 50
for (i in 1:res)
  grid.rect(y=1 - (i-1)/res, just="top",
            gp=gpar(col=NULL, fill=grey(0.5*i/res)))

moon <- function(x, y, size) {
  angle <- seq(-90, 90, length=50)/180*pi
  x1 <- x + size*cos(angle)
  y1 <- y + size*sin(angle)
  mod <- 0.6
  x2 <- x + mod*(x1 - x)
  grid.polygon(c(x1, rev(x2)), c(y1, rev(y1)),
               default.unit="native",
               gp=gpar(col="blue", fill="white"))
}

x = .1 
y = .9 
size = .1

moon(x,y,size)

grid.circle(runif(20, .2, 1), runif(20, .6, 1), r=.005,
            default.unit="native",
            gp=gpar(col=NULL, fill="white"))

star <- function(x, y, size) {
  x1 <- c(x,           x + size*.1, x + size*.5, x + size*.1,
          x,           x - size*.1, x - size*.5, x - size*.1) + .05
  y1 <- c(y - size,    y - size*.1, y,           y + size*.1,
          y + size*.7, y + size*.1, y,           y - size*.1) + .05
  grid.polygon(x1, y1, 
               default.unit="native",
               gp=gpar(col=NULL, fill="white"))
}

x = .5
y = .7
size = .02

star(.5, .7, .02)
star(.8, .9, .02)
star(.72, .74, .02)
star(.62, .88, .02)

hill <- function(height=0.1, col="black") {
  n <- 100
  x <- seq(0, 1, length=n)
  y1 <- sin(runif(1) + x*2*pi)
  y2 <- sin(runif(1) + x*4*pi)
  y3 <- sin(runif(1) + x*8*pi)
  y <- 0.6 + height*((y1 + y2 + y3)/3)
  grid.polygon(c(x, rev(x)), c(y, rep(0, n)),
               default.unit="native",
               gp=gpar(col=NULL, fill=col))
}

hill()

rdir <- function(n) {
  sample(seq(-45, 45, length=10), n)/180*pi
}

grid.text("Once upon a time ...",
          x=.15, y=.51, just="bottom",
          default.unit="native",
          gp=gpar(col="white", fontface="italic", fontsize=10))

rdir(10)


grid.roundrect(width=.5, height=.5, name="rr")
theta <- seq(0, 360, length=50)
for (i in 1:50)
  grid.circle(x=grobX("rr", theta[i]),
              y=grobY("rr", theta[i]),
              r=unit(1, "mm"),
              gp=gpar(fill="black"))


par(mfrow=c(2, 2), cex=0.6, mar=c(4, 4, 4, 2), mex=0.8)
plot(lm.SR <- lm(sr ~ pop15 + pop75 + dpi + ddpi, data = LifeCycleSavings),
     id.n=1)


par(mfrow=c(2, 2))
z <- 2 * volcano        # Exaggerate the relief
x <- 10 * (1:nrow(z))   # 10 meter spacing (S to N)
y <- 10 * (1:ncol(z))   # 10 meter spacing (E to W)
# Don't draw the grid lines :  border = NA
par(mar=rep(0, 4))
persp(x, y, z, theta = 135, phi = 30, col = "light grey", scale = FALSE,
      ltheta = -120, shade = 0.75, border = NA, box = FALSE)
mtext("persp()", side=3, line=-2)
contour(x, y, z, asp=1, labcex=0.35, axes=FALSE)
image(x, y, z, asp=1, col=grey(0.5 + 1:10/24), xlab="", ylab="", axes=FALSE)


################   grid lines ##############################

grid.newpage()
pushViewport(viewport(xscale=c(0, 1), yscale=c(0.5, 1),
                      clip=TRUE))

x <- c(.3, .7, .3)
y <- c(.2, .5, .8)
grid.rect(gp=gpar(col="grey"))
grid.lines(x, y, gp=gpar(lwd=10, lineend="square",
                         linejoin="mitre", col="black"))

grid.lines(x, y, gp=gpar(lwd=10, col="grey50"))
# lineend="round", linejoin="round"

grid.lines(x, y, gp=gpar(lwd=10, lineend="butt",
                         linejoin="bevel", col="grey80"))

grid.points(x, y, default="npc", pch=16, gp=gpar(cex=0.5))


plot.new()
par(new=TRUE)
grid.newpage()
pushViewport(viewport(xscale=c(0, 1), yscale=c(0.5, 1),
                      clip=TRUE))
grid.rect(gp=gpar(col="grey"))

ncol <- 4
nrow <- 4
xadj <- c(1, 0.5, NA, 0)
yadj <- c(1, 0.5, NA, 0)
size <- unit(3, "mm")
for (i in 1:nrow) {
  for (j in 1:ncol) {
    x <- i/(nrow + 1)
    y <- j/(ncol + 1)
    xu <- unit(x, "npc")
    yu <- unit(y, "npc")
    grid.segments(unit.c(xu - size, xu),
                  unit.c(yu, yu - size),
                  unit.c(xu + size, xu),
                  unit.c(yu, yu + size),
                  gp=gpar(col="grey"))
    text(x, y, paste("c(", xadj[j], ", ", yadj[i], ")", sep=""),
         adj=c(xadj[j], yadj[i]))
  }
}

##################  Trellis  ###############################
plot.new()
library(lattice)
tplot <- xyplot(lat ~ long, data=quakes, pch=".")

trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.set(list(dot.symbol=list(pch=1)))
print(
  update(tplot, 
         main="Earthquakes in the Pacific Ocean\n(since 1964)")
  
)


trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.set(list(dot.symbol=list(pch=1)))
print(
  xyplot(lat ~ long, data=quakes, pch=".")
  
)

grid.newpage()
grid.rect(gp=gpar(col="grey"))
grid.circle(x=seq(0.1, 0.9, length=100), 
            y=0.5 + 0.4*sin(seq(0, 2*pi, length=100)),
            r=abs(0.1*cos(seq(0, 2*pi, length=100))))

grid.lines()
# Using id (NOTE: locations are not in consecutive blocks)
grid.newpage()
grid.polygon(x=c((0:4)/10, rep(.5, 5), (10:6)/10, rep(.5, 5)),
              y=c(rep(.5, 5), (10:6/10), rep(.5, 5), (0:4)/10),
              id=rep(1:5, 4),
              gp=gpar(col=1:5, lwd=3))

grid.lines()
# Using id (NOTE: locations are not in consecutive blocks)
grid.newpage()
grid.polyline(x=c((0:4)/10, rep(.5, 5), (10:6)/10, rep(.5, 5)),
             y=c(rep(.5, 5), (10:6/10), rep(.5, 5), (0:4)/10),
             id=rep(1:5, 4),
             gp=gpar(col=1:5, lwd=3))

grid.newpage()
grid.polyline(x=outer(c(0, .5, 1, .5), 5:1/5),
              y=outer(c(.5, 1, .5, 0), 5:1/5),
              id.lengths=rep(4, 5),
              gp=gpar(col=1:5, lwd=3))

grid.newpage()
grid.rect(gp=gpar(col="grey"))
angle <- seq(0, 2*pi, length=10)[-10]
grid.polygon(x=0.25 + 0.15*cos(angle), y=0.5 + 0.3*sin(angle), 
             gp=gpar(fill="grey"))
grid.polygon(x=0.75 + 0.15*cos(angle), y=0.5 + 0.3*sin(angle), 
             id=rep(1:3, each=3),
             gp=gpar(fill="grey"))

grid.newpage()
grid.rect(gp=gpar(col="grey"))
pushViewport(viewport(gp=gpar(fontsize=10)))
grid.rect(x=0.33, height=0.7, width=0.2, gp=gpar(fill="black"))



######################  Lines   ####################################

require(RGraphics)

label <- function(label, row, col, title=FALSE, box=TRUE, tcol="black",
                  fill=NA, lwd=1) {
  if (title) {
    face <- "bold"
    cex <- 0.8
  } else {
    face <- "plain"
    cex <- 0.8
  }
  pushViewport(viewport(layout.pos.row=row,
                        layout.pos.col=col))
  if (box)
    grid.roundrect(w=unit(0.82, "inches"), h=unit(1.2, "lines"),
                   r=unit(0.2, "snpc"), 
                   gp=gpar(col=tcol, lwd=lwd, fill=fill))
  grid.text(label, gp=gpar(fontface=face, cex=cex, col=tcol))
  popViewport()
}

arrow <- function(row1, col1, row2, col2, col="black") {
  pushViewport(viewport(layout.pos.row=row1,
                        layout.pos.col=col1))
  grid.move.to(x=0.5, y=unit(0.5, "npc") - unit(0.6, "lines"))
  popViewport()
  pushViewport(viewport(layout.pos.row=row2,
                        layout.pos.col=col2))
  arrow(grob=lineToGrob(x=0.5,
                        y=unit(0.5, "npc") + unit(0.6, "lines")),
                      angle=10, type="closed", 
                     length=unit(0.15, "inches"), gp=gpar(col=col, fill=col))
  popViewport()
} 

pushViewport(viewport(width=unit(4.9, "inches"),
                      layout=grid.layout(7, 5,
                                         heights=unit(c(3,.7,3,.7,2,.7,3),
                                                      c("lines", "inches",
                                                        "lines", "inches",
                                                        "lines", "inches",
                                                        "lines")))
))

pushViewport(viewport(layout.pos.row=1,
                      layout.pos.col=2:5))
grid.roundrect(gp=gpar(lty="dashed"), r=unit(0.1, "snpc"))
popViewport()

label("Graphics\nPackages", 1, 1, title=TRUE, box=FALSE)
label("lattice", 1, 4, lwd=3, fill="grey80")
label("...", 1, 5)
label("maps", 1, 2)
label("...", 1, 3)
arrow(1, 2, 3, 3)
arrow(1, 3, 3, 3)
arrow(1, 4, 3, 4)
arrow(1, 5, 3, 4)

pushViewport(viewport(layout.pos.row=3,
                      layout.pos.col=2:5))
grid.roundRect(gp=gpar(lty="dashed"), r=unit(0.1, "snpc"))
popViewport()

label("Graphics\nSystems", 3, 1, title=TRUE, box=FALSE)
label("graphics", 3, 3, lwd=3, fill="grey80")
label("grid", 3, 4, lwd=3, fill="grey80")
arrow(3, 3, 5, 3:4)
arrow(3, 4, 5, 3:4)

label("Graphics\nEngine\n&\nDevices", 5, 1, title=TRUE, box=FALSE)
label("grDevices", 5, 3:4, lwd=3, fill="grey80")
# arrow(5, 3:4, 7, 2)
# arrow(5, 3:4, 7, 3)
# arrow(5, 3:4, 7, 4)
# arrow(5, 3:4, 7, 5)

pushViewport(viewport(layout.pos.row=7,
                      layout.pos.col=2:5))
grid.roundRect(gp=gpar(lty="dashed"), r=unit(0.1, "snpc"))
popViewport()

label("Graphics\nDevice\nPackages", 7, 1, title=TRUE, box=FALSE)
label("gtkDevice", 7, 3)
label("...", 7, 4)

popViewport()