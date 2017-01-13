

########  Pie charts for the strand breaks  ######################

library(plotrix)
percentage <- c(c(0.3, 0.3, 0.2, 0.2), c(0.2, 0.2, 0.3, 0.3))
bases <- c(c("A", "G", "C", "T"), c("A", "G", "C", "T"))
strand <- c(rep("left", 4), rep("right", 4))
df <- data.frame("percentage"=percentage,
                 "base" = bases,
                 "strand" = strand)

library(ggplot2)
p = ggplot(data=df,
           aes(x=strand,
               y=percentage,
               fill = factor(bases, levels = c("A", "G", "C", "T"))
           ),
)

p = p + coord_polar(theta="y")
p

ggplot(mtcars,aes(x = factor(1),fill=factor(cyl))) +
  facet_wrap(~gear) +
  geom_bar(width = 1,position = "fill") +
  coord_polar(theta="y")




library(reshape2)
y  = data.frame(category=c("left","left", "left","left", "right", "right", "right", "right"), value=c("A", "G", "C", "T", "A", "G", "C", "T"),
                percentage = c(0.4, 0.4, 0.1, 0.1, 0.1, 0.2, 0.3, 0.4))

# get counts and melt it
data.m = y
names(data.m)[3] = "percentage"

# calculate percentage:
m1 = ddply(data.m, .(category), summarize, ratio=percentage/sum(percentage))

#order data frame (needed to comply with percentage column):
#m2 = data.m[order(data.m$category),]
m2 = data.m

# combine them:
mydf = data.frame(m2,ratio=m1$ratio)

# get positions of percentage labels:
mydf = ddply(mydf, .(category), transform, position = 1 - cumsum(percentage) + 0.5*percentage)
bases = factor(mydf$value, levels=c("A", "G", "C", "T"))
# create bar plot
pie = ggplot(mydf, aes(x = factor(1), y = percentage, fill = bases)) +
  geom_bar(stat = "identity", width = 1) +
  facet_wrap(~category) +
  labs(x = "")+
  labs(y = "strand breaks")+
  labs(title = "Strand breaks composition")

# make a pie
pie = pie + coord_polar(theta = "y") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# add labels
pie + geom_text(aes(label = sprintf("%1.2f%%", 100*ratio), y = position))
