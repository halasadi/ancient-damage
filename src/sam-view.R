
###  SAM view

library(Rsamtools)
ncol <- 15 # adjust to number of cols in your file
# adapt the what part to fit your file:
mysam <- scan("125_all_chr.sam", what = list(rep(character(), ncol)), fill=TRUE)

lines <- readLines("125_all_chr.sam")
lines[1]
