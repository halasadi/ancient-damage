

#########  John Lindo moderns + ancients mixed data  ######################

lindo_moderns_metadata <- read.table("../utilities/Lindo_metadata.txt")

lindoancient_data <- get(load("../processed_data/lindo2016ancients-counts-table.rda"))
lindomoderns_data <- get(load("../processed_data/lindo2016moderns-counts-table.rda"))

lindoancient_data_clubbed <- club_signature_counts(lindoancient_data)
lindomoderns_data_clubbed <- club_signature_counts(lindomoderns_data)


pooled_names <- intersect(colnames(lindoancient_data_clubbed), colnames(lindomoderns_data_clubbed))
filtered_ancients <- lindoancient_data_clubbed[, match(pooled_names, colnames(lindoancient_data_clubbed))]
filtered_moderns <- lindomoderns_data_clubbed[, match(pooled_names, colnames(lindomoderns_data_clubbed))]

lindo_pooled <- rbind(filtered_ancients, filtered_moderns)

lindo_filtered <- filtered_moderns

par(mfrow=c(2,2))

temp <- filter_signatures_only_location(lindo_filtered, max_pos = 20)


modern_names <- as.character(sapply(rownames(lindo_filtered), function(x) strsplit(x, "[_]")[[1]][1]))

par(mfrow=c(5,5))
for(k in 1:25){
  plot(1:20, temp[k,]/sum(temp[k,], na.rm=TRUE), col="red", pch=20, main=paste0(modern_names[k], "_", lindo_moderns_metadata[k,2]), xlab="read position", ylab="prop pattern")
  axis(1, at=c(2, 5, 10, 15, 20), c(2, 5, 10, 15, 20))
}



sum_temp <- colSums(temp)
sum_temp = sum_temp/sum(sum_temp, na.rm=TRUE)

plot(1:20, sum_temp, col="red", pch=20, main="all patterns", xlab="read position", ylab="prop of pattern")
axis(1, at=c(2, 5, 10, 15, 20), c(2, 5, 10, 15, 20))



temp <- filter_signatures_per_substitution(lindo_filtered, pattern="T->C", use_prop = FALSE)
temp <- filter_signatures_only_location(temp, max_pos = 20)
sum_temp <- colSums(temp)
sum_temp = sum_temp/sum(sum_temp, na.rm=TRUE)

plot(1:20, sum_temp, col="red", pch=20, main="pattern T->C", xlab="read position", ylab="prop of pattern")
axis(1, at=c(2, 5, 10, 15, 20), c(2, 5, 10, 15, 20))


temp <- filter_signatures_per_substitution(lindo_filtered, pattern="C->T", use_prop = FALSE)
temp <- filter_signatures_only_location(temp, max_pos = 20)
sum_temp <- colSums(temp)
sum_temp = sum_temp/sum(sum_temp, na.rm=TRUE)

plot(1:20, sum_temp, col="red", pch=20, main="pattern C->T", xlab="read position", ylab="prop of pattern")
axis(1, at=c(2, 5, 10, 15, 20), c(2, 5, 10, 15, 20))

temp <- filter_signatures_per_substitution(lindo_filtered, pattern="C->A", use_prop = FALSE)
temp <- filter_signatures_only_location(temp, max_pos = 20)
sum_temp <- colSums(temp)
sum_temp = sum_temp/sum(sum_temp, na.rm=TRUE)
plot(1:20, sum_temp, col="red", pch=20, main="pattern C->A", xlab="read position", ylab="prop of pattern")
axis(1, at=c(2, 5, 10, 15, 20), c(2, 5, 10, 15, 20))



signature_set <- colnames(lindo_pooled)
sig_split <- t(sapply(1:length(signature_set), function(x) return(strsplit(signature_set[x], "")[[1]][1:8])))
new_sig_split <- matrix(0, dim(sig_split)[1], 5);
new_sig_split[,1] <- sig_split[,1]
new_sig_split[,2] <- sig_split[,2]
new_sig_split[,3] <- sapply(1:length(signature_set), function(x) return(paste(sig_split[x,3:6], collapse="")))
new_sig_split[,4] <- sig_split[,7]
new_sig_split[,5] <- sig_split[,8]

levels(new_sig_split[,1]) <- c("0", "1", "2", "3", "4")

pos <- t(sapply(1:length(signature_set), function(x)
{
  y = strsplit(signature_set[x], "")[[1]]
  return(paste(y[10:length(y)], collapse=""))
}))



mat <- matrix(0, dim(new_sig_split)[1], dim(new_sig_split)[2])
for(k in 1:dim(new_sig_split)[2]){
  temp <- as.factor(new_sig_split[,k])
  mat[,k] <- as.numeric(as.matrix(plyr::mapvalues(temp, from = levels(temp), to = 0:(length(levels(temp))-1))))
}

pos <- as.numeric(pos)
pos <- pos - min(pos)
pos <- factor(pos, levels = 0:22)

signatures <- mat;
signature_pos <- cbind.data.frame(signatures, pos)


out <- topics(lindo_pooled, K=5, tol=100, model="independent", signatures = signature_pos)
save(out, file="../processed_data/maptpx-runs/lindo2016-maptpx-independent-K-5.rda")

out <- get(load("../processed_data/maptpx-runs/lindo2016-maptpx-independent-K-3.rda"))

labs <- c(rep("Ancient", 25), rep("Moderns", 25))
omega <- out$omega

cols1 <- c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
                 "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
                 "brown4","darkorchid","magenta","yellow", "azure1","azure4")

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs)
)

CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = cols1,
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = FALSE,
                            figure_title = paste0("StructurePlot: K=", dim(omega)[2],""),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))


damageLogo_pos(out$theta, ic.scale = TRUE, renyi_alpha = 10)

x <- do.call(rbind, lapply(rownames(out$theta), function(x) return(strsplit(as.character(x), "")[[1]][1:8])))

theta <- out$theta
sig_names = NULL
ic.scale=TRUE
max_pos = 15
flanking_bases=2
yscale_change = TRUE
xaxis=TRUE
yaxis=TRUE
xaxis_fontsize=10
xlab_fontsize=15
y_fontsize=15
mut_width=2
start=0.0001
pop_names=paste0("Cluster ",1:dim(theta)[2])
logoport_x = 0.35
logoport_y= 0.5
logoport_width= 0.35
logoport_height= 0.5
lineport_x = 0.95
lineport_y=0.75
lineport_width=0.25
lineport_height=0.45
