

#####   number of damages per read added as a feature   ###########################

par(mfrow=c(1,1))
signature_counts <- get(load("../processed_data/annagosling2016-counts-table.rda"))
validation_check <- club_signature_validation_plot(signature_counts, log=TRUE)
clubbed_counts <- club_signature_counts(signature_counts)
clubbed_counts <- clubbed_counts[-28,];
filtered_counts <- filter_signatures_wo_location(clubbed_counts);

# library(readxl)
# mapped_reads <- read_excel("../data/AnnaGosling2016data/metadata.xlsx")
#
# mapped_reads <- read.csv("../data/AnnaGosling2016data/metadata.csv")
#
# damages_per_read <- rowSums(clubbed_counts)/ mapped_reads[,3]
#
# mapped_reads <- mapped_reads[-c(32,56),]
#
# grep(as.character(mapped_reads[1,1]), rownames(clubbed_counts))
#
# names <- as.character(sapply(rownames(clubbed_counts), function(x) return(strsplit(x, "[-,.]")[[1]][3])))
#
# mapped_reads_ordered <- mapped_reads[match(names, as.character(mapped_reads[,1])),];
#
# save(mapped_reads_ordered, file="../data/AnnaGosling2016data/metadata_ordered.rda")


metadata <- get(load(file="../data/AnnaGosling2016data/metadata_ordered.rda"))
prop_damages <- rowSums(clubbed_counts)/metadata[,3]

extra_col <- round(100*prop_damages);

filtered_counts_xx <- cbind(filtered_counts, extra_col);

names <- rownames(clubbed_counts_xx);
control_indices <- c(grep("EXN", names), grep("Libneg", names), grep("PCRneg", names))

labs <- character();
labs <- rep("ancient", dim(clubbed_counts)[1])
labs[control_indices] <- "controls"

topics_clus <- maptpx::topics(clubbed_counts_xx, K=3, tol=1);
save(topics_clus, file="../rda/annagosling2016_maptpx_with_damagecounts_perread_added_K_3.rda")

topics_clus <- get(load(file="../rda/annagosling2016_maptpx_with_damagecounts_perread_added_K_3.rda"))
theta <- topics_clus$theta
omega <- topics_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs)
)

CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Moderns vs Ancients",
                            order_sample = FALSE,
                            figure_title = paste0("StructurePlot: K=", dim(omega)[2],": pmsignature: withput position information"),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))
signature_set <- colnames(clubbed_counts_xx)
apply(CountClust::ExtractTopFeatures(topics_clus$theta, top_features = 10, method="poisson", options="min"), c(1,2), function(x) signature_set[x])

color <- c("black", "red")
plot(log(metadata[,3]), log(as.numeric(rowSums(clubbed_counts_xx))), xlab="log number of reads mapped",
     ylab="log total no. of damages", pch=20, col=color[factor(labs)], cex=1)
legend("topleft", legend = c("ancients", "controls"), fill=c("black", "red"), cex=0.5)

##############   slope check   #######################################

file = "../data/AnnaGosling2016data/ADR-T2-EXN1.dup.q30.csv"
out <- damage_slope_calc(file)

dir <- "../data/AnnaGosling2016data/";
files <- list.files(dir, pattern = ".csv");

slope_mat <- numeric();
out <- list();
for(l in 1:length(files)){
  out[[l]] <- damage_slope_calc(paste0(dir, files[l]), type="slope");
  slope_mat <- cbind(slope_mat, out[[l]]$left$`slope-counts`[-1,1]);
}

save(out, file="../rda/annagosling2016/slope_counts_data_gosling2016.rda")

colnames(slope_mat) <- files;
slope_data <-  t(slope_mat)
slope_data <-slope_data[-28,];

pr <- prcomp(log(slope_data))
par(mfrow=c(2,2))
plot(pr$x[,1], pr$x[,2], col=color[factor(labs)], pch=20, cex=1, xlab="PC1", ylab="PC2")
plot(pr$x[,1], pr$x[,3], col=color[factor(labs)], pch=20, cex=1, xlab="PC1", ylab="PC3")
plot(pr$x[,2], pr$x[,3], col=color[factor(labs)], pch=20, cex=1, xlab="PC2", ylab="PC3")


