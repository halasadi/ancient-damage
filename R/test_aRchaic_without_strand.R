

#############  test aRchaic cluster beta without strand ##################


mat1 <- get(load("../data/Fu_2016/Fu_2016.rda"))
mat2 <- get(load("../data/moderns_lite/moderns_lite.rda"))

mat_pooled <- rbind(mat1, mat2)

mat_reduced <- filter_out_strand(mat_pooled)


signature_set <- colnames(mat_reduced)
sig_split <- t(sapply(1:length(signature_set), function(x) return(strsplit(signature_set[x], "")[[1]][1:(flanking_bases+5)])))
new_sig_split <- matrix(0, dim(sig_split)[1], 3);
new_sig_split[,1] <- sig_split[,1]
new_sig_split[,2] <- sapply(1:length(signature_set), function(x) return(paste(sig_split[x,(flanking_bases+1):(flanking_bases+4)], collapse="")))
new_sig_split[,3] <- sig_split[,(flanking_bases+5)]

pos <- sapply(signature_set, function(x) return(strsplit(x, "_")[[1]][3]))

pos <- as.numeric(pos)
pos <- pos - min(pos)
pos <- factor(pos, levels = 0:21)

strand_break <- sapply(signature_set, function(x) return(strsplit(x, "_")[[1]][2]))
strand_break <- factor(strand_break)
signature <- new_sig_split

signature_pos_strand_break <- cbind.data.frame(new_sig_split, strand_break, pos)


K = 3
tol = 10

topics.control <- list()
topics.control.default <- list(bf = FALSE, kill = 2, ord = TRUE, verb = 1, admix = TRUE,
                               nbundles = 1, use_squarem = FALSE, init.adapt = FALSE, type = "full",
                               light = 1, method_admix = 1, sample_init = TRUE, tmax = 10000)
topics.control <- modifyList(topics.control.default, topics.control)


suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = mat_reduced, K=K, tol=tol, model = "independent", signatures = signature_pos_strand_break), topics.control)))

suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = mat_reduced, K=K, tol=tol, model = "full", signatures = NULL), topics.control)))

theta_pool <- topic_clus$theta

logo.control <- list()
logo.control.default <- list(sig_names = NULL, ic.scale=TRUE,
                             max_pos = 20, flanking_bases=1,
                             yscale_change = TRUE, xaxis=TRUE,
                             yaxis=TRUE, xlab = " ", xaxis_fontsize=20,
                             xlab_fontsize=10, title_aligner = 11,
                             y_fontsize=20, title_fontsize = 35,
                             mut_width=2, start=0.0001,
                             renyi_alpha = 5, inflation_factor = c(3,1,3),
                             pop_names = paste0("Cluster : ", 1:K),
                             logoport_x = 0.25,
                             logoport_y= 0.50,
                             logoport_width= 0.28,
                             logoport_height= 0.40,
                             lineport_x = 0.9,
                             lineport_y=0.40,
                             lineport_width=0.32,
                             lineport_height=0.28,
                             breaklogoport_x = 0.90,
                             breaklogoport_y = 0.40,
                             breaklogoport_width=0.35,
                             breaklogoport_height=0.50,
                             output_width = 1200, output_height = 700)

logo.control <- modifyList(logo.control.default, logo.control)

output_dir <- "../utilities/structure_mutation_wo_strand/moderns_Fu/"
if(is.null(output_dir)){ output_dir <- paste0(getwd(),"/")}
plot.new()
do.call(damageLogo_three, append(list(theta_pool = topic_clus$theta,
                                    output_dir = output_dir),
                               logo.control))
graphics.off()

