
###  Exploration Different signatures

signature_counts <- get(load("../summary_data/signature-counts-clubbed-Lindo2016.rda"))

sig_split <- do.call(rbind, lapply(colnames(signature_counts), function(x) strsplit(as.character(x), split="")[[1]]))

sig_split_red <- sig_split[,3:6];

substitute_sig <- sapply(1:dim(sig_split_red)[1], function(x) paste0(sig_split_red[x,], collapse=""))

normalize <- function(x) return(x/sum(x))
prop_substitute_sig <- do.call(rbind, 
          lapply(1:dim(signature_counts)[1], function(x) 
          normalize(tapply(signature_counts[x,], substitute_sig, sum))))

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(prop_substitute_sig))),
  tissue_label = factor(c(rep("Ancient",25), rep("Modern",25)))
)

rownames(prop_substitute_sig) <- annotation$sample_id;

StructureGGplot(omega = prop_substitute_sig,
                annotation = annotation,
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Development Phase",
                order_sample = FALSE,
                figure_title = paste0("StructurePlot: 6 substituion patterns (with C -> T)"),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))


##  Now we remove the C -> T substitutions 

prop_substitute_sig_filtered <- prop_substitute_sig[,-3]

prop_substitute_sig_filtered_norm <- t(apply(prop_substitute_sig_filtered, 1, normalize))

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(prop_substitute_sig_filtered_norm))),
  tissue_label = factor(c(rep("Ancient",25), rep("Modern",25)))
)

rownames(prop_substitute_sig_filtered_norm) <- annotation$sample_id;

StructureGGplot(omega = prop_substitute_sig_filtered_norm,
                annotation = annotation,
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Development Phase",
                order_sample = FALSE,
                figure_title = paste0("StructurePlot: 5 substituion patterns (without C -> T)"),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))


######  T -> AAA pattern counts ##########################



sig_T_AAA_counts <- rowSums(signature_counts[,grep("T->AAA", colnames(signature_counts))])
sig_counts <- rowSums(signature_counts)

T_AAA_mat <- cbind(sig_T_AAA_counts, sig_counts - sig_T_AAA_counts);
T_AAA_mat <- t(apply(T_AAA_mat, 1, normalize))
colnames(T_AAA_mat) <- c("T_AAA", "Other")

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(T_AAA_mat))),
  tissue_label = factor(c(rep("Ancient",25), rep("Modern",25)))
)

rownames(T_AAA_mat) <- annotation$sample_id;

barplot(sig_T_AAA_counts/sig_counts, main="Proportion of T->AAA (each sample)")
