

########  T -> AAA  effect  ####################


signature_counts_2 <- get(load("../summary_data/153-ancient-signature-counts-clubbed.rda"))

signature_counts_1 <- get(load("../summary_data/signature-counts-clubbed-Lindo2016.rda"))

signature_counts <- rbind(signature_counts_1, signature_counts_2)

signature_set <- colnames(signature_counts)

RISE_labels <- grep("RISE", rownames(signature_counts_2))
I_labels <- setdiff(1:152,RISE_labels)


sig_T_AAA_counts <- rowSums(signature_counts[,grep("T->AAA", colnames(signature_counts))])
sig_counts <- rowSums(signature_counts)

T_AAA_mat <- cbind(sig_T_AAA_counts, sig_counts - sig_T_AAA_counts);
T_AAA_mat <- t(apply(T_AAA_mat, 1, normalize))
colnames(T_AAA_mat) <- c("T_AAA", "Other")

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(T_AAA_mat))),
  tissue_label = factor(c(rep("Ancient_Lindo",25),
                          rep("Modern_Lindo",25),
                          rep("I",length(I_labels)), rep("RISE",length(RISE_labels))))
)


plot(sig_T_AAA_counts, main="Number of T->AAA (each sample)",
     col=annotation$tissue_label, pch=20, ylab="no. of T->AAA")
legend("topright", legend=c("Lindo_ancient", "Lindo_modern", "IO", "RISE"),
       fill=unique(annotation$tissue_label))


plot(sig_T_AAA_counts/sig_counts, main="Prop of T->AAA (each sample)",
     col=annotation$tissue_label, pch=20, ylab="prop. of T->AAA")
legend("topright", legend=c("Lindo_ancient", "Lindo_modern", "IO", "RISE"),
       fill=unique(annotation$tissue_label))