

############  R script damage counts processing (from csv)  #########################

df1 <- damage.build.counts(file="../summary_data/I0026.q30.csv")

df2 <- cbind.data.frame(paste0(df1[,1], "_", df1[,2]), df1[,3])
colnames(df2) <- c("pattern-bins", "counts")

C_to_A_data <- df1[grep("C->A", df1$pattern),]

bin_counts_C_to_A <- tapply(C_to_A_data[,3], C_to_A_data[,2], sum)

barplot(bin_counts_C_to_A, main="C to A damage numbers")

C_to_T_data <- df1[grep("C->T", df1$pattern),]

bin_counts_C_to_T <- tapply(C_to_T_data[,3], C_to_T_data[,2], sum)

barplot(bin_counts_C_to_T, main="C to T damage numbers")



