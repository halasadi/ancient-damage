


##########################   test   aRchaic_pool()   #####################################

folders <- c("../data/Jones_2015/", "../data/Sardinia2017/")
out <- aRchaic_pool(folders, pattern = c("Bichon.sort.rmdup.IR.q30.mapDamage.subsampled.q30.csv", 
                                         "MA85snp_sorted_deduped_filtered.q30"))



datalist <- vector("list", length(folders))
for(i in 1:2){
  datalist[[i]] <- get(load(paste0(folders[i], tail(strsplit(folders[i], "/")[[1]],1), ".rda")))
}

sig_names <- colnames(datalist[[1]])
row_names_pool <- rownames(datalist[[1]])
if(length(datalist) >= 2){
  for(num in 2:length(datalist)){
    sig_names <- union(sig_names, colnames(datalist[[num]]))
    row_names_pool <- c(row_names_pool, rownames(datalist[[num]]))
  }
}

pooled_data <- matrix(0, length(row_names_pool), length(sig_names))
rownames(pooled_data) <- row_names_pool
colnames(pooled_data) <- sig_names

pattern = c("Vestonice15.q30.csv")

temp <- pooled_data[grep(pattern = pattern, paste0(rownames(pooled_data), ".csv")),, drop=FALSE]
