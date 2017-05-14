

#########################   Building cluster_aRchaic   ########################################################



tol=100
labs = NULL
levels = NULL
run_from = "start"
run_index = 1:length(folders)
gom_method = "independent"
topic_cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
               "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
               "brown4","darkorchid","magenta","yellow", "azure1","azure4")
structure.control = list()
logo.control = list()
topics.control = list()
output_dir = NULL
save_plot = TRUE

labs <- c("ISB", "LON", rep("MA", 21), rep("SUA", 3), rep("Siberia",1), rep("Iceman", 2))
levels = rev(c("Siberia", "Iceman", "MA", "SUA", "ISB", "LON"))

folders <- c("../data/sardinia/", "../data/Siberian_MA-1_1stExtraction/", "../data/Tyrolean_Iceman/")

aRchaic_cluster(folders, K=4,
                tol=100,
                labs = labs,
                levels = levels,
                run_from = "start",
                run_index = 1:length(folders),
                gom_method = "independent",
                topic_cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
                               "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
                               "brown4","darkorchid","magenta","yellow", "azure1","azure4"),
                structure.control = list(),
                logo.control = list(),
                topics.control = list(),
                output_dir = "../utilities/model_output/",
                save_plot = TRUE)

for(i in 1:length(folders)){
  if(!file.exists(folders[i]))
    stop("A folder in the folder list does not exist:  aborting")
}

datalist <- vector("list", length(folders))

run_from_scratch = FALSE
run_index = 1:length(folders)

if(run_from_scratch){
  if(sum(run_index - 1:length(folders))^2 == 0){
    for(i in 1:length(folders)){
      file.remove(paste0(folders[i], tail(strsplit(folders[i], "/")[[1]],1), ".rda"))
    }
  }else{
      folders1 <- folders[run_index]
      for(i in 1:length(folders1)){
        file.remove(paste0(folders1[i], tail(strsplit(folders1[i], "/")[[1]],1), ".rda"))
      }
  }
}

for(i in 1:length(folders)){
  if(!file.exists(paste0(folders[i], tail(strsplit(folders[i], "/")[[1]],1), ".rda"))){
    out <- aggregate_signature_counts(dir = paste0(folders[i]),
                                      pattern = NULL,
                                      breaks = c(-1, seq(1,20,1)),
                                      flanking_bases = 1)
    clubbed_data <- club_signature_counts(out, flanking_bases = 1)
    save(clubbed_data, file = paste0(folders[i], tail(strsplit(folders[i], "/")[[1]],1), ".rda"))
    datalist[[i]] <- clubbed_data
  }else{
    datalist[[i]] <- get(load(paste0(folders[i], tail(strsplit(folders[i], "/")[[1]],1), ".rda")))
  }
}


###########################################   Pooling the data from multiple sources  ##########################################


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

for(num in 1:length(datalist)){
  pooled_data[match(rownames(datalist[[num]]), rownames(pooled_data)), match(colnames(datalist[[num]]), sig_names)] <- datalist[[num]]
}


signature_set <- colnames(pooled_data)
sig_split <- t(sapply(1:length(signature_set), function(x) return(strsplit(signature_set[x], "")[[1]][1:8])))
new_sig_split <- matrix(0, dim(sig_split)[1], 3);
new_sig_split[,1] <- sig_split[,1]
new_sig_split[,2] <- sapply(1:length(signature_set), function(x) return(paste(sig_split[x,2:5], collapse="")))
new_sig_split[,3] <- sig_split[,6]

levels(new_sig_split[,1]) <- c("0", "1", "2", "3", "4")

pos <- t(sapply(1:length(signature_set), function(x)
{
  y = strsplit(signature_set[x], "")[[1]]
  return(paste(y[12:length(y)], collapse=""))
}))



mat <- matrix(0, dim(new_sig_split)[1], dim(new_sig_split)[2])
for(k in 1:dim(new_sig_split)[2]){
  temp <- as.factor(new_sig_split[,k])
  mat[,k] <- as.numeric(as.matrix(plyr::mapvalues(temp, from = levels(temp), to = 0:(length(levels(temp))-1))))
}

pos <- as.numeric(pos)
pos <- pos - min(pos)
pos <- factor(pos, levels = 0:21)

signatures <- mat;
signature_pos <- cbind.data.frame(signatures, pos)


##################################  Grade of Membership  Model fit  ########################################################

output_dir <- "../utilities/"
K= 3
tol = 100
model = "independent"

if(is.null(output_dir)){
  output_dir <- getwd()
}

if(!file.exists(paste0(output_dir, "model.rda"))){
  topic_clus <- maptpx::topics(pooled_data, K=K, tol=tol, model="independent", signatures = signature_pos)
  save(topic_clus, file = paste0(output_dir, "model.rda"))
}else{
  topic_clus <- get(load(paste0(output_dir, "model.rda")))
}


###################################  STRUCTURE  plot representation  ##########################################################


omega <- topic_clus$omega
labs <- c("ISB", "LON", rep("MA", 21), rep("SUA", 3), rep("Siberia",1), rep("Iceman", 2))
annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs, levels = rev(c("Siberia", "Iceman", "MA", "SUA", "ISB", "LON")))
)

cols1 <- c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
           "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
           "brown4","darkorchid","magenta","yellow", "azure1","azure4")



save_plot = TRUE
if(save_plot){
  png(paste0(output_dir, "structure.png"), height = 700, width = 300)
  StructureGGplot(omega = omega,
                  annotation = annotation,
                  palette = cols1,
                  yaxis_label = "Moderns vs Ancients",
                  order_sample = FALSE,
                  figure_title = paste0("  StructurePlot: K=", dim(omega)[2],""),
                  axis_tick = list(axis_ticks_length = .1,
                                   axis_ticks_lwd_y = .1,
                                   axis_ticks_lwd_x = .1,
                                   axis_label_size = 7,
                                   axis_label_face = "bold"),
                  legend_title_size = 10, 
                  legend_key_size = 0.8, 
                  legend_text_size = 7)
  dev.off()
  
  damageLogo_five(topic_clus$theta, renyi_alpha=100, xaxis_fontsize=15,
                  y_fontsize=15, title_fontsize = 25, pop_names=paste0("Cluster ",1:dim(topic_clus$theta)[2]),
                  title_aligner = 14, output_dir = output_dir, output_width = 1000, output_height = 700, save_plot=TRUE)
  
  
}else{
  
  do.call(CountClust::StructureGGplot, structure_control)
  StructureGGplot(omega = omega,
                  annotation = annotation,
                  palette = cols1,
                  yaxis_label = "Moderns vs Ancients",
                  order_sample = FALSE,
                  figure_title = paste0("  StructurePlot: K=", dim(omega)[2],""),
                  axis_tick = list(axis_ticks_length = .1,
                                   axis_ticks_lwd_y = .1,
                                   axis_ticks_lwd_x = .1,
                                   axis_label_size = 7,
                                   axis_label_face = "bold"),
                  legend_title_size = 8, 
                  legend_key_size = 0.4, 
                  legend_text_size = 5)
  
  damageLogo_five(topic_clus$theta, renyi_alpha=100, xaxis_fontsize=15,
                  y_fontsize=15, title_fontsize = 25, pop_names=paste0("Cluster ",1:dim(topic_clus$theta)[2]),
                  title_aligner = 14, output_dir = output_dir, output_width = 1000, output_height = 700, save_plot = FALSE)
  
}





library(dplyr)
library(Logolas)
library(plyr)
library(grid)
library(gridBase)








i <- 3
dir = paste0(folders[i])
pattern = NULL
breaks = c(-1, seq(1,20,1))
flanking_bases = 1
