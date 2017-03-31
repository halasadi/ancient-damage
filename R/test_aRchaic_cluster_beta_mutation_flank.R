

#############  test aRchaic_cluster_beta_mutation_flank ######################


mat1 <- get(load("../data/Fu_2016/Fu_2016.rda"))
mat2 <- get(load("../data/moderns_lite/moderns_lite.rda"))

mat_pooled <- rbind(mat1, mat2)

mat_reduced <- filter_by_mutation_flank(mat_pooled)
labs <- c(rep("fu", 44), rep("moderns", 50))
flanking_bases = 1

signature_set <- colnames(mat_reduced)
sig_split <- t(sapply(1:length(signature_set), function(x) return(strsplit(signature_set[x], "")[[1]][1:(flanking_bases+5)])))
new_sig_split <- matrix(0, dim(sig_split)[1], 3);
new_sig_split[,1] <- sig_split[,1]
new_sig_split[,2] <- sapply(1:length(signature_set), function(x) return(paste(sig_split[x,(flanking_bases+1):(flanking_bases+4)], collapse="")))
new_sig_split[,3] <- sig_split[,(flanking_bases+5)]


topics.control <- list()
topics.control.default <- list(bf = FALSE, kill = 2, ord = TRUE, verb = 1, admix = TRUE,
                               nbundles = 1, use_squarem = FALSE, init.adapt = FALSE, type = "full",
                               light = 1, method_admix = 1, sample_init = TRUE, tmax = 10000)
topics.control <- modifyList(topics.control.default, topics.control)

suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = mat_reduced, K=K, tol=tol, model = "independent", signatures = new_sig_split), topics.control)))

suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = mat_reduced, K=K, tol=tol, model = "full", signatures = NULL), topics.control)))


theta_pool <- topic_clus$theta
signature_set <- rownames(theta_pool)
sig_split <- do.call(rbind,
                     lapply(signature_set,
                            function(x) strsplit(as.character(x), split="")[[1]][1:(4+2*flanking_bases)]))
ncol_sig <- (4+2*flanking_bases)

if(flanking_bases%%1 != 0){
  stop("flanking bases not evenly distributed")
}


sub_pattern <- sapply(1:dim(sig_split)[1],
                      function(x) paste(sig_split[x,(flanking_bases+1):(flanking_bases+4)], collapse=""))
new_sig_split <- cbind(sig_split[,1:flanking_bases], sub_pattern, sig_split[,((ncol_sig - flanking_bases +1):ncol_sig)])
colnames(new_sig_split) = NULL

prop_patterns_list <- list()

theta <- dplyr::tbl_df(data.frame(theta_pool)) %>% dplyr::mutate(sig = signature_set) %>% dplyr::group_by(sig) %>% dplyr::summarise_each(funs(sum)) %>% as.data.frame()
rownames(theta) <-  theta[,1]
theta <- theta[,-1]


for(l in 1:dim(theta)[2]){
  prop_patterns_list[[l]] <- numeric();
  for(j in 1:ncol(new_sig_split)){
    temp <- tapply(theta[,l], factor(new_sig_split[,j], levels=c("A", "C", "G", "T",
                                                                 "C->T", "C->A", "C->G",
                                                                 "T->A", "T->C", "T->G")), sum)

    temp[is.na(temp)]=0
    prop_patterns_list[[l]] <- cbind(prop_patterns_list[[l]], temp)
    rownames(prop_patterns_list[[l]]) <- gsub("->", ">", rownames(prop_patterns_list[[l]]))

  }
}

inflation_factor = c(2,1,2)
ic <- damage.ic(prop_patterns_list, alpha=renyi_alpha, inflation_factor = inflation_factor)


tab <- prop_patterns_list[[1]]
colnames(tab) <- c("left \n flank", "mutation", "right \n flank")

logo.control <- list()
logo.control.default <- list(hist = TRUE, frame_width =1,
                             ic.scale = TRUE,
                             xaxis_fontsize = 10, xlab_fontsize = 15,
                             y_fontsize = 15, main_fontsize = 16, start = 0.001,
                             yscale_change = TRUE, pop_name = paste0("logo plot mutation + flanking bases"), xlab = "",
                             ylab = "information content", col_line_split = "grey80", scale0 = 0.01,
                             scale1 = 0.99, newpage = TRUE)

logo.control <- modifyList(logo.control.default, logo.control)


cols = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(RColorBrewer::brewer.pal, cols$maxcolors, rownames(cols)))

total_chars = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O",
                "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "zero", "one", "two",
                "three", "four", "five", "six", "seven", "eight", "nine", "dot", "comma",
                "dash", "colon", "semicolon", "leftarrow", "rightarrow")

cols = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(RColorBrewer::brewer.pal, cols$maxcolors, rownames(cols)))
color_code = sample(col_vector, length(total_chars), replace=FALSE)

color_code[c(1,3,7,20, 43)] <- c("green", "blue", "orange", "red", "gray")

set.seed(20)
color_profile <- list("type" = "per_symbol",
                      "col" = color_code)

png(paste0(output_dir, "logo_clus_", l, ".png"), width=output_width, height = output_height)
do.call(Logolas::logomaker, c(list(table=tab, ic = ic[,1], color_profile = color_profile),
                              logo.control))
dev.off()

output_dir <- "../utilities/structure_mutation_flank/"
