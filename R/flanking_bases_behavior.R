

###########  Understanding flanking bases behavior  ########################

folders <- "../data/moderns_lite/"
files <- list.files(folders, pattern = ".csv")
ll <- list()
ll[["left"]] <- list()
ll[["right"]] <- list()

for(k in 1:length(files)){
  file <- read.csv(paste0(folders, files[k]), header=FALSE)
  file1 <- file[which(file[,6] == "+"),]
  file2 <- file1[union(which(file1[,2] < 20 ), which(file1[,3] < 20)),]
  sig_split <- do.call(rbind, lapply(file2[,1], function(x) return(strsplit(as.character(x), "")[[1]])))

  flanking_bases <- 1
  subs <- apply(sig_split, 1, function(x) return(paste0(x[(flanking_bases+1):(flanking_bases+4)], collapse = "")))
  left_flank <- apply(sig_split, 1, function(x) return(paste0(x[flanking_bases], collapse = "")))
  right_flank <- apply(sig_split, 1, function(x) return(paste0(x[(flanking_bases+5)], collapse = "")))

  file3 <- dplyr::tbl_df(file2) %>% dplyr::mutate(subs = subs) %>%
    dplyr::mutate(left_flank = left_flank) %>%
    dplyr::mutate(right_flank = right_flank) %>% filter(subs == "C->A") %>%
    select(V7, right_flank) %>% group_by(right_flank)  %>% as.data.frame()

  ll[["right"]][[k]] <- tapply(file3[,1], file3[,2], sum)

  file4 <- dplyr::tbl_df(file2) %>% dplyr::mutate(subs = subs) %>%
    dplyr::mutate(left_flank = left_flank) %>%
    dplyr::mutate(right_flank = right_flank) %>% filter(subs == "C->A") %>%
    select(V7, left_flank) %>% group_by(left_flank)  %>% as.data.frame()

  ll[["left"]][[k]] <- tapply(file4[,1], file4[,2], sum)

  cat("we are at file", k, "\n")
}

ll1<- list()
ll1[["left"]] <- do.call(rbind, ll[["left"]])
ll1[["right"]] <- do.call(rbind, ll[["right"]])
rownames(ll1[["left"]]) <- paste0(files)
rownames(ll1[["right"]]) <- paste0(files)

save(ll1, file = "../utilities/left_right_base_composition_moderns_CtoA.rda")

ll1 <- get(load(file="../utilities/left_right_base_composition_moderns_CtoA.rda"))

left_flank_prop <- t(apply(ll1["left"][[1]], 1, function(x) return(x/sum(x))))
right_flank_prop <- t(apply(ll1["right"][[1]], 1, function(x) return(x/sum(x))))


labs <- rep(c("GBR", "FIN", "PUR", "CLM", "MXL"), each = 10)
levels <- unique(labs)

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(right_flank_prop))),
  tissue_label = factor(labs, levels = levels)
)

topic_cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
               "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
               "brown4","darkorchid","magenta","yellow", "azure1","azure4")

rownames(right_flank_prop) <- annotation$sample_id
omega <- as.matrix(right_flank_prop)
StructureGGplot(omega = omega,
                annotation = annotation,
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Moderns",
                order_sample = FALSE,
                figure_title = paste0("right flanking base composition"),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))



annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(left_flank_prop))),
  tissue_label = factor(labs, levels = levels)
)

topic_cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
               "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
               "brown4","darkorchid","magenta","yellow", "azure1","azure4")

rownames(left_flank_prop) <- annotation$sample_id
omega <- as.matrix(left_flank_prop)
StructureGGplot(omega = omega,
                annotation = annotation,
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Moderns",
                order_sample = FALSE,
                figure_title = paste0("left flanking base composition"),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))


##########################   ancient files   ###################################



folders <- c("../data/Fu_2016/", "../data/AnnaGosling2016data/", "../data/Sherpa_data/",
             "../data/Iceman/", "../data/Skoglund/", "../data/Sardinia2017/",
             "../data/Reich/", "../data/Pinhasi/", "../data/neanderthal/",
             "../data/Nepal/", "../data/Lazaridis/", "../data/Jones_2015/",
             "../data/Allentoft/")

folders_pool <- c()
files_pool <- c()

for(k in 1:length(folders)){
  files <- list.files(folders[k], pattern = ".csv")
  folders_pool <- c(folders_pool, rep(folders[k], length(files)))
  files_pool <- c(files_pool, files)
}

files_folders_pool <- paste0(folders_pool, files_pool)
files_folders_pool <- files_folders_pool[-72]

ll <- list()
ll[["left"]] <- list()
ll[["right"]] <- list()

for(k in 1:length(files_folders_pool)){
  file <- read.csv(paste0(files_folders_pool[k]), header=FALSE)
  file1 <- file[which(file[,6] == "+"),]
  file2 <- file1[union(which(file1[,2] < 20 ), which(file1[,3] < 20)),]
  sig_split <- do.call(rbind, lapply(file2[,1], function(x) return(strsplit(as.character(x), "")[[1]])))

  flanking_bases <- 1
  subs <- apply(sig_split, 1, function(x) return(paste0(x[(flanking_bases+1):(flanking_bases+4)], collapse = "")))
  left_flank <- apply(sig_split, 1, function(x) return(paste0(x[flanking_bases], collapse = "")))
  right_flank <- apply(sig_split, 1, function(x) return(paste0(x[(flanking_bases+5)], collapse = "")))

  file3 <- dplyr::tbl_df(file2) %>% dplyr::mutate(subs = subs) %>%
    dplyr::mutate(left_flank = left_flank) %>%
    dplyr::mutate(right_flank = right_flank) %>% filter(subs == "C->A") %>%
    select(V7, right_flank) %>% group_by(right_flank)  %>% as.data.frame()

  ll[["right"]][[k]] <- tapply(file3[,1], file3[,2], sum)

  file4 <- dplyr::tbl_df(file2) %>% dplyr::mutate(subs = subs) %>%
    dplyr::mutate(left_flank = left_flank) %>%
    dplyr::mutate(right_flank = right_flank) %>% filter(subs == "C->A") %>%
    select(V7, left_flank) %>% group_by(left_flank)  %>% as.data.frame()

  ll[["left"]][[k]] <- tapply(file4[,1], file4[,2], sum)

  cat("we are at file", k, "\n")
}

#ll2 <- list()
#ll2[["right"]] <- ll[["left"]]
#ll2[["left"]] <- ll[["right"]]


ll1<- list()
ll1[["left"]] <- do.call(rbind, lapply(ll[["left"]], function(x) return(x[match(c("A", "C", "G", "T"), names(x))])))
ll1[["right"]] <- do.call(rbind, lapply(ll[["right"]], function(x) return(x[match(c("A", "C", "G", "T"), names(x))])))
rownames(ll1[["left"]]) <- paste0(files_folders_pool)
rownames(ll1[["right"]]) <- paste0(files_folders_pool)


save(ll1, file = "../utilities/left_right_base_composition_ancient_CtoA.rda")

ll1 <- get(load(file = "../utilities/left_right_base_composition_ancient_CtoA.rda"))
dim(ll1[["left"]][[1]])

left_flank_prop <- t(apply(ll1["left"][[1]], 1, function(x) return(x/sum(x))))
right_flank_prop <- t(apply(ll1["right"][[1]], 1, function(x) return(x/sum(x))))

folders_pool <- folders_pool[-72]
labs <- as.vector(sapply(folders_pool, function(x) return(strsplit(x, "/")[[1]][3])))
levels <- unique(labs)

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(right_flank_prop))),
  tissue_label = factor(labs, levels = levels)
)

topic_cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
               "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
               "brown4","darkorchid","magenta","yellow", "azure1","azure4")

rownames(right_flank_prop) <- annotation$sample_id
omega <- as.matrix(right_flank_prop)
StructureGGplot(omega = omega,
                annotation = annotation,
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Ancients",
                order_sample = FALSE,
                figure_title = paste0("right flanking base composition"),
                axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))



annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(left_flank_prop))),
  tissue_label = factor(labs, levels = levels)
)

topic_cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
               "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
               "brown4","darkorchid","magenta","yellow", "azure1","azure4")

rownames(left_flank_prop) <- annotation$sample_id
omega <- as.matrix(left_flank_prop)
StructureGGplot(omega = omega,
                annotation = annotation,
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Ancients",
                order_sample = FALSE,
                figure_title = paste0("left flanking base composition"),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))
