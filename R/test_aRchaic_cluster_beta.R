

#############  test aRchaic_cluster_beta_mutation()  #############################


mat1 <- get(load("../data/Fu_2016/Fu_2016.rda"))
mat2 <- get(load("../data/moderns_lite/moderns_lite.rda"))

mat_pooled <- rbind(mat1, mat2)
labs <- c(rep("fu", 44), rep("moderns", 50))

aRchaic_cluster_beta_mutation(mat_pooled,
                              K=3,
                              tol=0.01,
                              labs = labs,
                              output_dir = "../utilities/structure-mutation/moderns_Fu/")


mat1 <- get(load("../data/Pinhasi/Pinhasi.rda"))
mat2 <- get(load("../data/moderns_lite/moderns_lite.rda"))
labs <- c(rep("Pinhasi", 12), rep("moderns", 50))

mat_pooled <- rbind(mat1, mat2)

aRchaic_cluster_beta_mutation(mat_pooled,
                              K=2,
                              tol=0.01,
                              labs = labs,
                              output_dir = "../utilities/structure-mutation/moderns_Pinhasi/")




mat_reduced <- filter_by_mutation(mat_pooled)

signature_set <- colnames(mat_reduced)

K=3
tol = 0.01

topics.control <- list()
topics.control.default <- list(bf = FALSE, kill = 2, ord = TRUE, verb = 1, admix = TRUE,
                               nbundles = 1, use_squarem = FALSE, init.adapt = FALSE, type = "full",
                               light = 1, method_admix = 1, sample_init = TRUE, tmax = 10000)
topics.control <- modifyList(topics.control.default, topics.control)

suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = mat_reduced, K=K, tol=tol, model = "independent", signatures = signature_set), topics.control)))

suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = mat_reduced, K=K, tol=tol, model = "full", signatures = NULL), topics.control)))




message ("Structure plot and Logo plot representations of clusters")

omega <- topic_clus$omega
annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs, levels = levels)
)

rownames(omega) <- annotation$sample_id


structure.control.default <- list(yaxis_label = "aRchaic pops",
                                  order_sample = FALSE,
                                  figure_title = paste0("  StructurePlot: K=", K,""),
                                  axis_tick = list(axis_ticks_length = .1,
                                                   axis_ticks_lwd_y = .1,
                                                   axis_ticks_lwd_x = .1,
                                                   axis_label_size = 10,
                                                   axis_label_face = "bold"),
                                  legend_title_size = 10,
                                  legend_key_size = 0.7,
                                  legend_text_size = 8)

structure.control <- list()
structure.control <- modifyList(structure.control.default, structure.control)

if(is.null(output_dir)){ output_dir <- paste0(getwd(),"/")}

structure_width = 5
structure_height = 8


omega <- topic_clus$omega
annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(labs, levels = levels)
)

output_dir <- "../utilities/structure-mutation/"
if(is.null(output_dir)){ output_dir <- paste0(getwd(),"/")}
plot.new()
grid.newpage()
p <- do.call(StructureGGplot, append(list(omega= omega,
                                     annotation = annotation,
                                     palette = topic_cols),
                                structure.control))
print(p)
ggplot2::ggsave(paste0(output_dir, "structure.png"), width = structure_width,
                height = structure_height)

graphics.off()


theta_pool <- topic_clus$theta
signature_set <- rownames(theta_pool)
theta <- dplyr::tbl_df(data.frame(theta_pool)) %>% dplyr::mutate(sig = signature_set) %>% dplyr::group_by(sig) %>% dplyr::summarise_each(funs(sum)) %>% as.data.frame()
rownames(theta) <-  theta[,1]
theta <- theta[,-1]

mutations <- as.character(sapply(rownames(theta), function(x) return(paste0(strsplit(x,"")[[1]][-2], collapse = ""))))
rownames(theta) <- mutations

tab <- matrix(theta[,1], nrow=length(theta[,1]))
rownames(tab) <- rownames(theta)
colnames(tab) <- "mutation"

logo.control <- list()
logo.control.default <- list(ic = NULL, hist = TRUE, frame_width =1,
                             ic.scale = FALSE,
                             xaxis_fontsize = 10, xlab_fontsize = 15,
                             y_fontsize = 15, main_fontsize = 16, start = 0.001,
                             yscale_change = TRUE, pop_name = paste0("logo plot mutation"), xlab = "",
                             ylab = "composition", col_line_split = "grey80", scale0 = 0.01,
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
do.call(Logolas::logomaker, c(list(table=tab, color_profile = color_profile),
                              logo.control))
dev.off()


