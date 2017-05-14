


#############  test aRchaic_cluster_beta_mutation_flank ######################


mat1 <- get(load("../data/Fu_2016/Fu_2016.rda"))
mat2 <- get(load("../data/moderns_lite/moderns_lite.rda"))

mat_pooled <- rbind(mat1, mat2)

mat_reduced <- filter_by_pos_pattern(mat_pooled, max_pos = 20, pattern = "C->T")

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


theta_pool <- topic_clus$theta
pos <- sapply(rownames(theta_pool), function(x) return(strsplit(x, "_")[[1]][2]))
rownames(theta_pool) <- pos
max_prob <- max(theta_pool)

theta_pool_mod <- theta_pool[match(1:max_pos, as.numeric(rownames(theta_pool))),]

for(l in 1:dim(theta_pool_mod)[2]){
    plot_graph_2(theta_pool_mod[,l], max_pos, max_prob,
                 main = paste0("mutation trend:", pattern),
                 plot.control)
}


plot_graph_2 <- function(probs, max_pos, max_prob, col="red",
                       cex=unit(1, "npc"), pch=unit(16,"npc"),
                       xlab="position", ylab="prob. of mutation",
                       main="",
                       cex.axis=unit(1, "npc"),
                       cex.main=unit(1, "npc")){
  # if (length(probs) != max_pos){
  #   stop(cat('probability vector must be of length ', max_pos))
  # }
  par(font.axis = 2)
  plot(as.numeric(names(probs)), probs/max_prob, xlim = c(0, max_pos), ylim=c(0,1),
       type = "b", xaxt = "n", yaxt = "n", cex = cex, pch=pch, col=col, main=main,
       cex.main=cex.main, ylab="", xlab="")
  axis(side = 1, at = floor(seq(1, max_pos, length.out=5)), cex.axis = cex.axis, lwd.ticks = 1, tck=-0.05,
       cex.lab=2, mgp=c(2.5, 0.5, 0))
  title(xlab = xlab, mgp=c(2.5,1.5,0), cex.lab=1.8)
  ylimit <- c(0.0, 0.5, 1.0)*max_prob
  axis(side = 2, at = c(0.0, 0.5, 1.0), labels = round(ylimit,2), cex.axis = cex.axis, lwd.ticks=1, tck=-0.05,
       cex.lab=2, mgp=c(2.5, 0.5, 0))
  title(ylab = ylab, mgp=c(2.5,1,0), cex.lab=1.8)
}
